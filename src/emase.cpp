/*
 * Copyright (c) 2015 The Jackson Laboratory
 *
 * This software was developed by Gary Churchill's Lab at The Jackson
 * Laboratory (see http://research.jax.org/faculty/churchill).
 *
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software. If not, see <http://www.gnu.org/licenses/>.
 */


//
//  emase.cpp
//
//
//  Created by Glen Beane on 8/20/14.
//
//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <getopt.h>

#include "alignment_incidence_matrix.h"
#include "sample_allelic_expression.h"
#include "python_interface.h"
#include "kallisto_import.h"

#define VERSION "0.2.1"

void print_help();


int main(int argc, char **argv)
{
    AlignmentIncidenceMatrix *aim;

    bool binary_input = false;
    int verbose = 0;

    int num_iterations;
    int max_iterations = 999;
    int read_length = 100;

    SampleAllelicExpression::model model = SampleAllelicExpression::MODEL_1;
    int m;

    clock_t t1, t2;
    float diff;

    std::string input_filename;
    std::string output_filename;
    std::string transcript_length_file;
    std::string extension = ".pcl.bz2";
    std::string gene_file;

    //getopt related variables
    int c;
    int option_index = 0;
    bool bad_args = false;

    double tolerance = 0.0;

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"model", required_argument, 0, 'm'},
        {"output", required_argument, 0, 'o'},
        {"read-length", required_argument, 0, 'k'},
        {"transcript-lengths", required_argument, 0, 'l'},
        {"max-iterations", required_argument, 0, 'i'},
        {"bin", no_argument, 0, 'b'},
        {"gene-mappings", required_argument, 0, 'g'},
        {"verbose", no_argument, &verbose, 1},
        {"version", no_argument, 0, 'V'},
        {"tolerance", required_argument, 0, 't'},
        {0, 0, 0, 0}
    };

    while ((c = getopt_long(argc, argv, "hm:o:k:l:i:bg:vVc:", long_options,
                            &option_index)) != -1) {
        switch (c) {
            case 'h':
                print_help();
                return 0;

            case 'm':
                m = std::stoi(optarg);
                if (m < 1 || m > 4) {
                    std::cerr << "[ERROR] Invalid model number specified. Valid options: 1, 2, 3, 4\n";
                }

                switch (m) {
                    case 1:
                        model = SampleAllelicExpression::MODEL_1;
                        break;
                    case 2:
                        model = SampleAllelicExpression::MODEL_2;
                        break;
                    case 3:
                        model = SampleAllelicExpression::MODEL_3;
                        break;
                    case 4:
                        model = SampleAllelicExpression::MODEL_4;
                        break;
                }

                break;

            case 'o':
                output_filename = std::string(optarg);
                break;

            case 'k':
                read_length = std::stoi(optarg);
                break;

            case 'l':
                transcript_length_file = std::string(optarg);
                break;

            case 'i':
                max_iterations = std::stoi(optarg);
                break;

            case 'b':
                binary_input = true;
                break;

            case 'g':
                gene_file = std::string(optarg);
                break;

            case 'v':
                verbose = 1;
                break;

            case 'V':
               std::cout << VERSION << std::endl;
               return 0;

            case 't':
                tolerance = std::stod(optarg);
                break;

            case '?':
                bad_args = true;

        }
    }

    if (bad_args) {
        print_help();
        return 1;
    }


    if (argc - optind == 1) {
        input_filename = argv[optind];
    } else {
        std::cerr << "\n[ERROR] Missing required argument (input file name)\n";
        print_help();
        return 1;
    }

    if (!binary_input && (extension.size() >= input_filename.size()
                          || !std::equal(extension.rbegin(),
                                         extension.rend(),
                                         input_filename.rbegin()))) {
        std::cerr << "\n[ERROR] Expected file with .pcl.bz2 extension. Input file should be prepared with bam_to_pcl.py script.\n";
        return 1;
    }

    if (output_filename.empty()) {
        //use default, based on the input file name but placed in current working directdory

        if (!binary_input) {
            output_filename = input_filename.substr(0, input_filename.size() - extension.size()).append(".stacksum.tsv");
        } else {
           output_filename = input_filename;
           output_filename.append(".stacksum.tsv");
        }
        //check to see if there was a path in the input file name.  If so, trim it off
        std::size_t found = output_filename.rfind('/');
        if (found != std::string::npos) {
            output_filename = output_filename.substr(found+1);
        }

    }

    std::cout << "\nemase-zero Version " << VERSION << std::endl <<std::endl;

    std::cout << "Alignment File: " << input_filename << std::endl;

    if (!gene_file.empty()) {
        std::cout << "Grouping File: " << gene_file << std::endl;
    }
    else {
        std::cout << "Grouping File: None\n";
    }

    if (!transcript_length_file.empty()) {
        std::cout << "Transcript Length File: " << transcript_length_file
                  << std::endl;
    }
    else {
        std::cout << "Transcript Length File: None\n";
    }

    std::cout << "EM Model: " << model << std::endl
              << "Output File: " << output_filename << std::endl
              << "----------------------------------------------------\n\n\n";


    if (binary_input) {
        std::cout << "Loading " << input_filename << "..." << std::endl;
        aim = loadFromBin(input_filename);

        if (!aim) {
            std::cerr << "Error loading binary input file\n";
            return 1;
        }
    } else {
        PythonInterface pi = PythonInterface();
        if (pi.init()){
            std::cerr << "Error importing TranscriptHits Python module.\n";
            std::cerr << '\t' << pi.getErrorString() << std::endl;
            return 1;
        }

        std::cout << "Loading " << input_filename
                  << ". This may take a while..." << std::endl;
        aim = pi.load(input_filename);

        if (!aim) {
            std::cerr << "Error loading pcl file\n";
            std::cerr << '\t' << pi.getErrorString() << std::endl;
            return 1;
        }
    }


    std::vector<std::string> hap_names = aim->get_haplotype_names();

    if (verbose) {
        std::cout << "File had the following haplotype names:\n";
        for (auto it = hap_names.begin(); it != hap_names.end(); ++it) {
            std::cout << *it << "\t";
        }
        std::cout << std::endl;

        if (aim->has_equivalence_classes()) {
            std::cout << aim->num_alignment_classes()
                      << " alignment classes loaded ("
                      << aim->total_reads() << " total reads)\n";
        }
        else {
            std::cout << aim->total_reads() << " reads loaded\n";
        }
        std::cout << aim->num_transcripts() << " transcripts\n";

        std::cout << std::endl;
    }

    if (!transcript_length_file.empty()) {
        if (verbose) {
            std::cout << "Loading Transcript Length File "
                      << transcript_length_file << std::endl;
        }
        aim->loadTranscriptLengths(transcript_length_file);
    }

    if (!gene_file.empty()) {
        if (verbose) {
            std::cout << "Loading Gene Mapping File " << gene_file << std::endl;
        }
        aim->loadGeneMappings(gene_file);
    }

    if (model != SampleAllelicExpression::MODEL_4 && !aim->has_gene_mappings()) {
        if (!binary_input) {
            std::cerr << "[ERROR] File does not contain transcript to gene mapping information.  Only normalization Model 4 can be used.\n";
        }
        else {
            std::cerr << "[ERROR] Only model 4 is possible without gene to transcript information (--gene-mappings, -g). You specified another model.\n";
        }
        return 1;
    }

    t1 = clock();
    SampleAllelicExpression sae(aim, read_length, tolerance);
    t2 = clock();

    if (verbose) {
        diff = ((float)t2-(float)t1)/CLOCKS_PER_SEC;
        std::cout << "Time for initializing stack sum = " << diff << "s"
                  << std::endl;
    }

    if (max_iterations > 0) {
        num_iterations = 0;
        std::cout << "Beginning EM Iterations..." << std::endl;
        double change;
        double iteration_time;
        bool  converged;
        clock_t t1_inner, t2_inner;

        if (verbose) {
            std::cout << std::setw(7) << "Iter No" << "    "
                      << std::setw(8) << "Time(s)" << "    "
                      << std::setw(16) << "Change" << std::endl
                      << "-------    --------    ----------------\n";
        }

        t1 = clock();
        do {

            t1_inner = clock();
            sae.update(model);
            t2_inner = clock();
            converged = sae.converged(change);

            if (verbose) {
                std::cout << std::setw(7) << num_iterations + 1 << "    "
                          << std::setw(8) << std::setprecision(3)
                          << ((float)t2_inner - (float)t1_inner)/CLOCKS_PER_SEC
                          << "    " << std::setw(16) << std::setprecision(1)
                          << std::fixed << change << std::endl;
            }
        } while (++num_iterations < max_iterations && !converged);
        t2 = clock();

        diff = ((float)t2-(float)t1)/CLOCKS_PER_SEC;
        std::cout << "Time for " << num_iterations << " iterations = " << diff
                  << "s\n";
        std::cout << "Time per iteration " << std::setprecision(2) << diff/num_iterations << "s\n";
    }

    std::cout << "Saving results to " << output_filename << std::endl;
    sae.saveStackSums(output_filename);
    std::cout << "Done.\n";

    return 0;
}


void print_help()
{
    std::string title = "EMASE-Zero Version " VERSION " Help";

    std::cout << std::endl << std::endl
              << title << std::endl
              << std::string(title.length(), '-') << std::endl << std::endl
              << "USAGE: emase-zero [options] <alignment_incidence_file>\n\n"
              << "INPUT: Alignment Incidence file prepared with bam_to_pcl.py script or\n"
              << "       with kallisto-export (https://github.com/churchill-lab/kallisto-export)\n\n"
              << "OPTIONS\n"
              << "  --model (-m) <int>:\n"
              << "      Specify normalization model (can be 1-4, default=1)\n\n"
              << "  --bin (-b):\n"
              << "      Binary input mode, with this option emase2 will expect a binary\n"
              << "      file exported from Kallisto as input\n\n"
              << "  --output (-o) <filename>:\n"
              << "      Specify filename for output file (default is input_basename.stacksum.tsv)\n\n"
              << "  --transcript-lengths (-l) <filename>:\n"
              << "      Filename for transcript length file. Format of each line is\n"
              << "     \"TranscriptName_HaplotypeName\\tlength\"\n\n"
              << "  --gene-mappings (-g) <filename>:\n"
              << "      Filename containing transcript to gene mapping. Tab delimited file where\n"
              << "      the first field is the gene ID, all other fields are the transcript IDs\n"
              << "      that belong to that gene\n\n"
              << "  --read-length (-k) <int>:\n"
              << "      Specify read length for use when applying transcript length adjustment.\n"
              << "      Ignored unless combined with --transcript-lengths. (Default 100)\n\n"
              << "  --max-iterations (-i) <int>:\n"
              << "      Specify the maximum number of EM iterations. (Default 999)\n\n"
              << "  --tolerance (-t) <double>:\n"
              << "      Specify the convergence threshold. emase2 will terminate when\n"
              << "      the sum of the aboslute value of differences in the stack sum from\n"
              << "      one iteration to the next is lower than this value. (Default = 100\n"
              << "      when adjusting for transcript lengths, 0.001 * num_reads otherwise)\n\n"
              << "  --verbose (-v):\n"
              <<"       Run in verbose mode\n\n"
              << "  --version:\n"
              << "      Print the version and exit\n\n"
              << "  --help (-h):\n"
              << "      Print this message and exit\n\n";
}
