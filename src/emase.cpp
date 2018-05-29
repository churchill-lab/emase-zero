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

#include <cmath>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <stdio.h>

#include "alignment_incidence_matrix.h"
#include "sample_allelic_expression.h"
#include "alignment_import.h"

#define VERSION "0.3.0"

void print_help();


int main(int argc, char **argv) {
    AlignmentIncidenceMatrix *aim;
    SampleAllelicExpression::model model = SampleAllelicExpression::MODEL_1;

    clock_t t1, t2;
    float diff;

    int m;
    int max_iterations = 999;
    int num_iterations;
    int sample_start = 0;
    int sample_end = 0;
    int verbose = 0;

    double tolerance = 0.0001;

    std::string extension = ".pcl.bz2";
    std::string gene_file;
    std::string input_filename;
    std::string output_filename;
    std::string output_filename_counts;
    std::string samples_str;

    bool bad_args = false;

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"model", required_argument, 0, 'm'},
        {"outbase", required_argument, 0, 'o'},
        {"max-iterations", required_argument, 0, 'i'},
        {"gene-mappings", required_argument, 0, 'g'},
        {"samples", required_argument, 0, 's'},
        {"verbose", no_argument, &verbose, 1},
        {"version", no_argument, 0, 'V'},
        {"tolerance", required_argument, 0, 't'},
        {0, 0, 0, 0}
    };

    int c;
    int option_index = 0;

    while ((c = getopt_long(argc, argv, "hm:o:i:g:s:vVc:t:", long_options,
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
                output_filename_counts = std::string(optarg);
                output_filename.append(".tpm");
                output_filename_counts.append(".counts");
                break;

            case 'i':
                max_iterations = std::stoi(optarg);
                break;

            case 'g':
                gene_file = std::string(optarg);
                break;

            case 's':
                samples_str = std::string(optarg);
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

    if (output_filename.empty()) {
        //use default, based on the input file name but placed in current working directdory
        output_filename = "emase-zero.quantified.tpm";
        output_filename_counts = "emase-zero.quantified.counts";
    }

    std::cout << "\nemase-zero Version " << VERSION << std::endl <<std::endl;
    std::cout << "Alignment File: " << input_filename << std::endl;

    if (!gene_file.empty()) {
        std::cout << "Grouping File: " << gene_file << std::endl;
    } else {
        std::cout << "Grouping File: None\n";
    }

    int gzipped = isGZipped(input_filename);
    int format = getBinFormat(input_filename);

    std::cout << "Compressed: " << gzipped << "\n";
    std::cout << "Format: " << format << "\n";


    if (samples_str.length() > 0) {
        if (format != 2) {
            std::cerr << "\n[ERROR] Samples are not supported in format 0 or 1\n";
            print_help();
            return 1;
        }

        std::size_t location = samples_str.find(':');
        if (location > 0) {
            std::string s_start = samples_str.substr(0, location);
            std::string s_end = samples_str.substr(location + 1);

            sample_start = std::stoi(s_start);
            sample_end = std::stoi(s_end);
        } else {
            sample_start = std::stoi(samples_str);
            sample_end = sample_start;
        }
        
        std::cout << "Sample start: " << sample_start << std::endl;
        std::cout << "Sample end: " << sample_end << std::endl;
    }

    std::cout << "EM Model: " << model << std::endl
              << "Output File (TPM): " << output_filename << std::endl
              << "Output File (Counts): " << output_filename_counts << std::endl
              << "----------------------------------------------------\n\n\n";

    std::cout << "Loading " << input_filename << "..." << std::endl;
    aim = loadFromBin(input_filename);

    if (!aim) {
        std::cerr << "Error loading binary input file\n";
        return 1;
    }

    std::vector<std::string> hap_names = aim->get_haplotype_names();
    std::vector<std::string> sample_names = aim->get_sample_names();

    if (verbose) {
        std::cout << "File had the following haplotype names:\n";
        for (auto it = hap_names.begin(); it != hap_names.end(); ++it) {
            std::cout << "[" << *it << "]" << std::endl;
        }
        std::cout << std::endl;

        if (format != 2) {
            if (aim->has_equivalence_classes()) {
                std::cout << aim->num_alignment_classes()
                          << " alignment classes loaded ("
                          << aim->total_reads() << " total reads)\n";
            }
        }
        std::cout << aim->num_transcripts() << " transcripts\n";

        std::cout << std::endl;
    }

    if (!gene_file.empty()) {
        if (verbose) {
            std::cout << "Loading Gene Mapping File " << gene_file << std::endl;
        }
        aim->loadGeneMappings(gene_file);
    }

    if (model != SampleAllelicExpression::MODEL_4 && !aim->has_gene_mappings()) {
        std::cerr << "[ERROR] File does not contain transcript to gene mapping information.  Only normalization Model 4 can be used.\n";
        return 1;
    }

    // Loop through all the samples specified
    for (int i = sample_start; i < sample_end + 1; ++i) {
        std::cout << "SAMPLE = " << i << std::endl;

        if (i != sample_start) {
            //std::vector<int> hap_names = aim->get_haplotype_names();
        }

        if (format == 2) {
            loadNFromBin(input_filename, *aim, i);
        }

        if (verbose) {
            std::cout << aim->total_reads() << " reads loaded\n";
        }

        t1 = clock();
        std::cout<<"SampleAllelicExpression"<<std::endl;
        SampleAllelicExpression sae(aim, tolerance);
        std::cout<<"SampleAllelicExpression done"<<std::endl;
        t2 = clock();






        if (verbose) {
            diff = ((float)t2-(float)t1)/CLOCKS_PER_SEC;
            std::cout << "Time for initializing stack sum = " << diff << "s"
                      << std::endl;
        }

        if (max_iterations > 0) {
            num_iterations = 0;
            std::cout << "Running EM..." << std::endl;
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

            diff = ((float)t2 - (float)t1) / CLOCKS_PER_SEC;
            std::cout << "Time for " << num_iterations << " iterations = " << diff << "s\n";
            std::cout << "Time per iteration " << std::setprecision(2) << diff/num_iterations << "s\n";
        }

        std::cout << "Saving results..." << std::endl;

        // remove the file at the start, append each sample
        if (i == sample_start) {
            remove(output_filename.c_str());
            remove(output_filename_counts.c_str());
        }

        if (format == 2) {
            sae.saveStackSums(output_filename, sample_names[i]);
            sae.updateNoApplyTL(model);
            sae.saveStackSums(output_filename_counts, sample_names[i]);
            std::cout << "Done.\n";
        } else {
            sae.saveStackSums(output_filename);
            sae.updateNoApplyTL(model);
            sae.saveStackSums(output_filename_counts);
        }
    }

    std::cout << "Program finished.\n";

    return 0;
}


void print_help() {
    std::string title = "EMASE-Zero Version " VERSION " Help";

    std::cout << std::endl << std::endl
              << title << std::endl
              << std::string(title.length(), '-') << std::endl << std::endl
              << "USAGE: emase-zero [options] <alignment_incidence_file>\n\n"
              << "INPUT: Binary Alignment Incidence file prepared with alntools\n"
              << "       (https://churchill-lab.github.io/alntools/)\n\n"
              << "OPTIONS\n"
              << "  --model (-m) <int>:\n"
              << "      Specify normalization model (can be 1-4, default=1)\n\n"
              << "  --output (-o) <filename>:\n"
              << "      Specify filename for output file (default is input_basename.stacksum.tsv)\n\n"
              << "  --gene-mappings (-g) <filename>:\n"
              << "      Filename containing transcript to gene mapping. Tab delimited file where\n"
              << "      the first field is the gene ID, all other fields are the transcript IDs\n"
              << "      that belong to that gene\n\n"
              << "  --max-iterations (-i) <int>:\n"
              << "      Specify the maximum number of EM iterations. (Default 999)\n\n"
              << "  --samples (-s) <string>:\n"
              << "      Specify the sample indices.  Either one number or in the format\n"
              << "      of <sample_start>:<sample_end>\n\n"
              << "  --tolerance (-t) <double>:\n"
              << "      Specify the convergence threshold. emase will terminate when the\n"
              << "      sum of the aboslute value of differences in the stack sum from one\n"
              << "      iteration to the next is lower than this value multiplied by\n"
              << "      1,000,000 if performing length adjustment and multipled by #reads if\n"
              << "      not (Default = 0.0001)\n\n"
              << "  --verbose (-v):\n"
              <<"       Run in verbose mode\n\n"
              << "  --version:\n"
              << "      Print the version and exit\n\n"
              << "  --help (-h):\n"
              << "      Print this message and exit\n\n";
}
