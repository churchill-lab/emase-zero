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
//  populase.cpp
//
//
//  Created by Glen Beane on 8/20/14.
//
//

#include <iostream>

#include <getopt.h>

#include "alignment_incidence_matrix.h"
#include "sample_allelic_expression.h"
#include "python_interface.h"

void print_help();


int main(int argc, char **argv)
{
    clock_t t1, t2;
    float diff;

    int num_iterations;
    int max_iterations = 200;
    int read_length = 100;

    SampleAllelicExpression::model model = SampleAllelicExpression::MODEL_2;
    int m;

    std::string input_filename;
    std::string output_filename;
    std::string transcript_length_file;
    std::string extension = ".pcl.bz2";

    //getopt related variables
    opterr = 0;
    int c;
    int option_index = 0;

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"model", required_argument, 0, 'm'},
        {"output", required_argument, 0, 'o'},
        {"read-length", required_argument, 0, 'k'},
        {"transcript-lengths", required_argument, 0, 'l'},
        {0, 0, 0, 0}
    };

    while ((c = getopt_long(argc, argv, "hm:o:k:l:", long_options, &option_index)) != -1) {
        switch (c) {
            case 'h':
                print_help();
                return 0;

            case 'm':
                m = std::stoi(optarg);
                if (m < 1 || m > 4) {
                    std::cerr << "Invalid model number specified. Valid options: 1, 2, 3, 4\n";
                }

                switch (m) {
                    case 1:
                        model = SampleAllelicExpression::MODEL_1;
                        break;
                    case 2:
                        model = SampleAllelicExpression::MODEL_2;
                        break;
                    case 3:
                        std::cerr << "Model 3 is currently unimplemented, please specify a different model\n";
                        return 1;
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

        }
    }

    input_filename = argv[optind];

    if (extension.size() >= input_filename.size() || !std::equal(extension.rbegin(), extension.rend(), input_filename.rbegin())) {
        std::cerr << "Error, expected file with .pcl.bz2 extension. Input file should be prepared with bam_to_pcl.py script.\n";
        return 1;
    }

    if (output_filename.empty()) {
        //use default, based on the input file name but placed in current working directdory
        output_filename = input_filename.substr(0, input_filename.size() - extension.size()).append(".stacksum.tsv");

        //check to see if there was a path in the input file name.  If so, trim it off
        std::size_t found = output_filename.rfind('/');
        if (found != std::string::npos) {
            output_filename = output_filename.substr(found+1);
        }

    }

    PythonInterface pi = PythonInterface();
    if (pi.init()){
        std::cerr << "Error importing TranscriptHits Python module.\n";
        std::cerr << '\t' << pi.getErrorString() << std::endl;
        return 1;
    }

    std::cout << "Loading " << input_filename << ". This may take a while..." << std::endl;
    AlignmentIncidenceMatrix *aim = pi.load(input_filename);

    if (!aim) {
        std::cerr << "Error loading pcl file\n";
        std::cerr << '\t' << pi.getErrorString() << std::endl;
        return 1;
    }

    std::cout << "Loaded Pickled Alignment Incidence file " << input_filename << std::endl;
    std::vector<std::string> hap_names = aim->get_haplotype_names();

    std::cout << "File had the following haplotype names:\n";
    for (std::vector<std::string>::iterator it = hap_names.begin(); it != hap_names.end(); ++it) {
        std::cout << *it << "\t";
    }
    std::cout << std::endl;

    std::cout << aim->num_reads() << " reads loaded\n";
    std::cout << aim->num_transcripts() << " transcripts\n";

    std::cout << std::endl;

    if (!transcript_length_file.empty()) {
        std::cout << "Loading Transcript Lenght File " << transcript_length_file << std::endl;
        aim->loadTranscriptLengths(transcript_length_file);
    }

    if (model != SampleAllelicExpression::MODEL_4 && !aim->has_gene_mappings()) {
        std::cerr << "File does not contain transcript to gene mapping information.  Only normalization Model 4 can be used.\n";
        return 1;
    }

    t1 = clock();
    SampleAllelicExpression sae(aim, read_length);
    t2 = clock();
    diff = ((float)t2-(float)t1)/CLOCKS_PER_SEC;
    std::cout << "Time for initializing stack sum = " << diff << "s" << std::endl;


    num_iterations = 0;
    std::cout << "Beginning EM Iterations" << std::endl;
    t1 = clock();
    do {
        sae.update(model);
    } while (++num_iterations < max_iterations && !sae.converged());
    t2 = clock();

    diff = ((float)t2-(float)t1)/CLOCKS_PER_SEC;
    std::cout << "Time for " << num_iterations << " iterations = " << diff << "s\n";
    std::cout << "Time per iteration " << diff/num_iterations << "s\n";

    std::cout << "Saving results to " << output_filename << std::endl;
    sae.saveStackSums(output_filename);
    std::cout << "Done.\n";

    return 0;
}


void print_help()
{
    std::cout << "EMASE Help\n\n"
              << "USAGE: emase [options] <alignment_incidence>\n\n"
              << "INPUT: Alignment Incidence file prepared with bam_to_pcl.py script\n\n"
              << "OPTIONS\n"
              << "  --model (-m) : Specify normalization model (can be 1-4, default=2)\n"
              << "  --output (-o) : Specify filename for output file (default is input_basename.stacksum.tsv)\n"
              << "  --transcript-lengths (-l) : Filename for transcript length file. Format is \"transcript_id\tlength\"\n"
              << "  --read-length (-k) : Specify read length for use when apply transcript lengths\n"
              << "  --help (-h) : Print this message and exit\n\n";
}
