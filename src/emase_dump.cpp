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
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <zlib.h>

#include "alignment_incidence_matrix.h"
//#include "sample_allelic_expression.h"
#include "alignment_import.h"

#define VERSION "0.3.0"

void print_help();

// check the file magic's number to see if it is gzipped
bool isGZipped(std::string filename) {
    std::ifstream infile(filename, std::ios::in|std::ios::binary);

    if (!infile.is_open()) {
        // something went wrong reading from stream for now return NULL
        std::cerr << "ERROR LOADING FILE " << filename << std::endl;
        return -1;
    }

    char bytes[2];

    infile.read(bytes, 2);

    if ((bytes[0] == '\x1F') && (bytes[1] == '\x8B')) {
        return true;
    }

    return false;
}


// get the format version number
int getBinFormat(std::string filename) {
    bool gz = isGZipped(filename);
    int format = -1;

    if (gz) {
        gzFile infile = (gzFile) gzopen(filename.c_str(), "rb");
        gzrewind(infile);

        int len = gzread(infile, &format, sizeof(format));

        gzclose(infile);
    } else {
        std::ifstream infile(filename, std::ios::in|std::ios::binary);

        if (!infile.is_open()) {
            // something went wrong reading from stream for now return NULL
            std::cerr << "ERROR LOADING FILE " << filename << std::endl;
            return -1;
        }

        infile.read((char *) &format, sizeof(int));
    }

    if (format < 0 || format > 2) {
        std::cerr << "Unable to determine file format, may be corrupt!" << filename << std::endl;
        return -1;
    }

    return format;
}


int readIntFromFile(gzFile gzinfile, std::ifstream &infile) {
    int i;

    if (infile.is_open()) {
        infile.read((char*)&i, sizeof(i));
    } else {
        gzread(gzinfile, &i, sizeof(i));
    }

    return i;
}

char readCharFromFile(gzFile gzinfile, std::ifstream &infile) {
    char c;

    if (infile.is_open()) {
        infile.read(&c, sizeof(c));
    } else {
        gzread(gzinfile, &c, sizeof(c));
    }

    return c;
}

void fileSeek(gzFile gzinfile, std::ifstream &infile, long pos) {
    if (infile.is_open()) {
        infile.seekg(pos);
    } else {
        gzrewind(gzinfile);
        gzseek(gzinfile, pos, SEEK_CUR);
    }
}

long fileTell(gzFile gzinfile, std::ifstream &infile) {
    if (infile.is_open()) {
        return (long)infile.tellg();
    } else {
        return (long)gztell(gzinfile);
    }
}




int main(int argc, char **argv) {
    int verbose = 0;
    int sample_start = -1;
    int sample_end = -1;

    std::string input_filename;
    std::string samples_str;
    std::string samples_names_str;

    bool bad_args = false;
    bool find_sample_idx = false;

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"names", required_argument, 0, 'n'},
        {"samples", required_argument, 0, 's'},
        {"verbose", no_argument, &verbose, 1},
        {"version", no_argument, 0, 'V'},
        {0, 0, 0, 0}
    };

    int c;
    int option_index = 0;

    while ((c = getopt_long(argc, argv, "hn:s:vV", long_options,
                            &option_index)) != -1) {
        switch (c) {
            case 'h':
                print_help();
                return 0;

            case 'n':
                samples_names_str = std::string(optarg);
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

    std::cout << "\nemase-dump Version " << VERSION << std::endl <<std::endl;
    std::cout << "Input File: " << input_filename << std::endl;

    int gzipped = isGZipped(input_filename);
    int format = getBinFormat(input_filename);

    if (gzipped) {
        std::cout << "Compressed: Yes" << std::endl;
    } else {
        std::cout << "Compressed: No" << std::endl;
    }

    std::cout << "Format: " << format << std::endl;

    if (samples_names_str.length() > 0) {
        if (format != 2) {
            std::cerr << "\n[ERROR] Samples are not supported in format 0 or 1\n";
            print_help();
            return 1;
        }

        if (samples_str.length() > 0) {
            std::cerr << "\n[ERROR] Sample names and sample indices are not supported at the same time\n";
            print_help();
            return 1;
        }

        find_sample_idx = true;
        std::cout << "Searching for sample: " << samples_names_str << std::endl;
    }


    if (samples_str.length() > 0) {
        if (format != 2) {
            std::cerr << "\n[ERROR] Samples are not supported in format 0 or 1\n";
            print_help();
            return 1;
        }

        if (samples_names_str.length() > 0) {
            std::cerr << "\n[ERROR] Sample names and sample indices are not supported at the same time\n";
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

    std::vector<std::string> haplotypes;
    std::vector<std::string> transcripts;
    std::vector<std::string> reads;
    std::vector<std::string> samples;

    std::vector<int> indptr;
    std::vector<int> indices;
    std::vector<int> data;

    std::vector<int> ecs;

    std::ifstream infile;
    gzFile gzinfile;

    if (gzipped) {
        gzinfile = (gzFile) gzopen(input_filename.c_str(), "rb");
    } else {
        infile.open(input_filename, std::ios::binary);
    }

    int num_haplotypes;
    int num_transcripts;
    int size;

    std::vector<char> buffer;

    if (format != 0 && format != 1 && format != 2) {
        std::cerr << "Binary input file is unknown format\n";
        return -1;
    }

    // skip format version since we already know
    readIntFromFile(gzinfile, infile);

    //
    //load list of haplotype names
    //
    num_haplotypes = readIntFromFile(gzinfile, infile);
    haplotypes.reserve(num_haplotypes);

    std::cout << "Haplotypes: " << num_haplotypes << std::endl;

    for (int i = 0; i < num_haplotypes; ++i) {
        size = readIntFromFile(gzinfile, infile);
        buffer.clear();

        for (int j = 0; j < size; ++j) {
            char c = readCharFromFile(gzinfile, infile);
            buffer.push_back(c);
        }
        buffer.push_back('\0');
        haplotypes.push_back(std::string(buffer.data()));

        if (verbose) {
            std::cout << std::string(buffer.data()) << std::endl;
        }
    }

    //
    // load list of transcripts
    //
    num_transcripts = readIntFromFile(gzinfile, infile);
    transcripts.reserve(num_transcripts);

    std::cout << "Transcripts: " << num_transcripts << std::endl;

    int num_transcript_lengths = num_transcripts * num_haplotypes;
    std::vector<double> transcript_lengths(num_transcript_lengths);

    for (int i = 0; i < num_transcripts; ++i) {
        size = readIntFromFile(gzinfile, infile);
        buffer.clear();

        for (int j = 0; j < size; ++j) {
            char c = readCharFromFile(gzinfile, infile);
            buffer.push_back(c);
        }
        buffer.push_back('\0');
        transcripts.push_back(std::string(buffer.data()));

        if (verbose) {
            std::cout << std::string(buffer.data());
        }

        for (int h = 0; h < num_haplotypes; ++h) {
            int length = readIntFromFile(gzinfile, infile);
            transcript_lengths[(i * num_haplotypes) + h] = std::max((double)length, 1.0);

            if (verbose) {
                std::cout  << "\t" << length;
            }
        }

        if (verbose) {
            std::cout << std::endl;
        }
    }


    if (format == 0) {
        //
        // load list of reads
        //

        int num_reads = readIntFromFile(gzinfile, infile);
        reads.reserve(num_reads);

        std::cout << "Reads: " << num_reads << std::endl;

        for (int i = 0; i < num_reads; i++) {
            size = readIntFromFile(gzinfile, infile);
            buffer.clear();

            for (int j = 0; j < size; j++) {
                char c = readCharFromFile(gzinfile, infile);
                buffer.push_back(c);
            }
            buffer.push_back('\0');
            reads.push_back(std::string(buffer.data()));

            if (verbose) {
                std::cout << std::string(buffer.data()) << std::endl;
            }
        }

        //
        // load alignment tuples
        //

        int read_id;
        int transcript_id;
        int value;

        int num_alignments = readIntFromFile(gzinfile, infile);

        std::cout << "Alignments: " << num_alignments << std::endl;

        for (int i = 0; i < num_alignments; ++i) {
            read_id = readIntFromFile(gzinfile, infile);
            transcript_id = readIntFromFile(gzinfile, infile);
            value = readIntFromFile(gzinfile, infile);

            if (verbose) {
                std::cout << read_id << "\t" << transcript_id << "\t" << value << std::endl;
            }
        }
    } else if (format == 1) {

        /*
         * OLD WAY
         *
         *  infile.read((char*)&row_offsets[0], (num_ecs + 1) * sizeof(int));
         *  infile.read((char*)&columns[0], nnz * sizeof(int));
         *  infile.read((char*)&data[0], nnz * sizeof(int));
         *
         */

        //
        // load "A" matrix (csr format)
        //
        int num_indptr = readIntFromFile(gzinfile, infile);
        int nnz = readIntFromFile(gzinfile, infile);

        std::cout << "A MATRIX LENGTH INDPTR: " << num_indptr << std::endl;
        std::cout << "A MATRIX NNZ: " << nnz << std::endl;

        std::vector<int> indptr(num_indptr);
        std::vector<int> indices(nnz);
        std::vector<int> data(nnz);

        if (verbose) {
            std::cout << "A MATRIX INDPTR" << std::endl;
        }

        // TODO: better way to do this, but need to make utility function
        for (int j = 0; j < num_indptr; j++) {
            int r = readIntFromFile(gzinfile, infile);
            indptr[j] = r;

            if (verbose) {
                std::cout << "#" << j << ": " << r << std::endl;
            }
        }

        if (verbose) {
            std::cout << "A MATRIX INDICES" << std::endl;
        }

        for (int j = 0; j < nnz; j++) {
            int c = readIntFromFile(gzinfile, infile);
            indices[j] = c;

            if (verbose) {
                std::cout << "#" << j << ": " << c << std::endl;
            }
        }

        if (verbose) {
            std::cout << "A MATRIX DATA" << std::endl;
        }

        for (int j = 0; j < nnz; j++) {
            int d = readIntFromFile(gzinfile, infile);
            data[j] = d;

            if (verbose) {
                std::cout << "#" << j << ": " << d << std::endl;
            }
        }

        //
        // load list of ec (equivalence classes)
        //

        int num_ecs = readIntFromFile(gzinfile, infile);
        ecs.reserve(num_ecs);

        std::cout << "EC: " << num_ecs << std::endl;

        for (int j = 0; j < num_ecs; j++) {
            int c = readIntFromFile(gzinfile, infile);
            ecs.push_back(c);

            if (verbose) {
                std::cout << "#" << j << ": " << c << std::endl;
            }
        }

    } else if (format == 2) {

        //
        // load list of CRS (samples)
        //

        int num_crs = readIntFromFile(gzinfile, infile);
        samples.reserve(num_crs);

        std::cout << "CRS: " << num_crs << std::endl;

        std::string sample_name_tmp;
        int found_sample_index = -1;


        for (int i = 0; i < num_crs; ++i) {
            size = readIntFromFile(gzinfile, infile);
            buffer.clear();

            for (int j = 0; j < size; j++) {
                char c = readCharFromFile(gzinfile, infile);
                buffer.push_back(c);
            }
            buffer.push_back('\0');
            sample_name_tmp = std::string(buffer.data());

            samples.push_back(sample_name_tmp);

            if (find_sample_idx && (samples_names_str == sample_name_tmp)) {
                found_sample_index = i;
                sample_start = i;
                sample_end = sample_start;
            }

            if (verbose) {
                std::cout << sample_name_tmp << std::endl;
            }
        }

        if (find_sample_idx) {
            if (found_sample_index >= 0) {
                std::cout << "Sample Index of '" << samples_names_str << "':" << found_sample_index << std::endl;
            } else {
                std::cout << "Sample Index of '" << samples_names_str << "' NOT FOUND" << std::endl;
            }
        }

        //
        // load "A" matrix (csr format)
        //
        int num_indptr = readIntFromFile(gzinfile, infile);
        int nnz = readIntFromFile(gzinfile, infile);

        std::cout << "A MATRIX LENGTH INDPTR: " << num_indptr << std::endl;
        std::cout << "A MATRIX NNZ: " << nnz << std::endl;

        std::vector<int> indptr(num_indptr);
        std::vector<int> indices(nnz);
        std::vector<int> data(nnz);

        if (verbose) {
            std::cout << "A MATRIX INDPTR" << std::endl;
        }

        // TODO: better way to do this, but need to make utility function
        for (int j = 0; j < num_indptr; j++) {
            int r = readIntFromFile(gzinfile, infile);
            indptr[j] = r;

            if (verbose) {
                std::cout << "#" << j << ": " << r << std::endl;
            }
        }

        if (verbose) {
            std::cout << "A MATRIX INDICES" << std::endl;
        }

        for (int j = 0; j < nnz; j++) {
            int c = readIntFromFile(gzinfile, infile);
            indices[j] = c;

            if (verbose) {
                std::cout << "#" << j << ": " << c << std::endl;
            }
        }

        if (verbose) {
            std::cout << "A MATRIX DATA" << std::endl;
        }

        for (int j = 0; j < nnz; j++) {
            int d = readIntFromFile(gzinfile, infile);
            data[j] = d;

            if (verbose) {
                std::cout << "#" << j << ": " << d << std::endl;
            }
        }

        //
        // read "N" matrix (csc format)
        //
        int num_indptr_csc = readIntFromFile(gzinfile, infile);
        int nnz_csc = readIntFromFile(gzinfile, infile);

        std::cout << "N MATRIX IND PTR: " << num_indptr_csc << std::endl;
        std::cout << "N MATRIX NNZ: " << nnz_csc << std::endl;

        long n_tell = fileTell(gzinfile, infile);


        /**
         * To load ALL data from binary file into CSC representation.
         */
        if (verbose) {
            std::vector<int> indptr_csc(num_indptr_csc);
            std::vector<int> indices_csc(nnz_csc);
            std::vector<int> data_csc(nnz_csc);

            if (verbose) {
                std::cout << "N MATRIX INDPTR" << std::endl;
            }

            for (int j = 0; j < num_indptr_csc; j++) {
                int r = readIntFromFile(gzinfile, infile);
                indptr_csc.push_back(r);

                if (verbose) {
                    std::cout << "#" << j << ": " << r << std::endl;
                }
            }

            if (verbose) {
                std::cout << "N MATRIX INDICES" << std::endl;
            }

            for (int j = 0; j < nnz_csc; j++) {
                int c = readIntFromFile(gzinfile, infile);
                indices_csc.push_back(c);

                if (verbose) {
                    std::cout << "#" << j << ": " << c << std::endl;
                }
            }

            if (verbose) {
                std::cout << "N MATRIX DATA" << std::endl;
            }

            for (int j = 0; j < nnz_csc; j++) {
                int d = readIntFromFile(gzinfile, infile);
                data_csc.push_back(c);

                if (verbose) {
                    std::cout << "#" << j << ": " << d << std::endl;
                }
            }
        }


        //std::cout << "n_tell = " << n_tell << std::endl;
        if (sample_start <= -1) {
            return 0;
        }

        /**
         * To load only the information we need from the CSC matrix.
         */
        for (int sample_idx = sample_start; sample_idx < sample_end + 1; ++sample_idx) {
            // seek to start of N matrix
            fileSeek(gzinfile, infile, n_tell);

            /**
             * Determine the position in the file where we need to be
             * with regards to the indptr.
             *
             * n_tell + (sample_idx * sizeof(int))
             */

            int seek_pos = n_tell + (sample_idx * sizeof(int));
            //std::cout << "seek_pos=" << seek_pos << std::endl;
            fileSeek(gzinfile, infile, seek_pos);

            /**
             * Read indptr[sample_idx] and indptr[sample_idx + 1]
             */
            int idx_1 = readIntFromFile(gzinfile, infile);
            int idx_2 = readIntFromFile(gzinfile, infile);
            //std::cout << "idx_1:idx_2=" << idx_1 << ":" << idx_2 << std::endl;

            /**
             * Find the start on indices.
             *
             * n_tell + (len(indptr) * sizeof(int))
             */
            seek_pos = n_tell + (sizeof(int) * num_indptr_csc);

            /**
             * Determine the position in the file where we need to be based
             * upon:
             *     idx_1 = indptr[sample_idx]
             *     idx_2 = indptr[sample_idx + 1]
             *
             * n_tell + (idx_1 * sizeof(int))
             */
            seek_pos += (sizeof(int) * idx_1);
            //std::cout << "indices_pos ..." << seek_pos << std::endl;
            fileSeek(gzinfile, infile, seek_pos);

            /**
             * The values of indices[idx_1] through indices[idx_2] are
             * the row numbers that contain data in the CSC matrix.
             */
            std::vector<int> row_values(idx_2 - idx_1);
            //std::cout << "row_values.size()=" << row_values.size() << std::endl;
            for (int x = 0; x < (idx_2 - idx_1); x++) {
                int row = readIntFromFile(gzinfile, infile);
                row_values[x] = row;
                //std::cout << "x = " << x << ", row = " << row << std::endl;
            }

            /**
             * Find the start on data.
             *
             * n_tell + (len(indptr) * sizeof(int)) + (len(nnz) * sizeof(int))
             */
            seek_pos = n_tell + sizeof(int) * num_indptr_csc + sizeof(int) * nnz_csc;

            /**
             * Determine the position in the file where we need to be based
             * upon:
             *     idx_1 = indptr[sample_idx]
             *     idx_2 = indptr[sample_idx + 1]
             *
             * seek_pos + (idx_1 * sizeof(int))
             */
            seek_pos += (sizeof(int) * idx_1);
            //std::cout << "data_pos ..." << seek_pos << std::endl;
            fileSeek(gzinfile, infile, seek_pos);

            /**
             * The values of data[idx_1] through data[idx_2] are
             * the actual values in the CSC matrix.
             */
            std::vector<int> data(idx_2 - idx_1);
            //std::cout << "data.size()=" << data.size() << std::endl;
            for (int x = 0; x < (idx_2 - idx_1); x++) {
                int value = readIntFromFile(gzinfile, infile);
                data[x] = value;
                //std::cout << "x = " << x << ", value = " << value << std::endl;
            }


            std::cout << "Sample: " << sample_idx << ", " << samples[sample_idx] << std::endl;
            if (verbose) {
                std::cout << "Row\tValue" << std::endl;
            }
            int counter = 0;
            long sum = 0;
            for (int i = 0; i < row_values.size(); ++i) {
                if (verbose) {
                    std::cout << row_values[i] << "\t" << data[i] << std::endl;
                }
                sum += data[i];
                counter++;
            }

            std::cout << "Counter: " << counter << std::endl;
            std::cout << "Sum: " << sum << std::endl;
        }
    }

    return 0;
}


void print_help() {
    std::string title = "EMASE-Zero Version " VERSION " Help";

    std::cout << std::endl << std::endl
              << title << std::endl
              << std::string(title.length(), '-') << std::endl << std::endl
              << "USAGE: emase-dump [options] <alignment_incidence_file>\n\n"
              << "INPUT: Binary Alignment Incidence file prepared with alntools\n"
              << "       (https://churchill-lab.github.io/alntools/)\n\n"
              << "OPTIONS:\n"
              << "   --name (-s) <string>:\n"
              << "      Specify the sample name.\n\n"
              << "  --samples (-s) <string>:\n"
              << "      Specify the sample indices.  Either one number or in the format\n"
              << "      of <sample_start>:<sample_end>\n\n"
              << "  --verbose (-v):\n"
              <<"       Run in verbose mode\n\n"
              << "  --version:\n"
              << "      Print the version and exit\n\n"
              << "  --help (-h):\n"
              << "      Print this message and exit\n\n";
}
