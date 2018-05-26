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




int main(int argc, char **argv) {
    int sample_idx = -1;
    int verbose = 0;

    std::string gene_file;
    std::string input_filename;

    bool bad_args = false;

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"samples", required_argument, 0, 's'},
        {"verbose", no_argument, &verbose, 1},
        {"version", no_argument, 0, 'V'},
        {0, 0, 0, 0}
    };

    int c;
    int option_index = 0;

    while ((c = getopt_long(argc, argv, "hs:vV", long_options,
                            &option_index)) != -1) {
        switch (c) {
            case 'h':
                print_help();
                return 0;

            case 's':
                sample_idx = std::stoi(optarg);
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

    if ((sample_idx > 0) && (format != 2)) {
        std::cerr << "\n[ERROR] Samples are not supported in format 0 or 1\n";
        print_help();
        return 1;
    }

    if (gzipped) {
        std::cout << "Compressed: Yes" << std::endl;
    } else {
        std::cout << "Compressed: No" << std::endl;
    }

    std::cout << "Format: " << format << std::endl;

    std::vector<std::string> haplotypes;
    std::vector<std::string> reads;
    std::vector<std::string> transcripts;
    std::vector<std::string> samples;
    std::vector<int> values;
    std::vector<int> col_ind;
    std::vector<int> row_ptr;
    std::vector<int> counts;

    std::ifstream infile;
    gzFile gzinfile;

    if (gzipped) {
        gzinfile = (gzFile) gzopen(input_filename.c_str(), "rb");
    } else {
        infile.open(input_filename, std::ios::binary);
    }

    int num_transcripts;
    int num_haplotypes;
    int num_alignments;
    int size;

    std::vector<char> buffer;

    if (format != 0 && format != 1 && format != 2) {
        std::cerr << "Binary input file is unknown format\n";
        return -1;
    }

    // skip format version
    readIntFromFile(gzinfile, infile);

    //load list of haplotype names
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

    //load list of transcript names
    num_transcripts = readIntFromFile(gzinfile, infile);
    transcripts.reserve(num_transcripts);

    int total_elements_lengths = num_transcripts * num_haplotypes;
    std::vector<double> transcript_lengths(total_elements_lengths);

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
        // load alignments, use default counts of 1 per alignment
        // format is "read_id transcript_id value"

        int read_id;
        int transcript_id;
        int value;

        // load list of read names
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

        num_alignments = readIntFromFile(gzinfile, infile);

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
        // load equivalence class, also includes counts

        int equivalence_id;
        int transcript_id;
        int value;
        int num_classes = readIntFromFile(gzinfile, infile);
        counts.resize(num_classes);

        std::cout << "EC: " << num_classes << std::endl;

        for (int j = 0; j < num_classes; j++) {
            int c = readIntFromFile(gzinfile, infile);
            counts.push_back(c);

            if (verbose) {
                std::cout << c << std::endl;
            }
        }


        //infile.read((char*)&row_offsets[0], (num_ecs + 1) * sizeof(int));
        //infile.read((char*)&columns[0], nnz * sizeof(int));
        //infile.read((char*)&data[0], nnz * sizeof(int));
        /*

        num_alignments = readIntFromFile(gzinfile, infile);

        std::cout << "Alignments: " << num_alignments << std::endl;

        for (int i = 0; i < num_alignments; ++i) {
            equivalence_id = readIntFromFile(gzinfile, infile);
            transcript_id = readIntFromFile(gzinfile, infile);
            value = readIntFromFile(gzinfile, infile);

            if (verbose) {
                std::cout << equivalence_id << "\t" << transcript_id << "\t" << value << std::endl;
            }
        }
         */


        /*
        int idx = 0;

        std::vector<int> row_data(num_transcripts, 0);

        for (int row = 0; row < num_rowptr - 1; ++row) {
            //std::cout << "ROW=" << row << std::endl;
            row_data.clear();
            std::fill(row_data.begin(), row_data.begin() + num_transcripts, 0.0);

            bool found = false;

            for (int j = 0; j < row_offsets[row+1] - row_offsets[row]; ++j) {
                row_data[columns[idx]] = data[idx];
                //std::cout << "j=" << j << std::endl;
                //std::cout << "idx=" << idx << std::endl;
                //if (n_columns[idx]  == sample_idx) {
                //    column_data[row] = n_data[idx];
                //}

                if (data[idx] != 0) {
                    int _column = columns[idx];
                    std::cout << counts[row_offsets[row] << "\t" << transcripts[_column] << "\t" << data[idx] << std::endl;
                }

                ++idx;
            }

            //for (int i = 0; i < num_transcripts; i++) {
            //    std::cout << row_data[i] << "\t";
            //}


            //std::cout << std::endl;
        }
         */

        //read "A" matrix
        int num_rowptr = readIntFromFile(gzinfile, infile);
        std::cout << "A MATRIX ROW PTR: " << num_rowptr << std::endl;

        // read nnz
        int nnz = readIntFromFile(gzinfile, infile);
        std::cout << "A MATRIX NNZ: " << nnz << std::endl;

        std::vector<int> row_offsets(num_rowptr);
        std::vector<int> columns(nnz);
        std::vector<int> data(nnz);


        if (verbose) {
            std::cout << "A MATRIX INDPTR" << std::endl;
        }

        // TODO: better way to do this, but need to make utility function
        for (int j = 0; j < num_rowptr; j++) {
            int r = readIntFromFile(gzinfile, infile);
            row_offsets[j] = r;

            if (verbose) {
                std::cout << "#" << j << ": " << r << std::endl;
            }
        }

        if (verbose) {
            std::cout << "A MATRIX INDICES" << std::endl;
        }

        for (int j = 0; j < nnz; j++) {
            int c = readIntFromFile(gzinfile, infile);
            columns[j] = c;

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
    } else if (format == 2) {
        // read crs
        int num_crs = readIntFromFile(gzinfile, infile);
        samples.reserve(num_crs);

        std::cout << "CRS: " << num_crs << std::endl;

        for (int i = 0; i < num_crs; ++i) {
            size = readIntFromFile(gzinfile, infile);
            buffer.clear();

            for (int j = 0; j < size; j++) {
                char c = readCharFromFile(gzinfile, infile);
                buffer.push_back(c);
            }
            buffer.push_back('\0');
            samples.push_back(std::string(buffer.data()));

            if (verbose) {
                std::cout << std::string(buffer.data()) << std::endl;
            }
        }

        if (sample_idx > num_crs) {
            std::cerr << "ERROR: sample index must be between 0 and " << num_crs - 1 << "\n";
            return -1;
        }

        //read "N" matrix

        // read ec
        int num_rowptr = readIntFromFile(gzinfile, infile);
        std::cout << "N MATRIX ROW PTR: " << num_rowptr << std::endl;

        // read nnz
        int nnz = readIntFromFile(gzinfile, infile);
        std::cout << "N MATRIX NNZ: " << nnz << std::endl;

        std::vector<int> row_offsets;
        std::vector<int> columns;
        std::vector<int> data;

        row_offsets.reserve(num_rowptr);
        columns.reserve(nnz);
        data.reserve(nnz);

        //infile.read((char*)&row_offsets[0], (num_ecs + 1) * sizeof(int));
        //infile.read((char*)&columns[0], nnz * sizeof(int));
        //infile.read((char*)&data[0], nnz * sizeof(int));

        if (verbose) {
            std::cout << "N MATRIX INDPTR" << std::endl;
        }

        // TODO: better way to do this, but need to make utility function
        for (int j = 0; j < num_rowptr; j++) {
            int r = readIntFromFile(gzinfile, infile);
            row_offsets.push_back(r);

            if (verbose) {
                std::cout << "#" << j << ": " << r << std::endl;
            }
        }

        if (verbose) {
            std::cout << "N MATRIX INDICES" << std::endl;
        }

        for (int j = 0; j < nnz; j++) {
            int c = readIntFromFile(gzinfile, infile);
            columns.push_back(c);

            if (verbose) {
                std::cout << "#" << j << ": " << c << std::endl;
            }
        }

        if (verbose) {
            std::cout << "N MATRIX DATA" << std::endl;
        }

        for (int j = 0; j < nnz; j++) {
            int d = readIntFromFile(gzinfile, infile);
            data.push_back(c);

            if (verbose) {
                std::cout << "#" << j << ": " << d << std::endl;
            }
        }

        //read "A" matrix
        num_rowptr = readIntFromFile(gzinfile, infile);
        std::cout << "A MATRIX ROW PTR: " << num_rowptr << std::endl;

        // read nnz
        nnz = readIntFromFile(gzinfile, infile);
        std::cout << "A MATRIX NNZ: " << nnz << std::endl;

        row_offsets = std::vector<int>(num_rowptr);
        columns = std::vector<int>(nnz);
        data = std::vector<int>(nnz);


        //infile.read((char*)&row_offsets[0], (num_ecs + 1) * sizeof(int));
        //infile.read((char*)&columns[0], nnz * sizeof(int));
        //infile.read((char*)&data[0], nnz * sizeof(int));

        if (verbose) {
            std::cout << "A MATRIX INDPTR" << std::endl;
        }

        // TODO: better way to do this, but need to make utility function
        for (int j = 0; j < num_rowptr; j++) {
            int r = readIntFromFile(gzinfile, infile);
            row_offsets[j] = r;

            if (verbose) {
                std::cout << "#" << j << ": " << r << std::endl;
            }
        }

        if (verbose) {
            std::cout << "A MATRIX INDICES" << std::endl;
        }

        for (int j = 0; j < nnz; j++) {
            int c = readIntFromFile(gzinfile, infile);
            columns[j] = c;

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
              << "OPTIONS\n"
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
