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

#include <fstream>
#include <iostream>
#include <string>
#include <zlib.h>

#include "alignment_incidence_matrix.h"

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


/* This function will read in a binary file produced by Matt Vincent's
   Kallisto exporter (https://github.com/churchill-lab/kallisto-export) and
   create an AlignmentIncidenceMatrix instance.  A pointer to the new aim
   object will be returned.  If there is a problem loading the file, NULL will
   be returned.

   This code assumes that the reads/equivalence classes are stored in order.
 */
AlignmentIncidenceMatrix *loadFromBin(std::string filename, int sample_idx = -1) {
    AlignmentIncidenceMatrix *aim = NULL;

    // the following vectors will hold the data we read in from the file,
    // the values, col_ind, and row_ptr vectors are the compressed row format
    // sparse matrix representation of the alignments.  counts are used for
    // equivalence class data, it is not used for read data
    std::vector<std::string> haplotypes;
    std::vector<std::string> reads;
    std::vector<std::string> transcripts;
    std::vector<std::string> samples;
    std::vector<int> values;
    std::vector<int> col_ind;
    std::vector<int> row_ptr;
    std::vector<int> counts;

    bool gzipped = isGZipped(filename);
    int format = getBinFormat(filename);

    std::ifstream infile;
    gzFile gzinfile;

    if (gzipped) {
        std::cout << "GZIPPPPPPED" << std::endl;
        gzinfile = (gzFile) gzopen(filename.c_str(), "rb");
    } else {
        std::cout << "NORMAL" << std::endl;
        infile.open(filename, std::ios::binary);
    }

    int num_transcripts;
    int num_haplotypes;
    int num_reads;
    int num_alignments;
    int size;

    std::vector<char> buffer;

    if (format != 0 && format != 1 && format != 2) {
        std::cerr << "Binary input file is unknown format\n";
        return NULL;
    }

    readIntFromFile(gzinfile, infile);

    //load list of transcript names
    num_transcripts = readIntFromFile(gzinfile, infile);
    transcripts.reserve(num_transcripts);

    for (int i = 0; i < num_transcripts; ++i) {
        size = readIntFromFile(gzinfile, infile);
        buffer.clear();

        for (int j = 0; j < size; ++j) {
            char c = readCharFromFile(gzinfile, infile);
            buffer.push_back(c);
        }
        buffer.push_back('\0');
        transcripts.push_back(std::string(buffer.data()));
    }

    //load list of haplotype names
    num_haplotypes = readIntFromFile(gzinfile, infile);
    haplotypes.reserve(num_haplotypes);

    for (int i = 0; i < num_haplotypes; ++i) {
        size = readIntFromFile(gzinfile, infile);
        buffer.clear();

        for (int j = 0; j < size; ++j) {
            char c = readCharFromFile(gzinfile, infile);
            buffer.push_back(c);
        }
        buffer.push_back('\0');
        haplotypes.push_back(std::string(buffer.data()));
    }

    if (format == 0) {
        // load alignments, use default counts of 1 per alignment
        // format is "read_id transcript_id value"

        int read_id;
        int transcript_id;
        int value;

        // load list of read names
        num_reads = readIntFromFile(gzinfile, infile);
        reads.reserve(num_reads);

        for (int i = 0; i < num_reads; i++) {
            size = readIntFromFile(gzinfile, infile);
            buffer.clear();

            for (int j = 0; j < size; j++) {
                char c = readCharFromFile(gzinfile, infile);
                buffer.push_back(c);
            }
            buffer.push_back('\0');
            reads.push_back(std::string(buffer.data()));
        }

        num_alignments = readIntFromFile(gzinfile, infile);

        values.reserve(num_alignments);
        col_ind.reserve(num_alignments);

        // first read values start at index 0
        row_ptr.push_back(0);
        int last_read = 0;

        for (int i = 0; i < num_alignments; ++i) {
            read_id = readIntFromFile(gzinfile, infile);
            transcript_id = readIntFromFile(gzinfile, infile);
            value = readIntFromFile(gzinfile, infile);

            // sanity check that read_id is not less than last_read
            if (read_id < last_read) {
                // this is a problem with the file
                std::cerr << "ERROR: binary input file must be sorted\n";
                return NULL;
            }

            values.push_back(value);
            col_ind.push_back(transcript_id);

            if (read_id != last_read) {
                // we've just transitioned to a new read, so we need to
                // record the starting index in row_ptr;
                row_ptr.push_back(i);
                last_read = read_id;
            }
        }
        row_ptr.push_back(num_alignments);
        aim = new AlignmentIncidenceMatrix(haplotypes, reads, transcripts, samples, col_ind, row_ptr, values);

    } else if (format == 1) {
        // load equivalence class, also includes counts

        int equivalence_id;
        int transcript_id;
        int value;
        int num_classes = readIntFromFile(gzinfile, infile);

        counts.resize(num_classes);
        infile.read((char*)&counts[0], num_classes*sizeof(int));

        for (int j = 0; j < num_classes; j++) {
            counts.push_back(readIntFromFile(gzinfile, infile));
        }

        num_alignments = readIntFromFile(gzinfile, infile);

        values.reserve(num_alignments);
        col_ind.reserve(num_alignments);

        int last_ec = 0;
        row_ptr.push_back(0);

        for (int i = 0; i < num_alignments; ++i) {
            equivalence_id = readIntFromFile(gzinfile, infile);
            transcript_id = readIntFromFile(gzinfile, infile);
            value = readIntFromFile(gzinfile, infile);

            // sanity check that read_id is not less than last_read
            if (equivalence_id < last_ec) {
                // this is a problem with the file
                // XXX version 1 files are small enough that we could
                // probably read this all in and then sort it, and then
                // iterate over it and build the CRS representation.
                std::cerr << "ERROR: binary input file must be sorted\n";
                return NULL;
            }

            values.push_back(value);
            col_ind.push_back(transcript_id);

            if (equivalence_id != last_ec) {
                // record the start index when we transition to a new read
                row_ptr.push_back(i);
                last_ec = equivalence_id;
            }
        }
        row_ptr.push_back(num_alignments);
        aim = new AlignmentIncidenceMatrix(haplotypes, reads, transcripts, samples,
                                            col_ind, row_ptr, values,
                                            counts);
    } else if (format == 2) {
        // multisample
        std::cout << "Sample idx: " << sample_idx << std::endl;

        if (sample_idx == -1) {
            std::cerr << "ERROR: sample index must be greater than or equal to 0\n";
            return NULL;            
        }

        int equivalence_id;
        int transcript_id;
        int value;

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
        }

        if (sample_idx > num_crs) {
            std::cerr << "ERROR: sample index must be between 0 and " << num_crs - 1 << "\n";
            return NULL;
        }

        // read ec
        int num_ecs = readIntFromFile(gzinfile, infile);
        std::cout << "EC: " << num_ecs << std::endl;

        // read nnz
        int nnz = readIntFromFile(gzinfile, infile);
        std::cout << "NNZ: " << nnz << std::endl;

        //read "N" matrix 
        std::vector<int> row_offsets;
        std::vector<int> columns;
        std::vector<int> data;

        row_offsets.reserve(num_ecs + 1);
        columns.reserve(nnz);
        data.reserve(nnz);


        //infile.read((char*)&row_offsets[0], (num_ecs + 1) * sizeof(int));
        //infile.read((char*)&columns[0], nnz * sizeof(int));
        //infile.read((char*)&data[0], nnz * sizeof(int));

        // TODO: better way to do this, but need to make utility function
        for (int j = 0; j < num_ecs + 1; j++) {
            row_offsets.push_back(readIntFromFile(gzinfile, infile));
        }

        for (int j = 0; j < nnz; j++) {
            columns.push_back(readIntFromFile(gzinfile, infile));
        }

        for (int j = 0; j < nnz; j++) {
            data.push_back(readIntFromFile(gzinfile, infile));
        }

        // decipher csr to just column data
        std::vector<int> column_data(num_ecs, 0);

        int idx = 0;

        for (int row = 0; row < num_ecs; ++row) {
            //std::cout << "ROW=" << row << std::endl;
            for (int j = 0; j < row_offsets[row+1] - row_offsets[row]; ++j) {
                //std::cout << "j=" << j << std::endl;
                //std::cout << "idx=" << idx << std::endl;
                if (columns[idx]  == sample_idx) {
                    column_data[row] = data[idx];
                }
                //std::cout << column_data[row] << std::endl;
                ++idx;
            }
        }

        //read "A" matrix 
        int num_alignments = readIntFromFile(gzinfile, infile);
        std::cout << "NUM ALIGNMENTS: " << num_alignments << std::endl;

        values.reserve(num_alignments);
        col_ind.reserve(num_alignments);

        int last_ec = 0;
        row_ptr.push_back(0);

        for (int i = 0; i < num_alignments; ++i) {
            equivalence_id = readIntFromFile(gzinfile, infile);
            transcript_id = readIntFromFile(gzinfile, infile);
            value = readIntFromFile(gzinfile, infile);

            // sanity check that read_id is not less than last_read
            if (equivalence_id < last_ec) {
                // this is a problem with the file
                // XXX version 1 files are small enough that we could
                // probably read this all in and then sort it, and then
                // iterate over it and build the CRS representation.
                std::cerr << "ERROR: binary input file must be sorted\n";
                return NULL;
            }

            values.push_back(value);
            col_ind.push_back(transcript_id);

            if (equivalence_id != last_ec) {
                // record the start index when we transition to a new read
                row_ptr.push_back(i);
                last_ec = equivalence_id;
            }
        }
        row_ptr.push_back(num_alignments);

        std::cout << "Creating AlignmentIncidenceMatrix" << std::endl;
        
        aim = new AlignmentIncidenceMatrix(haplotypes, reads, transcripts, samples,
                                            col_ind, row_ptr, values,
                                            column_data);
    }

    return aim;
}
