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




/* This function will read in a binary file produced by Matt Vincent's
   Kallisto exporter (https://github.com/churchill-lab/kallisto-export) and
   create an AlignmentIncidenceMatrix instance.  A pointer to the new aim
   object will be returned.  If there is a problem loading the file, NULL will
   be returned.

   This code assumes that the reads/equivalence classes are stored in order.
 */
AlignmentIncidenceMatrix *loadFromBin(std::string filename) {
    AlignmentIncidenceMatrix *aim = NULL;

    // the following vectors will hold the data we read in from the file,
    // the values, col_ind, and row_ptr vectors are the compressed row format
    // sparse matrix representation of the alignments.  counts are used for
    // equivalence class data, it is not used for read data
    std::vector<int> values;
    std::vector<int> col_ind;
    std::vector<int> row_ptr;
    std::vector<int> counts;
    std::vector<std::string> reads;
    std::vector<std::string> samples;

    bool gzipped = isGZipped(filename);
    int format = getBinFormat(filename);

    std::ifstream infile;
    gzFile gzinfile;

    if (gzipped) {
        gzinfile = (gzFile) gzopen(filename.c_str(), "rb");
    } else {
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

    // skip format version
    readIntFromFile(gzinfile, infile);

    //load list of haplotype names
    num_haplotypes = readIntFromFile(gzinfile, infile);
    std::vector<std::string> haplotypes(num_haplotypes);

    /*
    std::cout << "==========" << std::endl;
    std::cout << "HAPLOTYPES" << std::endl;
    std::cout << "==========" << std::endl;
    std::cout << "#\tHAPLOTYPE" << std::endl;
     */

    for (int i = 0; i < num_haplotypes; ++i) {
        size = readIntFromFile(gzinfile, infile);
        buffer.clear();

        for (int j = 0; j < size; ++j) {
            char c = readCharFromFile(gzinfile, infile);
            buffer.push_back(c);
        }
        buffer.push_back('\0');
        haplotypes[i] = std::string(buffer.data());
        //std::cout << "[" << i << "]\t" << haplotypes[i] << std::endl;
    }

    //load list of transcript names
    num_transcripts = readIntFromFile(gzinfile, infile);
    std::vector<std::string> transcripts(num_transcripts);

    int total_elements_lengths = num_transcripts * num_haplotypes;
    std::vector<double> transcript_lengths(total_elements_lengths);

    std::cout << "===========" << std::endl;
    std::cout << "TRANSCRIPTS" << std::endl;
    std::cout << "===========" << std::endl;
    std::cout << "#\tTRANSCRIPT\tLENGTHS" << std::endl;

    for (int i = 0; i < num_transcripts; ++i) {
        size = readIntFromFile(gzinfile, infile);
        buffer.clear();

        for (int j = 0; j < size; ++j) {
            char c = readCharFromFile(gzinfile, infile);
            buffer.push_back(c);
        }
        buffer.push_back('\0');
        transcripts[i] = std::string(buffer.data());

        if (i < 10) {
            std::cout << "[" << i << "]\t" << transcripts[i];
        }

        for (int h = 0; h < num_haplotypes; ++h) {
            double length = (double)readIntFromFile(gzinfile, infile);
            if (i < 10) {
                std::cout << "\t" << length;
            }

            transcript_lengths[(i * num_haplotypes) + h] = std::max(length, 1.0);
        }

        if (i < 10) {
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

        std::cout << "Creating AlignmentIncidenceMatrix" << std::endl;

        aim = new AlignmentIncidenceMatrix(haplotypes, transcripts, reads,
                                           col_ind, row_ptr, values);

    } else if (format == 1) {
        //
        // read "A" matrix (CSR format)
        //
        int num_rowptr = readIntFromFile(gzinfile, infile);
        int nnz = readIntFromFile(gzinfile, infile);

        //std::cout << "A MATRIX ROW PTR: " << num_rowptr << std::endl;
        //std::cout << "A MATRIX NNZ: " << nnz << std::endl;

        row_ptr.reserve(num_rowptr);
        col_ind.reserve(nnz);
        values.reserve(nnz);

        //infile.read((char*)&row_offsets[0], (num_ecs + 1) * sizeof(int));
        //infile.read((char*)&columns[0], nnz * sizeof(int));
        //infile.read((char*)&data[0], nnz * sizeof(int));

        // TODO: better way to do this, but need to make utility function
        for (int j = 0; j < num_rowptr; j++) {
            row_ptr.push_back(readIntFromFile(gzinfile, infile));
        }

        for (int j = 0; j < nnz; j++) {
            col_ind.push_back(readIntFromFile(gzinfile, infile));
        }

        for (int j = 0; j < nnz; j++) {
            values.push_back(readIntFromFile(gzinfile, infile));
        }

        /*
        std::cout << "=============" << std::endl;
        std::cout << "HUMAN READOUT" << values.size() << std::endl;;
        std::cout << "=============" << std::endl;
        std::cout << "#\tEID\tTARGET\tVALUE" << std::endl;
        int idx = 0;
        for (int i = 0; i < num_rowptr - 1; ++i) {
            for (int j = 0; j < row_ptr[i+1] - row_ptr[i]; ++j) {
                std::cout << "[" << idx << "]\t";
                std::cout << counts[i] << "\t";
                std::cout << transcripts[col_ind[idx]] << "\t";
                std::cout << values[idx] << std::endl;
                ++idx;
            }
        }
        */

        int num_ec = readIntFromFile(gzinfile, infile);
        counts.reserve(num_ec);


        //infile.read((char*)&counts[0], num_classes*sizeof(int));
        std::cout << "EC\n";
        for (int j = 0; j < num_ec; j++) {
            int count = readIntFromFile(gzinfile, infile);
            counts.push_back(count);

            std::cout << "[" << j << "]\t" << counts[j] << std::endl;
        }


        std::cout << "Creating AlignmentIncidenceMatrix" << std::endl;

        aim = new AlignmentIncidenceMatrix(haplotypes, transcripts, reads,
                                           col_ind, row_ptr, values,
                                           counts, transcript_lengths);

    } else if (format == 2) {
        //
        // read samples (crs)
        //
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

        //
        // read "A" matrix (CSR matrix)
        //
        int num_rowptr = readIntFromFile(gzinfile, infile);
        int nnz = readIntFromFile(gzinfile, infile);

        std::cout << "A MATRIX ROW PTR: " << num_rowptr << std::endl;
        std::cout << "A MATRIX NNZ: " << nnz << std::endl;

        row_ptr.reserve(num_rowptr);
        col_ind.reserve(nnz);
        values.reserve(nnz);

        //infile.read((char*)&row_offsets[0], (num_ecs + 1) * sizeof(int));
        //infile.read((char*)&columns[0], nnz * sizeof(int));
        //infile.read((char*)&data[0], nnz * sizeof(int));

        // TODO: better way to do this, but need to make utility function
        for (int j = 0; j < num_rowptr; j++) {
            row_ptr.push_back(readIntFromFile(gzinfile, infile));
        }

        for (int j = 0; j < nnz; j++) {
            col_ind.push_back(readIntFromFile(gzinfile, infile));
        }

        for (int j = 0; j < nnz; j++) {
            values.push_back(readIntFromFile(gzinfile, infile));
        }

        /*
        std::cout << "row_ptr" << std::endl;
        for (int i = 0; i < 10; i++) {
            std::cout << "[" << i << "]\t" << row_ptr[i] << std::endl;
        }
        for (int i = row_ptr.size() - 10; i < row_ptr.size(); i++) {
            std::cout << "[" << i << "]\t" << row_ptr[i] << std::endl;
        }

        std::cout << "col_ind" << std::endl;
        for (int i = 0; i < 10; i++) {
            std::cout << "[" << i << "]\t" << col_ind[i] << std::endl;
        }
        for (int i = col_ind.size() - 10; i < col_ind.size(); i++) {
            std::cout << "[" << i << "]\t" << col_ind[i] << std::endl;
        }

        std::cout << "values" << std::endl;
        for (int i = 0; i < 10; i++) {
            std::cout << "[" << i << "]\t" << values[i] << std::endl;
        }
        for (int i = values.size() - 10; i < values.size(); i++) {
            std::cout << "[" << i << "]\t" << values[i] << std::endl;
        }
        */


        std::cout << "Creating AlignmentIncidenceMatrix" << std::endl;

        aim = new AlignmentIncidenceMatrix(haplotypes, transcripts, reads, samples,
                                           col_ind, row_ptr, values,
                                           transcript_lengths);

        aim->setNTell(fileTell(gzinfile, infile));
    }

    return aim;
}


void loadNFromBin(std::string filename, AlignmentIncidenceMatrix &aim, int sample_idx) {
    bool gzipped = isGZipped(filename);
    int format = getBinFormat(filename);

    std::ifstream infile;
    gzFile gzinfile;

    if (gzipped) {
        gzinfile = (gzFile) gzopen(filename.c_str(), "rb");
    } else {
        infile.open(filename, std::ios::binary);
    }


    if (format == 0) {
        std::cout << "N matrix is not in Format 0" << std::endl;
    } else if (format == 1) {
        std::cout << "N matrix is not in Format 0" << std::endl;
    } else if (format == 2) {
        //
        // read "N" matrix (csc format)
        //

        long n_tell = aim.getNTell();
        //std::cout << "n_tell = " << n_tell << std::endl;

        /**
         * To load only the information we need from the CSC matrix.
         */
        // seek to start of N matrix
        fileSeek(gzinfile, infile, n_tell);

        int num_indptr_csc = readIntFromFile(gzinfile, infile);
        int nnz_csc = readIntFromFile(gzinfile, infile);

        std::cout << "N MATRIX ROW PTR: " << num_indptr_csc << std::endl;
        std::cout << "N MATRIX NNZ: " << nnz_csc << std::endl;

        n_tell = fileTell(gzinfile, infile);

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

        aim.setSampleFilter(sample_idx, row_values);
        aim.setCounts(data);
    }

}