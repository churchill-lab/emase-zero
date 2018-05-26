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

        aim = new AlignmentIncidenceMatrix(haplotypes, reads, transcripts, samples,
                                           col_ind, row_ptr, values, counts,
                                           transcript_lengths);

    } else if (format == 1) {
        int num_ec = readIntFromFile(gzinfile, infile);
        counts.reserve(num_ec);

        /*
        std::cout << "==" << std::endl;
        std::cout << "EC" << std::endl;
        std::cout << "==" << std::endl;
        std::cout << "#\tEC" << std::endl;
        */

        //infile.read((char*)&counts[0], num_classes*sizeof(int));
        std::cout << "EC\n";
        for (int j = 0; j < num_ec; j++) {
            int count = readIntFromFile(gzinfile, infile);
            counts.push_back(count);

            std::cout << "[" << j << "]\t" << counts[j] << std::endl;
        }

        // NEW WAY


        //read "A" matrix
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

        // OLD WAY
        int read_id;
        int transcript_id;
        int value;
        int num_alignments = readIntFromFile(gzinfile, infile);
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
         */

/*
//////////////
        std::cout << "===============" << std::endl;
        std::cout << "row_ptr.size()=" << row_ptr.size() << std::endl;;
        std::cout << "===============" << std::endl;
        std::cout << "#\tVALUE\tEC" << std::endl;;
        for (int i = 0; i < row_ptr.size(); i++) {
            std::cout << "[" << i << "]\t" << row_ptr[i] << "\t" << counts[row_ptr[i]] << std::endl;
        }
//////////////


//////////////
        std::cout << "===============" << std::endl;
        std::cout << "col_ind.size()=" << col_ind.size() << std::endl;;
        std::cout << "===============" << std::endl;
        std::cout << "#\tVALUE\tTRANSCRIPT" << std::endl;;
        for (int i = 0; i < col_ind.size(); i++) {
            std::cout << "[" << i << "]\t" << col_ind[i] << "\t" << transcripts[col_ind[i]] << std::endl;
        }
//////////////

//////////////
        std::cout << "===============" << std::endl;
        std::cout << "values.size()=" << values.size() << std::endl;;
        std::cout << "===============" << std::endl;
        std::cout << "#\tVALUE" << std::endl;;
        for (int i = 0; i < values.size(); i++) {
            std::cout << "[" << i << "]\t" << values[i] << std::endl;
        }
//////////////
*/

/*
 * for i in xrange(len(indptr) - 1):
    for j in xrange(indptr[i+1] - indptr[i]):
        if data[idx] != 0:
            print ec[i], targets[indices[idx]], data[idx]
        idx += 1

 */
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


        std::cout << "Creating AlignmentIncidenceMatrix" << std::endl;

        aim = new AlignmentIncidenceMatrix(haplotypes, reads, transcripts, samples,
                                           col_ind, row_ptr, values, counts,
                                           transcript_lengths);

    } else if (format == 2) {
        // multisample
        std::cout << "Sample idx: " << sample_idx << std::endl;

        if (sample_idx == -1) {
            std::cerr << "ERROR: sample index must be greater than or equal to 0\n";
            return NULL;            
        }

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

        //read "N" matrix

        int num_rowptr = readIntFromFile(gzinfile, infile);
        int nnz = readIntFromFile(gzinfile, infile);

        std::cout << "N MATRIX ROW PTR: " << num_rowptr << std::endl;
        std::cout << "N MATRIX NNZ: " << nnz << std::endl;

        std::vector<int> n_row_offsets(num_rowptr);
        std::vector<int> n_columns(nnz);
        std::vector<int> n_data(nnz);

        //infile.read((char*)&row_offsets[0], (num_ecs + 1) * sizeof(int));
        //infile.read((char*)&columns[0], nnz * sizeof(int));
        //infile.read((char*)&data[0], nnz * sizeof(int));

        // TODO: better way to do this, but need to make utility function
        for (int j = 0; j < num_rowptr; j++) {
            n_row_offsets[j] = readIntFromFile(gzinfile, infile);
        }

        for (int j = 0; j < nnz; j++) {
            n_columns[j] = readIntFromFile(gzinfile, infile);
        }

        for (int j = 0; j < nnz; j++) {
            n_data[j] = readIntFromFile(gzinfile, infile);
        }


        std::cout << "n_row_offsets" << std::endl;
        for (int i = 0; i < 10; i++) {
            std::cout << "[" << i << "]\t" << n_row_offsets[i] << std::endl;
        }
        for (int i = n_row_offsets.size() - 10; i < n_row_offsets.size(); i++) {
            std::cout << "[" << i << "]\t" << n_row_offsets[i] << std::endl;
        }

        std::cout << "n_columns" << std::endl;
        for (int i = 0; i < 10; i++) {
            std::cout << "[" << i << "]\t" << n_columns[i] << std::endl;
        }
        for (int i = n_columns.size() - 10; i < n_columns.size(); i++) {
            std::cout << "[" << i << "]\t" << n_columns[i] << std::endl;
        }

        std::cout << "n_data" << std::endl;
        for (int i = 0; i < 10; i++) {
            std::cout << "[" << i << "]\t" << n_data[i] << std::endl;
        }
        for (int i = n_data.size() - 10; i < n_data.size(); i++) {
            std::cout << "[" << i << "]\t" << n_data[i] << std::endl;
        }


        // decipher csr to just column data
        std::vector<int> column_data(num_rowptr - 1, 0);

        std::cout << "DEBUG\n";

        int idx = 0;
        for (int i = 0; i < num_rowptr - 1; ++i) {
            for (int j = 0; j < n_row_offsets[i+1] - n_row_offsets[i]; ++j) {

                if (n_columns[idx] == sample_idx) {
                    column_data[i] = n_data[idx];
                    //std::cout << column_data[i] << std::endl;
                }

                ++idx;
            }
        }



        std::cout << "end DEBUG\n";

        std::cout << "SIZE OF COLUMN_DATA=" << column_data.size() << std::endl;

        std::cout << "column_data" << std::endl;
        for (int i = 0; i < 20; i++) {
            std::cout << "[" << i << "]\t" << column_data[i] << std::endl;
        }
        for (int i = column_data.size() - 20; i < column_data.size(); i++) {
            std::cout << "[" << i << "]\t" << column_data[i] << std::endl;
        }

        //read "A" matrix
        num_rowptr = readIntFromFile(gzinfile, infile);
        nnz = readIntFromFile(gzinfile, infile);

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

        std::cout << "Creating AlignmentIncidenceMatrix" << std::endl;
        
        aim = new AlignmentIncidenceMatrix(haplotypes, reads, transcripts, samples,
                                           col_ind, row_ptr, values, column_data,
                                           transcript_lengths);
    }

    return aim;
}
