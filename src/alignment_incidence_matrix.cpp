
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
//  alignment_incidence_matrix.cpp
//
//
//  Created by Glen Beane on 8/20/14.
//

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

#include "alignment_incidence_matrix.h"

AlignmentIncidenceMatrix::AlignmentIncidenceMatrix(std::vector<std::string> haplotypes,
                                                   std::vector<std::string> reads,
                                                   std::vector<std::string> transcripts,
                                                   std::vector<std::string> samples,
                                                   std::vector<int> col_ind,
                                                   std::vector<int> row_ptr,
                                                   std::vector<int> val) :
    haplotype_names(haplotypes),
    transcript_names(transcripts),
    read_names(reads),
    sample_names(samples),
    col_ind(col_ind),
    row_ptr(row_ptr),
    val(val)
{
    has_gene_mappings_ = false;
    has_equivalence_classes_ = false;
    total_reads_ = row_ptr.size() - 1;
    counts = std::vector<int>(total_reads_, 1);
}

AlignmentIncidenceMatrix::AlignmentIncidenceMatrix(std::vector<std::string> haplotypes,
                                                   std::vector<std::string> reads,
                                                   std::vector<std::string> transcripts,
                                                   std::vector<std::string> samples,
                                                   std::vector<int> col_ind,
                                                   std::vector<int> row_ptr,
                                                   std::vector<int> val,
                                                   std::vector<int> counts) :
    haplotype_names(haplotypes),
    transcript_names(transcripts),
    read_names(reads),
    sample_names(samples),
    col_ind(col_ind),
    row_ptr(row_ptr),
    val(val),
    counts(counts)
{
    has_gene_mappings_ = false;
    has_equivalence_classes_ = true;

    total_reads_ = 0;
    for (auto it = counts.begin(); it != counts.end(); ++it) {
        total_reads_ += *it;
    }
}


void AlignmentIncidenceMatrix::setGeneMappings(std::vector<int> tx_to_gene) {
    this->tx_to_gene = tx_to_gene;
    has_gene_mappings_ = true;
}

void AlignmentIncidenceMatrix::setGeneNames(std::vector<std::string> gene_names) {
    this->gene_names = gene_names;
}

/*
    This function reads in a file with the format described below and
    initializes tx_to_gene (transcript ID to gene ID lookup table) and
    gene_names (gene ID to gene name table)

    format:
    A tab-delimited file where the first field is the gene name, and the
    remaining fields are the transcript names that are mapped to this gene.
    gene_name\ttranscript_name_0\t...transcript_name_n

    NOTE, this currently does not work with whitespace separators other than tab

    Returns true if file is loaded sucessfully, false otherwise.

    NOTE: for some error conditions we just print an error message and exit the
    program. Future versions may return false for those errors, and the caller
    will have to exit.
 */
bool AlignmentIncidenceMatrix::loadGeneMappings(std::string filename) {
    std::map<std::string, int> transcript_name_to_id;
    std::ifstream input(filename);

    if (!input.is_open()) {
        // something went wrong reading from stream for now just bail out
        std::cerr << "ERROR LOADING GENE MAPPING FILE " << filename << std::endl;
        exit(1);
    }

    // in most cases tx_to_gene and gene_names should already have size zero,
    // but it is possible that a user would pass a mapping file while they are
    // using a pcl.bz file that contained the mappings. In that case we'll load
    // the mappings from the file and ignore the ones that were present in the
    // pcl.bz file.
    // We initialize all transcripts in tx_to_gene to map to -1, after loading
    // the file, we will make sure there are no -1s in the table, which means
    // that that transcript name was not included in the mapping file.
    tx_to_gene.resize(num_transcripts(), -1);
    gene_names.clear();


    // build a map so we can lookup transcript IDs by transcript names
    int t_index = 0;
    for (auto it = transcript_names.begin(); it != transcript_names.end(); ++it) {
        transcript_name_to_id[*it] = t_index++;
    }

    // read file
    std::string line;
    int gene_id = 0;
    while (std::getline(input, line)) {
        // first split the line on tabs
        std::stringstream ss(line);
        std::string item;
        std::vector<std::string> tokens;

        while (std::getline(ss, item, '\t')) {
            if (!item.empty()) {
                tokens.push_back(item);
            }
        }

        // now we have a vector of tokens from this line.  The first token is
        // the gene name,  the other tokens are the trascripts that belong to
        // that gene
        gene_names.push_back(tokens[0]);
        for (int i = 1; i < tokens.size(); ++i) {
            // make sure the transcript name is one we know about:
            auto transcript_search = transcript_name_to_id.find(tokens[i]);

            if (transcript_search == transcript_name_to_id.end()) {
                // handle error condition.
                // For now we will just print an error message and exit.

                std::cerr << "ERROR LOADING GENE MAPPING FILE " << filename << std::endl
                          << "UNKNOWN TRANSCRIPT ID: " << tokens[i] << std::endl;

                exit(1);
            }

            // transcript is known, set its value to the gene_id
            tx_to_gene[transcript_search->second] = gene_id;

        }

        ++gene_id;
    }

    bool bail = false;
    for (int i = 0; i < tx_to_gene.size(); ++i) {
        if (tx_to_gene[i] == -1) {
            std::cerr << "No mapping information for " << transcript_names[i] << std::endl;
            bail = true;
        }
    }

    if (bail) {
        exit(1);
    }

    has_gene_mappings_ = true;
    return true;
}

void AlignmentIncidenceMatrix::loadTranscriptLengths(std::string filename) {
    std::ifstream input(filename);


    // since we don't want to assume that the order the transcripts appear in
    // the length file is the same as the order they are obtained from the input
    // file we go through the pain of
    std::map<std::string, int> transcript_name_to_id;
    std::map<std::string, int> haplotype_name_to_id;

    if (!input.is_open()) {
        // something went wrong reading from stream for now just bail out
        std::cerr << "ERROR LOADING TRANSCRIPT LENGTH FILE " << filename << std::endl;
        exit(1);
    }

    for (int i = 0; i < num_transcripts(); ++i) {
        transcript_name_to_id[transcript_names[i]] = i;
    }

    for (int i = 0; i < num_haplotypes(); ++i) {
        haplotype_name_to_id[haplotype_names[i]] = i;
    }

    int lengths_loaded = 0;

    transcript_lengths_.resize(num_transcripts() * num_haplotypes());

    int total_elements = num_transcripts() * num_haplotypes();

    double length;
    std::string t_name;
    std::string hap_name;
    std::string id;

    while (input >> id >> length) {

        // need to split t_name:  form read from file is transcriptName_haplotypeName
        std::stringstream id_stringstream(id);
        std::getline(id_stringstream, t_name, '_');
        std::getline(id_stringstream, hap_name, '_');

        if (++lengths_loaded > total_elements) {
            // lengths file contained more transcripts than we expected.
            // for now, just bail out

            std::cerr << "ERROR LOADING TRANSCRIPT LENGTH FILE " << filename << std::endl
                      << "EXPECTED " << total_elements << " LENGTH VALUES BUT FILE CONTAINS MORE\n";
            exit(1);

        }

        auto transcript_search = transcript_name_to_id.find(t_name);

        if (transcript_search == transcript_name_to_id.end()) {
            // handle error condition.  Should we abort?  Should we just
            // continue without transcript lengths? For now we will just print
            // an error message and exit.

            std::cerr << "ERROR LOADING TRANSCRIPT LENGTH FILE " << filename << std::endl
                      << "UNKNOWN TRANSCRIPT ID: " << t_name << std::endl;

            exit(1);
        }

        auto hap_search = haplotype_name_to_id.find(hap_name);

        if (hap_search == haplotype_name_to_id.end()) {
            std::cerr << "ERROR LOADING TRANSCRIPT LENGTH FILE " << filename << std::endl
                      << "UNKNOWN HAPLOTYPE NAME: " << hap_name << std::endl;
        }

        // it exists,  add its lenght to our list of transcript lengths
        transcript_lengths_[transcript_search->second * num_haplotypes() + hap_search->second] = std::max(length, 1.0);

    }

    if (lengths_loaded != total_elements) {
        // didn't have enough transcripts in file.  for now, just bail out.
        std::cerr << "ERROR LOADING TRANSCRIPT LENGHT FILE " << filename << std::endl
                  << "EXPECTED " << total_elements << " VALUES BUT FILE CONTAINS FEWER\n";
        exit(1);
    }

    if (input.bad()) {
        // something went wrong reading from stream for now just bail out
        std::cerr << "ERROR LOADING TRANSCRIPT LENGHT FILE " << filename << std::endl;
        exit(1);
    }
}
