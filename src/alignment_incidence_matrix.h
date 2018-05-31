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

#ifndef ALIGNMENT_INCIDENCE_MATRIX_H
#define ALIGNMENT_INCIDENCE_MATRIX_H

#include <string>
#include <vector>

#include "utilities.h"


// need for forward declare SampleAllelicExpression so we can declare friend
// relationship
class SampleAllelicExpression;

class AlignmentIncidenceMatrix {

public:
    AlignmentIncidenceMatrix(std::vector<std::string> haplotypes,
                             std::vector<std::string> transcripts,
                             std::vector<std::string> reads,
                             std::vector<int> col_ind,
                             std::vector<int> row_ptr,
                             std::vector<int> val);

    AlignmentIncidenceMatrix(std::vector<std::string> haplotypes,
                             std::vector<std::string> transcripts,
                             std::vector<std::string> reads,
                             std::vector<int> col_ind,
                             std::vector<int> row_ptr,
                             std::vector<int> val,
                             std::vector<int> counts,
                             std::vector<double> transcript_lengths);

    AlignmentIncidenceMatrix(std::vector<std::string> haplotypes,
                             std::vector<std::string> transcripts,
                             std::vector<std::string> reads,
                             std::vector<std::string> samples,
                             std::vector<int> col_ind,
                             std::vector<int> row_ptr,
                             std::vector<int> val,
                             std::vector<double> transcript_lengths);

    inline std::vector<std::string>::size_type num_haplotypes() {
        return haplotype_names.size();
    }

    inline std::vector<std::string>::size_type num_transcripts() {
        return transcript_names.size();
    }

    inline std::vector<std::string>::size_type num_samples() {
        return sample_names.size();
    }

    inline std::vector<std::string>::size_type num_alignment_classes() {
        return row_ptr.size() - 1;
    }

    inline std::vector<std::string>::size_type total_reads() {
        return total_reads_;
    }

    inline std::vector<std::string> get_haplotype_names() {
        return haplotype_names;
    }

    inline std::vector<std::string> get_sample_names() {
        return sample_names;
    }

    inline std::vector<std::string> get_gene_names() {
        return gene_names;
    }

    inline void setTranscriptLengths(int length) {
        std::fill(transcript_lengths_.begin(), transcript_lengths_.end(), 1);
    }


    void setGeneMappings(std::vector<int> tx_to_gene);

    void setGeneNames(std::vector<std::string> gene_names);

    bool loadGeneMappings(std::string filename);

    inline int num_genes() {
        return (int)gene_names.size();
    }

    inline bool has_gene_mappings() {
        return has_gene_mappings_;
    }

    inline bool has_equivalence_classes() {
        return has_equivalence_classes_;
    }

    void setCounts(std::vector<int> counts_);

    std::vector<int> getCounts() {
        return counts;
    }

    void setSampleFilter(int sample_idx_, std::vector<int> rows);

    inline void setNTell(long tell) {
        n_tell = tell;
    }

    inline long getNTell() {
        return n_tell;
    }


private:

    // SampleAllelicExpression needs access to the sparse matrix may just merge these classes...
    friend class SampleAllelicExpression;

    std::vector<std::string> haplotype_names;
    std::vector<std::string> transcript_names;
    std::vector<std::string> read_names;
    std::vector<std::string> sample_names;
    std::vector<std::string> gene_names;

    std::vector<int> tx_to_gene;
    std::vector<int> col_ind_orig;
    std::vector<int> row_ptr_orig;
    std::vector<int> val_orig;
    std::vector<int> col_ind;
    std::vector<int> row_ptr;
    std::vector<int> val;
    std::vector<int> gene_mapping;
    std::vector<int> counts;
    std::vector<int> row_filters;

    std::vector<double> transcript_lengths_;

    bool has_gene_mappings_;
    bool has_equivalence_classes_;

    long total_reads_;

    long n_tell;
    int sample_index;



    DISALLOW_COPY_AND_ASSIGN(AlignmentIncidenceMatrix);
};


#endif /* defined(ALIGNMENT_INCIDENCE_MATRIX_H) */
