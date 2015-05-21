
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

#include <iostream>
#include "alignment_incidence_matrix.h"

AlignmentIncidenceMatrix::AlignmentIncidenceMatrix(std::vector<std::string> haplotypes,
                                                   std::vector<std::string> reads,
                                                   std::vector<std::string> transcripts,
                                                   std::vector<int> col_ind,
                                                   std::vector<int> row_ptr,
                                                   std::vector<int> val) : haplotype_names(haplotypes), transcript_names(transcripts), read_names(reads), col_ind(col_ind), row_ptr(row_ptr), val(val) {

    has_gene_mappings_ = false;

}


void AlignmentIncidenceMatrix::setGeneMappings(std::vector<int> tx_to_gene) {
    this->tx_to_gene = tx_to_gene;
    has_gene_mappings_ = true;
}

void AlignmentIncidenceMatrix::setGeneNames(std::vector<std::string> gene_names) {
    this->gene_names = gene_names;
    num_genes_ = (int)gene_names.size();
}