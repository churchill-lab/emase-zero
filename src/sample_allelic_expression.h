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
//  sample_allelic_expression.h
//
//  Created by Glen Beane on 8/22/14.
//

#ifndef SAMPLE_ALLELIC_EXPRESSION_H
#define SAMPLE_ALLELIC_EXPRESSION_H


#include "alignment_incidence_matrix.h"
#include "utilities.h"


class SampleAllelicExpression {

public:

    enum model {MODEL_1 = 1, MODEL_2 = 2, MODEL_3 = 3, MODEL_4 = 4};

    SampleAllelicExpression(AlignmentIncidenceMatrix *alignment_incidence, double tolerance);
    ~SampleAllelicExpression();
    void update(model m = MODEL_1);
    void updateNoApplyTL(model m = MODEL_1);
    bool converged(double &change);
    inline int size() {return num_haplotypes * num_transcripts;}
    void saveStackSums(std::string filename);
    void saveStackSums(std::string filename, std::string sample_name);
    void applyTranscriptLength();

private:


    AlignmentIncidenceMatrix *alignment_incidence_;

    double *current_;
    double *previous_;
    double *working_;

    double *gene_sums_;
    double *gene_sum_hits_only_;
    double *transcript_sums_;
    double *gene_sum_by_strain_;
    double *gene_sum_by_strain_hits_only_;
    double *read_sum_by_gene_;

    int *gene_masks_;

    int num_haplotypes;
    int num_transcripts;

    double threshold_;

    void init();
    void init_normalize_read();

    void updateModel1();
    void updateModel2();
    void updateModel3();
    void updateModel4();

    DISALLOW_COPY_AND_ASSIGN(SampleAllelicExpression);
};

#endif
