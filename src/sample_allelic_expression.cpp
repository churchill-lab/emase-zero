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

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "sample_allelic_expression.h"


SampleAllelicExpression::SampleAllelicExpression(AlignmentIncidenceMatrix *alignment_incidence, double tolerance) :
    alignment_incidence_(alignment_incidence)
{
    num_haplotypes = (int)alignment_incidence->num_haplotypes();
    num_transcripts = (int)alignment_incidence->num_transcripts();

    current_ = new double[size()];
    previous_ = new double[size()];

    // This is allocated during init(), since it's size depends on the
    // maximum number of transcripts with hits in a single read
    working_ = NULL;

    // only used by some normalization modes, will allocate if required.
    transcript_sums_ = NULL;
    gene_sums_ = NULL;
    gene_sum_by_strain_ = NULL;
    gene_sum_by_strain_hits_only_ = NULL;
    gene_sum_hits_only_ = NULL;
    gene_masks_ = NULL;
    read_sum_by_gene_ = NULL;

    /*
    if (alignment_incidence->transcript_lengths_.size() / num_haplotypes == alignment_incidence->num_transcripts()) {
        // we are doing length adjustment, mutltiply threshold by 1M
        threshold_ = tolerance * 1000000;
    }
    else {
        // no length adjustment, multiply threshold by number of reads
        threshold_ = tolerance * alignment_incidence->total_reads();
    }
    */

    threshold_ = tolerance * alignment_incidence->total_reads();


    init();
}


SampleAllelicExpression::~SampleAllelicExpression() {
    delete [] current_;
    delete [] previous_;

    if (transcript_sums_) {
        delete [] transcript_sums_;
    }

    if (gene_sum_hits_only_) {
        delete [] gene_sum_hits_only_;
    }

    if (working_) {
        delete [] working_;
    }

    if (gene_sums_) {
        delete [] gene_sums_;
    }

    if (gene_sum_by_strain_) {
        delete [] gene_sum_by_strain_;
    }

    if (gene_sum_by_strain_hits_only_) {
        delete [] gene_sum_by_strain_hits_only_;
    }

    if (gene_masks_) {
        delete [] gene_masks_;
    }

    if (read_sum_by_gene_) {
        delete [] read_sum_by_gene_;
    }
}


void SampleAllelicExpression::init() {
    init_normalize_read();
    applyTranscriptLength();
}


void SampleAllelicExpression::init_normalize_read() {
    std::fill(current_, current_ + size(), 0.0);


    int working_size;
    int start_index;
    int free_slots;
    int work_index;
    double read_sum;
    auto end = alignment_incidence_->row_ptr.size() - 1;

    working_size = num_haplotypes * 4;
    free_slots = working_size;
    working_ = new double[working_size];

    //iterate over each read
    for (int i = 0; i != end; ++i) {
        work_index = 0;
        read_sum = 0.0;
        free_slots = working_size;

        // alignment_incidence_->row_ptr[i] gives us the start index for the
        // nozero values that make up this row
        // and there will be alignment_incidence_->row_ptr[i+1] - alignment_incidence_->row_ptr[i] values
        for (long j = alignment_incidence_->row_ptr[i]; j < alignment_incidence_->row_ptr[i+1]; ++j) {
            // make sure we haven't run out of memory in our working_ buffer
            // if it is full, double the capacity.  The size will always be a
            // multiple of num_haplotypes.
            if (free_slots == 0) {
                double *tmp = working_;
                working_ = new double[working_size * 2];
                std::copy(tmp, tmp + working_size, working_);

                // increased capacity by working_size, so we have working_size
                // free slots now (was full).
                free_slots = working_size;
                working_size *= 2;
                delete tmp;
            }

            // going to consume num_haplotypes more elements in working_ so
            // decrement the amount of free space
            free_slots -= num_haplotypes;

            // okay, we can do the actual work for this read
            for (int k = 0; k < num_haplotypes; ++k) {
                // initialize to 1.0 if they had a hit, 0 otherwise
                // each haplotype is represented by a single bit in the value
                if (alignment_incidence_->val[j] & (1 << k)) {
                    working_[work_index++] = 1.0;
                    read_sum += 1.0;
                }
                else {
                    working_[work_index++] = 0.0;
                }
            }
        }

        work_index = 0;
        for (long j = alignment_incidence_->row_ptr[i]; j < alignment_incidence_->row_ptr[i+1]; ++j) {
            start_index = alignment_incidence_->col_ind[j] * num_haplotypes;
            for (int k = 0; k < num_haplotypes; ++k) {
                current_[start_index + k] += (working_[work_index++] / read_sum) * (double)alignment_incidence_->counts[i];
            }
        }
    }
}


void SampleAllelicExpression::update(SampleAllelicExpression::model m) {
    switch (m) {

        case MODEL_1:
            updateModel1();
            break;

        case MODEL_2:
            updateModel2();
            break;

        case MODEL_3:
            updateModel3();
            break;

        case MODEL_4:
            updateModel4();
            break;
    }

    applyTranscriptLength();

}


void SampleAllelicExpression::updateNoApplyTL(SampleAllelicExpression::model m) {
    switch (m) {

        case MODEL_1:
            updateModel1();
            break;

        case MODEL_2:
            updateModel2();
            break;

        case MODEL_3:
            updateModel3();
            break;

        case MODEL_4:
            updateModel4();
            break;
    }
}


void SampleAllelicExpression::updateModel1() {
    int start_index;
    int work_index;
    double read_sum;
    auto end = alignment_incidence_->row_ptr.size() - 1;

    if (!gene_sums_) {
        gene_sums_ =  new double[alignment_incidence_->num_genes()];
    }

    // allocate some extra working memory needed by this model
    if (!gene_sum_by_strain_) {
        gene_sum_by_strain_ = new double[num_haplotypes * alignment_incidence_->num_genes()];
    }

    if (!gene_sum_hits_only_) {
        gene_sum_hits_only_ = new double[alignment_incidence_->num_genes()];
    }

    if (!gene_sum_by_strain_hits_only_) {
        gene_sum_by_strain_hits_only_ = new double[num_haplotypes * alignment_incidence_->num_genes()];
    }

    if (!gene_masks_) {
        gene_masks_ = new int[alignment_incidence_->num_genes()];
    }

    if (!transcript_sums_) {
        transcript_sums_ = new double[num_transcripts];
    }

    std::swap(current_, previous_);

    // clear current_ so we can use it to accumulate sums
    std::fill(current_, current_ + size(), 0.0);

    // generate gene_sum_by_strain from previous value
    std::fill(gene_sum_by_strain_, gene_sum_by_strain_ + alignment_incidence_->num_genes() * num_haplotypes, 0.0);
    for (int i = 0; i < num_transcripts; ++i) {
        for (int j = 0; j < num_haplotypes; ++j) {
            gene_sum_by_strain_[alignment_incidence_->tx_to_gene[i] * num_haplotypes + j] += previous_[i * num_haplotypes + j];
        }
    }

    // generate transcript sums from previous value
    for (int i = 0; i < num_transcripts; ++i) {
        transcript_sums_[i] = 0.0;
        for (int j = 0; j < num_haplotypes; ++j) {
            transcript_sums_[i] += previous_[i*num_haplotypes + j];
        }
    }

    // generate gene sums from previous value
    std::fill(gene_sums_, gene_sums_ + alignment_incidence_->num_genes(), 0.0);
    double genes_total = 0.0;
    for (int i = 0; i < num_transcripts; ++i) {
        gene_sums_[alignment_incidence_->tx_to_gene[i]] += transcript_sums_[i];
        genes_total += transcript_sums_[i];
    }

    //std::cout << "END: " << end << std::endl;
    //iterate over each read
    for (int i = 0; i != end; ++i) {
        read_sum = 0.0;

        /* unfortunately, given the information we have and the way it is stored in memory,  we need to make three passes over
           this read computing some values we need for the actual update, which is yet another pass over this read, followed by
           a read-normalization step giving us a total of 5 iterations over the read.  It would be
           possible to come up with another way to store the data that might reduce the number of passes we need.  Initial
           benchmarking gave us around a 4s update time for a test file with 8 strains, ~30K genes, ~95K transcripts, and 5.9M reads
           on an Intel Xeon E5-2670.  This is still fast enough for us to converge in a few minutes in most cases, and this is
           before multi-threading is implemented.
         */

        for (long j = alignment_incidence_->row_ptr[i]; j < alignment_incidence_->row_ptr[i+1]; ++j) {
            // make sure gene_sum_by_strain_hits_only_ is initialized to zero for any gene that has a hit for this read
            for (int k = 0; k < num_haplotypes; ++k) {
                gene_sum_by_strain_hits_only_[alignment_incidence_->tx_to_gene[alignment_incidence_->col_ind[j]] * num_haplotypes + k] = 0.0;
            }
            // also need to clear out gene_sum_hits_only_ and gene_masks_ for any gene with a hit
            gene_sum_hits_only_[alignment_incidence_->tx_to_gene[alignment_incidence_->col_ind[j]]] = 0.0;
            gene_masks_[alignment_incidence_->tx_to_gene[alignment_incidence_->col_ind[j]]] = 0;
        }

        // compute gene_sum_by_strain_hits_only_, gene_masks_
        for (long j = alignment_incidence_->row_ptr[i]; j < alignment_incidence_->row_ptr[i+1]; ++j) {
            start_index = alignment_incidence_->col_ind[j] * num_haplotypes;

            for (int k = 0; k < num_haplotypes; ++k) {
                if (alignment_incidence_->val[j] & (1 << k)) {
                    gene_sum_by_strain_hits_only_[alignment_incidence_->tx_to_gene[alignment_incidence_->col_ind[j]] * num_haplotypes + k] += previous_[start_index + k];
                }
                gene_masks_[alignment_incidence_->tx_to_gene[alignment_incidence_->col_ind[j]]] |= alignment_incidence_->val[j];
            }
        }

        // compute gene_sum_hits_only_ using gene_sum_by_strain_ and gene_masks_
        for (long j = alignment_incidence_->row_ptr[i]; j < alignment_incidence_->row_ptr[i+1]; ++j) {
            start_index = alignment_incidence_->col_ind[j] * num_haplotypes;

            if (gene_masks_[alignment_incidence_->tx_to_gene[alignment_incidence_->col_ind[j]]]) {
                for (int k = 0; k < num_haplotypes; ++k) {
                    if (gene_masks_[alignment_incidence_->tx_to_gene[alignment_incidence_->col_ind[j]]] & (1 << k)) {
                        gene_sum_hits_only_[alignment_incidence_->tx_to_gene[alignment_incidence_->col_ind[j]]] += gene_sum_by_strain_[alignment_incidence_->tx_to_gene[alignment_incidence_->col_ind[j]] * num_haplotypes + k];
                    }
                }

                // don't contribute this gene to the sum more than once even if it appears multiple times in read
                // zeroing out the mask will ensure we skip this gene if we see it again
                gene_masks_[alignment_incidence_->tx_to_gene[alignment_incidence_->col_ind[j]]] = 0;
            }
        }

        // finally we can compute the actual values to add to current_
        work_index = 0;
        for (long j = alignment_incidence_->row_ptr[i]; j < alignment_incidence_->row_ptr[i+1]; ++j) {

            start_index = alignment_incidence_->col_ind[j] * num_haplotypes;

            for (int k = 0; k < num_haplotypes; ++k) {
                if (alignment_incidence_->val[j] & (1 << k)) {
                    working_[work_index] = (previous_[start_index + k] / gene_sum_by_strain_hits_only_[alignment_incidence_->tx_to_gene[alignment_incidence_->col_ind[j]] * num_haplotypes + k]) *
                                           (gene_sum_by_strain_[alignment_incidence_->tx_to_gene[alignment_incidence_->col_ind[j]] * num_haplotypes + k] / gene_sum_hits_only_[alignment_incidence_->tx_to_gene[alignment_incidence_->col_ind[j]]]) *
                                           (gene_sums_[alignment_incidence_->tx_to_gene[alignment_incidence_->col_ind[j]]] / genes_total);

                    if (working_[work_index] != working_[work_index]) {
                        // nan! underflow!

                        // set to zero
                        working_[work_index] = 0.0;

                        // skip this location next time
                        alignment_incidence_->val[j] &= ~(1 << k);
                    }
                    read_sum += working_[work_index++];
                }
                else {
                    working_[work_index++] = 0.0;
                }
            }
        }

        //std::cout << "read_sum: " << read_sum << std::endl;
        // now we can normalize by read and sum into current_
        if (read_sum > 0) {
            work_index = 0;
            for (long j = alignment_incidence_->row_ptr[i]; j < alignment_incidence_->row_ptr[i+1]; ++j) {

                start_index = alignment_incidence_->col_ind[j] * num_haplotypes;

                for (int k = 0; k < num_haplotypes; ++k) {
                    current_[start_index + k] += (working_[work_index] / read_sum) * (double)alignment_incidence_->counts[i];
                    work_index++;
                }
            }
        }
    }
}


void SampleAllelicExpression::updateModel2() {
    std::swap(current_, previous_);

    // clear current_ so we can use it to accumulate sums
    std::fill(current_, current_ + size(), 0.0);

    int start_index;
    int work_index;
    int transcript_hits;
    double transcript_sum;
    double read_sum;
    auto end = alignment_incidence_->row_ptr.size() - 1;

    if (!gene_sums_) {
        gene_sums_ =  new double[alignment_incidence_->num_genes()];
    }

    if (!gene_sum_hits_only_) {
        gene_sum_hits_only_ = new double[alignment_incidence_->num_genes()];
    }

    if (!transcript_sums_) {
        transcript_sums_ = new double[num_transcripts];
    }

    // generate transcript sums from previous value
    for (int i = 0; i < num_transcripts; ++i) {
        transcript_sums_[i] = 0.0;
        for (int j = 0; j < num_haplotypes; ++j) {
            transcript_sums_[i] += previous_[i*num_haplotypes + j];
        }
    }

    // generate gene sums from previous value
    std::fill(gene_sums_, gene_sums_ + alignment_incidence_->num_genes(), 0.0);
    double genes_total = 0.0;
    for (long i = 0; i < num_transcripts; ++i) {
        gene_sums_[alignment_incidence_->tx_to_gene[i]] += transcript_sums_[i];
        genes_total += transcript_sums_[i];
    }

    // normalize gene_sums_
    for (int i = 0; i < alignment_incidence_->num_genes(); i++) {
        gene_sums_[i] /= genes_total;
    }

    //iterate over each read
    for (int i = 0; i != end; ++i) {
        work_index = 0;
        read_sum = 0.0;

        for (long j = alignment_incidence_->row_ptr[i]; j < alignment_incidence_->row_ptr[i+1]; ++j) {
            gene_sum_hits_only_[alignment_incidence_->tx_to_gene[alignment_incidence_->col_ind[j]]] = 0.0;

            start_index = alignment_incidence_->col_ind[j] * num_haplotypes;
            for (int k = 0; k < num_haplotypes; ++k) {
                // initialize all strains at this locus
                if (alignment_incidence_->val[j] & (1 << k)) {
                    working_[work_index++] = previous_[start_index + k];
                }
                else {
                    working_[work_index++] = 0.0;
                }
            }
        }

        transcript_hits = work_index / num_haplotypes;

        // normalize working_ by transcript
        for (long j = 0; j < transcript_hits; j++) {
            transcript_sum = 0.0;
            for (int k = 0; k < num_haplotypes; ++k) {
                transcript_sum += working_[j * num_haplotypes + k];
            }
            for (int k = 0; k < num_haplotypes; ++k) {
                 working_[j * num_haplotypes + k] /= transcript_sum;
            }
        }

        // isoform_specificity
        for (long j = alignment_incidence_->row_ptr[i]; j < alignment_incidence_->row_ptr[i+1]; ++j) {
            long gene_index = alignment_incidence_->tx_to_gene[alignment_incidence_->col_ind[j]];
            gene_sum_hits_only_[gene_index] += transcript_sums_[alignment_incidence_->col_ind[j]];
        }

        work_index = 0;
        for (long j = alignment_incidence_->row_ptr[i]; j < alignment_incidence_->row_ptr[i+1]; ++j) {
            long gene_index = alignment_incidence_->tx_to_gene[alignment_incidence_->col_ind[j]];
            double isoform_specificity = transcript_sums_[alignment_incidence_->col_ind[j]] / gene_sum_hits_only_[gene_index];
            for (int k = 0; k < num_haplotypes; ++k) {
                working_[work_index] *= isoform_specificity * gene_sums_[gene_index];

                if (working_[work_index] != working_[work_index]) {
                    // nan! underflow!

                    // set to zero
                    working_[work_index] = 0.0;

                    // skip this location next time
                    alignment_incidence_->val[j] &= ~(1 << k);
                }
                read_sum += working_[work_index++];

            }
        }

        work_index = 0;
        for (long j = alignment_incidence_->row_ptr[i]; j < alignment_incidence_->row_ptr[i+1]; ++j) {
            start_index = alignment_incidence_->col_ind[j] * num_haplotypes;
            for (int k = 0; k < num_haplotypes; ++k) {
                current_[start_index + k] += (working_[work_index++] / read_sum) * (double)alignment_incidence_->counts[i];
            }
        }
    }
}


void SampleAllelicExpression::updateModel3() {
    long start_index;
    int work_index;
    double read_sum;
    auto end = alignment_incidence_->row_ptr.size() - 1;

    if (!read_sum_by_gene_) {
        read_sum_by_gene_ = new double[alignment_incidence_->num_genes()];
    }

    if (!gene_sums_) {
        gene_sums_ =  new double[alignment_incidence_->num_genes()];
    }

    if (!gene_sum_hits_only_) {
        gene_sum_hits_only_ = new double[alignment_incidence_->num_genes()];
    }

    if (!transcript_sums_) {
        transcript_sums_ = new double[num_transcripts];
    }

    std::swap(current_, previous_);
    std::fill(current_, current_ + size(), 0.0);

    // generate transcript sums from previous value
    for (int i = 0; i < num_transcripts; ++i) {
        transcript_sums_[i] = 0.0;
        for (int j = 0; j < num_haplotypes; ++j) {
            transcript_sums_[i] += previous_[i*num_haplotypes + j];
        }
    }

    // generate gene sums from previous value
    std::fill(gene_sums_, gene_sums_ + alignment_incidence_->num_genes(), 0.0);
    double genes_total = 0.0;
    for (int i = 0; i < num_transcripts; ++i) {
        gene_sums_[alignment_incidence_->tx_to_gene[i]] += transcript_sums_[i];
        genes_total += transcript_sums_[i];
    }

    // normalize gene_sums_
    for (int i = 0; i < alignment_incidence_->num_genes(); i++) {
        gene_sums_[i] /= genes_total;
    }

    //iterate over each read
    for (int i = 0; i != end; ++i) {
        read_sum = 0.0;

        for (long j = alignment_incidence_->row_ptr[i]; j < alignment_incidence_->row_ptr[i+1]; ++j) {
            // make sure read_sum_by_gene_ is initialized to zero for any gene that has a hit for this read
            read_sum_by_gene_[alignment_incidence_->tx_to_gene[alignment_incidence_->col_ind[j]]] = 0.0;
        }

        work_index = 0;

        for (long j = alignment_incidence_->row_ptr[i]; j < alignment_incidence_->row_ptr[i+1]; ++j) {
            gene_sum_hits_only_[alignment_incidence_->tx_to_gene[alignment_incidence_->col_ind[j]]] = 0.0;

            start_index = alignment_incidence_->col_ind[j] * num_haplotypes;
            for (int k = 0; k < num_haplotypes; ++k) {
                // initialize all strains at this locus
                if (alignment_incidence_->val[j] & (1 << k)) {
                    working_[work_index++] = previous_[start_index + k];
                    read_sum_by_gene_[alignment_incidence_->tx_to_gene[alignment_incidence_->col_ind[j]]] += previous_[start_index + k];
                }
                else {
                    working_[work_index++] = 0.0;
                }
            }
        }

        work_index = 0;

        for (long j = alignment_incidence_->row_ptr[i]; j < alignment_incidence_->row_ptr[i+1]; ++j) {
            for (int k = 0; k < num_haplotypes; ++k) {
                if (alignment_incidence_->val[j] & (1 << k)) {
                    working_[work_index] /= read_sum_by_gene_[alignment_incidence_->tx_to_gene[alignment_incidence_->col_ind[j]]];
                    working_[work_index] *= gene_sums_[alignment_incidence_->tx_to_gene[alignment_incidence_->col_ind[j]]];

                    if (working_[work_index] != working_[work_index]) {
                        // nan! underflow!

                        // set to zero
                        working_[work_index] = 0.0;

                        // skip this location next time
                        alignment_incidence_->val[j] &= ~(1 << k);
                    }
                }
                read_sum += working_[work_index++];
            }
        }

        work_index = 0;

        for (long j = alignment_incidence_->row_ptr[i]; j < alignment_incidence_->row_ptr[i+1]; ++j) {
            start_index = alignment_incidence_->col_ind[j] * num_haplotypes;
            for (int k = 0; k < num_haplotypes; ++k) {
                current_[start_index + k] += (working_[work_index++] / read_sum)  * (double)alignment_incidence_->counts[i];
            }
        }

    }
}


void SampleAllelicExpression::updateModel4() {
    std::swap(current_, previous_);

    // clear current_ so we can use it to accumulate sums
    std::fill(current_, current_ + size(), 0.0);


    int start_index;
    int work_index;
    double read_sum;
    auto end = alignment_incidence_->row_ptr.size() - 1;

    //iterate over each read
    for (int i = 0; i != end; ++i) {
        work_index = 0;
        read_sum = 0.0;

        for (long j = alignment_incidence_->row_ptr[i]; j < alignment_incidence_->row_ptr[i+1]; ++j) {
            start_index = alignment_incidence_->col_ind[j] * num_haplotypes;

            for (int k = 0; k < num_haplotypes; ++k) {
                // initialize all strains at this locus
                if (alignment_incidence_->val[j] & (1 << k)) {
                    working_[work_index++] = previous_[start_index + k];
                    read_sum += previous_[start_index + k];
                }
                else {
                    working_[work_index++] = 0.0;
                }
            }
        }

        work_index = 0;

        for (long j = alignment_incidence_->row_ptr[i]; j < alignment_incidence_->row_ptr[i+1]; ++j) {
            start_index = alignment_incidence_->col_ind[j] * num_haplotypes;
            for (int k = 0; k < num_haplotypes; ++k) {
                current_[start_index + k] += (working_[work_index++] / read_sum)  * (double)alignment_incidence_->counts[i];
            }
        }
    }
}


bool SampleAllelicExpression::converged(double &change) {
    change = 0.0;
    double current_sum = 0.0;
    double previous_sum = 0.0;

    for (int i = 0; i < size(); ++i) {
        // change += std::abs(current_[i] - previous_[i]);
        current_sum += current_[i] * alignment_incidence_->transcript_lengths_[i];
        previous_sum += previous_[i] * alignment_incidence_->transcript_lengths_[i];
        change += std::abs(current_[i] - previous_[i]) * alignment_incidence_->transcript_lengths_[i];
    }

    std::cout << "current sum: " << current_sum << ", previous sum: " << previous_sum << std::endl;

    if (change < threshold_) {
        return true;
    }

    return false;
}


void SampleAllelicExpression::saveStackSums(std::string filename) {
    std::ofstream outfile;
    outfile.open(filename, std::fstream::app);

    outfile << "#Transcript";

    for (int i = 0; i < num_haplotypes; i++) {
        outfile << '\t' << alignment_incidence_->haplotype_names[i];
    }
    outfile << '\t' << "sum" << std::endl;

    // use 4 fixed decimal places for output
    outfile << std::fixed;
    outfile << std::setprecision(6);

    double sum;

    for (int i = 0; i < num_transcripts; i++) {
        sum = 0.0;
        outfile << alignment_incidence_->transcript_names[i] << "\t";

        for (int j = 0; j < num_haplotypes;  j++) {
            outfile << current_[i * num_haplotypes + j];
            sum += current_[i * num_haplotypes + j];
            outfile << '\t';
        }
        outfile << sum << std::endl;
    }
}


void SampleAllelicExpression::saveStackSums(std::string filename, std::string sample_name) {
    std::ofstream outfile;
    outfile.open(filename, std::fstream::app);

    outfile << "##Sample: " << sample_name << std::endl;
    outfile << "#Transcript";

    for (int i = 0; i < num_haplotypes; i++) {
        outfile << '\t' << alignment_incidence_->haplotype_names[i];
    }
    outfile << '\t' << "sum" << std::endl;

    // use 4 fixed decimal places for output
    outfile << std::fixed;
    outfile << std::setprecision(6);

    double sum;

    for (int i = 0; i < num_transcripts; i++) {
        sum = 0.0;
        outfile << alignment_incidence_->transcript_names[i] << "\t";
        
        for (int j = 0; j < num_haplotypes;  j++) {
            outfile << current_[i * num_haplotypes + j];
            sum += current_[i * num_haplotypes + j];
            outfile << '\t';
        }
        outfile << sum << std::endl;
    }
}


void SampleAllelicExpression::applyTranscriptLength(bool in_tpm) {

    if (alignment_incidence_->transcript_lengths_.size() == 0) {
        // didn't load transcript lengths, don't do anything
        return;
    }

    if (in_tpm) {
        double sum = 0.0;
        for (int i = 0; i < num_transcripts; i++) {
            for (int j = 0; j < num_haplotypes;  j++) {
                current_[i * num_haplotypes + j] /= alignment_incidence_->transcript_lengths_[i * num_haplotypes + j];
                sum += current_[i * num_haplotypes + j];
            }
        }
        for (int i = 0; i < num_transcripts; i++) {
            for (int j = 0; j < num_haplotypes; j++) {
                current_[i * num_haplotypes + j] *= 1000000.0 / sum;
            }
        }
    } else {
        for (int i = 0; i < num_transcripts; i++) {
            for (int j = 0; j < num_haplotypes;  j++) {
                current_[i * num_haplotypes + j] /= alignment_incidence_->transcript_lengths_[i * num_haplotypes + j];
            }
        }
    }

}
