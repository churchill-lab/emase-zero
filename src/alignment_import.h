/*
 * Copyright (c) 2020 The Jackson Laboratory
 *
 * This software was developed by Gary Churchill's Lab at The Jackson
 * Laboratory (see http://research.jax.org/faculty/churchill).
 *
 */

#ifndef ALIGNMENT_IMPORT_H
#define ALIGNMENT_IMPORT_H


AlignmentIncidenceMatrix *loadFromBin(std::string filename);
void loadNFromBin(std::string filename, AlignmentIncidenceMatrix &aim, int sample_idx);
bool isGZipped(std::string filename);
int getBinFormat(std::string filename);

#endif
