/*
 * Copyright (c) 2020 The Jackson Laboratory
 *
 * This software was developed by Gary Churchill's Lab at The Jackson
 * Laboratory (see http://research.jax.org/faculty/churchill).
 *
 */

//
//  utilities.h
//
//  Created by Glen Beane on 9/4/14.
//

#ifndef UTILITIES_H
#define UTILITIES_H


// A macro to disallow the copy constructor and operator= functions
// This should be used in the private: declarations for a class
#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
  TypeName(const TypeName&);               \
  void operator=(const TypeName&)


#endif
