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
//  python_interface.h
//
//  Created by Glen Beane on 10/21/14.
//

#ifndef PYTHON_INTERFACE
#define PYTHON_INTERFACE


#include "alignment_incidence_matrix.h"

class PythonInterface {
    
public:
    PythonInterface();
    int init();
    AlignmentIncidenceMatrix *load(std::string filename);
    
    inline std::string getErrorString() {
        return err_string;
    }

private:
    
    void *module_;
    void *module_dict_;
    void *transcript_hits_;
    
    std::string err_string;
    
};


#endif 
