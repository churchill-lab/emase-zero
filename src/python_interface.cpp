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
//  python_interface.cpp
//
//  Created by Glen Beane on 10/21/14.
//

#include "Python.h"
#include "python_interface.h"

#include <string>

static const char *module_name = "TranscriptHits";

std::vector<std::string> pythonStringListToVector(PyObject *list);
std::vector<int> pythonIntListToVector(PyObject *list);

PythonInterface::PythonInterface()
{
    module_ = NULL;
    module_dict_ = NULL;
    transcript_hits_ = NULL;
    err_string = "";
}

int PythonInterface::init()
{
    Py_Initialize();
    module_ = PyImport_ImportModule(module_name);

    if (!module_) {
        if (PyErr_Occurred()) {
            setErrorStringFromPythonError();
        }
        else {
            err_string = "Unable to load TranscriptHits Python Module. Check PYTHONPATH.";
        }
        return 1;
    }

    module_dict_ = PyModule_GetDict((PyObject*)module_);

    PyObject *constructor = PyDict_GetItemString((PyObject*)module_dict_, "TranscriptHits");
    transcript_hits_ = PyObject_CallObject(constructor, NULL);
    Py_DECREF(constructor);

    return 0;
}

AlignmentIncidenceMatrix *PythonInterface::load(std::string filename) {
    AlignmentIncidenceMatrix *aim = NULL;

    // some variables we will reuse as we interact with the Python object
    PyObject *results;

    //read in file
    results = PyObject_CallMethod((PyObject*)transcript_hits_, (char*)"load", (char*)"s", filename.c_str());

    if (!results) {
        if (PyErr_Occurred()) {
            setErrorStringFromPythonError();
        }
        return NULL;
    }

    Py_DECREF(results);

    //TODO need to do some error checking and make sure these three calls don't return empty vectors


    PyObject *vectors = PyObject_CallMethod((PyObject*)transcript_hits_, (char*)"to_populase_array", NULL);
    PyObject *indptr = PyDict_GetItemString(vectors, (char*)"indptr");
    PyObject *indices = PyDict_GetItemString(vectors, (char*)"indices");
    PyObject *values = PyDict_GetItemString(vectors, (char*)"data");
    PyObject *transcripts = PyDict_GetItemString(vectors, (char*)"transcripts");
    PyObject *reads = PyDict_GetItemString(vectors, (char*)"reads");
    PyObject *genes = PyDict_GetItemString(vectors, (char*)"genes");
    PyObject *tx_to_gene = PyDict_GetItemString(vectors, (char*)"tx_to_gene");
    PyObject *strains = PyDict_GetItemString(vectors, (char*)"strains");


    // These are all cadidates for return value optimization, which should avoid unnecessary copies
    std::vector<std::string> read_names = pythonStringListToVector(reads);
    std::vector<std::string> transcript_names = pythonStringListToVector(transcripts);
    std::vector<std::string> gene_names = pythonStringListToVector(genes);
    std::vector<int> transcript_to_gene = pythonIntListToVector(tx_to_gene);
    std::vector<int> row_ptr = pythonIntListToVector(indptr);
    std::vector<int> col_ind = pythonIntListToVector(indices);
    std::vector<int> val = pythonIntListToVector(values);
    std::vector<std::string> haplotype_names = pythonStringListToVector(strains);


    // do some limited error checking

    if (col_ind.size() != val.size()) {
        err_string = "Sparse matrix column index vector does not match value vector. File is corrupt.";
        return NULL;
    }

    if (transcript_to_gene.size() > 0 && transcript_to_gene.size() < transcript_names.size()) {
        // file contained transcript to gene mappings, but lenght of the mapping vector does not match the number of transcripts
        err_string = "Input file contained transcript to gene mapping information, but size of mapping vector does not match number of transcripts.";
        return NULL;
    }


    aim = new AlignmentIncidenceMatrix(haplotype_names, read_names, transcript_names, col_ind, row_ptr, val);

    if (transcript_to_gene.size()) {
        aim->setGeneMappings(transcript_to_gene);
    }

    aim->setGeneNames(gene_names);
    return aim;
}

std::vector<std::string> pythonStringListToVector(PyObject *list) {
    std::vector<std::string> strings;

    Py_ssize_t size = PySequence_Size(list);
    strings.reserve(size);

    for (Py_ssize_t i = 0; i < size; ++i) {
        PyObject *s = PySequence_ITEM(list, i);
        strings.push_back(std::string(PyString_AsString(s)));
        Py_DECREF(s);
    }
    return strings;
}

std::vector<int> pythonIntListToVector(PyObject *list) {
    std::vector<int> v;

    Py_ssize_t size = PySequence_Size(list);
    v.reserve(size);

    for (Py_ssize_t i = 0; i < size; ++i) {
        PyObject *s = PySequence_ITEM(list, i);
        v.push_back((int)PyInt_AsLong(s));
        Py_DECREF(s);
    }
    return v;
}

void PythonInterface::setErrorStringFromPythonError() {
    PyObject *excType, *excValue, *excTraceback;
    PyErr_Fetch(&excType, &excValue, &excTraceback);
    PyErr_NormalizeException(&excType, &excValue, &excTraceback);

    if (excValue) {
        PyObject *pystr = PyObject_Str(excValue);
        err_string = std::string(PyString_AsString(pystr));

        Py_DECREF(pystr);
    }

    Py_DECREF(excType);
    Py_XDECREF(excValue);
    Py_XDECREF(excTraceback);
}
