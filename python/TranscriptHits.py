#!/usr/bin/env python

"""
    TranscriptHits.py
    Al Simons
    August, 2014

    Works with a compressed version of read data that tracks only
    what contig a read hits against in which haplotype.  Typically,
    the contigs are transcripts, and the haplotypes are CC founder
    strains, but this is not required.

    While the storage format is most efficient for eight or fewer
    haplotypes, it is intended to be able to handle more.  This has
    not been tested--but any problems with more than eight haplotypes
    will be considered a bug.
"""
__author__ = 'asimons'

import sys
import inspect
import pysam

import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse import csc_matrix
from scipy.sparse import csr_matrix
from scipy import array

try:
    import cPickle as pickle
except:
    print >> sys.stderr, 'Falling back to using normal pickle'
    import pickle

import bz2
import time


class AlignmentPropertyMatrix(object):
    """
    The object we expose to clients for passing Transcript names,
    Read names, Strain names and the actual transcript hit arrays.

    properties:
        shape: shape of the 3d sparse matrix in the following
            order (num_of_transcripts, num_of_strains, num_of_reads)
        hname <== a list of haplotype (strain) names.
            Must be the same length as the 'data' property, below,
            and must not exceed length 8.  The order of the strain
            names in this list matches the order of the associated
            2D sparse matrices in the 'data' property.
        lname: a list of locus (transcript) names.  The order
            of the transcript names in this list matches the columns
            in all the sparse arrays.
        rname: a list of read names.  The order of the read names
            in this list matches the rows in all the sparse arrays.
        data: a list of sparse matrices.  The length of this list must
            match the length of the 'hname' property, above,
            and must not exceed length 8.
    """
    def __init__(self,
                 shape = None,
                 transcripts = None,
                 reads = None,
                 strains = None,
                 arrays = None):
        self._shape = shape
        self._hname = strains
        self._lname = transcripts
        self._rname = reads
        self._data  = arrays

    @property
    def shape(self):
        return self._shape

    @shape.setter
    def shape(self, value):
        assert len(value) == 3
        self._shape = value

    @property
    def hname(self):
        return self._hname

    @hname.setter
    def hname(self, value):
        if self._data is not None:
            assert len(value) == len(self._data)
        assert len(value) <= 8
        self._hname = value

    # Transcripts are "loci," thus, lname, not tname
    @property
    def lname(self):
        return self._lname

    @lname.setter
    def lname(self, value):
        self._lname = value

    @property
    def rname(self):
        return self._rname

    @rname.setter
    def rname(self, value):
        self._rname = value

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, value):
        self._data = value



    # DEBUGGING AND DUMPING ROUTINES
    def dump_base_arrays(self, ofn = None):
        if ofn:
            f = open(ofn, 'w')
        else:
            f = sys.stderr

        print >> f, '\nDumping base arrays'
        for i, strain in enumerate(self.hname):
            print >> f, strain
            print >> f, '    shape:  ', self.data[i].shape
            print >> f, '    data:   ', self.data[i].data
            print >> f, '    indices:', self.data[i].indices
            print >> f, '    indptrs:', self.data[i].indptr
        print >> f

        if ofn:
            f.close()

    def dump_KB_form(self, ofn = None):
        # For a pretty output of sets of test data, KB uses a form
        # Read:
        #    Transcript:   [[ strain list ]]
        # Creating that naively incurs about the worst possible
        # traversal of these data.  So we'll use the underlying lists
        # to get at the data.
        #
        map = {}
        for s_idx, s_data in enumerate(self.data):
            # Get the array's three vectors
            # indptrs = columns = transcripts
            # indices = read name indices
            # data = data points.
            indptrs = s_data.indptr
            indices = s_data.indices
            data = s_data.data
            # THere is a 1-1 correspondence between the indptrs vector
            # and the list of transcript names.

            # Loop over the transcripts
            for n in range(len(indptrs) - 1):
                if indptrs[n] == indptrs[n + 1]:
                    # There is no data in this column.
                    continue

                # For comparison with previous software, include the
                # transcript's index with its name.
                t_name = '{0} ({1})'.format(self.lname[n], n)

                # Loop over the reads for this transcript
                for m in range(indptrs[n],indptrs[n+1]):
                    # Record the Read index, the strain, and the value.
                    # Our map is read -> transcript -> strain -> value
                    r_idx = indices[m]
                    r_name = self.rname[r_idx]

                    # Insert the read in the map if this is the first
                    #  time we've seen it.
                    if r_name not in map:
                        map[r_name] = {}
                    r_map = map[r_name]

                    # We're converting a transcript(column) major
                    # matrix to a read-major map.
                    # Each read in the main map will have a sub-map
                    # for each transcript the read aligns to.  The
                    # values in the transcript is a string built up
                    # to show the strains.

                    # Insert this transcript in the read sub-map if this
                    # is the first time we've seen this transcript for
                    # this read.
                    if t_name not in r_map:
                        r_map[t_name] = ''
                    strains = r_map[t_name]

                    # Account for previous strains that didn't align
                    # this read to this transcript.
                    n = s_idx - len(strains)
                    strains += ' ' * n

                    # Record this one
                    strains += self.hname[s_idx]
                    # And save the modified string in the map
                    r_map[t_name] = strains

        # We've built the map that is able to efficiently output our
        # data in the required form.  Let's do it.
        ret = []

        for r in sorted(map.keys()):
            read = map[r]
            ret.append(r)
            for t in sorted(read.keys()):
                transcript = read[t]
                # Pad so that the closing brackets match for all
                # transcript lines.
                n = len(self.hname) - len(transcript)
                transcript += ' ' * n
                ret.append('    {0}: [{1}]'.format(t, transcript))

        ret.append('Shape:{0}'.format(self.shape))
        ret.append('Transcript length: {0}'.format(len(self.lname)))

        if ofn:
            f = open(ofn, 'w')
        else:
            f = sys.stderr
        print >> f, '\n'.join(ret)
        if ofn:
            f.close()



class TranscriptHits(object):
    """
        A class to manage read data that has been reduced to transcript hits
        per strain.

        The data is stored in a pickled file containing:
            1: a map of read hits structured as:
                {transcript_id: {read_idx: strain_mask, ...}, ...}
                where the transcript_id is the transcript name with the
                leading ENSMUST removed, the read_idx is the read name index (
                see below), and the mask is a bit-or'd mask of all the strains
                in which this read hit somewhere in this transcript.

                The haplotype bitmasks are assigned in the same order as
                haplotype names appear in the haplotype list.  For
                instance when used to track CC founder strains,
                with a strain name list of
                ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'],
                the mask values for the strains are:
                   'A' : 1,
                   'B' : 2,
                   'C' : 4,
                   'D' : 8,
                   'E' : 16,
                   'F' : 32,
                   'G' : 64,
                   'H' : 128

            2: a list of read names, used in conjunction with the read_idx
            saved in the map.

            3: A list of haplotype (strain) names

            4: A list of all contig names from which the transcripts in
            the data are taken.
    """

    def __init__(self,
                 strains=('A','B','C','D','E','F','G','H'),
                 transcripts=None,
                 transcript_file=None,
                 genes_file=None,
                 debug=[]):
        """
        Create a new TranscriptHits object.

        :param strains: list of strain names to be used

        :param transcripts: List of transcript names, without the
        strain name decorations.

        :param transcript_file: Path to a file of transcript names
        without the strain name decorations, with one name per line.

        If the 'transcripts' parameter is also specified,
        this parameter is ignored. If neither parameter is specified,
        the transcript names will be harvested from the bam file.

        If this TranscriptHits object is populated from an
        AlignmentPropertyMatrix object, the transcript names from
        that object overwrite any specified here.

        :param genes_file: Path to a TSV file containing gene to
        transcript mappings. First column is the gene name; remaining
        columns are the contained transcripts.  If using genes_file,
        either transcripts or transcript_file must be provided and
        it is an error to have a transcript without a gene mapping.

        :param debug: A list of string debug feature names, or the
        value True.  If true, all currently defined debug features
        are turned on.

        :return: None
        """
        if debug is True:
            # We were passed the value True.  Turn on all debugging:
            # Order in approximate frequency of use.
            self._debug = [
                           'Times',
                           # 'Details',
                           #'Files',
                           #'Sanity'
                          ]
        elif debug:
            # We were passed a list of debug items. Use it.
            self._debug = debug
        else:
            # No debugging.  Set an empty list
            self.debug = []

        self.debug_time(True)

        if genes_file and not (transcript_file or transcripts):
            raise ValueError('Transcripts must be provided to use'
                             'the genes-to-transcripts file')

        self._read_names = []
        self._map = {}
        self._array_list = []
        self._arrays_done = False
        self._c_aln = 0
        self._current_read = ''
        self._c_read_n = -1
        self._strains = strains
        self._num_strains = len(strains)
        self._strain_bits = {}
        self._transcript_to_idx = None

        # If we are provided with transcript names, we will use those,
        # and it will be an error to have one not in the list.  If,
        # however, we were not provided with the list of transcript
        # names a priori, we will harvest the names from the bam
        # file(s).
        #
        # We can take in either a list of transcript names or the path
        # to a file of transcript names, one per line.
        if transcripts:
            self._transcript_names = transcripts
        elif transcript_file:
            self._transcript_names = []
            for line in open(transcript_file):
                self._transcript_names.append(line.rstrip())
        else:
            self._transcript_names = []
        self._transcripts_provided = self.transcript_names
        self._harvest_transcripts = not self._transcript_names

        self._t_idx_to_g_idx = []
        self._gene_names = []
        if genes_file:
            # Create the gene names list
            # and stow the gene list index against the transcript name
            # in a dict.
            t_to_g_idx_map = {}
            # Track the current index in the gene names list
            gene_idx = 0
            for line in open(genes_file, 'rU'):
                parts = [x.strip() for x in line.split('\t')]
                self._gene_names.append(parts[0])
                for t in parts[1:]:
                    if t in t_to_g_idx_map:
                        raise ValueError('transcript already in map')
                    t_to_g_idx_map[t] = gene_idx
                gene_idx += 1

            # Now create the transcript to gene_idx list
            # This is indexed the same as the transcript names, and
            # contains the index of the name of the gene containing the
            # transcript in the gene_names list.
            for n in range(len(self._transcript_names)):
                tn = self._transcript_names[n]
                # Don't know what else to do if there is no gene
                # mapping for a transcript
                if tn not in t_to_g_idx_map:
                    raise ValueError('Transcript {0} has no gene '
                                     'mapping'.format(tn))
                self.t_idx_to_g_idx.append(
                    t_to_g_idx_map[tn])

        # Set up the strain mask bits. Look up by strain name.
        for n in range(self._num_strains):
            self._strain_bits[strains[n]] = 2**n

        self.debug_time(False)



    # Getters and setters
    @property
    def debug(self):
        return self._debug

    @debug.setter
    def debug(self, value):
        self._debug = value

    @property
    def read_names(self):
        return self._read_names

    @read_names.setter
    def read_names(self, value):
        self._read_names = value

    @property
    def transcript_names(self):
        return self._transcript_names

    @transcript_names.setter
    def transcript_names(self, value):
        self._transcript_names = value

    @property
    def map(self):
        return self._map

    @map.setter
    def map(self, value):
        self._map = value

    @property
    def array_list(self):
        return self._array_list

    @array_list.setter
    def array_list(self, value):
        self._array_list = value

    @property
    def strains(self):
        return self._strains

    @strains.setter
    def strains(self, value):
        self._strains = value

    @property
    def gene_names(self):
        return self._gene_names

    # No setter; only set in constructor

    @property
    def t_idx_to_g_idx(self):
        return self._t_idx_to_g_idx

    # No setter; only set in constructor


    # The remaining attributes are considered private
    # implementation artifacts.


    def write_file(self, ofn_prefix):
        """
        Write out the compressed file. The output filename prefix
        will have the string ".pcl.bz2" appended to it to form
        the full output file name.
        """
        self.debug_time(True)

        pickle_fn = '{0}.pcl.bz2'.format(ofn_prefix)
        pick = bz2.BZ2File(pickle_fn, 'wb')
        p = pickle.Pickler(pick, 2)

        wrapper = [self.map, self.read_names, self.transcript_names,
                   self.gene_names, self.t_idx_to_g_idx, self.strains]
        p.dump(wrapper)
        pick.close()
        self.debug_time(False)

    def load (self,pfn):
        """
        Read a compressed file and from it build the internal
        map-based representation.
        """
        self.debug_time(True)
        pick = bz2.BZ2File(pfn, 'rb')
        p = pickle.Unpickler(pick)
        parts = p.load()
        # Old files don't have the gene names and the transcript to
        # gene mapping.  So save the data away in two steps.
        self._map, self._read_names, self._transcript_names = parts[0:3]
        if len(parts) > 3:
            # We have the gene names and transcript to gene mapping
            self._gene_names, self._t_idx_to_g_idx = parts[3:5]
        else:
            self._gene_names = []
            self._t_idx_to_g_idx = []

        if len(parts) > 5:
            self.strains = parts[5]

        pick.close()
        self.debug_time(False)

    def transcript_to_idx(self, name=None):

        # Handle our initialization
        if not self._transcript_to_idx:
            self._transcript_to_idx = {}
            for i, nam in enumerate(self.transcript_names):
                # We only need one index for a name. Ignore
                # duplicates
                if nam in self._transcript_to_idx:
                    continue
                self._transcript_to_idx[nam] = i

        if not name:
            return

        if name not in self._transcript_to_idx:
            if not self._harvest_transcripts:
                print >> sys.stderr, 'Transcript name {0} not found ' \
                                     'in supplied transcript ' \
                                     'list.\nAdding anyway.'.format(
                    name)
            self._transcript_to_idx[name] = len(self.transcript_names)
            self.transcript_names.append(name)
        return self._transcript_to_idx[name]

    def to_populase_array(self):
        """
        The populase program written in C++ / boost needs these data
        in the form of three CSR vectors: indptrs, indices,
        strain-bitmaps.  Since the data are stored in the map with the
        primary key the transcripts (loci), it is naturally in CSC
        format.  So this routine uses roughly the same code flow as
        to_array_list(), below.  The difference is that this produces
        one array with the the strain bitmaps, instead of one array per
        strain.

        It will:
            - create the underlying CSC lists from the map
            - use them to create a csc_matrix
            - convert it to a csr_matrix,
            - strip out the numpy arrays that make up the three
              underlying vectors,
            - turn them back into python lists with python ints,
            - return.

        :return: A dict of lists representing:
            - indptrs
            - indexes
            - data
            - transcript names
            - read names
        """
        self.debug_time(True)

        # Make sure that we've been loaded successfully.
        if not self._map:
            raise Exception('TranscriptHits not initialized properly.')

        self.debug_time('Beginning csc_matrix vectors construction')

        indptrs = []
        indices = []
        data = []

        cur_trans_idx = -1
        if 'Details' in self._debug:
            print >> sys.stderr, 'Building arrays, trans len', \
                len(self.transcript_names), len(self.map)

        for transcript in sorted(self.map):
            # Transcripts are columns. Record where this column
            # starts in each of the eight lists

            # The transcript as stored in the map is the integer
            # part of the transcript name "ENSMUST0000XXXXXXXXX"
            # Need to convert that to the index of the transcript
            # name in the tx names list.
            t_name = self.int_to_transcript(transcript)
            t_idx = self.transcript_to_idx(t_name)

            for i in range(cur_trans_idx, t_idx):
                indptrs.append(len(indices))
            cur_trans_idx = t_idx

            # Process the reads that hit against this transcript
            reads = self._map[transcript]
            for read in reads:
                # The read is the index into the read_names list.
                strains = reads[read]
                indices.append(int(read))
                data.append(strains)


        # Done processing the map.  Put on the last column indptrs
        for i in range(cur_trans_idx, len(self.transcript_names)):
            indptrs.append(len(indices))

        # Some debugging code.
        if 'Details' in self._debug:
            print >> sys.stderr, 'Map length (transcripts)', \
                len(self._map)
            print >> sys.stderr, 'Lengths:'
            print >> sys.stderr, '    data:', len(data)
            print >> sys.stderr, '    indexes:', len(indices)
            print >> sys.stderr, '    indptrs:', len(indptrs)

        csc = csc_matrix((array(data),
                          array(indices),
                          array(indptrs)),
                         shape=(len(self._read_names),
                                len(self._transcript_names)))

        self.debug_time('Starting conversion to csr_matrix')
        csr = csc.tocsr()
        self.debug_time('csr conversion complete')

        indptr = csr.indptr.tolist()
        indices = csr.indices.tolist()
        data = csr.data.tolist()
        for n in range(len(indptr)):
            indptr[n] = int(indptr[n])
        for n in range(len(indices)):
            indices[n] = int(indices[n])
        for n in range(len(data)):
            data[n] = int(data[n])

        # DEBUGGING CODE
        if 'Details' in self._debug:
            s = csr.shape
            print 'Details of csr array with shape {0}\n'.format(s)
            if len(self._read_names) <= 100:
                for r in range(s[0]):
                    print self.read_names[r]
                    for c in range(s[1]):
                        if csr[r, c] != 0:
                            print '    {0}: {1}'.format(
                                self.transcript_names[c],
                                csr[r, c])
            print >> sys.stderr, 'Lengths: indptr: {0} ' \
                                 'indices: {1} data: {2}\n' \
                                 'transcripts: {3}, reads: {4}'.\
                format(
                    len(indptr),
                    len(indices),
                    len(data),
                    len(self.transcript_names),
                    len(self.read_names)
                )

        self.debug_time(False)
        return {'indptr': indptr,
                'indices': indices,
                'data': data,
                'transcripts': self.transcript_names,
                'reads': self.read_names,
                'genes': self.gene_names,
                'tx_to_gene': self.t_idx_to_g_idx,
                'strains': self.strains
        }

    def to_array_list(self):
        """
        Converts the internal map-based representation to a list of
        #strains sparse arrays.

        We are creating a list of up to eight
        scipy.sparse.CSC_matrices, one per strain. When using the
        default strains (the eight CC founder strains), the array for
        strain A will be in index 0, H will be in 7.

        Within each array, reads (~ 10^7) will be on rows and
        transcripts (~ 10^5) in the columns.

        We will use CSC_matrix constructor
            csc_matrix((data, indices, indptr), [shape=(M, N)])

        data and indices will both be the length of all of the hits;
        indptr will be the length of the number of transcripts, since
        all of the transcripts we found in the original BAM file by
        definition have at least one hit.
        """
        self.debug_time(True)

        # If we've already converted from the map to the arrays,
        # simply return the list.
        if self._arrays_done:
            self.debug_time(False)
            return self._array_list

        #print >> sys.stderr, self.map
        try:
            # Make sure that we've been loaded successfully.
            if not self._map:
                raise Exception(
                    'TranscriptHits not initialized properly.')

            # We're building multiple arrays at a time, one per strain
            # so each of these lists contain num_strains lists
            indptrs = []
            indices = []
            data = []
            for n in range(self._num_strains):
                indptrs.append([])
                indices.append([])
                data.append([])

            cur_trans_idx = -1
            #print >> sys.stderr, 'Building arrays, trans len', \
            #    len(self.transcript_names), len(self.map)

            for transcript in sorted(self.map):
                # Transcripts are columns. Record where this column
                # starts in each of the eight lists

                # The transcript as stored in the map is the integer
                # part of the transcript name "ENSMUST0000XXXXXXXXX"
                # Need to convert that to the index of the transcript
                # name in the tx names list.
                t_name = self.int_to_transcript(transcript)
                t_idx = self.transcript_to_idx(t_name)

                for i in range(cur_trans_idx, t_idx):
                    for n in range(self._num_strains):
                        indptrs[n].append(len(indices[n]))
                cur_trans_idx = t_idx

                # Process the reads that hit against this transcript
                reads = self._map[transcript]
                for read in reads:
                    # The read is the index into the read_names list.
                    strains = reads[read]
                    for n in range(self._num_strains):
                        if strains & 2 ** n:
                            indices[n].append(int(read))
                            data[n].append(1.0)

            # Done processing the map.  Put on the last column indptrs
            for i in range(cur_trans_idx, len(self.transcript_names)):
                for n in range(self._num_strains):
                    indptrs[n].append(len(indices[n]))

            # Some debugging code.
            if 'Details' in self._debug:
                print >> sys.stderr, 'Map length (transcripts)', \
                    len(self._map)
                for n in range(self._num_strains):
                    print >> sys.stderr, 'Lengths for', n
                    print >> sys.stderr, '    data:', len(data[n])
                    print >> sys.stderr, '    indexes:', len(indices[n])
                    print >> sys.stderr, '    indptrs:', len(indptrs[n])

            # Now we can actually create the eight csc_matrices.
            for n in range(self._num_strains):
                #print >> sys.stderr, 'Creating array for', \
                #    self.strains[n]
                #print >> sys.stderr, 'Type of indices is', \
                # type(indices[n]), indices[n]
                self._array_list.append(
                    csc_matrix((array(data[n]),
                                array(indices[n]),
                                array(indptrs[n])),
                               shape=(len(self._read_names),
                                      len(self._transcript_names))))

            # Mark that we've done this step
            self._arrays_done = True

            self.debug_time(False)
            return self._array_list
        except:
            self.debug_time(False)
            raise

    def migrate_clean_numpy_types(self):
        """
        This is a migration aid.  KB's HDF5 files stored things as
        numpy types; turn everything to plain Python types.

        :return None
        """
        for i, name in enumerate(self.transcript_names):
            self.transcript_names[i] = str(name)
        for i, name in enumerate(self.read_names):
            self.read_names[i] = str(name)

    def cleanup_intermediate(self):
        # Save space when the user says we don't need the map any more.
        self._map = None

    def int_to_transcript(self, i):
        name = 'ENSMUST{0:011d}'.format(int(i))
        if 'Details' in self.debug and not name in \
                self.transcript_names:
            print >> sys.stderr, 'Reconstructed transcript name {0}' \
                'not found in self.transcript_names'.format(name)
            raise Exception("Oops!")
        return name

    @staticmethod
    def transcript_to_int(transcript):
                # The transcripts in the index are named
                #    ENSMUSTnnnnnnnnnn_<letter>,
                # where <letter> is a CC founder code letter, A-H.
                # Isolate the transcript from the strain
                if '_' in transcript:
                    t_name, strain = transcript.split('_')
                else:
                    t_name = transcript
                    # Receiver should never expect a strain returned
                    # if it didn't pass in a transcript with appended
                    # strain.
                    strain = 'ERROR STRAIN'

                # We don't actually store the transcript name.  They are
                # all of the form ENSMUSTnnnnnnnnnnn.  We check that the
                # prefix is indeed ENSMUST, and then just store the
                # integer form of the trailing numbers.
                if t_name[:7] != 'ENSMUST':
                    raise ValueError('Found non-mouse transcript{0}.'
                                     .format(t_name))

                try:
                    t_int = int(t_name[7:])
                except ValueError:
                    raise ValueError(
                        'Unable to convert suffix of transcript name '
                        '{0} to an int.'.format(t_name))

                return t_int, strain

    def from_bam(self, bam_n):
        """
        Process a bam file containing transcript alignments, where the
        transcripts have a suffix indicating the strain for
        the transcript.

        Creates an internal form of the bam, which can either be used
        to write out a compressed transcript hits file, or create
        sparse arrays for conducting population studies.

        Returns the number of alignments processed.
        """

        self.debug_time(True)

        self._harvest_transcripts = not self._transcript_names
        if (not self._harvest_transcripts) and \
           (not self._transcript_to_idx):
            self.transcript_to_idx()

        try:
            # Invalidate the existing arrays, if any.
            self._arrays_done = False

            # Open the bam
            bam_f = pysam.Samfile(bam_n, "rb")

            # Get an iterator for all the alignments in the bam file
            alignments = bam_f.fetch(until_eof=True)

            c_aln = 0

            if 'Times' in self._debug:
                print 'parsing {1}: {0}'.format(time.strftime(
                    "%I:%M:%S"), bam_n)

            for aln in alignments:
                c_aln += 1
                read_name = aln.qname
                # A single read may have many alignments
                if read_name != self._current_read:
                    # Save the new current read and bump our index
                    self._current_read = read_name
                    self._c_read_n += 1
                    self._read_names.append(read_name)
                    #print >> sys.stderr, 'New read', read_name

                # Done with the read name processing.  Now work with the
                # transcript.
                transcript = bam_f.getrname(aln.tid)

                t_int, strain = TranscriptHits.transcript_to_int(
                    transcript)

                # The top level map is the transcripts.  Will eventually
                # be the columns (with the reads being the rows).
                # We track transcripts by the integer portion of
                # their name.
                if t_int not in self.map:
                    # Create the second level map for reads
                    self.map[t_int] = {}
                    if self._harvest_transcripts:
                        self._transcript_names.append(
                            self.int_to_transcript(t_int)
                        )
                tx = self.map[t_int]
                if self._c_read_n not in tx:
                    tx[self._c_read_n] = 0
                try:
                    tx[self._c_read_n] |= self._strain_bits[strain]
                except ValueError:
                    if not strain in self._strain_bits:
                        print >> sys.stderr, \
                            'BAM file {0} contains a transcript with ' \
                            'strain name {1}, which is not in the strain ' \
                            'list: {2}'.format(bam_n, strain,
                                               self.strains)
                    else:
                        # Let any other ValueError flow out; it isn't
                        # expected.
                        raise

            # To ensure the transcript names will align correctly,
            # they need to be sorted.
            self._transcript_names.sort()

            self.debug_time(False)

            # Done processing this bam. Return the number of alignments.
            return c_aln
        except:
            self.debug_time(False)
            raise

    def to_AlignmentPropertyMatrix(self):
        """
        Create an AlignmentPropertyMatrix object.

        :return: an initialized AlignmentPropertyMatrix object
        """
        self.debug_time(True)

        if not (self._arrays_done or self._map):
                raise ValueError('TranscriptHits object not '
                                 'initialized')

        return AlignmentPropertyMatrix(
            shape=(len(self.transcript_names),
                   len(self.strains),
                   len(self.read_names)),
            transcripts=self.transcript_names,
            reads=self.read_names,
            strains=self.strains,
            arrays= self.to_array_list())

        self.debug_time(False)

    def from_AlignmentPropertyMatrix(self, aln_prop_mat):
        """
        Takes in an aks_AlignmentPropertyMatrix, and converts its data
        to the internal format, suitable to write out to the
        compressed file format.

        This is primarily to convert KB's old format files to the new
        compressed file format.

        The aks_AlignmentPropertyMatrix has four lists:
          lname: the list of transcript names
          rname: the list of read names
          strains: the list of strains (must be <= 8)
          data: a list of eight scipy csc sparse matrices, one
                per strain, in the same order as strains.

        The arrays contain a row per read in the list of read names,
        where the array row index is the index into the list of read
        names. The arrays contain a column per transcript name in the
        list of transcript names, again with indices matching.

        We're building a doubly nested map, where the keys in the
        outer map is the integer portion of the transcript name (
        omitting 'ENSMUST'). The values of in the outer map are maps
        with the keys being the index of read names in the read name
        list.  The value of each read name map is a bit map of the
        strains in which that read mapped in that transcript.
        """

        self.debug_time(True)
        # Don't do this if we're not starting out with a new instance.
        assert not self.map

        self.transcript_names = aln_prop_mat.lname
        self.read_names = aln_prop_mat.rname
        #print >> sys.stderr, 'len transcript_names',
        #    len(self.transcript_names)
        #print >> sys.stderr, 'len read_names', len(self.read_names)

        data = aln_prop_mat.data

        # The strains list and arrays list (in data) must be the same.
        assert len(data) == len(self.strains)

        for n in range(len(data)):
            self.array_to_map(data[n], n)

        self.debug_time(False)

    def array_to_map(self, arr, idx):
        """
        Used by from
        :param arr:
        :param idx:
        :return:
        """
        # Some sanity checking.
        if 'Sanity' in self.debug:
            if not arr.shape[1] == len(self.transcript_names):
                raise ValueError('array.shape[1] ({0}) does not match '
                                 'length of transcript_names ({1})'
                    .format(arr.shape[1], len(self.transcript_names)))

        # Compute the bit mask that will represent this strain
        strain_bit = 2 ** idx

        #print >> sys.stderr, '\n\nProcessing strain', idx,
        #    self.strains[idx]

        # The columns are the transcripts. The indptr has num columns
        # + 1 items. The indices into indptr are the column indices.
        # The values in indptr[i] and indptr[i+1] are the open range
        # of indexes in the indices array, containing the row index
        # in the sparse array.  Confusing?  Yeah, me too.  See
        # http://docs.scipy.org/doc/scipy/reference/generated/\
        # scipy.sparse.csc_matrix.html
        # (Link visited 20140617)

        cols = arr.indptr
        if isinstance(cols, np.ndarray):
            cols = cols.tolist()
        if not isinstance(cols, list):
            raise Exception('Unexpected type for indptr', type(cols))

        idxs = arr.indices
        if isinstance(idxs, np.ndarray):
            idxs = idxs.tolist()
        if not isinstance(idxs, list):
            raise Exception('Unexpected type for indices',
                            type(idxs))

        for c in range(len(cols) - 1):
            # Were there any reads that mapped to this transcript?
            if cols[c] == cols[c + 1]:
                # No; a completely empty column.
                continue

            #print >> sys.stderr, 'Processing transcript idx', c
            # Get the integer representing this transcript
            transcript = self.transcript_names[c]
            t_int, _ = TranscriptHits.transcript_to_int(transcript)
            #print >> sys.stderr, 'Tx name', transcript, t_int

            if not t_int in self.map:
                self.map[t_int] = {}
                #print >> sys.stderr, 'added to map; new len', \
                #    len(self.map)
            t_map = self.map[t_int]


            # We're interested in just the reads for this transcript
            # Extract the range of the overall indices list that
            # represent the reads of this transcript, ie the non-zero
            # rows (reads) in this column (transcript).
            # Remember that the row indices are also the indices into
            # the read_names list.  But we don't care about the read
            # names.  We're just going to store the row / read
            # name indices in the map.
            reads = arr.indices[cols[c]:cols[c+1]]
            for read_idx in range(len(reads)):
                # Indirect through the indices:
                read = reads[read_idx]
                # Ensure that this read index in in the transcript's
                # map. If not add it with an initial value of
                # integer zero (we'll use it as a bit map
                # representing the strains).
                #print >> sys.stderr, 'Processing read', \
                #    self.read_names[read], read
                if not read in t_map:
                    t_map[read] = 0
                    #print >> sys.stderr, "added read to map"

                # Record this strain
                t_map[read] |= strain_bit
                # print >> sys.stderr, 'strain_bits now', t_map[read]

                # Note that we never referenced the array's data list.
                # We don't care about the data (every entered value
                # is a 1).  The only thing we care about is that
                # there is an entry for this transcript/read pair.

    # DEBUGGING AND DUMPING FUNCTIONS
    def debug_time(self, is_start):
        # Do we want to do this at all?
        if not 'Times' in self._debug:
            return

        if (is_start is True) or (is_start is False):
            which = 'start' if is_start else 'end'
        else:
            # Passed in an arbitrary message
            which = is_start
        caller_name = inspect.stack()[1][3]
        print >> sys.stderr, '{0}: {1}() {2}' \
            .format(time.strftime("%I:%M:%S"), caller_name, which)

    def dump_compressed_pickle_file(self, pfn, ofn=None):
        if ofn:
            of = open(ofn, 'w')
        else:
            of = sys.stdout

        self.load(pfn)
        print >> of, 'Length of read names:', len(self._read_names)
        print >> of, 'Length of transcript names:', \
            len(self._transcript_names)
        print >> of, 'Length of map:', len(self._map)
        print >> of
        print >> of, 'Read names'
        for n in self._read_names:
            print >> of, '   ', n
        print >> of
        print >> of, 'Transcript names'
        for n in self._transcript_names:
            print >> of, '   ', n
        print >> of
        print >> of, 'Map'
        for k in sorted(self._map.keys()):
            print >> of, str(k) + ':', self._map[k]

        if ofn:
            of.close()

    def dump_map(self, ofn=None):
        self.debug_time(True)
        if ofn:
            f = open(ofn, 'w')
        else:
            f = sys.stderr
        print >> f, self.map
        if ofn:
            f.close()
        self.debug_time(False)

    def dump_map_pretty(self, ofn=None):
        self.debug_time(True)
        if ofn:
            f = open(ofn, 'w')
        else:
            f = sys.stderr
        tot_num_xscripts = len(self.transcript_names)
        num_xscripts = len(self.map)
        num_alignments = 0
        num_reads = len(self.read_names)
        for transcript in sorted(self.map):
            print >> f, self.int_to_transcript(transcript)
            t = self.map[transcript]
            num_alignments += len(t)
            for read in sorted(t):
                print >> f, '\t{0}\t{1}'.format(
                    self.read_names[read], t[read])
        print >> f, '\nTotal Transcripts:\t{0}'.format(
            tot_num_xscripts)
        print >> f, 'Transcripts in file:\t{0}'.format(num_xscripts)
        print >> f, 'Number of reads:\t{0}'.format(num_reads)
        print >> f, 'Number of alignments\t{0}'.format(num_alignments)
        if ofn:
            f.close()
        self.debug_time(False)


# A main() for unit testing only.  Assumes the presence of an input
# BAM file, whose name is passed in as sys.arg[1].
def main():
    debug = True

    if debug:
        print >> sys.stderr, '*** Unit test begins: {0} ***\n' \
            .format(time.strftime("%I:%M:%S"))

    if len(sys.argv) == 3:
        transcript_fn = sys.argv[2]
    else:
        transcript_fn = None

    fn = sys.argv[1]
    if fn.endswith('.bam'):
        comp_fn = fn[:-4]
    else:
        comp_fn = fn

    th = TranscriptHits(debug=debug,
                         strains=['A','B','C','D','E','F','G','H'],
                         transcript_file=transcript_fn)

    # Create the internal representation
    if debug:
        print >> sys.stderr, '*** Create internal rep from bam: {0} ***' \
            .format(time.strftime("%I:%M:%S"))

    th.from_bam(fn)

    th.dump_map('map_from_{0}_1.txt'.format(comp_fn))
    th.dump_map_pretty('map_from_{0}_pretty_1.txt'.format(comp_fn))
    # Write out the old style to compare sizes
    if debug:
        print >> sys.stderr, '*** Writing compressed map: {0} ***' \
            .format(time.strftime("%I:%M:%S"))

    th.write_file(comp_fn)

    # Turn it into an alignmentPropertyMatrix and save to disk.
    if debug:
        print >> sys.stderr, '*** Building APM: {0} ***' \
            .format(time.strftime("%I:%M:%S"))

    apm = th.to_AlignmentPropertyMatrix()

    if debug:
        print >> sys.stderr, '*** Writing APM: {0} ***' \
            .format(time.strftime("%I:%M:%S"))

    if debug:
        apm.dump_KB_form('kb1.txt')
        th.dump_compressed_pickle_file(comp_fn + '.pcl.bz2',
                                       'raw_dump_of_pickle.txt')

    print >> sys.stderr, '*** Reloading from the pcl file'
    th = TranscriptHits(debug=debug,
                         strains=['A','B','C','D','E','F','G','H'],
                         transcript_file=transcript_fn)

    th.load(comp_fn + '.pcl.bz2')

    th.dump_map('map_from_{0}_2.txt'.format(comp_fn))
    th.dump_map_pretty('map_from_{0}_pretty_2.txt'.format(comp_fn))

    apm = th.to_AlignmentPropertyMatrix()

    if debug:
        apm.dump_KB_form('kb2.txt')

    th.write_file(comp_fn + '2')

    if debug:
        print >> sys.stderr, '*** Unit test ends: {0} ***' \
            .format(time.strftime("%I:%M:%S"))

if __name__ == '__main__':
    main()
