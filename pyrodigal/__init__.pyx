# coding: utf-8
# cython: language_level=3, linetrace=True

"""Bindings to Prodigal, an ORF finder for genomes, progenomes and metagenomes.
"""

# ----------------------------------------------------------------------------

cimport libc.errno
from libc.stdlib cimport free, qsort
from libc.string cimport memchr, memcmp, memcpy, memset, strcpy, strstr
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
from cpython.unicode cimport PyUnicode_FromUnicode

from pyrodigal.prodigal cimport bitmap, dprog, gene, node, sequence
from pyrodigal.prodigal.bitmap cimport bitmap_t
from pyrodigal.prodigal.gene cimport MAX_GENES, _gene
from pyrodigal.prodigal.metagenomic cimport NUM_META, _metagenomic_bin, initialize_metagenomic_bins
from pyrodigal.prodigal.node cimport _node
from pyrodigal.prodigal.sequence cimport calc_most_gc_frame, gc_content
from pyrodigal.prodigal.training cimport _training

# ----------------------------------------------------------------------------

import warnings
import threading

# ----------------------------------------------------------------------------

cdef size_t MIN_SINGLE_GENOME = 20000
cdef size_t IDEAL_SINGLE_GENOME = 100000
cdef size_t MIN_NODES = 98

# ----------------------------------------------------------------------------

cdef void sequence_to_bitmap(
    const Py_UNICODE* text,
    size_t slen,
    bitmap_t* seq,
    bitmap_t* rseq,
    bitmap_t* useq
):
    """Create bitmaps from a textual sequence in ``text``.

    Arguments:
        text (const Py_UNICODE*): A pointer to the raw unicode buffer storing
            the sequence. Characters other than 'ATGC' or 'atgc' will be
            ignored and added to the ``useq`` bitmap.
        slen (size_t): The length of the input sequence.
        seq (bitmap_t*): An address where to put the bitmap storing the
          sequence.
        rseq (bitmap_t*): An address where to put the bitmap storing the
          reverse complement.
        useq (bitmap_t*): An address where to put the bitmap storing the
          *unknown character* bitmap.

    """
    # allocate memory for the bitmaps
    cdef size_t blen = slen//4 + (slen%4 != 0)
    cdef size_t ulen = slen//8 + (slen%8 != 0)
    seq[0] = <bitmap_t> PyMem_Malloc(blen * sizeof(unsigned char))
    rseq[0] = <bitmap_t> PyMem_Malloc(blen * sizeof(unsigned char))
    useq[0] = <bitmap_t> PyMem_Malloc(ulen * sizeof(unsigned char))
    if not seq[0] or not useq[0] or not rseq[0]:
        PyMem_Free(seq[0])
        PyMem_Free(rseq[0])
        PyMem_Free(useq[0])
        raise MemoryError()

    cdef size_t i, j
    with nogil:
        # clear memory
        memset(seq[0], 0, blen * sizeof(unsigned char))
        memset(rseq[0], 0, blen * sizeof(unsigned char))
        memset(useq[0], 0, ulen * sizeof(unsigned char))

        # fill the bitmaps depending on the sequence
        for i,j in enumerate(range(0, slen*2, 2)):
            letter = text[i]
            if letter == u'A' or letter == u'a':
                pass
            elif letter == u'T' or letter == u't':
                bitmap.set(seq[0], j)
                bitmap.set(seq[0], j+1)
            elif letter == u'G' or letter == u'g':
                bitmap.set(seq[0], j)
            elif letter == u'C' or letter == u'c':
                bitmap.set(seq[0], j+1)
            else:
                bitmap.set(useq[0], i)
            j += 2
            i += 1

        # compute reverse complement
        sequence.rcom_seq(seq[0], rseq[0], useq[0], slen)


cdef size_t count_genes(_node* nodes, int path) nogil:
    """Count the number of genes found in the node list.

    Arguments:
        nodes (_node*): An array of dynamic programming nodes.
        path (int): An index found by `dprog.dprog`.

    Returns:
        size_t: The number of genes that can be extracted from the nodes list.

    """
    cdef size_t ctr = 0

    if path == -1:
        return 0

    while nodes[path].traceb != -1:
        path = nodes[path].traceb

    while path != -1 and ctr < gene.MAX_GENES:
        if nodes[path].elim != 1:
          if nodes[path].strand == 1 and nodes[path].type == sequence.STOP:
              ctr += 1
          elif nodes[path].strand == -1 and nodes[path].type != sequence.STOP:
              ctr += 1
        path = nodes[path].tracef

    return ctr


# ----------------------------------------------------------------------------

# Initializing the metagenomic bins is fast enough that we can afford to do it
# only once, when the module is imported; storing them in a global means that:
#   1. we don't need to copy them in the results
#   2. we don't have to worry about their lifetime within the results

cdef _metagenomic_bin META_BINS[NUM_META]
cdef size_t _i
for _i in range(NUM_META):
    memset(&META_BINS[_i], 0, sizeof(_metagenomic_bin))
    strcpy(<char*> &META_BINS[_i].desc, "None")
    META_BINS[_i].tinf = <_training*> PyMem_Malloc(sizeof(_training))
    if not META_BINS[_i].tinf:
        raise MemoryError()
    memset(META_BINS[_i].tinf, 0, sizeof(_training))
initialize_metagenomic_bins(META_BINS)


# ----------------------------------------------------------------------------

cdef class TrainingInfo:
    cdef _training* raw

    def __dealloc__(self):
        PyMem_Free(self.raw)


# ----------------------------------------------------------------------------

cdef class Genes:
    """A collection of all the genes found by Prodigal in a single sequence.

    It implements all the methods of `collection.abc.Sequence`, so genes can
    be accessed individually using an index, or iterated upon in forward or
    reverse mode. Genes are sorted by their leftmost coordinate, independently
    of the strand.
    """

    # the list of nodes in the input sequence
    cdef _node* nodes
    cdef size_t nn
    # the array of genes in the input sequence
    cdef _gene* genes
    cdef size_t ng
    # the training information
    cdef TrainingInfo tinf_rc
    cdef _training* tinf
    # the sequence length and bitmaps
    cdef size_t slen
    cdef bitmap_t seq
    cdef bitmap_t rseq
    cdef bitmap_t useq

    def __dealloc__(self):
        PyMem_Free(self.nodes)
        PyMem_Free(self.genes)
        PyMem_Free(self.seq)
        PyMem_Free(self.rseq)
        PyMem_Free(self.useq)

    def __len__(self):
        return self.ng

    def __getitem__(self, index):
        cdef size_t index_
        if index < 0:
            index += self.ng
        if index >= self.ng or index < 0:
            raise IndexError("list index out of range")
        return self._gene(index)

    def __iter__(self):
        return (self._gene(i) for i in range(self.ng))

    def __reversed__(self):
        return (self._gene(self.ng-i) for i in range(1, self.ng+1))

    cdef _gene(self, index):
        return Gene.__new__(Gene, self, index)


# ----------------------------------------------------------------------------

cdef class Gene:
    """A single gene found by Prodigal within a DNA sequence.
    """

    # a hard reference to the Genes instance that created this object
    # to avoid the data referenced by other pointers to be deallocated.
    cdef Genes genes
    #
    cdef _node* nodes
    cdef _gene* gene
    #
    cdef size_t slen
    cdef bitmap_t* seq
    cdef bitmap_t* rseq
    cdef bitmap_t* useq
    #
    cdef _training* tinf

    def __cinit__(self, Genes genes, size_t index):
        if index > genes.ng:
            raise IndexError(index)
        self.genes = genes
        self.nodes = genes.nodes
        self.gene = &genes.genes[index]
        self.slen = genes.slen
        self.seq = &genes.seq
        self.useq = &genes.useq
        self.rseq = &genes.rseq
        self.tinf = genes.tinf

    @property
    def _data(self):
        return (<bytes> self.gene.gene_data).decode('ascii')

    @property
    def begin(self):
        """`int`: The coordinate at which the gene begins.
        """
        return self.gene.begin

    @property
    def end(self):
        """`int`: The coordinate at which the gene ends.
        """
        return self.gene.end

    @property
    def strand(self):
        """`int`: *-1* if the gene is on the reverse strand, *+1* otherwise.
        """
        return self.nodes[self.gene.start_ndx].strand

    @property
    def partial_begin(self):
        """`bool`: whether the gene overlaps with the start of the sequence.
        """
        if self.strand == 1:
            return self.nodes[self.gene.start_ndx].edge == 1
        else:
            return self.nodes[self.gene.stop_ndx].edge == 1

    @property
    def partial_end(self):
        """`bool`: whether the gene overlaps with the end of the sequence.
        """
        if self.strand == 1:
            return self.nodes[self.gene.stop_ndx].edge == 1
        else:
            return self.nodes[self.gene.start_ndx].edge == 1

    @property
    def start_type(self):
        """`str`: The start codon of this gene.

        Can be one of ``ATG``, ``GTG`` or ``TTG``, or ``Edge`` if `Pyrodigal`
        has been initialized in open ends mode and the gene starts right at the
        beginning of the input sequence.
        """
        node = self.nodes[self.gene.start_ndx]
        start_type = 3 if node.edge else node.type
        return ["ATG", "GTG", "TTG" , "Edge"][start_type]

    @property
    def rbs_motif(self):
        """``str``, optional: The motif of the Ribosome Binding Site.

        Possible non-`None` values are ``GGA/GAG/AGG``, ``3Base/5BMM``,
        ``4Base/6BMM``, ``AGxAG``, ``GGxGG``, ``AGGAG(G)/GGAGG``, ``AGGA``,
        ``AGGA/GGAG/GAGG``, ``GGAG/GAGG``, ``AGGAG/GGAGG``, ``AGGAG``
        ``GGAGG`` or ``AGGAGG``.

        """
        cdef char* data = self.gene.gene_data
        cdef char* i = strstr(data, "rbs_motif")
        cdef char* j = <char*> memchr(i, b';', 30)
        cdef size_t length = j - i
        if i[10:length] == b"None":
            return None
        return i[10:length].decode("ascii")

    @property
    def rbs_spacer(self):
        """`str`, optional: The number of base pair between the RBS and the CDS.

        Possible non-`None` values are ``3-4bp``, ``5-10bp``, ``11-12bp`` or
        ``13-15bp``.

        """
        cdef char* data = self.gene.gene_data
        cdef char* i = strstr(data, "rbs_spacer")
        cdef char* j = <char*> memchr(i, b';', 30)
        cdef size_t length = j - i
        if i[11:length] == b"None":
            return None
        return i[11:length].decode("ascii")

    @property
    def gc_cont(self):
        """`float`: The GC content of the gene (between *0* and *1*).
        """
        cdef char* data = self.gene.gene_data
        cdef char* i = strstr(data, "gc_cont")
        cdef char* j = <char*> memchr(i, b'\0', 30)
        cdef size_t length = j - i
        return float(i[8:length])

    @property
    def translation_table(self):
        """`int`: The translation table used to translate the gene.
        """
        cdef int table = self.tinf[0].trans_table
        return table

    cpdef translate(self):
        """Translate the gene into a protein sequence.

        Returns:
          `str`: The proteins sequence as a string using the right translation
          table and the standard single letter alphabet for proteins.

        """
        # create a new PyUnicode string of the right length to hold the protein
        cdef size_t nucl_length = (<size_t> self.gene.end) - (<size_t> self.gene.begin)
        cdef size_t prot_length = nucl_length//3 + (nucl_length%3 != 0)
        cdef unicode string = PyUnicode_FromUnicode(NULL, prot_length)

        # extract the boundaries / bitmap depending on
        cdef bitmap_t* seq
        cdef size_t begin, end
        if self.nodes[self.gene.start_ndx].strand == 1:
            begin = self.gene.begin
            end = self.gene.end
            seq = self.seq
        else:
            begin = self.slen + 1 - self.gene.end
            end = self.slen + 1 - self.gene.begin
            seq = self.rseq

        # copy the aminoacids to the sequence buffer
        cdef size_t i = 0
        cdef size_t j = begin
        while j < end:
            (<Py_UNICODE*> string)[i] = sequence.amino(seq[0], j-1, self.tinf, i==0)
            j += 3
            i += 1

        # return the string containing the protein sequence
        return string


# ----------------------------------------------------------------------------

cdef class Pyrodigal:
    """An efficient ORF finder for genomes, progenomes and metagenomes.
    """
    #
    cdef public bint closed
    cdef public object lock
    cdef readonly bint meta
    cdef readonly size_t _num_seq

    #
    cdef size_t nn
    cdef _node* nodes
    cdef size_t max_slen

    #
    cdef size_t ng
    cdef _gene* genes
    cdef size_t max_genes

    #
    cdef TrainingInfo tinf

    def __init__(self, meta=False, closed=False):
        """Instantiate and configure a new ORF finder.

        Arguments:
            meta (`bool`): Set to `True` to run in metagenomic mode, using a
                pre-trained profiles for better results with metagenomic or
                progenomic inputs. Defaults to `False`.
            closed (`bool`): Set to `True` to consider sequences ends 'closed',
                which prevents proteins from running off edges. Defaults to
                `False`.

        """
        self.meta = meta
        self.closed = closed
        self.lock = threading.Lock()

    def __cinit__(self, meta=False, closed=False):
        self._num_seq = 1
        # node array, uninitialized on object creation to reduce memory usage
        self.max_slen = 0
        self.nn = 0
        self.nodes = NULL
        # gene array, uninitialized on object creation to reduce memory usage
        self.max_genes = 0
        self.ng = 0
        self.genes = NULL

    def __dealloc__(self):
        PyMem_Free(self.nodes)
        PyMem_Free(self.genes)

    def find_genes(self, sequence):
        """Find all the genes in the input DNA sequence.

        Arguments:
            sequence (`str`): A DNA sequence to process. Letters not
                corresponding to a usual nucleotide (not any of "ATGC") will
                be ignored.

        Returns:
            `Genes`: A collection of all the genes found in the input.

        Raises:
            `MemoryError`: when allocation of an internal buffers fails.
            `RuntimeError`: on calling this method without `train` in *single* mode.
            `TypeError`: when ``sequence`` is not a string.

        """
        if not self.meta and self.tinf is None:
            raise RuntimeError("cannot find genes without having trained in single mode")
        if not isinstance(sequence, str):
            raise TypeError(f"sequence must be a string, not {type(sequence).__name__}")

        cdef size_t slen = len(sequence)
        cdef bitmap_t seq = NULL
        cdef bitmap_t rseq = NULL
        cdef bitmap_t useq = NULL
        sequence_to_bitmap(sequence, slen, &seq, &rseq, &useq)

        with self.lock:
            if self.meta:
                return self._find_genes_meta(slen, seq, useq, rseq)
            else:
                return self._find_genes_single(slen, seq, useq, rseq)

    cdef _reallocate_nodes(self, size_t slen):
        cdef size_t new_length = slen//8 + (slen%8 != 0)
        if new_length < MIN_NODES:
            new_length = MIN_NODES
        self.nodes = <_node*> PyMem_Realloc(self.nodes, new_length*sizeof(_node))
        if not self.nodes:
            raise MemoryError()
        self.max_slen = new_length*8

    cdef _find_genes_meta(self, size_t slen, bitmap_t seq, bitmap_t useq, bitmap_t rseq):
        cdef size_t i
        cdef size_t gene_count
        cdef size_t new_length

        # reallocate memory for the nodes if this is the biggest sequence
        # processed by this object so far
        if slen > self.max_slen:
            self._reallocate_nodes(slen)

        cdef size_t gc_count = 0
        cdef double gc, low, high
        cdef double max_score = -100
        cdef size_t max_phase = 0

        with nogil:
          # compute the GC% of the sequence
            for i in range(slen):
                gc_count += sequence.is_gc(seq, i)
            gc = (<double> gc_count) / slen

            # compute the min/max acceptable gc for the sequence to only
            # use appropriate metagenomic bins
            low = 0.88495*gc - 0.0102337
            high = 0.86596*gc + .1131991
            if low > 0.65:
                low = 0.65
            if high < 0.35:
                high = 0.35

            # check which of the metagenomic bins gets the best results
            for i in range(NUM_META):
                # recreate the node list if the translation table changed
                if i == 0 or META_BINS[i].tinf.trans_table != META_BINS[i-1].tinf.trans_table:
                    memset(self.nodes, 0, self.nn*sizeof(_node))
                    self.nn = node.add_nodes(seq, rseq, slen, self.nodes, self.closed, NULL, 0, META_BINS[i].tinf)
                    qsort(self.nodes, self.nn, sizeof(_node), node.compare_nodes)

                # check the GC% is compatible with the current bin
                if META_BINS[i].tinf.gc < low or META_BINS[i].tinf.gc > high:
                    continue

                # compute the score for the current metagenomic bin
                node.reset_node_scores(self.nodes, self.nn)
                node.score_nodes(seq, rseq, slen, self.nodes, self.nn, META_BINS[i].tinf, self.closed, True)
                node.record_overlapping_starts(self.nodes, self.nn, META_BINS[i].tinf, 1)
                ipath = dprog.dprog(self.nodes, self.nn, META_BINS[i].tinf, 1)

                # update if the current bin gave a better score
                if self.nodes[ipath].score > max_score:
                    max_phase = i
                    max_score = self.nodes[ipath].score
                    dprog.eliminate_bad_genes(self.nodes, ipath, META_BINS[i].tinf)
                    # reallocate memory for the nodes if this is the largest amount
                    # of genes found so far
                    gene_count = count_genes(self.nodes, ipath)
                    if gene_count > self.max_genes:
                        with gil:
                            self.genes = <_gene*> PyMem_Realloc(self.genes, gene_count*sizeof(_gene))
                            if not self.genes:
                                raise MemoryError()
                        self.max_genes = gene_count
                    # extract the genes from the dynamic programming array
                    self.ng = gene.add_genes(self.genes, self.nodes, ipath)
                    gene.tweak_final_starts(self.genes, self.ng, self.nodes, self.nn, META_BINS[i].tinf)
                    gene.record_gene_data(self.genes, self.ng, self.nodes, META_BINS[i].tinf, self._num_seq)

            # recover the nodes corresponding to the best run
            memset(self.nodes, 0, self.nn*sizeof(_node))
            self.nn = node.add_nodes(seq, rseq, slen, self.nodes, self.closed, NULL, 0, META_BINS[max_phase].tinf)
            qsort(self.nodes, self.nn, sizeof(_node), node.compare_nodes)
            node.score_nodes(seq, rseq, slen, self.nodes, self.nn, META_BINS[max_phase].tinf, self.closed, True)

        # make a `Genes` instance to store the results
        cdef Genes genes = Genes.__new__(Genes)
        genes.tinf_rc = None
        # copy nodes
        genes.nn = self.nn
        genes.nodes = <_node*> PyMem_Malloc(self.nn*sizeof(_node))
        if not genes.nodes: raise MemoryError()
        memcpy(genes.nodes, self.nodes, self.nn*sizeof(_node))
        # copy genes
        genes.ng = self.ng
        genes.genes = <_gene*> PyMem_Malloc(self.ng*sizeof(_gene))
        if not genes.genes: raise MemoryError()
        memcpy(genes.genes, self.genes, self.ng*sizeof(_gene))
        # take ownership of bitmaps
        genes.slen = slen
        genes.seq = seq
        genes.rseq = rseq
        genes.useq = useq
        # keep reference to training information
        # (this reference should never be invalid since META_BINS has been
        # allocated constantly and therefore will not be unallocated while
        # `genes` exist )
        genes.tinf = META_BINS[max_phase].tinf

        # free resources
        memset(self.nodes, 0, self.nn*sizeof(_node))
        memset(self.genes, 0, self.ng*sizeof(_gene))
        self.ng = self.nn = 0
        self._num_seq += 1

        #
        return genes

    cdef _find_genes_single(self, size_t slen, bitmap_t seq, bitmap_t useq, bitmap_t rseq):
        # reallocate memory for the nodes if this is the biggest sequence
        # processed by this object so far
        if slen > self.max_slen:
            self._reallocate_nodes(slen)

        # find all the potential starts and stops, and sort them
        self.nn = node.add_nodes(seq, rseq, slen, self.nodes, self.closed, NULL, 0, self.tinf.raw)
        qsort(self.nodes, self.nn, sizeof(_node), node.compare_nodes)

        # second dynamic programming, using the dicodon statistics as the scoring
        # function
        node.score_nodes(seq, rseq, slen, self.nodes, self.nn, self.tinf.raw, self.closed, False)
        node.record_overlapping_starts(self.nodes, self.nn, self.tinf.raw, True)
        cdef int ipath = dprog.dprog(self.nodes, self.nn, self.tinf.raw, True)
        dprog.eliminate_bad_genes(self.nodes, self.nn, self.tinf.raw)

        # reallocate memory for the nodes if this is the largest amount
        # of genes found so far
        gene_count = count_genes(self.nodes, ipath)
        if gene_count > self.max_genes:
            self.genes = <_gene*> PyMem_Realloc(self.genes, gene_count*sizeof(_gene))
            if not self.genes: raise MemoryError()
            self.max_genes = gene_count

        # extract the genes from the dynamic programming array
        self.ng = gene.add_genes(self.genes, self.nodes, ipath)
        gene.tweak_final_starts(self.genes, self.ng, self.nodes, self.nn, self.tinf.raw)
        gene.record_gene_data(self.genes, self.ng, self.nodes, self.tinf.raw, self._num_seq)

        # make a `Genes` instance to store the results
        cdef Genes genes = Genes.__new__(Genes)
        genes.tinf_rc = self.tinf
        # copy nodes
        genes.nn = self.nn
        genes.nodes = <_node*> PyMem_Malloc(self.nn*sizeof(_node))
        if not genes.nodes: raise MemoryError()
        memcpy(genes.nodes, self.nodes, self.nn*sizeof(_node))
        # copy genes
        genes.ng = self.ng
        genes.genes = <_gene*> PyMem_Malloc(self.ng*sizeof(_gene))
        if not genes.genes: raise MemoryError()
        memcpy(genes.genes, self.genes, self.ng*sizeof(_gene))
        # take ownership of bitmaps
        genes.slen = slen
        genes.seq = seq
        genes.rseq = rseq
        genes.useq = useq
        # keep reference to training information
        # (this reference should never be invalid since `self.tinf` is a
        # reference-counted Python object, which has at least the same lifetime
        # as the `Genes` instance)
        genes.tinf = self.tinf.raw

        # free resources
        memset(self.nodes, 0, self.nn*sizeof(_node))
        memset(self.genes, 0, self.ng*sizeof(_gene))
        self.ng = self.nn = 0
        self._num_seq += 1

        #
        return genes

    def train(self, sequence, force_nonsd=False, st_wt=4.35, trans_table=11):
      """Search optimal parameters for the ORF finder using a training sequence.

      Arguments:
          sequence (`str`): The sequence to use, as a string of nucleotides.
          force_nonsd (`bool`, optional): Set to ``True`` to bypass the heuristic
              algorithm that tries to determine if the organism the training
              sequence belongs to uses a Shine-Dalgarno motif or not.
          st_wt (`float`, optional): The start score weight to use. The default
              value has been manually selected by the PRODIGAL authors as an
              appropriate value for 99% of genomes.
          trans_table (`int`, optional): The translation table to use.

      Raises:
          `MemoryError`: when allocation of an internal buffers fails.
          `RuntimeError`: when calling this method while in *metagenomic* mode.
          `TypeError`: when ``sequence`` is not a string.
          `ValueError`: when ``trans_table`` is not a valid translation table.

      """
      if self.meta:
          raise RuntimeError("cannot use training sequence in metagenomic mode")
      if not isinstance(sequence, str):
          raise TypeError(f"sequence must be a string, not {type(sequence).__name__}")
      if trans_table not in set((*range(1, 7), *range(9, 17), *range(21, 26))):
          raise ValueError(f"{trans_table} is not a valid translation table index")

      # check we have enough nucleotides to train
      cdef size_t slen = len(sequence)
      if slen < MIN_SINGLE_GENOME:
          raise ValueError(
            f"sequence must be at least {MIN_SINGLE_GENOME} characters ({slen} found)"
          )
      elif slen < IDEAL_SINGLE_GENOME:
          warnings.warn(
            f"sequence should be at least {IDEAL_SINGLE_GENOME} characters ({slen} found)"
          )

      # convert sequence to bitmap for dynamic programming
      cdef bitmap_t seq = NULL
      cdef bitmap_t rseq = NULL
      cdef bitmap_t useq = NULL
      sequence_to_bitmap(sequence, slen, &seq, &rseq, &useq)

      # create the training structure and compute GC content
      cdef _training* tinf = <_training*> PyMem_Malloc(sizeof(_training))
      memset(tinf, 0, sizeof(_training));
      tinf.gc = gc_content(seq, 0, slen-1)
      tinf.st_wt = st_wt
      tinf.trans_table = trans_table

      cdef int* gc_frame;
      cdef int ipath;
      with self.lock:
          # check if we need to reallocate the node array
          if slen > self.max_slen:
              self._reallocate_nodes(slen)

          # find all the potential starts and stops and sort them
          self.nn = node.add_nodes(seq, rseq, slen, self.nodes, self.closed, NULL, 0, tinf)
          qsort(self.nodes, self.nn, sizeof(_node), node.compare_nodes)

          # scan all the ORFs looking for a potential GC bias in a particular
          # codon position, in order to acquire a good initial set of genes
          gc_frame = calc_most_gc_frame(seq, slen)
          if not gc_frame:
              raise MemoryError()
          node.record_gc_bias(gc_frame, self.nodes, self.nn, tinf)
          free(gc_frame)

          # do an initial dynamic programming routine with just the GC frame bias
          # used as a scoring function.
          node.record_overlapping_starts(self.nodes, self.nn, tinf, 0)
          ipath = dprog.dprog(self.nodes, self.nn, tinf, 0)

          # gather dicodon statistics for the training set
          node.calc_dicodon_gene(tinf, seq, rseq, slen, self.nodes, ipath)
          node.raw_coding_score(seq, rseq, slen, self.nodes, self.nn, tinf)

          # determine if this organism uses Shine-Dalgarno and score the node
          node.rbs_score(seq, rseq, slen, self.nodes, self.nn, tinf)
          node.train_starts_sd(seq, rseq, slen, self.nodes, self.nn, tinf)
          if force_nonsd:
              tinf.uses_sd = False
          else:
              node.determine_sd_usage(tinf)
          if not tinf.uses_sd:
              node.train_starts_nonsd(seq, rseq, slen, self.nodes, self.nn, tinf)

          # reset internal buffers and free allocated memory
          PyMem_Free(seq)
          PyMem_Free(rseq)
          PyMem_Free(useq)
          memset(self.nodes, 0, self.nn*sizeof(_node))
          self.nn = 0

          # store the training information in a Python object so it can be
          # shared with reference counting in the later `find_genes` calls
          self.tinf = TrainingInfo.__new__(TrainingInfo)
          self.tinf.raw = tinf
