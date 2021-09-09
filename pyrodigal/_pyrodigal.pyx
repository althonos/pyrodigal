# coding: utf-8
# cython: language_level=3, linetrace=True

from . cimport _api
from ._api cimport Nodes, Genes, TrainingInfo, Sequence

import ctypes
import logging
import multiprocessing.pool
import queue
import threading

from ._api import METAGENOMIC_BINS

LOGGER = logging.getLogger("pyrodigal")


class _MetaThread(threading.Thread):

    def __init__(
        self,
        seq,
        closed,
        input_queue,
        max_phase,
        max_score,
        genes,
        genes_lock,
        kill_switch,
    ):
        super().__init__()
        self.seq = seq
        self.closed = closed
        self.nodes = _api.Nodes()
        self.translation_table = -1
        self.input_queue = input_queue
        self.max_phase = max_phase
        self.max_score = max_score
        self.genes = genes
        self.genes_lock = genes_lock
        self.kill_switch = kill_switch
        self.error = None

        self.gc = self.seq.gc
        self.low =  min(0.65, 0.88495*self.gc - 0.0102337)
        self.high = max(0.35, 0.86596*self.gc + 0.1131991)

    def run(self):
        while not self.kill_switch.is_set():
            # attempt to get the next argument, with a timeout
            # so that the thread can periodically check if it has
            # been killed, even when the query queue is empty
            try:
                args = self.input_queue.get(timeout=1)
            except queue.Empty:
                continue
            # check if arguments from the queue are a poison-pill (`None`),
            # in which case the thread will stop running
            if args is None:
                self.input_queue.task_done()
                return
            else:
                index, metagenomic_bin = args
            # process the arguments, making sure to capture any exception
            # raised while processing the query, and then mark the hits
            # as "found" using a `threading.Event` for each.
            try:
                self.score(index, metagenomic_bin)
                self.input_queue.task_done()
            except BaseException as exc:
                print(exc)
                self.error = exc
                self.kill()
                return

    def kill(self):
        self.kill_switch.set()

    def score(self, i, metagenomic_bin):

        cdef int          ipath
        cdef bint         closed = self.closed
        cdef Genes        genes  = self.genes
        cdef Nodes        nodes  = self.nodes
        cdef Sequence     seq    = self.seq
        cdef TrainingInfo tinf   = metagenomic_bin.training_info

        # check which of the metagenomic bins gets the best results
        if tinf.gc < self.low or tinf.gc > self.high:
            return
        # recreate the node list if the translation table changed
        if tinf.translation_table != self.translation_table:
            self.translation_table = tinf.translation_table
            nodes._clear()
            _api.add_nodes(nodes, seq, tinf, closed=closed)
            nodes._sort()
        # compute the score for the current bin
        # with nogil:
        _api.reset_node_scores(nodes)
        _api.score_nodes(nodes, seq, tinf, closed=closed, is_meta=True)
        _api.record_overlapping_starts(nodes, tinf, is_meta=True)
        ipath = _api.dynamic_programming(nodes, tinf, is_meta=True)



        # update genes if the current bin had a better score
        if self.nodes[ipath].score > self.max_score.value:
            # record best phase and score
            self.max_phase.value = i
            self.max_score.value = self.nodes[ipath].score
            # eliminate eventual bad genes in the nodes
            _api.eliminate_bad_genes(self.nodes, ipath, tinf)
            with self.genes_lock:
                # clear the gene array
                self.genes.clear()
                # extract the genes from the dynamic programming array
                with nogil:
                    _api.add_genes(genes, nodes, ipath)
                    _api.tweak_final_starts(genes, nodes, tinf)
                    _api.record_gene_data(genes, nodes, tinf, sequence_index=1)


class Pyrodigal:

    def __init__(self, meta=False, closed=False):
        self.training_info = None
        self.meta = meta
        self.closed = closed
        self.__num_seq = 1

    def find_genes(self, sequence):
        if not self.meta and self.training_info is None:
            raise RuntimeError("cannot find genes without having trained in single mode")

        if isinstance(sequence, str):
            seq = _api.Sequence.from_string(sequence)
        else:
            seq = _api.Sequence.from_bytes(sequence)

        if self.meta:
            return self._find_genes_meta(seq)
        else:
            raise NotImplementedError("TODO")

    def _find_genes_meta(self, seq):
        predictions = _api.find_genes_meta(seq, closed=self.closed, sequence_index=self.__num_seq)
        self.__num_seq += 1
        return predictions
        # compute the min/max acceptable gc for the sequence to only
        # use appropriate metagenomic bins

        # # initialize stuff
        # cdef Nodes  nodes
        # cdef Genes  genes       = _api.Genes()
        # cdef object genes_lock  = threading.Lock()
        # cdef object kill_switch = threading.Event()
        # cdef object input_queue = queue.Queue()
        # cdef object max_phase   = multiprocessing.Value(ctypes.c_int)
        # cdef object max_score   = multiprocessing.Value(ctypes.c_double)
        #
        # # create threads
        # threads = []
        # for i in range(8): # NUM_THREADS
        #     threads.append(_MetaThread(
        #         seq,
        #         self.closed,
        #         input_queue,
        #         max_phase,
        #         max_score,
        #         genes,
        #         genes_lock,
        #         kill_switch
        #     ))
        #     threads[i].start()
        #
        # # fill input queue
        # for i, metagenomic_bin in enumerate(METAGENOMIC_BINS):
        #     input_queue.put((i, metagenomic_bin))
        # for i in range(8):
        #     input_queue.put(None)
        #
        # # wait for the thread to be done
        # input_queue.join()
        #
        # # recover the nodes corresponding to the best run
        # nodes = threads[0].nodes
        # training_info = METAGENOMIC_BINS[max_phase.value].training_info
        # nodes._clear()
        # _api.add_nodes(nodes, seq, training_info, closed=self.closed)
        # nodes._sort()
        # _api.score_nodes(nodes, seq, training_info, closed=self.closed, is_meta=True)
        #
        # #
        # return _api.Predictions(genes, nodes, seq, training_info)
