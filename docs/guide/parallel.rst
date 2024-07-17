Parallelism
===========

.. currentmodule:: pyrodigal

Pyrodigal is entirely thread-safe, which allows processing different contigs 
in parallel with the same `GeneFinder` object. 


Command Line
------------

The command line supports processing input sequences in parallel, but this needs
to be enabled explicitly. Use the `-j` flag to specify any number of jobs to 
run in parallel, or `-j0` to run as many jobs as there are available cores on
the machine (as reported by `os.cpu_count`):

.. code:: console

    $ pyrodigal -j0 ...

Reentrancy
----------

The `GeneFinder.find_genes` method is re-entrant, so calling it across 
different threads doesn't cause any issue. The easiest way to call the 
`~GeneFinder.find_genes` method in parallel is with a 
`multiprocessing.pool.ThreadPool`, which can easily split the work into chunks
to be processed across different threads:

.. code:: python

    from multiprocess.pool import ThreadPool
    from pyrodigal import GeneFinder

    sequences = [ ... ]                  # a list of sequences to process
    gene_finder = GeneFinder(meta=True)  # a single gene finder object

    with ThreadPool() as pool:
        genes = pool.map(gene_finder.find_genes, sequences)

This is internally what the Pyrodigal CLI does when called with the `-j` flag
set to any number of jobs but 1.

Processes
---------

On some setups, such as workstations with virtualized CPUs, threads may not be 
as efficient because of the requirement for a shared memory space accessible by 
all threads which may cross physical boundaries. In that case, using processes
will be faster, despite the initial requirement to copy data in each worker process.

To use processes instead of threads in the command line, use the following flag:

.. code:: console

    $ pyrodigal --pool=process ...

In the API example from above, simply import `~multiprocessing.pool.Pool` from
`multiprocessing.pool` instead of `~multiprocessing.pool.ThreadPool`. 
