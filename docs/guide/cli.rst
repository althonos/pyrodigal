Command-Line interface
======================

Pyrodigal comes with a command-line interface (CLI) that can be used as a 
drop-in replacement to the original ``prodigal`` executable. While using the 
API is recommended to get the most out of Pyrodigal, the CLI can be useful
to benefit from performance enhancements of Pyrodigal inside workflows or 
shell scripts.

Changes from Prodigal
---------------------

The Pyrodigal CLI behaves mostly like the Prodigal CLI, with the following
exceptions:

- The input file for Pyrodigal can only be in FASTA format, GenBank is not 
  supported.
- Pyrodigal supports getting the input sequences through a pipe to *stdin*, 
  however the stream cannot be compressed.
- The GenBank output of Pyrodigal is a full GenBank record including the 
  input sequence, unlike Prodigal which only outputs the features section.

Flags
-----

The Pyrodigal CLI has all the flags of the original CLI:

.. code-block:: text

    -a trans_file         Write protein translations to the selected file.
    -c                    Closed ends. Do not allow genes to run off edges.
    -d nuc_file           Write nucleotide sequences of genes to the selected file.
    -f output_type        Select output format.
    -g tr_table           Specify a translation table to use.
    -i input_file         Specify FASTA input file.
    -m                    Treat runs of N as masked sequence and don't build genes across them.
    -n                    Bypass Shine-Dalgarno trainer and force a full motif scan.
    -o output_file        Specify output file.
    -p mode               Select procedure.
    -s start_file         Write all potential genes (with scores) to the selected file.
    -t training_file      Write a training file (if none exists); otherwise, read and use the specified training file.

In addition, the following *new* flags can be used to control the new features
of Pyrodigal, such as multi-threading or smaller gene prediction:

.. code-block:: text

    -j jobs, --jobs jobs           The number of threads to use if input contains multiple sequences.
    --min-gene MIN_GENE            The minimum gene length.
    --min-edge-gene MIN_EDGE_GENE  The minimum edge gene length.
    --max-overlap MAX_OVERLAP      The maximum number of nucleotides that can overlap between two genes on the same strand.
                                   This must be lower or equal to the minimum gene length.
    --no-stop-codon                Disables translation of stop codons into star characters (*) for complete genes.
    --pool {thread,process}        The sort of pool to use to process genomes in parallel. Processes may be faster than
                                   threads on some machines, refer to documentation. (default: thread)

Piping
------

Pyrodigal supports reading the input through a pipe:

.. code-block:: console

  $ wget 'https://example.com/genome.fna' -O- | pyrodigal -a proteins.faa

If an input is given with the `-i` flag, it takes priority over the pipe.