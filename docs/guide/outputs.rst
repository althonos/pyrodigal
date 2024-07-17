Output Formats
==============

.. currentmodule:: pyrodigal

Pyrodigal was created among other things to skip parsing the results
from Prodigal while building a larger pipeline. Therefore, it's advised
to manipulate the genes predicted by Pyrodigal through the object layer
instead of writing them to a file only to parse it later.

Nevertheless, Pyrodigal provides convenience methods to write the genes
from a `~pyrodigal.Genes` collection to a variety of formats, as found in the
Prodigal command line. Note that all these functions need to know the
name of the sequence they were obtained from.

.. note::

    The methods shown before are implemented using only the Python
    interface, so it should be possible to replicate them using the
    attributes of the individual `Gene` objects.


Genes
-----

Use `Genes.write_genes` to write the nucleotide sequences of all the genes
found by Pyrodigal to a file or a file-like object, similar to ``prodigal -d``:

.. code:: python

    genes = orf_finder.find_genes(sequence)
    with open("genes.fna", "w") as dst:
        genes.write_genes(dst, sequence_id="seqXYZ")

It's also possible to change the width of the FASTA columns,
if needed (use `math.inf` as the width to write 2-line FASTA):

.. code:: python

    genes = orf_finder.find_genes(sequence)
    with open("genes.fna", "w") as dst:
        genes.write_genes(dst, sequence_id="seqXYZ", width=80)


Translations
------------

Use `Genes.write_translations` to write the protein sequences of all the
genes found by Pyrodigal to a file or file-object, similar to ``prodigal -a``:

.. code:: python

    genes = orf_finder.find_genes(sequence)
    with open("proteins.faa", "w") as dst:
        genes.write_translations(dst, sequence_id="seqXYZ")


It's also possible to change the translation table and to control the
width of the FASTA columns, if needed:

.. code:: python

    genes = orf_finder.find_genes(sequence)
    with open("proteins.faa", "w") as dst:
        genes.write_translations(dst, sequence_id="seqXYZ", width=80, translation_table=11)

GFF
---

Use `Genes.write_gff` to write the genes to a file or file-like object in GFF
format, similar to ``prodigal -f gff -o``:

.. code:: python

    genes = orf_finder.find_genes(sequence)
    with open("genes.gff", "w") as dst:
        genes.write_gff(dst, sequence_id="seqXYZ")

Use ``header=False`` to skip writing the GFF3 header. This can be useful when
writing genes from different sequences to the same result file:

.. code:: python

    with open("genes.gff", "w") as dst:
        for i, record in enumerate(Bio.SeqIO.parse("contig.fna")):
            genes = orf_finder.find_genes(str(record.seq))
            genes.write_gff(dst, sequence_id=record.id, header=(i == 0))

Use ``include_translation_table=True`` to include the translation table
of the genes in the GFF attributes. This may be useful for genes predicted in
metagenomic mode, since they may have been predicted with non-standard
genetic codes.

Genbank
-------

Use `Genes.write_genbank` to write the genes to a file or file-like object
in `GenBank <https://www.insdc.org/submitting-standards/feature-table/#3.1>`_
format, similar to ``prodigal -f gbk -o``:

.. code:: python

    genes = orf_finder.find_genes(sequence)
    with open("genes.gbk", "w") as dst:
        genes.write_genbank(dst, sequence_id="seqXYZ")

.. note::

    The original Prodigal outputs incomplete GenBank files containing only
    the coordinates of the predicted genes inside
    `CDS <https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html#CDSB>`_ features,
    without including the translation or the original sequence. Since this is
    not the most useful output, and often requires additional post-processing,
    Pyrodigal outputs a complete GenBank record instead.


Scores
------

Use `Genes.write_scores` to write the scores of all potential genes in a format
proprietary to Prodigal, similar to ``prodigal -s``:

.. code:: python

    genes = orf_finder.find_genes(sequence)
    with open("scores.txt", "w") as dst:
        genes.write_scores(dst, sequence_id="seqXYZ")





