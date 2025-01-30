Reference Sequence Formats
=============================

This appendix explains the details of how |breseq| handles different reference sequence formats. Most importantly, this includes how different types of feature annotations are used to improve mutation predictions.

Each reference sequence file (the ``-r`` option to |breseq| and many |gdtools| subcommands) can contain **sequence** information (the nucleotide sequences of chromosomes or plasmids) and/or **annotations** (the locations and identities of features such as genes on those DNA sequences).

Three types of input files are accepted for reference sequences:

* `GenBank Format <https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html>`_ (sequences and/or annotations)
* `GFF3 Format <http://gmod.org/wiki/GFF3>`_ (sequences and/or annotations)
* `FASTA Format <https://www.ncbi.nlm.nih.gov/genbank/fastaformat/>`_ (sequences only)

Each loaded sequence is assigned a ``SEQ_ID`` as explained below. Sequences and their annotations can be input in different files as long as the ``SEQ_ID`` matches between files.

.. note::
  During a run, |breseq| merges and converts all input reference sequences into one annotated reference file that is output as ``data/reference.gff3``. If you are having trouble interpreting how |breseq| is loading your reference files, you should examine this file.

Sequences
------------------------
The header of each reference sequence and sometimes a special feature corresponding to the entire sequence are loaded to determine the **sequence id** (``SEQ_ID``). The **length** of the sequence is also provided in the header for some formats. Any length provided here will be checked against the actual nucleotide sequence that is loaded. Finally, the **topology** of the sequence can be set to linear (default) or circular as described below in certain formats.

FASTA
^^^^^^^^
The ``SEQ_ID`` is assigned as the first word on the sequence description line (i.e. all characters before encountering whitespace). Any later descriptive information on this line is ignored. Example:

``>SEQ_ID DESCRIPTON``

FASTA format does not support setting the sequence to have circular topology.

GenBank
^^^^^^^^

The ``SEQ_ID`` is assigned with this order of preference from the ``LOCUS`` > ``ACCESSION`` > ``VERSION`` lines. This behavior can be overriden with the ``--genbank-field-for-seq-id`` command-line option, which can have the values ``AUTOMATIC``, ``LOCUS``, ``VERSION``, ``ACCESSION``

The length provided in the ``LOCUS`` line  is used. If the ``source`` feature annotation has a different length, then a warning is shown. The ``LOCUS`` line will contain either ``LINEAR`` or ``CIRCULAR`` which sets the sequence topology.

GFF3
^^^^^

The line that begins with ``##sequence-region`` has this whitespace delimited format:

``##sequence-region SEQ_ID START END``

The ``SEQ_ID`` is taken from the first item. Then, the length is determined as ``END - START + 1``.

Sequences with a circular topology have the attribute ``Is_circular=true`` for the ``region`` feature that corresponds to the entire sequence .

Annotations
----------------------------

GenBank and GFF3 format support providing a list of feature annotations, which are sequence locations having start, end, and strand attributes that together define the bases constituting the feature.

Each feature may be composed of a list of one or more nucleotide segments which may be discontiguous (for example, in the case of an ORF defined by multiple exons on a spliced RNA). In some cases, the start or end position of a feature may be indeterminate (ambiguous) because the sequence fragment is truncated (for example, in a *de novo* assembly or draft genome sequence).

Both GenBank and GFF3 files define a type for each sequence feature and then have various information stored as key/value pairs.

The types that |breseq| recognizes are:

* ``CDS`` (protein-coding sequence)
* ``fCDS`` (fragmented CDS, for pseudogenes)
* ``rRNA`` (ribosomal RNA)
* ``tRNA`` (transfer RNA)
* ``ncRNA`` (noncoding RNA)
* ``RNA`` (generic RNA)
* ``mobile_element`` (transposon or other mobile DNA element)
* ``repeat_region`` (transposon or other mobile DNA element)

Features marked with a type that is only ``gene`` are not used on their own because it cannot be determined whether they encode a protein or are noncoding. If another identical feature of the other type exists, auxiliary information is loaded from the corresponding ```gene`` feature.

|breseq| will annotate the effects of mutations differently in features that are marked as pseudogenes rather than normal coding sequence (CDS) features. Pseudogenes can be marked as described below in each format. If a CDS is encountered that does not have a length that is a multiple of three, |breseq| fill print a warning that suggests adding the pseudo tag to that feature.

Internally, |breseq| tries to load three pieces of information describing each feature: ``name``, ``accession`` (like a unique ``locus_tag``), and ``description``.

|breseq| is able to more accurately predict the locations of **transposon insertions** if these elements are annotated in the reference genome. They must have a feature type of ``repeat_region`` or ``mobile_element`` to be recognized. The ends of these features should correspond to the entire unit that is inserted when the DNA "moves" (e.g., encompassing the inverted repeats on the end of an IS element and everything between them, not just the transposase gene). If there are multiple copies of an element in the genome, then all of them should have the exact same name (correct: ``IS150`` and ``IS150``; incorrect: ``IS150a`` and ``IS150b``). This is important for letting |breseq| match up junction evidence from different copies.

GenBank
^^^^^^^^

Genbank files can name features using many different tags. |breseq| uses this order of preference in deciding on the main name to use for a gene:

The ``name`` for a feature is determined by |breseq| by checking in this order for ``/name=``, ``/locus_tag=``, ``/label=``, and then ``/note=`` tags.

The ``accession`` is loaded from the ``/locus_tag=`` tag. (It may end up being the same as the ``name``.)

The ``product`` for a feature is assigned from the ``/product=`` tag if it exists, and then from the ``/note=`` tag as a backup.

Complex positions and indeterminate start/end positions are described in the line that gives the location of each feature according to the Genbank format specification.

Pseudogenes are CDS features marked by adding a line that consistes solely of the ``/pseudo`` tag.

GFF3
^^^^^

The ``name`` for a feature is determined by |breseq| by checking in this order for ``Name=``, ``gene=``, ``accession=`` attributes.

The ``accession`` is loaded in order of preference from the first attribute that exists from ``accession=``, ``locus_tag=``, ``ID=`` or ``Alias=``.

The ``product`` for a feature is assigned from the ``product=`` attribute if it exists, and then from the ``note=`` attribute as a backup.

If multiple feature lines have identical accessions and types, then the locations from each one are concatenated together in one feature. This is how you represent a programmed frameshift or exons in a spliced gene, for example. Indeterminate (ambiguous) start/end coordinates for a segment are specified by adding an ``indeterminate_coordinate=start`` or ``indeterminate_coordinate=end`` as an attribute to the semicolon-delimited list on the line for a location.

Pseudogenes are marked by adding ``Pseudo=true`` to the semicolon-delimited list of attributes at the end of the feature line line. Additionally, pseudogenes are reassigned a different feature type of ``fCDS``.

Adding IS Element Annotations
-------------------------------

Many sequence files don't have IS elements annotated. To have |breseq| automatically predict IS elements as single events versus two JC evidence items that you have to figure out, we highly recommend adding these annotations. You can accomplish it using these steps:

1. Install and run `ISEScan <https://github.com/xiezhq/ISEScan>`_ to generate a CSV file of IS predictions.

  .. code-block:: bash

    isescan.py --nthread 4 --seqfile reference.fasta --output output

  .. note::

    If you don't have a FASTA version of your reference, you can generate one using ```breseq CONVERT-REFERENCE``.

  .. note::

    The current version of ISEScan and its dependencies installed through Conda crashes on MacOSX due to a problem with FragGeneScan. If you really want to get it running on a Mac, you can install ISEScan via Conda and then install this fixed version inside of the same environment: `FragGeneScan with bug fix <https://github.com/barricklab/FragGeneScan>`_ .

2. Merge these predictions into your reference file using ```breseq CONVERT-REFERENCE``.

  .. code-block:: bash

    breseq CONVERT-REFERENCE -f GENBANK -s output/reference.fasta.csv -o reference_with_IS.gbk reference.gbk

  .. note::

    You can also output as a GFF3 (substitute ``-f GFF3`` and ``-o reference_with_IS.gff``).

3. Now run |breseq| with the updated reference file!

Illegal Characters
--------------------

For all sequence formats:

#. In nucleotide sequences, all characters are converted to uppercase and all non [``ATCG``] characters are converted to [``N``].
#. In gene names and locus tags, the characters [``,;/\|``] are replaced with [``_``].
#. In gene descriptions, the character [``|``] is replaced with [``;``].
