Reference Sequence Formats
=============================

This appendix explains the details of how |breseq| handles different reference sequence formats. Most importantly, this includes how different types of feature annotations are used.

Illegal Characters
--------------------

For all sequence formats:

#. In nucleotide sequences, all characters are converted to uppercase and all non [ATCG] characters are converted to [N].
#. In gene names and locus tags, the characters [,;/|] are replaced with [_].
#. In gene descriptions, the character [|] is replaced with [;].


Feature Annotations
----------------------------

|breseq| is able to more accurately predict the locations of transposon insertions if these elements are annotated in the reference genome. They must have a feature type of ``repeat_region`` or ``mobile_element``.
