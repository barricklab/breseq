Methods
==============

This section describes the algorithms used by :program:`breseq` in detail for those who want to know how it predicts mutations.


New junction (NJ) evidence
-----------------------------


Read alignment (RA) evidence
------------------------------

Read end trimming
*****************

Alignments of the ends of short reads.


1.	Examining portion of read that matches reference
2.	Examining portion of reference matched by read
3.	Examining entire read
4.	Examining extended region of reference

Note to self: we really need to do this once for the entire genome,
keeping track of how much to trim when end is there from both sides.
Then we can load this and use it.


Missing coverage (MC) evidence
------------------------------