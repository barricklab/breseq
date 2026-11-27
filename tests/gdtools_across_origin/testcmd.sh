#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

# Every type of mutation crossing the origin of a circular sequence: what gdtools APPLY makes of
# it (GFF3 has the features and the sequence), and the VCF and GVF that describe the same change.
# The IS element for the MOBs comes from the second reference file.

# The date and the breseq version are in the VCF and GVF headers
DIFF_IGNORE='^##(fileDate=|source=|file-date|source-method)'

CASES="sub del inv amp con mob mob_del"
REFERENCE_ARG="-r ${DATADIR}/tmv_plasmid/tmv-plasmid.gbk -r ${DATADIR}/REL606/REL606.fragment.gbk"

TESTCMD="true"
i=0
for CASE in ${CASES}; do
  for FORMAT in gff3 vcf gvf; do
    CURRENT_OUTPUTS[$i]="${SELF}/output.${CASE}.${FORMAT}"
    EXPECTED_OUTPUTS[$i]="${SELF}/expected.${CASE}.${FORMAT}"
    i=$((i+1))
  done
  TESTCMD="${TESTCMD} \
    && ${GDTOOLS} APPLY -f GFF3 -s TMV-plasmid -o ${SELF}/output.${CASE}.gff3 ${REFERENCE_ARG} ${SELF}/${CASE}.gd \
    && ${GDTOOLS} CONVERT -f VCF -o ${SELF}/output.${CASE}.vcf ${REFERENCE_ARG} ${SELF}/${CASE}.gd \
    && ${GDTOOLS} CONVERT -f GVF -o ${SELF}/output.${CASE}.gvf ${REFERENCE_ARG} ${SELF}/${CASE}.gd \
    "
done

# A junction across the origin of a circular sequence whose origin lies INSIDE a repeat.
#
# circular_origin_repeat.fna has three 2000-bp sequences that each begin with a 300-bp repeat whose
# second, identical copy is at 1001-1300. FASTA cannot say circular, so the first step makes them
# circular with '#=TOPOLOGY' header lines (origin_repeat.topology.gd) and writes GFF3 -- which is
# also the test that a sequence read from FASTA, with no 'region' feature of its own, gets one
# written carrying Is_circular=true.
#
# origin_repeat.evidence.gd is what breseq reports when the reads cross such an origin: the side
# that enters the repeat matches both copies, is flagged redundant, and comes back at whichever
# copy the aligner returned -- here the FAR one. gdtools MUTATIONS predicts from that evidence
# alone, so no reads are needed. The predictor must re-place that side at the origin:
#
#   origin-flush    JC 2000 (-1) / 1001 (+1), no overlap. Becomes 2000 / 1, which is just the
#                   origin: ignore=CIRCULAR_CHROMOSOME, and no mutation.
#   origin-shifted  JC 2000 (-1) / 1002 (+1), 1-bp overlap: the reads have one T fewer than the
#                   four that run across the origin. Becomes 2000 / 2, which is a 1-bp deletion
#                   sitting on the origin: DEL at position 1, size 1.
#
#   origin-inplace  JC 2000 (-1) / 1 (+1), no overlap, and the redundant side is ALREADY at the
#                   origin because that is the copy the aligner returned. Nothing to re-place,
#                   but it is still just the origin: ignore=CIRCULAR_CHROMOSOME.
#
# and two MC items, for the missing coverage such an origin produces. The piece of an
# origin-crossing read on a junction's redundant side is not written to the reference alignment, so
# position 1 has no coverage of any kind, which seeds an MC that then runs through the whole repeat
# (unique coverage is zero across it):
#
#   origin-inplace  MC 1-300, all but 2 bases in its redundant range. An artifact of the origin:
#                   ignore=CIRCULAR_CHROMOSOME.
#   origin-flush    MC 1-340, 40 bases NOT in its redundant range. That is real missing coverage
#                   next to the repeat and must be left alone.
#
# The last MUTATIONS step predicts again FROM THE OUTPUT of the first, where both sides are already
# at the origin. The breseq pipeline does exactly that -- junctions are prepared once before the
# evidence is saved and again when mutations are predicted -- and the deletion was lost on the
# second pass, because what marked the junction as re-placed was not saved with it.
#
# Both were left as unexplained junctions to the far copy (seen on a PacBio HiFi run of a
# chromosome whose position 1 starts an rRNA operon). side_1_placed_at_origin_from records where
# the side was reported.
CURRENT_OUTPUTS[$i]="${SELF}/output.origin_repeat.gff3"
EXPECTED_OUTPUTS[$i]="${SELF}/expected.origin_repeat.gff3"
i=$((i+1))
CURRENT_OUTPUTS[$i]="${SELF}/output.origin_repeat.gd"
EXPECTED_OUTPUTS[$i]="${SELF}/expected.origin_repeat.gd"
i=$((i+1))
TESTCMD="${TESTCMD} \
  && ${GDTOOLS} APPLY -f GFF3 -o ${SELF}/output.origin_repeat.gff3 \
       -r ${DATADIR}/circular_origin_repeat/circular_origin_repeat.fna ${SELF}/origin_repeat.topology.gd \
  && ${GDTOOLS} MUTATIONS -o ${SELF}/output.origin_repeat.gd \
       -r ${SELF}/output.origin_repeat.gff3 ${SELF}/origin_repeat.evidence.gd \
  && test \$(grep -c 'ignore=CIRCULAR_CHROMOSOME' ${SELF}/output.origin_repeat.gd) -eq 3 \
  && test \$(cut -f1,4,5,6 ${SELF}/output.origin_repeat.gd | grep -c '^DEL.origin-shifted.1.1\$') -eq 1 \
  && test \$(grep -v '^#' ${SELF}/output.origin_repeat.gd | cut -f1 | grep -cE '^(SNP|SUB|DEL|INS|MOB|AMP|INV|CON|INT)\$') -eq 1 \
  && ${GDTOOLS} MUTATIONS -o ${SELF}/output.origin_repeat.second_pass.gd \
       -r ${SELF}/output.origin_repeat.gff3 ${SELF}/output.origin_repeat.gd \
  && test \$(grep -c 'ignore=CIRCULAR_CHROMOSOME' ${SELF}/output.origin_repeat.second_pass.gd) -eq 3 \
  && test \$(grep '^MC' ${SELF}/output.origin_repeat.gd | grep -c 'origin-flush.*ignore=') -eq 0 \
  && test \$(cut -f1,4,5,6 ${SELF}/output.origin_repeat.second_pass.gd | grep -c '^DEL.origin-shifted.1.1\$') -eq 1 \
  "

do_test $1 ${SELF}
