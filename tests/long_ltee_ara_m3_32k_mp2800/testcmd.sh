#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=7
# Large real-world data, downloaded on demand and md5-verified rather than
# committed. See tests/long_ltee_clone/testcmd.sh for what TEST_DATA does and
# why it has to be declared here, before common.sh is sourced.
TEST_DATA=ena_SRR098033,ltee_REL606
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"

# LTEE clone ZDB30: population Ara-3, generation 32,000 (Blount et al. 2012).
# Paired-end Illumina Genome Analyzer II, 25,400,140 read pairs.
#
# READ LENGTHS ARE 36 AND 34, NOT 2x35, and the ENA record looks like an off-by-one:
# R1 is 36 bp, R2 is 34, and they sum to the 70 a 2x35 run would give. R1's 36th base
# is 95% T at an ordinary Q32.2 -- not a dead cycle, so --quality-score-trim will
# never remove it and breseq has no fixed-end trim. R2's last cycle is normal. Both
# 5' ends are T-skewed as well (R1 cycle 1 is 80% T, cycle 2 66%).
#
# It costs R1 about 1.1 points of alignment rate (88.62% against R2's 92.52%; a
# direct bowtie2 comparison recovers 224 of 20,000 reads by trimming that one base).
# Do NOT trim it here: a run on reads trimmed one base at each end gives FEWER
# junctions, not more (JC 10 -> 7, still no MOB), because the shorter reads cost more
# than the bad base does. It is also not a meaningful contributor to MP evidence --
# if it were, the ~58,000 R1 reads it costs would show up as mate-never-aligned, and
# the whole run carries 7,525.
#
# WHY THIS DATA SET: insert size, and nothing else comes close. ENA records
# nominal_length 2834 for this library, so with 2x35 reads each pair spans about
# 2.76 kb of sequence that was never read. Every other paired data set in the
# LTEE collection -- and every other long test -- is a few hundred bases at
# most, so this is the only place the pair-distance and discordant-pair machinery
# sees a kilobase-scale span. Alongside long_ltee_ara_m1_40k_pe36 (135 bp, also
# recorded) it brackets the range: the two libraries whose insert sizes we
# actually know, twenty-fold apart.
#
# It is deliberately a second Ara-3 clone (long_ltee_ara_m3_38k_se36 is ZDB107 at
# 38,000 generations). Duplicating the population is the price of the insert
# size; the two clones' curated genome diffs are unalike anyway -- this one is
# MOB-heavy with no AMP, that one is AMP-heavy.
#
# COVERAGE CAP: as served this run is ~384x, far above the 40-80x band the rest
# of these tests were chosen for, and at ~35 bp per read that is 50 million reads. '-l 80'
# caps what breseq PARSES at roughly 80-fold. The limit is applied inside
# normalize_fastq_paired (see read_file_base_limit in breseq_cmdline.cpp), so it
# truncates within the paired set -- it takes the head of the files in order and
# never splits a pair, which is what makes the golden reproducible. The full
# 1.39 GB is still downloaded and md5-verified: the manifest checksum has to
# cover exactly the bytes ENA serves, which is why the data is not subsampled on
# disk instead. This is also the largest download of any test in the tree.
#
# EXPERIMENTAL PAIR EVIDENCE: this test runs the three PAIR predictors, and it is
# the only place any of them see a kilobase-scale span. What they currently produce
# here is MP 0, PD 28 (25 accepted) and DP 0 -- the zeros are expected and worth
# pinning.
#
# THE PD ENTRIES ARE NOT OVER-PREDICTION, which is the natural first reading of 25
# calls against 2 to 7 on the short-insert paired tests. Every accepted one has
# independent corroboration:
#
#   * 21 have a NEGATIVE size_shift matching an annotated IS element's length, i.e.
#     an insertion. 12 carry a repeat_name (6 of those anchored on a junction whose
#     other side is in the element, the rest named by size where exactly one family
#     fits); the other 9 carry an honest repeat_size_candidates list because IS150
#     (1443) and IS186 (1343) are closer together than the estimator's scatter.
#     13 of the 21 match a MOB in long_ltee_ara_m3_38k_se36 -- the same Ara-3
#     population 6,000 generations later -- to within a few bases and the right
#     family.
#   * 4 have a POSITIVE shift and all four sit on a missing-coverage region, i.e.
#     they are deletions. Two are multi-kilobase (MC 1609173-1615472 and
#     MC 3894999-3901455) and PD did not see them at all before.
#
# WHY '--junction-alignment-pair-limit 2000000' IS PASSED HERE. At the stock 100000
# this run examines 4.7% of the junction evidence it has: the reference offers
# 2,139,024 passed alignment pairs and the cap stops at 100,000. The cost is severe
# and specific to this library -- JC 10 instead of 37 and MOB 0 instead of 4 -- and
# it is not a read-length problem, which is what it looks like at first. Mate-pair
# circularization junctions are the cause: soft clipping finds 2,610 clipped tails
# matching back 2.4-3.3 kb (a 1305x enrichment, see below), and each is a ONE-OFF
# molecule, so each yields a candidate with position-hash 1 that the >=2 threshold
# discards. They flood the budget and dilute the real IS junctions -- which do have
# many reads sharing one breakpoint -- below the threshold. Raising the cap costs
# nothing measurable in runtime (12 min against 14 at the default) and the curve is
# flat past 2M: unlimited gives only JC 37 / MOB 4 against 36 / 3.
#
# Even with every junction, only 4 of the 25 IS insertion sites that HAVE a junction
# get a second one, so predictJCplusJCtoMOB fires 4 times while PD sees 21. That gap
# is not going to be closed by pairing a lone junction with a PD call: a MOB has to
# be placed to base-pair resolution and PD cannot do that -- its position_range is
# tens of bases and its size estimates scatter by 4-72 against a known element
# length. So the remaining insertions stay PD evidence, annotated with the family
# they are consistent with, and are not promoted to mutations.
#
# PD was 6 here until 6a303834 loosened the seed in favour of a genome-wide e-value;
# the 24 it gained were overwhelmingly these real events. It then went 30 -> 28 when
# PD stopped counting pairs whose placement was CHOSEN rather than measured
# (kBreseqAmbiguousPairPlacementBAMTag). That removed 8 artifacts -- 5 of them in
# rRNA operons, where a mate's copy is picked by a per-locus cluster vote so the
# distance error is coherent across every pair there and looks exactly like a real
# collective shift -- and, by sharpening the statistic, revealed 5 real events
# including the two multi-kilobase deletions above. 3.66% of pairs are refused here,
# which is the worst case in the tree: a repeat-rich reference and a 2.8 kb insert.
#
# MP was 98 (73 accepted) until MP started counting only reads whose mate produced NO
# ALIGNMENT AT ALL, rather than every read flagged BAM_FMUNMAP. Those are very
# different populations here: this library is 2x35, so --require-match-fraction 0.9
# demands 32 of 35 bases align, and 1,171,383 of the reads in reference.bam carry
# BAM_FMUNMAP while only 7,525 -- 0.64% -- have a mate bowtie2 truly placed nowhere.
# The other 99.36% were mates the aligner DID place and breseq then rejected, i.e.
# partially-aligning reads, which is SC's signal and not evidence that a mate's
# sequence is missing from the reference. Shorter reads make this worse, which is why
# this test had 98 MP entries and the 2x101/2x150 paired tests had 1-2.
#
# 0 is the right answer, not a loss of sensitivity. The 3 calls here that used to look
# real (2121958, 3901476, 4561292; frequency 1.0000, spanning_pair_count=0) have only
# 2, 17 and 4 genuinely-unplaced mates within 3 kb, against ~10 expected at the
# genome-wide rate. A real clonal insertion costs essentially every spanning fragment
# its mate -- some 2,900 reads in a window that size -- so these were positions where
# every crossing mate happened to be rejected, which is exactly what a frequency of
# 1.0 with zero opposing pairs looks like when the numerator is an artifact. The seed
# now never fires at all here (0 candidate regions), so the score gate is not even
# reached. The fitted null is p0 = 0.00073 with Pearson phi 1.35 -- near-binomial.
#
# The DP zero has a different cause: the distance distribution is unimodal about
# 2837 with a discordance cutoff of ~5381, so essentially no pair is an outlier
# INDIVIDUALLY, which is precisely the regime PD exists for and DP cannot see. If a
# change ever makes DP fire in bulk here, it is calling library structure, not
# variants. The golden will move whenever the DP/MP/PD code changes; that is the
# point, so rebuild it and review the diff rather than treating the failure as a
# surprise.
#
# SOFT CLIPPING IS DELIBERATELY OFF HERE, unlike the other three paired long tests.
# On a mate-pair library SC does not measure the genome, it measures the library
# prep, and it buries everything else:
#
#   * 15,086 SC entries against 21-46 on the other paired tests, 4,351 of them
#     accepted, spread evenly at roughly one per 1,064 bp of the 4.6 Mb genome, and
#     97% with consensus_fraction 1.0000. The golden was 4.1 MB, about 90% SC.
#   * They are real, and they are circularization junctions. Matching each clipped
#     tail back to the reference puts 2,610 of them 2.4-3.3 kb from their own SC
#     position -- against 2 expected by chance, a 1305x enrichment -- and 91% land
#     within 2.4-5 kb, centred on this library's 2837 bp fragment size. A soft clip
#     here is the far end of the same fragment, joined during circularization.
#
# Trimming would not have helped: the 2011 Illumina mate-pair v2 junction is a
# direct genomic-to-genomic blunt ligation with no adapter in it. Scanning 400,000
# raw reads finds zero Nextera Mate Pair junction adapter and 0.01% TruSeq
# read-through, so there is nothing for fastp to cut, and NxTrim-style splitters key
# on exactly the junction adapter this protocol lacks. SC is behaving correctly; the
# artifact is simply indistinguishable from structural variation without pairing
# information, which is why it is the pair predictors that earn their keep here.
#
# Consequence to be aware of when comparing goldens: --predict-soft-clipping lowers
# require-match-fraction from 0.9 to 0.5 unless set explicitly (settings.cpp). This
# test does NOT pass it, so it maps at the stock 0.9 while the other three paired
# long tests map at 0.5.
#
# Also covers the manB/cpsG paralogous deletion at LOW coverage (~10x) and from a
# badly fragmented MC. The single large MC 2032710-2054821 implies a register of
# 22112, which is 1181 bp short of the true 23293 -- the register search recovers
# the right one anyway (identity 0.947 against a best decoy of 0.411), so
# DEL 2032711 23293 is predicted here with no JC and no DP.
#
# That deletion also REMOVES a SNP this test used to call: cpsG T423S at 2054830.
# It was never a point mutation. 2054830 sits in the right copy; its homologue in
# the left copy is 2031537, one of the columns where the two paralogs differ, and
# it lies before the crossover -- so the surviving hybrid carries the left copy's
# A there, and reads placed on the right copy read as T->A. Without the deletion
# call breseq promoted that to an annotated amino-acid change in a gene that is in
# fact deleted. Do not "restore" it.
#
# See tests/long_ltee_clone/testcmd.sh for the notes that apply to every long
# test: rebuild with BRESEQ_TEST_DATA_DIR unset so the committed expected.gd
# keeps repo-relative paths, and note that #=MAPPED-BASES / #=MAPPED-READS are
# compared, so a bowtie2 version change fails this test even with breseq
# unchanged.
REFERENCE_ARG="-r ${DOWNLOADDIR}/ltee_REL606/REL606.gbk"

TESTCMD="\
    ${BRESEQ} \
    ${BRESEQ_TEST_THREAD_ARG} \
    -o ${SELF} \
    ${REFERENCE_ARG} \
    -l 80 \
    --junction-alignment-pair-limit 2000000 \
    --predict-copy-number \
    --predict-discordant-pairs \
    --predict-missing-pairs \
    --predict-pair-distance \
    ${DOWNLOADDIR}/ena_SRR098033/SRR098033_1.fastq.gz \
    ${DOWNLOADDIR}/ena_SRR098033/SRR098033_2.fastq.gz \
    "

do_test $1 ${SELF}
