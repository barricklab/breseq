#!/bin/bash
#
# Read type detection and the <read type>:<path> prefix (see READ TYPES in src/breseq/fastq.h).
#
# Builds small read files from the lambda reads, with the read names and base qualities rewritten
# to look like each kind of read, and writes the type breseq detects for each to output.txt. The
# files whose names identify the platform must be detected from the "header"; the others, as an
# SRA download with the original names gone would be, from "quality" (lengths and base qualities).
#
# Then checks how the prefix is parsed: nothing here runs the pipeline, because a prefix that is
# not a read type, and a pair given two different types, must both stop breseq before it starts.
#
# $1 == this test's directory, $2 == breseq binary, $3 == tests/data
set -u

SELF="$1"
BRESEQ="$2"
DATADIR="$3"

WORK="${SELF}/data"            # gitignored, and removed by common.sh's do_clean
OUT="${SELF}/output.txt"
mkdir -p "${WORK}"
: > "${OUT}"

READS="${WORK}/reads.fastq"
gunzip -c "${DATADIR}/lambda/lambda_mixed_population.fastq.gz" | head -8000 > "${READS}"

# Rewrites the @NAME line of every read. $1 = output name, $2 = awk expression for the new line
# (n is the read number).
rename_reads() {
  awk "NR % 4 == 1 { n++; print $2; next } { print }" "${READS}" > "${WORK}/$1.fastq"
}

# Joins every 40 reads into one long read, named as an SRA download with the original name gone.
# $1 = output name, $2 = string of quality characters cycled along each read.
long_reads() {
  awk -v qchars="$2" '
    NR % 4 == 2 { seq = seq $0; k++ }
    k == 40 && NR % 4 == 0 {
      n++; q = ""; len = length(seq); nq = length(qchars)
      for (i = 1; i <= len; i++) q = q substr(qchars, (i + n) % nq + 1, 1)
      printf "@SRR0000001.%d %d length=%d\n%s\n+\n%s\n", n, n, len, seq, q
      seq = ""; k = 0
    }' "${READS}" > "${WORK}/$1.fastq"
}

UUID='sprintf("%08x-0000-4000-8000-%012x", n, n)'

# ---- Identified by the read names
cp "${READS}" "${WORK}/illumina_old_names.fastq"
rename_reads illumina_casava_names '"@HWI-ST278R:417:D29LHACXX:2:1101:" n ":2224 1:N:0:ATCACG"'
rename_reads illumina_sra_names    '"@SRR0000001." n " HWI-ST278R:417:D29LHACXX:2:1101:" n ":2224/1"'
rename_reads nanopore_fast_model   '"@" '"${UUID}"' " runid=0a1b ch=12 basecall_model_version_id=dna_r10.4.1_e8.2_400bps_fast@v4.3.0"'
rename_reads nanopore_hac_model    '"@" '"${UUID}"' " runid=0a1b ch=12 basecall_model_version_id=dna_r10.4.1_e8.2_400bps_hac@v4.3.0"'
# super accuracy is reported as nanopore-hac; this is also the TAB-separated form samtools fastq -T writes
rename_reads nanopore_sup_bam_tag  '"@" '"${UUID}"' "\tRG:Z:0a1b_dna_r10.4.1_e8.2_400bps_sup@v5.0.0"'
rename_reads pacbio_ccs_names      '"@m64011_190830_220126/" n "/ccs"'
rename_reads pacbio_subread_names  '"@SRR0000001." n " m54119U_210220_034433/" n "/0_35 length=35"'

# ---- Names gone: lengths and qualities
rename_reads short_no_names        '"@SRR0000001." n " " n " length=35"'
long_reads long_q0                 '!'
long_reads long_rs2_qualities      '!"#$%&()*+,-./0'
long_reads long_hifi_qualities     '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~5:?DINSX]bglqv{'
long_reads long_binned_qualities   '(/5<CI'
long_reads long_low_qualities      '#$%&'"'"'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQR'
long_reads long_high_qualities     '.0123456789:;<=>?@ABCDEFGHIJKLMNOPQRS/0123456789:;<=>?@ABCDEFGHI'
# an ONT read name but no basecaller model: the platform from the name, the accuracy from qualities
awk 'NR % 4 == 1 { n++; printf "@%08x-0000-4000-8000-%012x\n", n, n; next } { print }' "${WORK}/long_high_qualities.fastq" > "${WORK}/nanopore_names_no_model.fastq"

echo "# file / detected read type / decided from" >> "${OUT}"
for name in illumina_old_names illumina_casava_names illumina_sra_names \
            nanopore_fast_model nanopore_hac_model nanopore_sup_bam_tag \
            pacbio_ccs_names pacbio_subread_names \
            short_no_names long_q0 long_rs2_qualities long_hifi_qualities long_binned_qualities \
            long_low_qualities long_high_qualities nanopore_names_no_model; do
  "${BRESEQ}" SUMMARIZE-FASTQ "${WORK}/${name}.fastq" 2>/dev/null \
    | tr ',{}' '\n\n\n' | tr -d '"' | awk -F: -v name="${name}" '
        $1 == "read_type" { type = $2 } $1 == "read_type_source" { source = $2 }
        END { printf "%-24s %-14s %s\n", name, type, source }' >> "${OUT}"
done

# ---- The prefix
echo "# read type prefixes" >> "${OUT}"
REF="${DATADIR}/lambda/lambda.gbk"

# $1 = description, $2 = text that must be in the output, rest = read file arguments.
# --dry-run validates the command line and stops, so no pipeline is run and nothing is written.
check_prefix() {
  local description="$1" expect="$2"; shift 2
  "${BRESEQ}" --dry-run -o "${WORK}/dry_run" -r "${REF}" "$@" > "${WORK}/prefix.log" 2>&1
  local status=$?
  if grep -qF -- "${expect}" "${WORK}/prefix.log"; then
    echo "${description}: exit $(( status != 0 )), message found" >> "${OUT}"
  else
    echo "${description}: exit $(( status != 0 )), MESSAGE NOT FOUND" >> "${OUT}"
  fi
}

cp "${READS}" "${WORK}/pair_R1.fastq"
cp "${READS}" "${WORK}/pair_R2.fastq"

check_prefix "known type"            "fastq [input file]"                             "nanopore-hac:${WORK}/long_high_qualities.fastq"
check_prefix "synonym, any case"     "fastq [input file]"                             "ONT-SUP:${WORK}/long_high_qualities.fastq"
check_prefix "unknown type"          "Unrecognized read type prefix \"ont-hca\""      "ont-hca:${WORK}/long_high_qualities.fastq"
check_prefix "missing file"          "does not exist"                                 "illumina:${WORK}/no_such_file.fastq"
check_prefix "pair, one type given"  "fastq [input file]"                             "illumina:${WORK}/pair_R1.fastq" "${WORK}/pair_R2.fastq"
check_prefix "pair, two types"       "Paired read files were given different read types" "illumina:${WORK}/pair_R1.fastq" "pacbio-hifi:${WORK}/pair_R2.fastq"

exit 0
