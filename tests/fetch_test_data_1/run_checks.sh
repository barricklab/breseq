#!/bin/bash
#
# Offline self-test for tests/fetch_test_data.sh.
#
# Exercises the download/verify/cache logic end to end using curl's file://
# support against small files already in the repository, so it runs in well under
# a second, needs no network, and is safe to run in CI on every pull request --
# unlike the 'long' tests the fetcher actually serves.
#
# What it pins down (each of these is a way the fetcher could silently do the
# wrong thing, and the expensive real downloads would not catch any of them
# quickly):
#   * a good file is fetched and md5-verified                     -> exit 0
#   * a second run trusts the sidecar stamp                       -> [cached]
#   * a CHECKSUM MISMATCH is fatal AND leaves nothing in the cache -> exit 3
#   * an unreachable URL is fatal AND leaves nothing in the cache  -> exit 2
#   * tampering with a cached file invalidates the stamp          -> re-fetched
#   * a malformed manifest is rejected by --lint                  -> exit 1
#
# $1 == this test's directory
set -u

SELF="$1"
TESTS_DIR="${SELF}/.."
FETCH="${TESTS_DIR}/fetch_test_data.sh"

WORK="${SELF}/data"            # gitignored, and removed by common.sh's do_clean
CACHE="${WORK}/cache"
MANIFEST="${WORK}/fixture_manifest.tsv"
BAD_MANIFEST="${WORK}/bad_manifest.tsv"
OUT="${SELF}/output.txt"

rm -rf "${WORK}"
mkdir -p "${CACHE}"

# Absolute paths for the file:// URLs, computed here rather than committed,
# since they necessarily differ per machine and checkout.
DATA_ABS=$(cd "${TESTS_DIR}/data" && pwd)
GOOD_GBK="${DATA_ABS}/lambda/lambda.3.gbk"
GOOD_GZ="${DATA_ABS}/lambda/lambda_mixed_population.3.fastq.gz"

GOOD_GBK_MD5=7f0ebcf4ec51b462a3d11321c6e80d02
GOOD_GZ_MD5=c5fe95ecba6816bc516d9e4b47709a4a
WRONG_MD5=00000000000000000000000000000000

printf '#dataset\tfile\tmd5\turl\n'                                     > "${MANIFEST}"
printf 'good_plain\tlambda.3.gbk\t%s\tfile://%s\n' \
	"${GOOD_GBK_MD5}" "${GOOD_GBK}"                                >> "${MANIFEST}"
printf 'good_gz\treads.3.fastq.gz\t%s\tfile://%s\n' \
	"${GOOD_GZ_MD5}"  "${GOOD_GZ}"                                 >> "${MANIFEST}"
printf 'bad_md5\tlambda.3.gbk\t%s\tfile://%s\n' \
	"${WRONG_MD5}"    "${GOOD_GBK}"                                >> "${MANIFEST}"
printf 'missing\tnope.gbk\t%s\tfile://%s\n' \
	"${GOOD_GBK_MD5}" "${DATA_ABS}/lambda/this-file-does-not-exist" >> "${MANIFEST}"

# A manifest with a dataset name containing a '.', which tests/Snakefile's
# TEST_DATA regex could not express, and a truncated row.
printf '#dataset\tfile\tmd5\turl\n'                                      > "${BAD_MANIFEST}"
printf 'bad.name\tx.gbk\t%s\tfile://%s\n' \
	"${GOOD_GBK_MD5}" "${GOOD_GBK}"                                 >> "${BAD_MANIFEST}"
printf 'no_url\ty.gbk\t%s\thttp://insecure.example/y.gbk\n' \
	"${GOOD_GBK_MD5}"                                               >> "${BAD_MANIFEST}"

fetch() { # $1 = label, rest = args to fetch_test_data.sh
	local label="$1"; shift
	local output status
	output=$("${FETCH}" --manifest "${MANIFEST}" --dir "${CACHE}" "$@" 2>&1)
	status=$?
	# Keep only the bracketed status verbs and drop everything after them: the
	# '[fetching]' lines embed absolute file:// URLs and the '[fetched]' lines a
	# byte count, neither of which belongs in a committed expected output.
	echo "${label}: exit=${status}"
	echo "${output}" | sed -n 's/^\(\[[a-z]*\]\) *\([^ ]*\).*/  \1 \2/p' | grep -v '^  \[fetching\]'
}

exists() { # $1 = label, $2 = path
	if [[ -e "$2" ]]; then echo "$1: present"; else echo "$1: absent"; fi
}

{
	echo "=== a good file is downloaded and verified ==="
	fetch "fetch good_plain" good_plain
	exists "  cached file" "${CACHE}/good_plain/lambda.3.gbk"
	exists "  stamp" "${CACHE}/good_plain/lambda.3.gbk.md5ok"

	echo
	echo "=== re-running trusts the stamp instead of re-hashing ==="
	fetch "refetch good_plain" good_plain

	echo
	echo "=== a gzipped file is stored as served, not decompressed ==="
	fetch "fetch good_gz" good_gz
	exists "  cached file" "${CACHE}/good_gz/reads.3.fastq.gz"
	exists "  no decompressed sibling" "${CACHE}/good_gz/reads.3.fastq"

	echo
	echo "=== a checksum mismatch is fatal and poisons nothing ==="
	fetch "fetch bad_md5" bad_md5
	exists "  cached file" "${CACHE}/bad_md5/lambda.3.gbk"
	exists "  partial file" "${CACHE}/bad_md5/lambda.3.gbk.part.$$"

	echo
	echo "=== an unreachable URL is fatal and poisons nothing ==="
	fetch "fetch missing" missing
	exists "  cached file" "${CACHE}/missing/nope.gbk"

	echo
	echo "=== tampering with a cached file invalidates its stamp ==="
	echo "extra" >> "${CACHE}/good_plain/lambda.3.gbk"
	fetch "refetch after tampering" good_plain

	echo
	echo "=== --verify re-hashes even when the stamp is good ==="
	fetch "verify good_plain" --verify good_plain

	echo
	echo "=== --lint rejects a malformed manifest ==="
	"${FETCH}" --manifest "${BAD_MANIFEST}" --lint > /dev/null 2>&1
	echo "lint bad manifest: exit=$?"
	"${FETCH}" --manifest "${MANIFEST}" --lint > /dev/null 2>&1
	echo "lint good manifest: exit=$?"

	echo
	echo "=== an unknown dataset name is rejected ==="
	fetch "fetch nonexistent" no_such_dataset
} > "${OUT}" 2>&1

# The harness compares ${OUT} against expected.txt; always succeed here so that a
# behavior change shows up as a readable diff rather than an opaque failure.
exit 0
