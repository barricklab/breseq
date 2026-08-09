#!/bin/bash
#
# Download, md5-verify and cache the large real-world data sets used by breseq's
# 'long' tests (see tests/test_data_manifest.tsv).
#
# Only 'long'-named tests use this -- 'make test' never touches the network.
# A test declares what it needs with a TEST_DATA=<dataset>[,<dataset>...] line in
# its testcmd.sh. Two independent callers act on that declaration:
#
#   * tests/Snakefile makes each dataset a node in the job graph, so a data set
#     shared by several long tests is downloaded once, and distinct data sets
#     download in parallel ('make test-long').
#   * tests/common.sh calls this script from do_test, so running a single test
#     directly ('./tests/test.sh test <name>') -- which bypasses Snakemake
#     entirely -- fetches the same data. When Snakemake already did the work this
#     is a couple of stat() calls, not a re-hash of multi-gigabyte files.
#
# Design notes that are easy to get wrong if this is ever rewritten:
#
#   * Nothing is written to the cache until its md5 has been verified. The
#     download goes to a temp file IN THE DESTINATION DIRECTORY (so the rename is
#     atomic -- same filesystem) and is renamed into place only after the hash
#     matches. A failed, interrupted or corrupt download therefore cannot leave
#     a plausible-looking file behind for the next run to trust.
#   * There is deliberately no 'curl -C -' resume. Resuming into a file whose
#     provenance is unknown is how you end up with 90% of one version of a data
#     set and 10% of another, which would then fail the md5 check with no
#     explanation. An interrupted download costs the whole file again.
#   * A checksum mismatch is a HARD failure, never a warning. It means the remote
#     data changed, so any expected.gd built from it is no longer meaningful.
#
# Usage:
#   tests/fetch_test_data.sh [options] <dataset> [<dataset> ...]
#   tests/fetch_test_data.sh [options] --all
#
#   --all              act on every dataset in the manifest
#   --verify           ignore the cached-and-verified markers; re-hash everything
#   --lint             validate the manifest and exit (no network, no writes)
#   --list             print dataset names, one per line, and exit
#   --print-dir        print the resolved cache directory and exit
#   --clean            delete the named datasets' directories from the cache
#   --manifest <path>  use an alternate manifest (used by the offline self-test)
#   --dir <path>       override the cache directory
#
# Exit codes:
#   0  every requested file is present and verified
#   1  usage error, unknown dataset, malformed manifest, or missing tool
#   2  download failed (offline, DNS failure, 404, timeout, stalled transfer)
#   3  CHECKSUM MISMATCH
#   4  could not acquire the per-file lock before timing out
#
# Portability: must work with bash 3.2, which is what /bin/bash is on macOS (and
# on GitHub's macos runners). In particular that means no mapfile/readarray and no
# associative arrays, and note that under 'set -u' bash 3.2 treats "${arr[@]}" on
# an EMPTY array as an unbound variable -- so every such expansion below is
# guarded by a preceding emptiness check. ("${!arr[@]}" is safe.)

set -u

EXIT_OK=0
EXIT_USAGE=1
EXIT_DOWNLOAD=2
EXIT_CHECKSUM=3
EXIT_LOCK=4

# Resolved from this script's own location rather than the working directory, so
# the cache is the same whether we are invoked by the Snakefile (which runs from
# the top-level source directory) or by a testcmd.sh.
SELFDIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

MANIFEST="${SELFDIR}/test_data_manifest.tsv"
OPT_DIR=""
DO_ALL=0
DO_VERIFY=0
DO_LINT=0
DO_LIST=0
DO_PRINT_DIR=0
DO_CLEAN=0
DATASETS=()

# How long to wait for another process that is downloading the same file.
LOCK_TIMEOUT_SECONDS=$((3 * 60 * 60))
LOCK_POLL_SECONDS=5

banner_start() { echo "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" >&2; }
banner_end()   { echo "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" >&2; }

usage() {
	sed -n '/^# Usage:/,/^#   4  could not/p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//' >&2
	exit ${EXIT_USAGE}
}

##############################################################################
# Portable helpers.
#
# These fall back through several implementations on purpose: GNU coreutils and
# BSD/macOS differ here, and the conda environment supplies neither. A GitHub
# macos runner, for example, has 'md5' but not 'md5sum'.
##############################################################################

md5_of() { # $1 = path; prints 32 lowercase hex characters
	local out
	if   command -v md5sum  > /dev/null 2>&1; then out=$(md5sum "$1" | cut -d' ' -f1)
	elif command -v md5     > /dev/null 2>&1; then out=$(md5 -q "$1")
	elif command -v openssl > /dev/null 2>&1; then out=$(openssl md5 -r "$1" | cut -d' ' -f1)
	elif command -v python3 > /dev/null 2>&1; then
		out=$(python3 -c 'import hashlib,sys
h = hashlib.md5()
with open(sys.argv[1], "rb") as f:
    for block in iter(lambda: f.read(1 << 20), b""):
        h.update(block)
print(h.hexdigest())' "$1")
	else
		banner_start
		echo "No md5 tool found. Need one of: md5sum, md5, openssl, python3." >&2
		banner_end
		exit ${EXIT_USAGE}
	fi
	echo "${out}" | tr 'A-F' 'a-f'
}

size_and_mtime() { # $1 = path; prints "<size_bytes> <mtime_epoch>"
	stat -c '%s %Y' "$1" 2> /dev/null || stat -f '%z %m' "$1"
}

human_size() { # $1 = bytes
	awk -v b="$1" 'BEGIN {
		if (b >= 1073741824) printf "%.1f GB", b/1073741824;
		else if (b >= 1048576) printf "%.1f MB", b/1048576;
		else if (b >= 1024) printf "%.1f KB", b/1024;
		else printf "%d bytes", b;
	}'
}

##############################################################################
# Argument parsing
##############################################################################

while [[ $# -gt 0 ]]; do
	case "$1" in
		--all)        DO_ALL=1 ;;
		--verify)     DO_VERIFY=1 ;;
		--lint)       DO_LINT=1 ;;
		--list)       DO_LIST=1 ;;
		--print-dir)  DO_PRINT_DIR=1 ;;
		--clean)      DO_CLEAN=1 ;;
		--manifest)   shift; [[ $# -gt 0 ]] || usage; MANIFEST="$1" ;;
		--dir)        shift; [[ $# -gt 0 ]] || usage; OPT_DIR="$1" ;;
		-h|--help)    usage ;;
		-*)           echo "Unknown option: $1" >&2; usage ;;
		*)            DATASETS+=("$1") ;;
	esac
	shift
done

if [[ ! -f "${MANIFEST}" ]]; then
	banner_start
	echo "Manifest not found: ${MANIFEST}" >&2
	banner_end
	exit ${EXIT_USAGE}
fi

# BRESEQ_TEST_DATA_DIR lets several worktrees (or a CI job with a restored cache)
# share one copy of the data. The default is inside tests/data/ but deliberately
# NOT inside any test directory: common.sh's do_clean does 'rm -Rf <testdir>/data',
# which would delete a cache placed there. Note that no clean target removes this
# directory -- see 'make clean-test-data'.
CACHE_DIR="${OPT_DIR:-${BRESEQ_TEST_DATA_DIR:-${SELFDIR}/data/downloads}}"

if [[ ${DO_PRINT_DIR} -eq 1 ]]; then
	echo "${CACHE_DIR}"
	exit ${EXIT_OK}
fi

##############################################################################
# Manifest reading
#
# Rows are read into parallel arrays so the rest of the script can iterate
# without re-parsing, and so line numbers survive for error messages.
##############################################################################

M_DATASET=()
M_FILE=()
M_MD5=()
M_URL=()
M_LINE=()

read_manifest() {
	local lineno=0 dataset file md5 url extra
	while IFS=$'\t' read -r dataset file md5 url extra || [[ -n "${dataset:-}" ]]; do
		lineno=$((lineno + 1))
		# Skip comments and blank lines.
		[[ "${dataset}" =~ ^[[:space:]]*# ]] && continue
		[[ -z "${dataset// /}" ]] && continue

		if [[ -z "${file:-}" || -z "${md5:-}" || -z "${url:-}" || -n "${extra:-}" ]]; then
			banner_start
			echo "Malformed manifest row (expected exactly 4 tab-separated columns)" >&2
			echo "  ${MANIFEST} line ${lineno}" >&2
			echo "  dataset='${dataset}' file='${file:-}' md5='${md5:-}' url='${url:-}'" >&2
			banner_end
			exit ${EXIT_USAGE}
		fi

		M_DATASET+=("${dataset}")
		M_FILE+=("${file}")
		M_MD5+=("${md5}")
		M_URL+=("${url}")
		M_LINE+=("${lineno}")
	done < "${MANIFEST}"
}

read_manifest

if [[ ${#M_DATASET[@]} -eq 0 && ${DO_LINT} -eq 0 && ${DO_LIST} -eq 0 ]]; then
	banner_start
	echo "Manifest contains no data rows: ${MANIFEST}" >&2
	banner_end
	exit ${EXIT_USAGE}
fi

all_dataset_names() {
	local i
	for i in "${!M_DATASET[@]}"; do echo "${M_DATASET[$i]}"; done | sort -u
}

if [[ ${DO_LIST} -eq 1 ]]; then
	all_dataset_names
	exit ${EXIT_OK}
fi

##############################################################################
# --lint : offline validation, safe to run on every pull request
##############################################################################

if [[ ${DO_LINT} -eq 1 ]]; then
	errors=0
	seen_keys=""
	for i in "${!M_DATASET[@]}"; do
		d="${M_DATASET[$i]}"; f="${M_FILE[$i]}"; m="${M_MD5[$i]}"; u="${M_URL[$i]}"; l="${M_LINE[$i]}"

		if [[ ! "${d}" =~ ^[A-Za-z0-9_]+$ ]]; then
			echo "${MANIFEST}:${l}: dataset '${d}' must match [A-Za-z0-9_]+ (tests/Snakefile parses TEST_DATA as [\\w,]+)" >&2
			errors=$((errors + 1))
		fi
		if [[ "${f}" == */* || "${f}" == *..* ]]; then
			echo "${MANIFEST}:${l}: file '${f}' must be a bare file name" >&2
			errors=$((errors + 1))
		fi
		if [[ ! "${m}" =~ ^[0-9a-f]{32}$ ]]; then
			echo "${MANIFEST}:${l}: md5 '${m}' must be 32 lowercase hex characters" >&2
			errors=$((errors + 1))
		fi
		if [[ "${u}" != https://* && "${u}" != file://* ]]; then
			echo "${MANIFEST}:${l}: url must start with https:// (file:// allowed only for test fixtures)" >&2
			errors=$((errors + 1))
		fi
		key="${d}/${f}"
		case " ${seen_keys} " in
			*" ${key} "*)
				echo "${MANIFEST}:${l}: duplicate entry for ${key}" >&2
				errors=$((errors + 1)) ;;
			*) seen_keys="${seen_keys} ${key}" ;;
		esac
	done

	# Every dataset referenced by a test must actually exist in the manifest.
	# This is the check that catches a typo in a TEST_DATA= line at pull-request
	# time rather than an hour into a long-test run.
	#
	# Only meaningful for the real manifest: an alternate one (--manifest, used by
	# the offline self-test) describes its own fixtures and knows nothing about
	# what the committed tests declare.
	if [[ "${MANIFEST}" == "${SELFDIR}/test_data_manifest.tsv" ]]; then
		known=$(all_dataset_names)
		for testcmd in "${SELFDIR}"/*/testcmd.sh; do
			[[ -e "${testcmd}" ]] || continue
			# Same anchored pattern tests/Snakefile parses TEST_DATA with.
			decl=$(sed -n 's/^[[:space:]]*TEST_DATA=\([A-Za-z0-9_,]*\)[[:space:]]*$/\1/p' "${testcmd}")
			[[ -z "${decl}" ]] && continue
			for ds in ${decl//,/ }; do
				if ! echo "${known}" | grep -qx "${ds}"; then
					echo "${testcmd}: TEST_DATA names '${ds}', which is not in ${MANIFEST}" >&2
					errors=$((errors + 1))
				fi
			done
		done
	fi

	if [[ ${errors} -gt 0 ]]; then
		banner_start
		echo "Manifest validation FAILED with ${errors} problem(s)." >&2
		banner_end
		exit ${EXIT_USAGE}
	fi
	echo "[lint]     ${MANIFEST}: ${#M_DATASET[@]} file(s) in $(all_dataset_names | wc -l | tr -d ' ') dataset(s) -- OK"
	exit ${EXIT_OK}
fi

##############################################################################
# Work out which datasets were requested
##############################################################################

if [[ ${DO_ALL} -eq 1 ]]; then
	while IFS= read -r d; do DATASETS+=("${d}"); done < <(all_dataset_names)
fi

if [[ ${#DATASETS[@]} -eq 0 ]]; then
	echo "No datasets requested." >&2
	usage
fi

for ds in "${DATASETS[@]}"; do
	if ! all_dataset_names | grep -qx "${ds}"; then
		banner_start
		echo "Unknown dataset: ${ds}" >&2
		echo "Known datasets in ${MANIFEST}:" >&2
		all_dataset_names | sed 's/^/  /' >&2
		banner_end
		exit ${EXIT_USAGE}
	fi
done

##############################################################################
# --clean
##############################################################################

if [[ ${DO_CLEAN} -eq 1 ]]; then
	for ds in "${DATASETS[@]}"; do
		target="${CACHE_DIR}/${ds}"
		if [[ -d "${target}" ]]; then
			rm -rf "${target}"
			echo "[removed]  ${ds}"
		else
			echo "[absent]   ${ds}"
		fi
	done
	# Only the per-dataset directories are removed, never CACHE_DIR itself:
	# BRESEQ_TEST_DATA_DIR may point at a shared or user-owned location.
	exit ${EXIT_OK}
fi

##############################################################################
# Fetching
##############################################################################

if ! command -v curl > /dev/null 2>&1; then
	banner_start
	echo "curl is required to download test data but was not found on PATH." >&2
	banner_end
	exit ${EXIT_USAGE}
fi

# Cleaned up by the trap below, so an interrupted run leaves no partial file and
# no orphaned lock.
CUR_TMP=""
CUR_LOCK=""
cleanup() {
	[[ -n "${CUR_TMP}"  && -e "${CUR_TMP}"  ]] && rm -f "${CUR_TMP}"
	[[ -n "${CUR_LOCK}" && -d "${CUR_LOCK}" ]] && rmdir "${CUR_LOCK}" 2> /dev/null
	return 0
}
trap cleanup EXIT INT TERM

# Is the cached file present and known-good? Reads the sidecar stamp written by
# mark_verified() rather than re-hashing, which is what keeps a fully-cached run
# down to a few stat() calls instead of re-reading gigabytes.
is_cached() { # $1 = dest path, $2 = expected md5
	local dest="$1" want="$2" stamp="$1.md5ok"
	# --verify ignores the stamp; see verify_existing() for what it does instead.
	[[ ${DO_VERIFY} -eq 1 ]] && return 1
	[[ -f "${dest}" && -f "${stamp}" ]] || return 1

	local s_md5 s_size s_mtime cur
	read -r s_md5 s_size s_mtime < "${stamp}" || return 1
	[[ "${s_md5}" == "${want}" ]] || return 1

	# Size and mtime guard against the file being edited or truncated after it
	# was verified; anything unexpected forces a full re-download.
	cur=$(size_and_mtime "${dest}") || return 1
	[[ "${cur}" == "${s_size} ${s_mtime}" ]] || return 1
	return 0
}

mark_verified() { # $1 = dest path, $2 = md5
	local dest="$1" md5="$2" stamp="$1.md5ok" tmp="$1.md5ok.part.$$"
	echo "${md5} $(size_and_mtime "${dest}")" > "${tmp}"
	mv -f "${tmp}" "${stamp}"
}

acquire_lock() { # $1 = dest path
	local lock="$1.lock" waited=0
	while ! mkdir "${lock}" 2> /dev/null; do
		if [[ ${waited} -ge ${LOCK_TIMEOUT_SECONDS} ]]; then
			banner_start
			echo "Timed out waiting for another process to finish downloading" >&2
			echo "  $1" >&2
			echo "If no other test run is active, that process died holding the lock." >&2
			echo "Remove it by hand and retry:" >&2
			echo "  rmdir '${lock}'" >&2
			banner_end
			return 1
		fi
		if [[ $((waited % 60)) -eq 0 ]]; then
			echo "[waiting]  another process is downloading $(basename "$1") ..."
		fi
		sleep ${LOCK_POLL_SECONDS}
		waited=$((waited + LOCK_POLL_SECONDS))
	done
	CUR_LOCK="${lock}"
	return 0
}

release_lock() {
	[[ -n "${CUR_LOCK}" && -d "${CUR_LOCK}" ]] && rmdir "${CUR_LOCK}" 2> /dev/null
	CUR_LOCK=""
	return 0
}

fetch_one() { # $1 = index into the manifest arrays
	local i="$1"
	local dataset="${M_DATASET[$i]}" file="${M_FILE[$i]}"
	local md5="${M_MD5[$i]}" url="${M_URL[$i]}" line="${M_LINE[$i]}"
	local dir="${CACHE_DIR}/${dataset}"
	local dest="${dir}/${file}"

	if is_cached "${dest}" "${md5}"; then
		echo "[cached]   ${dataset}/${file}"
		return ${EXIT_OK}
	fi

	# Under --verify, re-hash what is already on disk rather than re-downloading
	# it: the point of 'make verify-test-data' is to detect local corruption, and
	# re-fetching hundreds of megabytes to answer that question would defeat it.
	# A file that fails this check falls through to a fresh download, which
	# repairs corruption and still hard-fails if the remote data itself drifted.
	if [[ ${DO_VERIFY} -eq 1 && -f "${dest}" ]]; then
		local on_disk
		on_disk=$(md5_of "${dest}")
		if [[ "${on_disk}" == "${md5}" ]]; then
			mark_verified "${dest}" "${md5}"
			echo "[verified] ${dataset}/${file}"
			return ${EXIT_OK}
		fi
		echo "[corrupt]  ${dataset}/${file} (md5 ${on_disk}, expected ${md5}) -- re-downloading" >&2
	fi

	mkdir -p "${dir}"
	acquire_lock "${dest}" || return ${EXIT_LOCK}

	# Re-check: whoever held the lock may have just finished this exact file.
	if is_cached "${dest}" "${md5}"; then
		echo "[cached]   ${dataset}/${file}"
		release_lock
		return ${EXIT_OK}
	fi

	echo "[fetching] ${dataset}/${file}  <-  ${url}"

	CUR_TMP="${dest}.part.$$"
	local curl_status=0
	curl --location --fail --show-error --silent \
	     --retry 3 --retry-delay 5 --retry-connrefused \
	     --connect-timeout 30 --speed-limit 1024 --speed-time 120 \
	     --output "${CUR_TMP}" "${url}" || curl_status=$?

	if [[ ${curl_status} -ne 0 ]]; then
		rm -f "${CUR_TMP}"; CUR_TMP=""
		release_lock
		banner_start
		echo "FAILED to download test data" >&2
		echo "  dataset : ${dataset}" >&2
		echo "  file    : ${file}" >&2
		echo "  url     : ${url}" >&2
		echo "  curl exited ${curl_status}" >&2
		echo "" >&2
		echo "Long tests require network access. If you are offline, run 'make test'," >&2
		echo "which never downloads anything, instead of 'make test-long'." >&2
		echo "Nothing was written to the cache." >&2
		banner_end
		return ${EXIT_DOWNLOAD}
	fi

	local actual
	actual=$(md5_of "${CUR_TMP}")
	if [[ "${actual}" != "${md5}" ]]; then
		rm -f "${CUR_TMP}"; CUR_TMP=""
		release_lock
		banner_start
		echo "CHECKSUM MISMATCH -- test data has changed or the download is corrupt" >&2
		echo "  dataset  : ${dataset}" >&2
		echo "  file     : ${file}" >&2
		echo "  url      : ${url}" >&2
		echo "  expected : ${md5}   (${MANIFEST} line ${line})" >&2
		echo "  actual   : ${actual}" >&2
		echo "" >&2
		echo "The downloaded file has been DELETED; the cache was not modified." >&2
		echo "If the remote data legitimately changed, update the md5 in the manifest" >&2
		echo "AND rebuild the expected.gd of every test that uses this dataset --" >&2
		echo "reviewing the resulting diff. Do not silently accept the new checksum." >&2
		banner_end
		return ${EXIT_CHECKSUM}
	fi

	mv -f "${CUR_TMP}" "${dest}"
	CUR_TMP=""
	mark_verified "${dest}" "${md5}"
	release_lock

	local sz
	sz=$(size_and_mtime "${dest}" | cut -d' ' -f1)
	echo "[fetched]  ${dataset}/${file}  ($(human_size "${sz}"), md5 ok)"
	return ${EXIT_OK}
}

status=${EXIT_OK}
for ds in "${DATASETS[@]}"; do
	for i in "${!M_DATASET[@]}"; do
		[[ "${M_DATASET[$i]}" == "${ds}" ]] || continue
		fetch_one "$i" || { status=$?; break 2; }
	done
done

exit ${status}
