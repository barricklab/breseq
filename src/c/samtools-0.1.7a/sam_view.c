#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include "sam_header.h"
#include "sam.h"
#include "faidx.h"

static int g_min_mapQ = 0, g_flag_on = 0, g_flag_off = 0;
static char *g_library, *g_rg;
static int g_sol2sanger_tbl[128];

static void sol2sanger(bam1_t *b)
{
	int l;
	uint8_t *qual = bam1_qual(b);
	if (g_sol2sanger_tbl[30] == 0) {
		for (l = 0; l != 128; ++l) {
			g_sol2sanger_tbl[l] = (int)(10.0 * log(1.0 + pow(10.0, (l - 64 + 33) / 10.0)) / log(10.0) + .499);
			if (g_sol2sanger_tbl[l] >= 93) g_sol2sanger_tbl[l] = 93;
		}
	}
	for (l = 0; l < b->core.l_qseq; ++l) {
		int q = qual[l];
		if (q > 127) q = 127;
		qual[l] = g_sol2sanger_tbl[q];
	}
}

static inline int __g_skip_aln(const bam_header_t *h, const bam1_t *b)
{
	if (b->core.qual < g_min_mapQ || ((b->core.flag & g_flag_on) != g_flag_on) || (b->core.flag & g_flag_off))
		return 1;
	if (g_rg) {
		uint8_t *s = bam_aux_get(b, "RG");
		if (s && strcmp(g_rg, (char*)(s + 1)) == 0) return 0;
	}
	if (g_library) {
		const char *p = bam_get_library((bam_header_t*)h, b);
		return (p && strcmp(p, g_library) == 0)? 0 : 1;
	}
	return 0;
}

// callback function for bam_fetch()
static int view_func(const bam1_t *b, void *data)
{
	if (!__g_skip_aln(((samfile_t*)data)->header, b))
		samwrite((samfile_t*)data, b);
	return 0;
}

static int usage(int is_long_help);

int main_samview(int argc, char *argv[])
{
	int c, is_header = 0, is_header_only = 0, is_bamin = 1, ret = 0, is_uncompressed = 0, is_bamout = 0, slx2sngr = 0;
	int of_type = BAM_OFDEC, is_long_help = 0;
	samfile_t *in = 0, *out = 0;
	char in_mode[5], out_mode[5], *fn_out = 0, *fn_list = 0, *fn_ref = 0;

	/* parse command-line options */
	strcpy(in_mode, "r"); strcpy(out_mode, "w");
	while ((c = getopt(argc, argv, "Sbt:hHo:q:f:F:ul:r:xX?T:C")) >= 0) {
		switch (c) {
		case 'C': slx2sngr = 1; break;
		case 'S': is_bamin = 0; break;
		case 'b': is_bamout = 1; break;
		case 't': fn_list = strdup(optarg); is_bamin = 0; break;
		case 'h': is_header = 1; break;
		case 'H': is_header_only = 1; break;
		case 'o': fn_out = strdup(optarg); break;
		case 'f': g_flag_on = strtol(optarg, 0, 0); break;
		case 'F': g_flag_off = strtol(optarg, 0, 0); break;
		case 'q': g_min_mapQ = atoi(optarg); break;
		case 'u': is_uncompressed = 1; break;
		case 'l': g_library = strdup(optarg); break;
		case 'r': g_rg = strdup(optarg); break;
		case 'x': of_type = BAM_OFHEX; break;
		case 'X': of_type = BAM_OFSTR; break;
		case '?': is_long_help = 1; break;
		case 'T': fn_ref = strdup(optarg); is_bamin = 0; break;
		default: return usage(is_long_help);
		}
	}
	if (is_uncompressed) is_bamout = 1;
	if (is_header_only) is_header = 1;
	if (is_bamout) strcat(out_mode, "b");
	else {
		if (of_type == BAM_OFHEX) strcat(out_mode, "x");
		else if (of_type == BAM_OFSTR) strcat(out_mode, "X");
	}
	if (is_bamin) strcat(in_mode, "b");
	if (is_header) strcat(out_mode, "h");
	if (is_uncompressed) strcat(out_mode, "u");
	if (argc == optind) return usage(is_long_help);

	// generate the fn_list if necessary
	if (fn_list == 0 && fn_ref) fn_list = samfaipath(fn_ref);
	// open file handlers
	if ((in = samopen(argv[optind], in_mode, fn_list)) == 0) {
		fprintf(stderr, "[main_samview] fail to open file for reading.\n");
		goto view_end;
	}
	if (in->header == 0) {
		fprintf(stderr, "[main_samview] fail to read the header.\n");
		goto view_end;
	}
	if ((out = samopen(fn_out? fn_out : "-", out_mode, in->header)) == 0) {
		fprintf(stderr, "[main_samview] fail to open file for writing.\n");
		goto view_end;
	}
	if (is_header_only) goto view_end; // no need to print alignments

	if (argc == optind + 1) { // convert/print the entire file
		bam1_t *b = bam_init1();
		int r;
		while ((r = samread(in, b)) >= 0) { // read one alignment from `in'
			if (!__g_skip_aln(in->header, b)) {
				if (slx2sngr) sol2sanger(b);
				samwrite(out, b); // write the alignment to `out'
			}
		}
		if (r < -1) fprintf(stderr, "[main_samview] truncated file.\n");
		bam_destroy1(b);
	} else { // retrieve alignments in specified regions
		int i;
		bam_index_t *idx = 0;
		if (is_bamin) idx = bam_index_load(argv[optind]); // load BAM index
		if (idx == 0) { // index is unavailable
			fprintf(stderr, "[main_samview] random alignment retrieval only works for indexed BAM files.\n");
			ret = 1;
			goto view_end;
		}
		for (i = optind + 1; i < argc; ++i) {
			int tid, beg, end;
			bam_parse_region(in->header, argv[i], &tid, &beg, &end); // parse a region in the format like `chr2:100-200'
			if (tid < 0) { // reference name is not found
				fprintf(stderr, "[main_samview] fail to get the reference name. Continue anyway.\n");
				continue;
			}
			bam_fetch(in->x.bam, idx, tid, beg, end, out, view_func); // fetch alignments
		}
		bam_index_destroy(idx); // destroy the BAM index
	}

view_end:
	// close files, free and return
	free(fn_list); free(fn_ref); free(fn_out); free(g_library); free(g_rg);
	samclose(in);
	samclose(out);
	return ret;
}

static int usage(int is_long_help)
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   samtools view [options] <in.bam>|<in.sam> [region1 [...]]\n\n");
	fprintf(stderr, "Options: -b       output BAM\n");
	fprintf(stderr, "         -h       print header for the SAM output\n");
	fprintf(stderr, "         -H       print header only (no alignments)\n");
	fprintf(stderr, "         -S       input is SAM\n");
	fprintf(stderr, "         -u       uncompressed BAM output (force -b)\n");
	fprintf(stderr, "         -x       output FLAG in HEX (samtools-C specific)\n");
	fprintf(stderr, "         -X       output FLAG in string (samtools-C specific)\n");
	fprintf(stderr, "         -t FILE  list of reference names and lengths (force -S) [null]\n");
	fprintf(stderr, "         -T FILE  reference sequence file (force -S) [null]\n");
	fprintf(stderr, "         -o FILE  output file name [stdout]\n");
	fprintf(stderr, "         -f INT   required flag, 0 for unset [0]\n");
	fprintf(stderr, "         -F INT   filtering flag, 0 for unset [0]\n");
	fprintf(stderr, "         -q INT   minimum mapping quality [0]\n");
	fprintf(stderr, "         -l STR   only output reads in library STR [null]\n");
	fprintf(stderr, "         -r STR   only output reads in read group STR [null]\n");
	fprintf(stderr, "         -?       longer help\n");
	fprintf(stderr, "\n");
	if (is_long_help)
		fprintf(stderr, "Notes:\n\
\n\
  1. By default, this command assumes the file on the command line is in\n\
     the BAM format and it prints the alignments in SAM. If `-t' is\n\
     applied, the input file is assumed to be in the SAM format. The\n\
     file supplied with `-t' is SPACE/TAB delimited with the first two\n\
     fields of each line consisting of the reference name and the\n\
     corresponding sequence length. The `.fai' file generated by `faidx'\n\
     can be used here. This file may be empty if reads are unaligned.\n\
\n\
  2. SAM->BAM conversion: `samtools view -bT ref.fa in.sam.gz'.\n\
\n\
  3. BAM->SAM conversion: `samtools view in.bam'.\n\
\n\
  4. A region should be presented in one of the following formats:\n\
     `chr1', `chr2:1,000' and `chr3:1000-2,000'. When a region is\n\
     specified, the input alignment file must be an indexed BAM file.\n\
\n\
  5. Option `-u' is preferred over `-b' when the output is piped to\n\
     another samtools command.\n\
\n\
  6. In a string FLAG, each character represents one bit with\n\
     p=0x1 (paired), P=0x2 (properly paired), u=0x4 (unmapped),\n\
     U=0x8 (mate unmapped), r=0x10 (reverse), R=0x20 (mate reverse)\n\
     1=0x40 (first), 2=0x80 (second), s=0x100 (not primary), \n\
     f=0x200 (failure) and d=0x400 (duplicate). Note that `-x' and\n\
     `-X' are samtools-C specific. Picard and older samtools do not\n\
     support HEX or string flags.\n\
\n");
	return 1;
}

int main_import(int argc, char *argv[])
{
	int argc2, ret;
	char **argv2;
	if (argc != 4) {
		fprintf(stderr, "Usage: bamtk import <in.ref_list> <in.sam> <out.bam>\n");
		return 1;
	}
	argc2 = 6;
	argv2 = calloc(6, sizeof(char*));
	argv2[0] = "import", argv2[1] = "-o", argv2[2] = argv[3], argv2[3] = "-bt", argv2[4] = argv[1], argv2[5] = argv[2];
	ret = main_samview(argc2, argv2);
	free(argv2);
	return ret;
}
