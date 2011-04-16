#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <zlib.h>
#include <errno.h>
#include "bcf.h"
#include "prob1.h"
#include "kstring.h"
#include "time.h"

#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 16384)

#define VC_NO_GENO 2
#define VC_BCFOUT  4
#define VC_CALL    8
#define VC_VARONLY 16
#define VC_VCFIN   32
#define VC_UNCOMP  64
#define VC_KEEPALT 256
#define VC_ACGT_ONLY 512
#define VC_QCALL   1024
#define VC_CALL_GT 2048
#define VC_ADJLD   4096
#define VC_NO_INDEL 8192
#define VC_ANNO_MAX 16384
#define VC_FIX_PL   32768

typedef struct {
	int flag, prior_type, n1, n_sub, *sublist, n_perm;
	char *prior_file, **subsam, *fn_dict;
	uint8_t *ploidy;
	double theta, pref, indel_frac, min_perm_p, min_smpl_frac;
	void *bed;
} viewconf_t;

void *bed_read(const char *fn);
void bed_destroy(void *_h);
int bed_overlap(const void *_h, const char *chr, int beg, int end);

static double test_hwe(const double g[3])
{
	extern double kf_gammaq(double p, double x);
	double fexp, chi2, f[3], n;
	int i;
	n = g[0] + g[1] + g[2];
	fexp = (2. * g[2] + g[1]) / (2. * n);
	if (fexp > 1. - 1e-10) fexp = 1. - 1e-10;
	if (fexp < 1e-10) fexp = 1e-10;
	f[0] = n * (1. - fexp) * (1. - fexp);
	f[1] = n * 2. * fexp * (1. - fexp);
	f[2] = n * fexp * fexp;
	for (i = 0, chi2 = 0.; i < 3; ++i)
		chi2 += (g[i] - f[i]) * (g[i] - f[i]) / f[i];
	return kf_gammaq(.5, chi2 / 2.);
}

typedef struct {
	double p[4];
	int mq, depth, is_tested, d[4];
} anno16_t;

static double ttest(int n1, int n2, int a[4])
{
	extern double kf_betai(double a, double b, double x);
	double t, v, u1, u2;
	if (n1 == 0 || n2 == 0 || n1 + n2 < 3) return 1.0;
	u1 = (double)a[0] / n1; u2 = (double)a[2] / n2;
	if (u1 <= u2) return 1.;
	t = (u1 - u2) / sqrt(((a[1] - n1 * u1 * u1) + (a[3] - n2 * u2 * u2)) / (n1 + n2 - 2) * (1./n1 + 1./n2));
	v = n1 + n2 - 2;
//	printf("%d,%d,%d,%d,%lf,%lf,%lf\n", a[0], a[1], a[2], a[3], t, u1, u2);
	return t < 0.? 1. : .5 * kf_betai(.5*v, .5, v/(v+t*t));
}

static int test16_core(int anno[16], anno16_t *a)
{
	extern double kt_fisher_exact(int n11, int n12, int n21, int n22, double *_left, double *_right, double *two);
	double left, right;
	int i;
	a->p[0] = a->p[1] = a->p[2] = a->p[3] = 1.;
	memcpy(a->d, anno, 4 * sizeof(int));
	a->depth = anno[0] + anno[1] + anno[2] + anno[3];
	a->is_tested = (anno[0] + anno[1] > 0 && anno[2] + anno[3] > 0);
	if (a->depth == 0) return -1;
	a->mq = (int)(sqrt((anno[9] + anno[11]) / a->depth) + .499);
	kt_fisher_exact(anno[0], anno[1], anno[2], anno[3], &left, &right, &a->p[0]);
	for (i = 1; i < 4; ++i)
		a->p[i] = ttest(anno[0] + anno[1], anno[2] + anno[3], anno+4*i);
	return 0;
}

static int test16(bcf1_t *b, anno16_t *a)
{
	char *p;
	int i, anno[16];
	a->p[0] = a->p[1] = a->p[2] = a->p[3] = 1.;
	a->d[0] = a->d[1] = a->d[2] = a->d[3] = 0.;
	a->mq = a->depth = a->is_tested = 0;
	if ((p = strstr(b->info, "I16=")) == 0) return -1;
	p += 4;
	for (i = 0; i < 16; ++i) {
		errno = 0; anno[i] = strtol(p, &p, 10);
		if (anno[i] == 0 && (errno == EINVAL || errno == ERANGE)) return -2;
		++p;
	}
	return test16_core(anno, a);
}

static void rm_info(bcf1_t *b, const char *key)
{
	char *p, *q;
	if ((p = strstr(b->info, key)) == 0) return;
	for (q = p; *q && *q != ';'; ++q);
	if (p > b->info && *(p-1) == ';') --p;
	memmove(p, q, b->l_str - (q - b->str));
	b->l_str -= q - p;
	bcf_sync(b);
}

static int update_bcf1(int n_smpl, bcf1_t *b, const bcf_p1aux_t *pa, const bcf_p1rst_t *pr, double pref, int flag)
{
	kstring_t s;
	int has_I16, is_var = (pr->p_ref < pref);
	double fq, r = is_var? pr->p_ref : pr->p_var;
	anno16_t a;

	has_I16 = test16(b, &a) >= 0? 1 : 0;
	rm_info(b, "I16=");

	memset(&s, 0, sizeof(kstring_t));
	kputc('\0', &s); kputs(b->ref, &s); kputc('\0', &s);
	kputs(b->alt, &s); kputc('\0', &s); kputc('\0', &s);
	kputs(b->info, &s);
	if (b->info[0]) kputc(';', &s);
//	ksprintf(&s, "AF1=%.4lg;AFE=%.4lg;CI95=%.4lg,%.4lg", 1.-pr->f_em, 1.-pr->f_exp, pr->cil, pr->cih);
	ksprintf(&s, "AF1=%.4g;CI95=%.4g,%.4g;G3=%.4g,%.4g,%.4g", 1.-pr->f_em, pr->cil, pr->cih, pr->g[2], pr->g[1], pr->g[0]);
	if (n_smpl > 5) {
		double hwe = test_hwe(pr->g);
		if (hwe < 0.1) ksprintf(&s, ";HWE=%.4g", hwe);
	}
	if (has_I16) ksprintf(&s, ";DP4=%d,%d,%d,%d;MQ=%d", a.d[0], a.d[1], a.d[2], a.d[3], a.mq);
	fq = pr->p_ref_folded < 0.5? -4.343 * log(pr->p_ref_folded) : 4.343 * log(pr->p_var_folded);
	if (fq < -999) fq = -999;
	if (fq > 999) fq = 999;
	ksprintf(&s, ";FQ=%.3g", fq);
	if (pr->cmp[0] >= 0.) { // two sample groups
		int i, q[3], pq;
		for (i = 1; i < 3; ++i) {
			double x = pr->cmp[i] + pr->cmp[0]/2.;
			q[i] = x == 0? 255 : (int)(-4.343 * log(x) + .499);
			if (q[i] > 255) q[i] = 255;
		}
		pq = (int)(-4.343 * log(pr->p_chi2) + .499);
		if (pr->perm_rank >= 0) ksprintf(&s, ";PR=%d", pr->perm_rank);
		ksprintf(&s, ";QCHI2=%d;PCHI2=%.3g;PC2=%d,%d", pq, q[1], q[2], pr->p_chi2);
		ksprintf(&s, ";AF2=%.4g,%.4g", 1.-pr->f_em2[0], 1.-pr->f_em2[1]);
//		ksprintf(&s, ",%g,%g,%g", pr->cmp[0], pr->cmp[1], pr->cmp[2]);
	}
	if (has_I16 && a.is_tested) ksprintf(&s, ";PV4=%.2g,%.2g,%.2g,%.2g", a.p[0], a.p[1], a.p[2], a.p[3]);
	kputc('\0', &s);
	kputs(b->fmt, &s); kputc('\0', &s);
	free(b->str);
	b->m_str = s.m; b->l_str = s.l; b->str = s.s;
	b->qual = r < 1e-100? 999 : -4.343 * log(r);
	if (b->qual > 999) b->qual = 999;
	bcf_sync(b);
	if (!is_var) bcf_shrink_alt(b, 1);
	else if (!(flag&VC_KEEPALT))
		bcf_shrink_alt(b, pr->rank0 < 2? 2 : pr->rank0+1);
	if (is_var && (flag&VC_CALL_GT)) { // call individual genotype
		int i, x, old_n_gi = b->n_gi;
		s.m = b->m_str; s.l = b->l_str - 1; s.s = b->str;
		kputs(":GT:GQ", &s); kputc('\0', &s);
		b->m_str = s.m; b->l_str = s.l; b->str = s.s;
		bcf_sync(b);
		for (i = 0; i < b->n_smpl; ++i) {
			x = bcf_p1_call_gt(pa, pr->f_em, i);
			((uint8_t*)b->gi[old_n_gi].data)[i] = (x&3) == 0? 1<<3|1 : (x&3) == 1? 1 : 0;
			((uint8_t*)b->gi[old_n_gi+1].data)[i] = x>>2;
		}
	}
	return is_var;
}

static char **read_samples(const char *fn, int *_n)
{
	gzFile fp;
	kstream_t *ks;
	kstring_t s;
	int dret, n = 0, max = 0;
	char **sam = 0;
	*_n = 0;
	s.l = s.m = 0; s.s = 0;
	fp = gzopen(fn, "r");
	if (fp == 0) return 0; // fail to open file
	ks = ks_init(fp);
	while (ks_getuntil(ks, 0, &s, &dret) >= 0) {
		int l;
		if (max == n) {
			max = max? max<<1 : 4;
			sam = realloc(sam, sizeof(void*)*max);
		}
		l = s.l;
		sam[n] = malloc(s.l + 2);
		strcpy(sam[n], s.s);
		sam[n][l+1] = 2; // by default, diploid
		if (dret != '\n') {
			if (ks_getuntil(ks, 0, &s, &dret) >= 0) { // read ploidy, 1 or 2
				int x = (int)s.s[0] - '0';
				if (x == 1 || x == 2) sam[n][l+1] = x;
				else fprintf(stderr, "(%s) ploidy can only be 1 or 2; assume diploid\n", __func__);
			}
			if (dret != '\n') ks_getuntil(ks, '\n', &s, &dret);
		}
		++n;
	}
	ks_destroy(ks);
	gzclose(fp);
	free(s.s);
	*_n = n;
	return sam;
}

static void write_header(bcf_hdr_t *h)
{
	kstring_t str;
	str.l = h->l_txt? h->l_txt - 1 : 0;
	str.m = str.l + 1; str.s = h->txt;
	if (!strstr(str.s, "##INFO=<ID=DP,"))
		kputs("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw read depth\">\n", &str);
	if (!strstr(str.s, "##INFO=<ID=DP4,"))
		kputs("##INFO=<ID=DP4,Number=4,Type=Integer,Description=\"# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases\">\n", &str);
	if (!strstr(str.s, "##INFO=<ID=MQ,"))
		kputs("##INFO=<ID=MQ,Number=1,Type=Integer,Description=\"Root-mean-square mapping quality of covering reads\">\n", &str);
	if (!strstr(str.s, "##INFO=<ID=FQ,"))
		kputs("##INFO=<ID=FQ,Number=1,Type=Float,Description=\"Phred probability of all samples being the same\">\n", &str);
	if (!strstr(str.s, "##INFO=<ID=AF1,"))
		kputs("##INFO=<ID=AF1,Number=1,Type=Float,Description=\"Max-likelihood estimate of the site allele frequency of the first ALT allele\">\n", &str);
	if (!strstr(str.s, "##INFO=<ID=G3,"))
		kputs("##INFO=<ID=G3,Number=3,Type=Float,Description=\"ML estimate of genotype frequencies\">\n", &str);
	if (!strstr(str.s, "##INFO=<ID=HWE,"))
		kputs("##INFO=<ID=HWE,Number=1,Type=Float,Description=\"Chi^2 based HWE test P-value based on G3\">\n", &str);
	if (!strstr(str.s, "##INFO=<ID=CI95,"))
		kputs("##INFO=<ID=CI95,Number=2,Type=Float,Description=\"Equal-tail Bayesian credible interval of the site allele frequency at the 95% level\">\n", &str);
	if (!strstr(str.s, "##INFO=<ID=PV4,"))
		kputs("##INFO=<ID=PV4,Number=4,Type=Float,Description=\"P-values for strand bias, baseQ bias, mapQ bias and tail distance bias\">\n", &str);
    if (!strstr(str.s, "##INFO=<ID=INDEL,"))
        kputs("##INFO=<ID=INDEL,Number=0,Type=Flag,Description=\"Indicates that the variant is an INDEL.\">\n", &str);
    if (!strstr(str.s, "##INFO=<ID=PC2,"))
        kputs("##INFO=<ID=PC2,Number=2,Type=Integer,Description=\"Phred probability of the nonRef allele frequency in group1 samples being larger (,smaller) than in group2.\">\n", &str);
    if (!strstr(str.s, "##INFO=<ID=PCHI2,"))
        kputs("##INFO=<ID=PCHI2,Number=1,Type=Float,Description=\"Posterior weighted chi^2 P-value for testing the association between group1 and group2 samples.\">\n", &str);
    if (!strstr(str.s, "##INFO=<ID=QCHI2,"))
        kputs("##INFO=<ID=QCHI2,Number=1,Type=Integer,Description=\"Phred scaled PCHI2.\">\n", &str);
    if (!strstr(str.s, "##INFO=<ID=RP,"))
        kputs("##INFO=<ID=PR,Number=1,Type=Integer,Description=\"# permutations yielding a smaller PCHI2.\">\n", &str);
    if (!strstr(str.s, "##FORMAT=<ID=GT,"))
        kputs("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n", &str);
    if (!strstr(str.s, "##FORMAT=<ID=GQ,"))
        kputs("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n", &str);
    if (!strstr(str.s, "##FORMAT=<ID=GL,"))
        kputs("##FORMAT=<ID=GL,Number=3,Type=Float,Description=\"Likelihoods for RR,RA,AA genotypes (R=ref,A=alt)\">\n", &str);
	if (!strstr(str.s, "##FORMAT=<ID=DP,"))
		kputs("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"# high-quality bases\">\n", &str);
	if (!strstr(str.s, "##FORMAT=<ID=SP,"))
		kputs("##FORMAT=<ID=SP,Number=1,Type=Integer,Description=\"Phred-scaled strand bias P-value\">\n", &str);
	if (!strstr(str.s, "##FORMAT=<ID=PL,"))
		kputs("##FORMAT=<ID=PL,Number=-1,Type=Integer,Description=\"List of Phred-scaled genotype likelihoods, number of values is (#ALT+1)*(#ALT+2)/2\">\n", &str);
	h->l_txt = str.l + 1; h->txt = str.s;
}

double bcf_ld_freq(const bcf1_t *b0, const bcf1_t *b1, double f[4]);

int bcfview(int argc, char *argv[])
{
	extern int bcf_2qcall(bcf_hdr_t *h, bcf1_t *b);
	extern void bcf_p1_indel_prior(bcf_p1aux_t *ma, double x);
	extern int bcf_fix_gt(bcf1_t *b);
	extern int bcf_anno_max(bcf1_t *b);
	extern int bcf_shuffle(bcf1_t *b, int seed);
	bcf_t *bp, *bout = 0;
	bcf1_t *b, *blast;
	int c, *seeds = 0;
	uint64_t n_processed = 0;
	viewconf_t vc;
	bcf_p1aux_t *p1 = 0;
	bcf_hdr_t *hin, *hout;
	int tid, begin, end;
	char moder[4], modew[4];

	tid = begin = end = -1;
	memset(&vc, 0, sizeof(viewconf_t));
	vc.prior_type = vc.n1 = -1; vc.theta = 1e-3; vc.pref = 0.5; vc.indel_frac = -1.; vc.n_perm = 0; vc.min_perm_p = 0.01; vc.min_smpl_frac = 0;
	while ((c = getopt(argc, argv, "FN1:l:cHAGvbSuP:t:p:QgLi:IMs:D:U:X:d:")) >= 0) {
		switch (c) {
		case '1': vc.n1 = atoi(optarg); break;
		case 'l': vc.bed = bed_read(optarg); break;
		case 'D': vc.fn_dict = strdup(optarg); break;
		case 'F': vc.flag |= VC_FIX_PL; break;
		case 'N': vc.flag |= VC_ACGT_ONLY; break;
		case 'G': vc.flag |= VC_NO_GENO; break;
		case 'A': vc.flag |= VC_KEEPALT; break;
		case 'b': vc.flag |= VC_BCFOUT; break;
		case 'S': vc.flag |= VC_VCFIN; break;
		case 'c': vc.flag |= VC_CALL; break;
		case 'v': vc.flag |= VC_VARONLY | VC_CALL; break;
		case 'u': vc.flag |= VC_UNCOMP | VC_BCFOUT; break;
		case 'g': vc.flag |= VC_CALL_GT | VC_CALL; break;
		case 'I': vc.flag |= VC_NO_INDEL; break;
		case 'M': vc.flag |= VC_ANNO_MAX; break;
		case 't': vc.theta = atof(optarg); break;
		case 'p': vc.pref = atof(optarg); break;
		case 'i': vc.indel_frac = atof(optarg); break;
		case 'Q': vc.flag |= VC_QCALL; break;
		case 'L': vc.flag |= VC_ADJLD; break;
		case 'U': vc.n_perm = atoi(optarg); break;
		case 'X': vc.min_perm_p = atof(optarg); break;
		case 'd': vc.min_smpl_frac = atof(optarg); break;
		case 's': vc.subsam = read_samples(optarg, &vc.n_sub);
			vc.ploidy = calloc(vc.n_sub + 1, 1);
			for (tid = 0; tid < vc.n_sub; ++tid) vc.ploidy[tid] = vc.subsam[tid][strlen(vc.subsam[tid]) + 1];
			tid = -1;
			break;
		case 'P':
			if (strcmp(optarg, "full") == 0) vc.prior_type = MC_PTYPE_FULL;
			else if (strcmp(optarg, "cond2") == 0) vc.prior_type = MC_PTYPE_COND2;
			else if (strcmp(optarg, "flat") == 0) vc.prior_type = MC_PTYPE_FLAT;
			else vc.prior_file = strdup(optarg);
			break;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage: bcftools view [options] <in.bcf> [reg]\n\n");
		fprintf(stderr, "Input/output options:\n\n");
		fprintf(stderr, "       -A        keep all possible alternate alleles at variant sites\n");
		fprintf(stderr, "       -b        output BCF instead of VCF\n");
		fprintf(stderr, "       -D FILE   sequence dictionary for VCF->BCF conversion [null]\n");
		fprintf(stderr, "       -F        PL generated by r921 or before (which generate old ordering)\n");
		fprintf(stderr, "       -G        suppress all individual genotype information\n");
		fprintf(stderr, "       -l FILE   list of sites (chr pos) or regions (BED) to output [all sites]\n");
		fprintf(stderr, "       -L        calculate LD for adjacent sites\n");
		fprintf(stderr, "       -N        skip sites where REF is not A/C/G/T\n");
		fprintf(stderr, "       -Q        output the QCALL likelihood format\n");
		fprintf(stderr, "       -s FILE   list of samples to use [all samples]\n");
		fprintf(stderr, "       -S        input is VCF\n");
		fprintf(stderr, "       -u        uncompressed BCF output (force -b)\n");
		fprintf(stderr, "\nConsensus/variant calling options:\n\n");
		fprintf(stderr, "       -c        SNP calling\n");
		fprintf(stderr, "       -d FLOAT  skip loci where less than FLOAT fraction of samples covered [0]\n");
		fprintf(stderr, "       -g        call genotypes at variant sites (force -c)\n");
		fprintf(stderr, "       -i FLOAT  indel-to-substitution ratio [%.4g]\n", vc.indel_frac);
		fprintf(stderr, "       -I        skip indels\n");
		fprintf(stderr, "       -p FLOAT  variant if P(ref|D)<FLOAT [%.3g]\n", vc.pref);
		fprintf(stderr, "       -P STR    type of prior: full, cond2, flat [full]\n");
		fprintf(stderr, "       -t FLOAT  scaled substitution mutation rate [%.4g]\n", vc.theta);
		fprintf(stderr, "       -v        output potential variant sites only (force -c)\n");
		fprintf(stderr, "\nContrast calling and association test options:\n\n");
		fprintf(stderr, "       -1 INT    number of group-1 samples [0]\n");
		fprintf(stderr, "       -U INT    number of permutations for association testing (effective with -1) [0]\n");
		fprintf(stderr, "       -X FLOAT  only perform permutations for P(chi^2)<FLOAT [%g]\n", vc.min_perm_p);
		fprintf(stderr, "\n");
		return 1;
	}

	if ((vc.flag & VC_VCFIN) && (vc.flag & VC_BCFOUT) && vc.fn_dict == 0) {
		fprintf(stderr, "[%s] For VCF->BCF conversion please specify the sequence dictionary with -D\n", __func__);
		return 1;
	}
	if (vc.n1 <= 0) vc.n_perm = 0; // TODO: give a warning here!
	if (vc.n_perm > 0) {
		seeds = malloc(vc.n_perm * sizeof(int));
		srand48(time(0));
		for (c = 0; c < vc.n_perm; ++c) seeds[c] = lrand48();
	}
	b = calloc(1, sizeof(bcf1_t));
	blast = calloc(1, sizeof(bcf1_t));
	strcpy(moder, "r");
	if (!(vc.flag & VC_VCFIN)) strcat(moder, "b");
	strcpy(modew, "w");
	if (vc.flag & VC_BCFOUT) strcat(modew, "b");
	if (vc.flag & VC_UNCOMP) strcat(modew, "u");
	bp = vcf_open(argv[optind], moder);
	hin = hout = vcf_hdr_read(bp);
	if (vc.fn_dict && (vc.flag & VC_VCFIN))
		vcf_dictread(bp, hin, vc.fn_dict);
	bout = vcf_open("-", modew);
	if (!(vc.flag & VC_QCALL)) {
		if (vc.n_sub) {
			vc.sublist = calloc(vc.n_sub, sizeof(int));
			hout = bcf_hdr_subsam(hin, vc.n_sub, vc.subsam, vc.sublist);
		}
		if (vc.flag & VC_CALL) write_header(hout);
		vcf_hdr_write(bout, hout);
	}
	if (vc.flag & VC_CALL) {
		p1 = bcf_p1_init(hout->n_smpl, vc.ploidy);
		if (vc.prior_file) {
			if (bcf_p1_read_prior(p1, vc.prior_file) < 0) {
				fprintf(stderr, "[%s] fail to read the prior AFS.\n", __func__);
				return 1;
			}
		} else bcf_p1_init_prior(p1, vc.prior_type, vc.theta);
		if (vc.n1 > 0) {
			bcf_p1_set_n1(p1, vc.n1);
			bcf_p1_init_subprior(p1, vc.prior_type, vc.theta);
		}
		if (vc.indel_frac > 0.) bcf_p1_indel_prior(p1, vc.indel_frac); // otherwise use the default indel_frac
	}
	if (optind + 1 < argc && !(vc.flag&VC_VCFIN)) {
		void *str2id = bcf_build_refhash(hout);
		if (bcf_parse_region(str2id, argv[optind+1], &tid, &begin, &end) >= 0) {
			bcf_idx_t *idx;
			idx = bcf_idx_load(argv[optind]);
			if (idx) {
				uint64_t off;
				off = bcf_idx_query(idx, tid, begin);
				if (off == 0) {
					fprintf(stderr, "[%s] no records in the query region.\n", __func__);
					return 1; // FIXME: a lot of memory leaks...
				}
				bgzf_seek(bp->fp, off, SEEK_SET);
				bcf_idx_destroy(idx);
			}
		}
	}
	while (vcf_read(bp, hin, b) > 0) {
		int is_indel;
		if ((vc.flag & VC_VARONLY) && strcmp(b->alt, "X") == 0) continue;
		if ((vc.flag & VC_VARONLY) && vc.min_smpl_frac > 0.) {
			extern int bcf_smpl_covered(const bcf1_t *b);
			int n = bcf_smpl_covered(b);
			if ((double)n / b->n_smpl < vc.min_smpl_frac) continue;
		}
		if (vc.n_sub) bcf_subsam(vc.n_sub, vc.sublist, b);
		if (vc.flag & VC_FIX_PL) bcf_fix_pl(b);
		is_indel = bcf_is_indel(b);
		if ((vc.flag & VC_NO_INDEL) && is_indel) continue;
		if ((vc.flag & VC_ACGT_ONLY) && !is_indel) {
			int x;
			if (b->ref[0] == 0 || b->ref[1] != 0) continue;
			x = toupper(b->ref[0]);
			if (x != 'A' && x != 'C' && x != 'G' && x != 'T') continue;
		}
		if (vc.bed && !bed_overlap(vc.bed, hin->ns[b->tid], b->pos, b->pos + strlen(b->ref))) continue;
		if (tid >= 0) {
			int l = strlen(b->ref);
			l = b->pos + (l > 0? l : 1);
			if (b->tid != tid || b->pos >= end) break;
			if (!(l > begin && end > b->pos)) continue;
		}
		++n_processed;
		if (vc.flag & VC_QCALL) { // output QCALL format; STOP here
			bcf_2qcall(hout, b);
			continue;
		}
		if (vc.flag & (VC_CALL|VC_ADJLD)) bcf_gl2pl(b);
		if (vc.flag & VC_CALL) { // call variants
			bcf_p1rst_t pr;
			bcf_p1_cal(b, p1, &pr); // pr.g[3] is not calculated here
			if (n_processed % 100000 == 0) {
				fprintf(stderr, "[%s] %ld sites processed.\n", __func__, (long)n_processed);
				bcf_p1_dump_afs(p1);
			}
			if (pr.p_ref >= vc.pref && (vc.flag & VC_VARONLY)) continue;
			if (vc.n_perm && vc.n1 > 0 && pr.p_chi2 < vc.min_perm_p) { // permutation test
				bcf_p1rst_t r;
				int i, n = 0;
				for (i = 0; i < vc.n_perm; ++i) {
					bcf_shuffle(b, seeds[i]);
					bcf_p1_cal(b, p1, &r);
					if (pr.p_chi2 >= r.p_chi2) ++n;
				}
				pr.perm_rank = n;
			}
			update_bcf1(hout->n_smpl, b, p1, &pr, vc.pref, vc.flag);
		}
		if (vc.flag & VC_ADJLD) { // compute LD
			double f[4], r2;
			if ((r2 = bcf_ld_freq(blast, b, f)) >= 0) {
				kstring_t s;
				s.m = s.l = 0; s.s = 0;
				if (*b->info) kputc(';', &s);
				ksprintf(&s, "NEIR=%.3f;NEIF4=%.3f,%.3f,%.3f,%.3f", r2, f[0], f[1], f[2], f[3]);
				bcf_append_info(b, s.s, s.l);
				free(s.s);
			}
			bcf_cpy(blast, b);
		}
		if (vc.flag & VC_ANNO_MAX) bcf_anno_max(b);
		if (vc.flag & VC_NO_GENO) { // do not output GENO fields
			b->n_gi = 0;
			b->fmt[0] = '\0';
			b->l_str = b->fmt - b->str + 1;
		} else bcf_fix_gt(b);
		vcf_write(bout, hout, b);
	}
	if (vc.prior_file) free(vc.prior_file);
	if (vc.flag & VC_CALL) bcf_p1_dump_afs(p1);
	if (hin != hout) bcf_hdr_destroy(hout);
	bcf_hdr_destroy(hin);
	bcf_destroy(b); bcf_destroy(blast);
	vcf_close(bp); vcf_close(bout);
	if (vc.fn_dict) free(vc.fn_dict);
	if (vc.ploidy) free(vc.ploidy);
	if (vc.n_sub) {
		int i;
		for (i = 0; i < vc.n_sub; ++i) free(vc.subsam[i]);
		free(vc.subsam); free(vc.sublist);
	}
	if (vc.bed) bed_destroy(vc.bed);
	if (seeds) free(seeds);
	if (p1) bcf_p1_destroy(p1);
	return 0;
}
