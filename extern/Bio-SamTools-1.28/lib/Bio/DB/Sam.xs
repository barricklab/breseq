#ifdef PERL_CAPI
#define WIN32IO_IS_STDIO
#endif
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#ifdef FCGI
 #include <fcgi_stdio.h>
#else
 #ifdef USE_SFIO
  #include <config.h>
 #else
  #include <stdio.h>
 #endif
 #include <perlio.h>
#endif

#include <unistd.h>
#include <math.h>
#include "bam.h"
#include "khash.h"
#include "faidx.h"

/* stolen from bam_aux.c */
#define MAX_REGION 1<<29

typedef tamFile         Bio__DB__Tam;
typedef faidx_t*        Bio__DB__Sam__Fai;
typedef bamFile         Bio__DB__Bam;
typedef bam_header_t*   Bio__DB__Bam__Header;
typedef bam1_t*         Bio__DB__Bam__Alignment;
typedef bam_index_t*    Bio__DB__Bam__Index;
typedef bam_pileup1_t*  Bio__DB__Bam__Pileup;
typedef struct {
  SV* callback;
  SV* data;
} fetch_callback_data;
typedef fetch_callback_data *fetch_callback_dataptr;
typedef struct {
  int    start;
  int    end;
  double width;
  int    reads;
  int*   bin;
} coverage_graph;
typedef coverage_graph *coverage_graph_ptr;

void XS_pack_charPtrPtr( SV * arg, char ** array, int count) {
  int i;
  AV * avref;
  avref = (AV*)sv_2mortal((SV*)newAV());
  for (i=0; i<count; i++) {
    av_push(avref, newSVpv(array[i], strlen(array[i])));
  }
  SvSetSV( arg, newRV((SV*)avref));
}

int bam_fetch_fun (const bam1_t *b, void *data) {
  dSP;
  int count;

  fetch_callback_dataptr fcp;
  SV* callback;
  SV* callbackdata;
  SV* alignment_obj;
  bam1_t *b2;

  fcp          = (fetch_callback_dataptr) data;
  callback     = fcp->callback;
  callbackdata = fcp->data;

  /* turn the bam1_t into an appropriate object */
  /* need to dup it here so that the C layer doesn't reuse the address under Perl */
  b2 = bam_dup1(b);

  alignment_obj = sv_setref_pv(newSV(sizeof(bam1_t)),"Bio::DB::Bam::Alignment",(void*) b2);

  /* set up subroutine stack for the call */
  ENTER;
  SAVETMPS;

  PUSHMARK(SP);
  XPUSHs(sv_2mortal(alignment_obj));
  XPUSHs(callbackdata);
  PUTBACK;

  /* execute the call */
  count = call_sv(callback,G_SCALAR|G_DISCARD);

  FREETMPS;
  LEAVE;

  return 1;
}

int invoke_pileup_callback_fun(uint32_t tid, 
			       uint32_t pos, 
			       int n, 
			       const bam_pileup1_t *pl,
			       void *data) {
  dSP;
  int count,i;
  fetch_callback_dataptr fcp;
  SV*  callback;
  SV*  callbackdata;
  SV*  pileup_obj;
  SV* p;
  SV** pileups;
  AV*  pileup;

  fcp          = (fetch_callback_dataptr) data;
  callback     = fcp->callback;
  callbackdata = fcp->data;

  /* turn the bam_pileup1_t into the appropriate object */
  /* this causes a compiler warning -- ignore it */
 if (0) {
  Newxz(pileups,n,SV*);
  for (i=0;i<n;i++)
        pileups[i] = sv_setref_pv(sv_2mortal(newSV(sizeof(bam_pileup1_t))),
			      "Bio::DB::Bam::Pileup",
			      (void*) &pl[i]);
  pileup = av_make(n,pileups);
  Safefree(pileups);
} else {
  pileup = newAV();
  av_extend(pileup,n);
  for (i=0;i<n;i++) {
    p = newSV(sizeof(bam_pileup1_t));	
    sv_setref_pv(p,"Bio::DB::Bam::Pileup",(void*) &pl[i]);
    av_push(pileup,p);
  } 
}
  
    /* set up subroutine stack for the call */
  ENTER;
  SAVETMPS;

  PUSHMARK(SP);
  XPUSHs(sv_2mortal(newSViv(tid)));
  XPUSHs(sv_2mortal(newSViv(pos)));
  XPUSHs(sv_2mortal(newRV_noinc((SV*)pileup)));
  XPUSHs(callbackdata);
  PUTBACK;

  /* execute the call */
  count = call_sv(callback,G_SCALAR|G_DISCARD);

  FREETMPS;
  LEAVE;
}

int add_pileup_line (const bam1_t *b, void *data) {
  bam_plbuf_t *pileup = (bam_plbuf_t*) data;
  bam_plbuf_push(b,pileup);
  return 0;
}

int add_lpileup_line (const bam1_t *b, void *data) {
  bam_lplbuf_t *pileup = (bam_lplbuf_t*) data;
  bam_lplbuf_push(b,pileup);
  return 0;
}

int coverage_from_pileup_fun (uint32_t tid, 
			      uint32_t pos, 
			      int n, 
			      const bam_pileup1_t *pl, 
			      void *data) {
  coverage_graph_ptr  cgp;
  int                 bin;
  int                 i;
  int                 valid;

  cgp = (coverage_graph_ptr) data;
  cgp->reads += n;

  valid = 0;
  for (i=0;i<n;i++) {
    if (!pl[i].is_del && !pl[i].is_refskip)
        valid++;
  }

  if (n > 0 && pos >= cgp->start && pos <= cgp->end) {
    bin = (pos-cgp->start)/cgp->width;
    cgp->bin[bin] += valid;
  }

  return 0;
}

MODULE = Bio::DB::Sam PACKAGE = Bio::DB::Tam PREFIX=tam_

Bio::DB::Tam
tam_open(packname="Bio::DB::Tam", filename)
   char * packname
   char * filename
 PROTOTYPE: $$
 CODE:
    RETVAL = sam_open(filename);
 OUTPUT:
    RETVAL

void
tam_DESTROY(tam)
  Bio::DB::Tam tam
  PROTOTYPE: $
  CODE:
     sam_close(tam);

Bio::DB::Bam::Header
tam_header_read2(packname="Bio::DB::Tam", filename)
    char * packname
    char * filename
    PROTOTYPE: $$
    CODE:
      RETVAL = sam_header_read2(filename);
    OUTPUT:
      RETVAL

Bio::DB::Bam::Header
tam_header_read(tam)
    Bio::DB::Tam            tam
    PROTOTYPE: $$
    CODE:
      RETVAL = sam_header_read(tam);
    OUTPUT:
      RETVAL

int
tam_read1(tam,header,alignment)
    Bio::DB::Tam            tam
    Bio::DB::Bam::Header    header
    Bio::DB::Bam::Alignment alignment
    CODE:
       RETVAL = sam_read1(tam,header,alignment);
    OUTPUT:
       RETVAL

MODULE = Bio::DB::Sam PACKAGE = Bio::DB::Sam::Fai PREFIX=fai_

Bio::DB::Sam::Fai
fai_load(packname="Bio::DB::Sam::Fai", filename)
  char * packname
  char * filename
 PROTOTYPE: $$
 CODE:
    RETVAL = fai_load(filename);
 OUTPUT:
    RETVAL

void
fai_destroy(fai)
  Bio::DB::Sam::Fai fai
  PROTOTYPE: $
  CODE:
    fai_destroy(fai);

SV*
fai_fetch(fai,reg)
  Bio::DB::Sam::Fai    fai
    const char *reg
  PROTOTYPE: $$$
  PREINIT:
    char     *seq;
    int       len;
  CODE:
    seq = fai_fetch(fai,reg,&len);
    if (seq == NULL)
       XSRETURN_EMPTY;
    RETVAL = newSVpv(seq,len);
    free((void*)seq);
  OUTPUT:
    RETVAL


MODULE = Bio::DB::Sam PACKAGE = Bio::DB::Bam PREFIX=bam_

Bio::DB::Bam
bam_open(packname, filename, mode="r")
      char * packname
      char * filename
      char * mode
      PROTOTYPE: $$$
      CODE:
        RETVAL = bam_open(filename,mode);
      OUTPUT:
      RETVAL

void
bam_DESTROY(bam)
   Bio::DB::Bam bam
PROTOTYPE: $
CODE:
   bam_close(bam);

int
bam_index_build(packname, filename)
   char *      packname
   const char * filename
  CODE:
     RETVAL = bam_index_build(filename);
  OUTPUT:
     RETVAL

void
bam_sort_core(packname, is_by_qname=0, filename, prefix, max_mem=500000000)
   char * packname
   int    is_by_qname
   char * filename
   char * prefix
   int    max_mem
 PROTOTYPE: $$$$$
 CODE:
   bam_sort_core(is_by_qname,filename,prefix,max_mem);

Bio::DB::Bam::Index
bam_index_open(packname="Bio::DB::Bam", filename)
      char * packname
      char * filename
    PROTOTYPE: $$
    CODE:
    RETVAL = bam_index_load(filename);
    OUTPUT:
    RETVAL

Bio::DB::Bam::Header
bam_header(bam)
    Bio::DB::Bam bam
    PROTOTYPE: $
    PREINIT:
      bam_header_t *bh;
      int64_t       result;
    CODE:
      result = bgzf_seek(bam,0,0);
      bh = bam_header_read(bam);
      RETVAL = bh;
    OUTPUT:
      RETVAL

int
bam_header_write(bam,header)
    Bio::DB::Bam         bam
    Bio::DB::Bam::Header header
    PROTOTYPE: $$
    CODE:
      bgzf_seek(bam,0,0);
      RETVAL= bam_header_write(bam,header);
    OUTPUT:
      RETVAL

char*
bam_tell(bam)
    Bio::DB::Bam bam
PROTOTYPE: $
CODE:
    int64_t t = bam_tell(bam);
    char    string[128];
    sprintf(string,"%llu",t);
    RETVAL = string;
OUTPUT:
    RETVAL

void
bam_seek(bam,pos,dir)
    Bio::DB::Bam bam
    int pos
    int dir
PROTOTYPE: $$$
CODE:
    bam_seek(bam,pos,dir);

Bio::DB::Bam::Alignment
bam_read1(bam)
    Bio::DB::Bam bam
  PROTOTYPE: $
  PREINIT:
    bam1_t *b;
  CODE:
    b = bam_init1();
    if (bam_read1(bam,b) >= 0) {
      RETVAL = b;
    }
    else
       XSRETURN_EMPTY;
  OUTPUT:
    RETVAL      

int
bam_write1(bam,align)
    Bio::DB::Bam            bam
    Bio::DB::Bam::Alignment align
   PROTOTYPE: $$
   CODE:
      RETVAL = bam_write1(bam,align);
   OUTPUT:
      RETVAL

MODULE = Bio::DB::Sam PACKAGE = Bio::DB::Bam::Alignment PREFIX=bama_

Bio::DB::Bam::Alignment
bama_new(package="Bio::DB::Bam::Alignment")
   char * package
   PROTOTYPE: $
   CODE:
      RETVAL = bam_init1();
   OUTPUT:
      RETVAL

void
bama_DESTROY(b)
  Bio::DB::Bam::Alignment b
PROTOTYPE: $
CODE:
    bam_destroy1(b);

int
bama_tid(b,...)
    Bio::DB::Bam::Alignment b
PROTOTYPE: $;$
CODE:
    if (items > 1)
      b->core.tid = SvIV(ST(1));
    RETVAL=b->core.tid;
OUTPUT:
    RETVAL

int
bama_pos(b,...)
    Bio::DB::Bam::Alignment b
PROTOTYPE: $;$
CODE:
    if (items > 1)
      b->core.pos = SvIV(ST(1));
    RETVAL=b->core.pos;
OUTPUT:
    RETVAL

int
bama_calend(b)
  Bio::DB::Bam::Alignment b
PROTOTYPE: $
CODE:
   RETVAL=bam_calend(&b->core,bam1_cigar(b));
OUTPUT:
   RETVAL    

int
bama_cigar2qlen(b)
  Bio::DB::Bam::Alignment b
PROTOTYPE: $
CODE:
   RETVAL=bam_cigar2qlen(&b->core,bam1_cigar(b));
OUTPUT:
   RETVAL    

int
bama_qual(b,...)
    Bio::DB::Bam::Alignment b
PROTOTYPE: $;$
CODE:
    if (items > 1)
      b->core.qual = SvIV(ST(1));
    RETVAL=b->core.qual;
OUTPUT:
    RETVAL

int
bama_flag(b,...)
    Bio::DB::Bam::Alignment b
PROTOTYPE: $;$
CODE:
    if (items > 1)
      b->core.flag = SvIV(ST(1));
    RETVAL=b->core.flag;
OUTPUT:
    RETVAL

int
bama_n_cigar(b,...)
    Bio::DB::Bam::Alignment b
PROTOTYPE: $;$
CODE:
  if (items > 1)
    b->core.n_cigar = SvIV(ST(1));
    RETVAL=b->core.n_cigar;
OUTPUT:
    RETVAL

int
bama_l_qseq(b,...)
    Bio::DB::Bam::Alignment b
PROTOTYPE: $;$
CODE:
    if (items > 1)
      b->core.l_qseq = SvIV(ST(1));
    RETVAL=b->core.l_qseq;
OUTPUT:
    RETVAL

SV*
bama_qseq(b)
Bio::DB::Bam::Alignment b
PROTOTYPE: $
PREINIT:
    char* seq;
    int   i;
CODE:
    seq = Newxz(seq,b->core.l_qseq+1,char);
    for (i=0;i<b->core.l_qseq;i++) {
      seq[i]=bam_nt16_rev_table[bam1_seqi(bam1_seq(b),i)];
    }
    RETVAL = newSVpv(seq,b->core.l_qseq);
    Safefree(seq);
OUTPUT:
    RETVAL

SV*
bama__qscore(b)
Bio::DB::Bam::Alignment b
PROTOTYPE: $
CODE:
    RETVAL = newSVpv(bam1_qual(b),b->core.l_qseq);
OUTPUT:
    RETVAL

int
bama_mtid(b,...)
    Bio::DB::Bam::Alignment b
PROTOTYPE: $;$
CODE:
    if (items > 1)
      b->core.mtid = SvIV(ST(1));
    RETVAL=b->core.mtid;
OUTPUT:
    RETVAL

int
bama_mpos(b,...)
    Bio::DB::Bam::Alignment b
PROTOTYPE: $;$
CODE:
    if (items > 1)
      b->core.mpos = SvIV(ST(1));
    if (b->core.pos <= 0)
      XSRETURN_UNDEF;
    RETVAL=b->core.mpos;
OUTPUT:
    RETVAL

int
bama_isize(b,...)
    Bio::DB::Bam::Alignment b
PROTOTYPE: $;$
CODE:
    if (items > 1)
      b->core.isize = SvIV(ST(1));
    RETVAL=b->core.isize;
OUTPUT:
    RETVAL

int
bama_l_aux(b,...)
    Bio::DB::Bam::Alignment b
PROTOTYPE: $;$
CODE:
    if (items > 1)
      b->l_aux = SvIV(ST(1));
    RETVAL=b->l_aux;
OUTPUT:
    RETVAL

char*
bama_aux(b)
   Bio::DB::Bam::Alignment b
PREINIT:
   uint8_t *s;
   uint8_t type, key[2];
   char    str[8192];
CODE:
   s = bam1_aux(b);
   str[0] = '\0';

   int left  = sizeof(str) - strlen(str);
   while (left > 0 && (s < b->data + b->data_len)) {
        char* d   = str+strlen(str); 

	key[0] = s[0]; 
	key[1] = s[1];
 	left -= snprintf(d, left, "%c%c:", key[0], key[1]);

	d    += 3;
	s    += 2;
	type = *s++; 

	if (left <= 0) continue;

	if (type == 'A')      { left -= snprintf(d, left, "A:%c", *s);           s++; }
	else if (type == 'C') { left -= snprintf(d, left, "i:%u", *s);           s++; }
	else if (type == 'c') { left -= snprintf(d, left, "i:%d", *s);           s++; }
	else if (type == 'S') { left -= snprintf(d, left, "i:%u", *(uint16_t*)s);s += 2; }
	else if (type == 's') { left -= snprintf(d, left, "i:%d", *(int16_t*)s); s += 2; }
	else if (type == 'I') { left -= snprintf(d, left, "i:%u", *(uint32_t*)s);s += 4; }
	else if (type == 'i') { left -= snprintf(d, left, "i:%d", *(int32_t*)s); s += 4; }
	else if (type == 'f') { left -= snprintf(d, left, "f:%g", *(float*)s);   s += 4; }
	else if (type == 'd') { left -= snprintf(d, left, "d:%lg", *(double*)s); s += 8; }
	else if (type == 'Z' || type == 'H') { left -= snprintf(d, left, "%c:", type); 
	                                       strncat(d,s,left);
					       while (*s++) {}
					       left = sizeof(str) - strlen(str);
	                                     }
	if (left <= 0) continue;
	strncat(d,"\t",left);
	left--;
   }	  
   str[strlen(str)-1] = '\0';
   RETVAL = str;
OUTPUT:
   RETVAL

SV*
bama_aux_get(b,tag)
   Bio::DB::Bam::Alignment b
   char*               tag
PROTOTYPE: $$
PREINIT:
   int           type;
   uint8_t       *s;
CODE:
   s    = bam_aux_get_core(b,tag);
   if (s==0)
      XSRETURN_EMPTY;
   type = *s++;
   switch (type) {
   case 'c':
     RETVAL = newSViv((int32_t)*(int8_t*)s);
     break;
   case 'C':
     RETVAL = newSViv((int32_t)*(uint8_t*)s);
     break;
   case 's':
     RETVAL = newSViv((int32_t)*(int16_t*)s);
     break;
   case 'S':
     RETVAL = newSViv((int32_t)*(uint16_t*)s);
     break;
   case 'i':
     RETVAL = newSViv(*(int32_t*)s);
     break;
   case 'I':
     RETVAL = newSViv((int32_t)*(uint32_t*)s);
     break;
   case 'f':
     RETVAL = newSVnv(*(float*)s);
     break;
   case 'Z':
   case 'H':
     RETVAL = newSVpv((char*)s,0);
     break;
   case 'A':
     RETVAL = newSVpv((char*)s,1);
     break;
   default:
     XSRETURN_EMPTY;
   }
OUTPUT:
   RETVAL

void
bama_aux_keys(b)
Bio::DB::Bam::Alignment b
PROTOTYPE: $
PREINIT:
   uint8_t *s;
   uint8_t type;
PPCODE:
   {
     s = bam1_aux(b);  /* s is a khash macro */
     while (s < b->data + b->data_len) {
       XPUSHs(sv_2mortal(newSVpv(s,2)));
       s   += 2; 
       type = *s++;
       if      (type == 'A') { ++s; }
       else if (type == 'C') { ++s; }
       else if (type == 'c') { ++s; }
       else if (type == 'S') { s += 2; }
       else if (type == 's') { s += 2; }
       else if (type == 'I') { s += 4; }
       else if (type == 'i') { s += 4; }
       else if (type == 'f') { s += 4; }
       else if (type == 'Z' || type == 'H') { while (*s) ++(s); ++(s); }
     }
   }

SV*
bama_data(b,...)
    Bio::DB::Bam::Alignment b
PROTOTYPE: $;$
PREINIT:
    STRLEN  len;
CODE:
    if (items > 1) {
      b->data     = SvPV(ST(1),len);
      b->data_len = len;
    }
    RETVAL=newSVpv(b->data,b->data_len);
OUTPUT:
    RETVAL

int
bama_data_len(b,...)
    Bio::DB::Bam::Alignment b
PROTOTYPE: $;$
CODE:
    if (items > 1)
      b->data_len = SvIV(ST(1));
    RETVAL=b->data_len;
OUTPUT:
    RETVAL

int
bama_m_data(b,...)
    Bio::DB::Bam::Alignment b
PROTOTYPE: $;$
CODE:
    if (items > 1) {
      b->m_data = SvIV(ST(1));
    }
    RETVAL=b->m_data;
OUTPUT:
    RETVAL

SV*
bama_qname(b)
  Bio::DB::Bam::Alignment b
PROTOTYPE: $
CODE:
    RETVAL=newSVpv(bam1_qname(b),0);
OUTPUT:
    RETVAL

int
bama_paired(b)
  Bio::DB::Bam::Alignment b
PROTOTYPE: $
CODE:
    RETVAL=(b->core.flag&BAM_FPAIRED) != 0;
OUTPUT:
  RETVAL

int
bama_proper_pair(b)
  Bio::DB::Bam::Alignment b
PROTOTYPE: $
CODE:
    RETVAL=(b->core.flag&BAM_FPROPER_PAIR) != 0;
OUTPUT:
  RETVAL

int
bama_unmapped(b)
  Bio::DB::Bam::Alignment b
PROTOTYPE: $
CODE:
    RETVAL=(b->core.flag&BAM_FUNMAP) != 0;
OUTPUT:
  RETVAL

int
bama_munmapped(b)
  Bio::DB::Bam::Alignment b
PROTOTYPE: $
CODE:
    RETVAL=(b->core.flag&BAM_FMUNMAP) != 0;
OUTPUT:
  RETVAL

int
bama_reversed(b)
  Bio::DB::Bam::Alignment b
PROTOTYPE: $
CODE:
  RETVAL=bam1_strand(b);
OUTPUT:
  RETVAL

int
bama_mreversed(b)
  Bio::DB::Bam::Alignment b
PROTOTYPE: $
CODE:
  RETVAL=bam1_mstrand(b);
OUTPUT:
  RETVAL

SV*
bama_cigar(b)
  Bio::DB::Bam::Alignment b
PROTOTYPE: $
PREINIT:
    int        i;
    uint32_t  *c;
    AV        *avref;
CODE:
    avref = (AV*) sv_2mortal((SV*)newAV());
    c     = bam1_cigar(b);
    for (i=0;i<b->core.n_cigar;i++)
      av_push(avref, newSViv(c[i]));
    RETVAL = (SV*) newRV((SV*)avref); 
OUTPUT:
  RETVAL

MODULE = Bio::DB::Sam PACKAGE = Bio::DB::Bam::Header PREFIX=bam_

Bio::DB::Bam::Header
bam_new(packname=Bio::DB::Bam::Header)
PROTOTYPE: $
CODE:
    RETVAL = bam_header_init();
OUTPUT:
    RETVAL

int
bam_n_targets(bamh)
  Bio::DB::Bam::Header bamh
  PROTOTYPE: $
  CODE:
    RETVAL = bamh->n_targets;
  OUTPUT:
    RETVAL

SV*
bam_target_name(bamh)
  Bio::DB::Bam::Header bamh
  PROTOTYPE: $
  PREINIT:
    int i;
    AV * avref;
  CODE:
    avref = (AV*) sv_2mortal((SV*)newAV());
    for (i=0;i<bamh->n_targets;i++)
      av_push(avref, newSVpv(bamh->target_name[i],0));
    RETVAL = (SV*) newRV((SV*)avref); 
  OUTPUT:
    RETVAL

SV*
bam_target_len(bamh)
    Bio::DB::Bam::Header bamh
  PROTOTYPE: $
  PREINIT:
    int i;
    AV * avref;
  CODE:
    avref = (AV*) sv_2mortal((SV*)newAV());
    for (i=0;i<bamh->n_targets;i++)
       av_push(avref, newSViv(bamh->target_len[i]));
    RETVAL = (SV*) newRV((SV*)avref); 
  OUTPUT:
    RETVAL

SV*
bam_text(bamh, ...)
  Bio::DB::Bam::Header bamh
  PREINIT:
    char   *newtext;
    STRLEN n;
  CODE:
    /* in case text is not null terminated, we copy it */
    RETVAL = newSVpv(bamh->text,bamh->l_text);
    if (items > 1) {
      newtext = (char*) SvPV(ST(1),n);
      bamh->text   = newtext;
      bamh->l_text = n;
    }
  OUTPUT:
    RETVAL

void
bam_parse_region(bamh,region)
    Bio::DB::Bam::Header bamh
    char*            region
    PROTOTYPE: $
    PREINIT:
       int seqid,start,end;
    PPCODE:
    {
      bam_parse_region(bamh,
		       region,
		       &seqid,
		       &start,
		       &end);
      if (seqid < 0)
	XSRETURN_EMPTY;
      else {
	EXTEND(sp,3);
	PUSHs(sv_2mortal(newSViv(seqid)));
	PUSHs(sv_2mortal(newSViv(start)));
	PUSHs(sv_2mortal(newSViv(end)));
      }
    }

void
bam_view1(bamh,alignment)
     Bio::DB::Bam::Header     bamh
     Bio::DB::Bam::Alignment  alignment
     PROTOTYPE: $$
     CODE:
       bam_view1(bamh,alignment);

void
bam_DESTROY(bamh)
  Bio::DB::Bam::Header bamh
  PROTOTYPE: $
  CODE:
    bam_header_destroy(bamh);

MODULE = Bio::DB::Sam PACKAGE = Bio::DB::Bam::Index PREFIX=bami_

int
bami_fetch(bai,bfp,ref,start,end,callback,callbackdata=&PL_sv_undef)
  Bio::DB::Bam::Index bai
  Bio::DB::Bam        bfp
  int   ref
  int   start
  int   end
  CV*   callback
  SV*   callbackdata
PREINIT:
  fetch_callback_data fcd;
CODE:
  {
    fcd.callback = (SV*) callback;
    fcd.data     = callbackdata;
    RETVAL = bam_fetch(bfp,bai,ref,start,end,&fcd,bam_fetch_fun);
  }
OUTPUT:
    RETVAL

void
bami_lpileup(bai,bfp,ref,start,end,callback,callbackdata=&PL_sv_undef)
  Bio::DB::Bam::Index bai
  Bio::DB::Bam        bfp
  int   ref
  int   start
  int   end
  CV*   callback
  SV*   callbackdata
PREINIT:  
  fetch_callback_data fcd;
  bam_lplbuf_t        *pileup;
CODE:
  fcd.callback = (SV*) callback;
  fcd.data     = callbackdata;
  pileup       = bam_lplbuf_init(invoke_pileup_callback_fun,(void*)&fcd);
  bam_fetch(bfp,bai,ref,start,end,(void*)pileup,add_lpileup_line);
  bam_lplbuf_push(NULL,pileup);
  bam_lplbuf_destroy(pileup);

void
bami_pileup(bai,bfp,ref,start,end,callback,callbackdata=&PL_sv_undef)
  Bio::DB::Bam::Index bai
  Bio::DB::Bam        bfp
  int   ref
  int   start
  int   end
  CV*   callback
  SV*   callbackdata
PREINIT:  
  fetch_callback_data fcd;
  bam_plbuf_t        *pileup;
CODE:
  fcd.callback = (SV*) callback;
  fcd.data     = callbackdata;
  pileup       = bam_plbuf_init(invoke_pileup_callback_fun,(void*)&fcd);
  bam_fetch(bfp,bai,ref,start,end,(void*)pileup,add_pileup_line);
  bam_plbuf_push(NULL,pileup);
  bam_plbuf_destroy(pileup);

AV*
bami_coverage(bai,bfp,ref,start,end,bins=0)
    Bio::DB::Bam::Index bai
    Bio::DB::Bam        bfp
    int             ref
    int             start
    int             end
    int             bins
PREINIT:
    coverage_graph  cg;
    bam_plbuf_t    *pileup;
    AV*             array;
    SV*             cov;
    int             i;
    bam_header_t   *bh;
CODE:
  {
      if (end >= MAX_REGION) {
          bgzf_seek(bfp,0,0);
          bh  = bam_header_read(bfp);
          end = bh->target_len[ref];
          bam_header_destroy(bh);
      }
      if ((bins==0) || (bins > (end-start)))
         bins = end-start;

      /* coverage graph used to communicate to our callback
	  the region we are sampling */
      cg.start = start;
      cg.end   = end;
      cg.reads = 0;
      cg.width = ((double)(end-start))/bins;
      Newxz(cg.bin,bins+1,int);

      /* accumulate coverage into the coverage graph */
      pileup   = bam_plbuf_init(coverage_from_pileup_fun,(void*)&cg);
      bam_fetch(bfp,bai,ref,start,end,(void*)pileup,add_pileup_line);
      bam_plbuf_push(NULL,pileup);
      bam_plbuf_destroy(pileup);

      /* now normalize to coverage/bp and convert into an array */
      array = newAV();
      av_extend(array,bins);
      if (cg.reads > 0) {
            for  (i=0;i<bins;i++)
	    	av_store(array,i,newSVnv(((float)cg.bin[i])/cg.width));
      }
      Safefree(cg.bin);
      RETVAL = array;
      sv_2mortal((SV*)RETVAL);  /* this fixes a documented bug in perl typemap */

  }
OUTPUT:
    RETVAL

void
bami_DESTROY(bai)
  Bio::DB::Bam::Index bai
  CODE:
    bam_index_destroy(bai);

MODULE = Bio::DB::Sam PACKAGE = Bio::DB::Bam::Pileup PREFIX=pl_

int
pl_qpos(pl)
  Bio::DB::Bam::Pileup pl
  CODE:
    RETVAL = pl->qpos;
  OUTPUT:
    RETVAL

int
pl_pos(pl)
  Bio::DB::Bam::Pileup pl
  CODE:
    RETVAL = pl->qpos+1;
  OUTPUT:
    RETVAL

int
pl_indel(pl)
  Bio::DB::Bam::Pileup pl
  CODE:
    RETVAL = pl->indel;
  OUTPUT:
    RETVAL

int
pl_level(pl)
  Bio::DB::Bam::Pileup pl
  CODE:
    RETVAL = pl->level;
  OUTPUT:
    RETVAL

int
pl_is_del(pl)
  Bio::DB::Bam::Pileup pl
  CODE:
    RETVAL = pl->is_del;
  OUTPUT:
    RETVAL

int
pl_is_refskip(pl)
  Bio::DB::Bam::Pileup pl
  CODE:
    RETVAL = pl->is_refskip;
  OUTPUT:
    RETVAL

int
pl_is_head(pl)
  Bio::DB::Bam::Pileup pl
  CODE:
    RETVAL = pl->is_head;
  OUTPUT:
    RETVAL

int
pl_is_tail(pl)
  Bio::DB::Bam::Pileup pl
  CODE:
    RETVAL = pl->is_tail;
  OUTPUT:
    RETVAL

Bio::DB::Bam::Alignment
pl_b(pl)
  Bio::DB::Bam::Pileup pl
  CODE:
    RETVAL = bam_dup1(pl->b);
  OUTPUT:
     RETVAL

Bio::DB::Bam::Alignment
pl_alignment(pl)
  Bio::DB::Bam::Pileup pl
  CODE:
    RETVAL = bam_dup1(pl->b);
  OUTPUT:
     RETVAL



