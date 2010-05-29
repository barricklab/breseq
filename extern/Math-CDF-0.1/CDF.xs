#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#include "cdflib/cdflib.h"

#ifndef PATCHLEVEL
#include <patchlevel.h>
#endif

#if defined(PATCHLEVEL) && (PATCHLEVEL < 5)
#define PL_sv_undef    sv_undef
#endif

static int
not_here(char *s)
{
    croak("%s not implemented on this architecture", s);
    return -1;
}

static double
constant(char *name, int arg)
{
    errno = 0;
    switch (*name) {
    }
    errno = EINVAL;
    return 0;

not_there:
    errno = ENOENT;
    return 0;
}


MODULE = Math::CDF		PACKAGE = Math::CDF		

PROTOTYPES: ENABLE


double
constant(name,arg)
	char *		name
	int		arg


SV *
pnorm(z)
	double z
	PREINIT:
	int which=1, status;
	double p, q, mean=0.0, sd=1.0, bound;
	CODE:
	ST(0) = sv_newmortal();
	(void)cdfnor(&which, &p, &q, &z, &mean, &sd, &status, &bound);
	if(status == 0) {
		sv_setnv( ST(0), (double)p );
	}
	else{
		ST(0) = &PL_sv_undef;
	}

SV *
qnorm(p)
	double p
	PREINIT:
	int which=2, status;
	double z, q, mean=0.0, sd=1.0, bound;
	CODE:
	q = 1.0 - p;
	ST(0) = sv_newmortal();
	(void)cdfnor(&which, &p, &q, &z, &mean, &sd, &status, &bound);
	if(status == 0) {
		sv_setnv( ST(0), (double)z );
	}
	else{
		ST(0) = &PL_sv_undef;
	}

SV *
pt(t, df, ncp = 0.0)
	double t
	double df
	double ncp
	PREINIT:
	int which=1, status;
	double p, q, bound;
	CODE:
	ST(0) = sv_newmortal();
	(void)cdftnc(&which, &p, &q, &t, &df, &ncp, &status, &bound);
	if(status == 0) {
		sv_setnv( ST(0), (double)p );
	}
	else{
		ST(0) = &PL_sv_undef;
	}

SV *
qt(p, df, ncp = 0.0)
	double p
	double df
	double ncp
	PREINIT:
	int which=2, status;
	double t, q, bound;
	CODE:
	q = 1.0 - p;
	ST(0) = sv_newmortal();
	(void)cdftnc(&which, &p, &q, &t, &df, &ncp, &status, &bound);
	if(status == 0) {
		sv_setnv( ST(0), (double)t );
	}
	else{
		ST(0) = &PL_sv_undef;
	}

SV *
pbeta(x, a, b)
	double x
	double a
	double b
	PREINIT:
	int which=1, status;
	double y, p, q, bound;
	CODE:
	y = 1.0 - x;
	ST(0) = sv_newmortal();
	(void)cdfbet(&which, &p, &q, &x, &y, &a, &b, &status, &bound);
	if(status == 0) {
		sv_setnv( ST(0), (double)p );
	}
	else{
		ST(0) = &PL_sv_undef;
	}

SV *
qbeta(p, a, b)
	double p
	double a
	double b
	PREINIT:
	int which=2, status;
	double x, y, q, bound;
	CODE:
	q = 1.0 - p;
	ST(0) = sv_newmortal();
	(void)cdfbet(&which, &p, &q, &x, &y, &a, &b, &status, &bound);
	if(status == 0) {
		sv_setnv( ST(0), (double)x );
	}
	else{
		ST(0) = &PL_sv_undef;
	}

SV *
pchisq(x, df, ncp = 0.0)
	double x
	double df
	double ncp
	PREINIT:
	int which=1, status;
	double p, q, bound;
	CODE:
	ST(0) = sv_newmortal();
	(void)cdfchn(&which, &p, &q, &x, &df, &ncp, &status, &bound);
	if(status == 0) {
		sv_setnv( ST(0), (double)p );
	}
	else{
		ST(0) = &PL_sv_undef;
	}

SV *
qchisq(p, df, ncp = 0.0)
	double p
	double df
	double ncp
	PREINIT:
	int which=2, status;
	double x, q, bound;
	CODE:
	q = 1.0 - p;
	ST(0) = sv_newmortal();
	(void)cdfchn(&which, &p, &q, &x, &df, &ncp, &status, &bound);
	if(status == 0) {
		sv_setnv( ST(0), (double)x );
	}
	else{
		ST(0) = &PL_sv_undef;
	}

SV *
pf(f, dfn, dfd, ncp = 0.0)
	double f
	double dfn
	double dfd
	double ncp
	PREINIT:
	int which=1, status;
	double p, q, bound;
	CODE:
	ST(0) = sv_newmortal();
	(void)cdffnc(&which, &p, &q, &f, &dfn, &dfd, &ncp, &status, &bound);
	if(status == 0) {
		sv_setnv( ST(0), (double)p );
	}
	else{
		ST(0) = &PL_sv_undef;
	}

SV *
qf(p, dfn, dfd, ncp = 0.0)
	double p
	double dfn
	double dfd
	double ncp
	PREINIT:
	int which=2, status;
	double f, q, bound;
	CODE:
	q = 1.0 - p;
	ST(0) = sv_newmortal();
	(void)cdffnc(&which, &p, &q, &f, &dfn, &dfd, &ncp, &status, &bound);
	if(status == 0) {
		sv_setnv( ST(0), (double)f );
	}
	else{
		ST(0) = &PL_sv_undef;
	}

SV *
pgamma(x, shape, scale)
	double x
	double shape
	double scale
	PREINIT:
	int which=1, status;
	double p, q, bound;
	CODE:
	ST(0) = sv_newmortal();
	(void)cdfgam(&which, &p, &q, &x, &shape, &scale, &status, &bound);
	if(status == 0) {
		sv_setnv( ST(0), (double)p );
	}
	else{
		ST(0) = &PL_sv_undef;
	}

SV *
qgamma(p, shape, scale)
	double p
	double shape
	double scale
	PREINIT:
	int which=2, status;
	double x, q, bound;
	CODE:
	q = 1.0 - p;
	ST(0) = sv_newmortal();
	(void)cdfgam(&which, &p, &q, &x, &shape, &scale, &status, &bound);
	if(status == 0) {
		sv_setnv( ST(0), (double)x );
	}
	else{
		ST(0) = &PL_sv_undef;
	}

SV *
ppois(x, lambda)
	double x
	double lambda
	PREINIT:
	int which=1, status;
	double p, q, bound;
	CODE:
	ST(0) = sv_newmortal();
	(void)cdfpoi(&which, &p, &q, &x, &lambda, &status, &bound);
	if(status == 0) {
		sv_setnv( ST(0), (double)p );
	}
	else{
		ST(0) = &PL_sv_undef;
	}

SV *
qpois(p, lambda)
	double p
	double lambda
	PREINIT:
	int which=2, status;
	double x, q, bound;
	CODE:
	q = 1.0 - p;
	ST(0) = sv_newmortal();
	(void)cdfpoi(&which, &p, &q, &x, &lambda, &status, &bound);
	if(status == 0) {
		sv_setnv( ST(0), (double)x );
	}
	else{
		ST(0) = &PL_sv_undef;
	}

SV *
pbinom(x, n, pr)
	double x
	double n
	double pr
	PREINIT:
	int which=1, status;
	double p, q, ompr, bound;
	CODE:
	ompr = 1.0 - pr;
	ST(0) = sv_newmortal();
	(void)cdfbin(&which, &p, &q, &x, &n, &pr, &ompr, &status, &bound);
	if(status == 0) {
		sv_setnv( ST(0), (double)p );
	}
	else{
		ST(0) = &PL_sv_undef;
	}

SV *
pnbinom(x, n, pr)
	double x
	double n
	double pr
	PREINIT:
	int which=1, status;
	double ompr, p, q, bound;
	CODE:
	ompr = 1.0 - pr;
	ST(0) = sv_newmortal();
	(void)cdfnbn(&which, &p, &q, &x, &n, &pr, &ompr, &status, &bound);
	if(status == 0) {
		sv_setnv( ST(0), (double)p );
	}
	else{
		ST(0) = &PL_sv_undef;
	}
