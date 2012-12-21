/*
   Copyright (c) 2008, 2010 by Attractive Chaos <attractor@live.co.uk>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

#ifndef KMIN_H
#define KMIN_H

#define KMIN_RADIUS  0.5
#define KMIN_EPS     1e-7
#define KMIN_MAXCALL 50000

typedef double (*kmin_f)(int, double*, void*);
typedef double (*kmin1_f)(double, void*);

#ifdef __cplusplus
extern "C" {
#endif

	double kmin_hj(kmin_f func, int n, double *x, void *data, double r, double eps, int max_calls);
	double kmin_brent(kmin1_f func, double a, double b, void *data, double tol, double *xmin);

#ifdef __cplusplus
}
#endif

#endif
