#include "pari.h"
#include "anal.h"

GEN   gen_0, gen_1, gen_m1, gen_2;
pari_sp top, bot, avma;
size_t memused = 0;
ulong  DEBUGLEVEL,DEBUGMEM = 0;

void specinit()
{
  long size = 100000L;
  bot = (pari_sp)malloc(size);
  top = avma = bot + size;
  gen_0 = (GEN)malloc(2);
  gen_0[0] = evaltyp(t_INT)|evallg(2);
  gen_0[1] = evallgefint(2);
  gen_1 = utoipos(1);
  gen_m1= utoineg(1);
  gen_2 = utoipos(2);
}

void gerepileall(pari_sp av, int n, ...){}
void pari_err(long x, ...) { exit (0); }
void pari_warn(long x, ...) { exit (0); }
char*gprealloc(void *pointer, size_t size) { return (char*)realloc(pointer,size); }
char*gpmalloc(size_t size) { return (char*)malloc(size); }
GEN gerepileuptoint(pari_sp av, GEN q){ return q; }
void gerepilemanysp(pari_sp av, pari_sp tetpil, GEN* gptr[], int n){}
void gerepilemany(pari_sp av, GEN* gptr[], int n){}
void fprintferr(const char *x, ...){}
GEN gcopy(GEN x){return x;}
GEN gtofp(GEN x, long prec){return x;}
void gaffect(GEN x, GEN y){}
GEN ceilr(GEN x){return x;}
GEN roundr(GEN x){return x;}
const long lontyp[] = {0};

void sorstring(ulong x)
{
#ifdef LONG_IS_64BIT
  printf("%016lx  ", x);
#else
  printf("%08lx  ", x);
#endif
}

void _voiri(GEN x)
{
  long i, lx = lgefint(x);
  GEN y = int_MSW(x);
  /* sorstring(x[0]); depends on the kernel and contains no useful info */
  sorstring(x[1]);
  for (i=2; i < lx; i++, y = int_precW(y)) sorstring(*y);
  printf("\n");
}
void _voirr(GEN x)
{
  long i, lx = lg(x);
  for (i=1; i < lx; i++) sorstring(x[i]);
  printf("\n");
}

int main()
{
  GEN x,y,r,z, xr,yr;
  specinit();
  x = utoipos(187654321UL);
  y = utoineg(12345678UL);
  printf("INT: %ld\n", itos(x));
  printf("conv:"); _voiri(x);
  printf("+:"); _voiri(addii(x,y));
  printf("-:"); _voiri(subii(x,y));
  printf("*:"); _voiri(mulii(x,y));
  printf("/:"); _voiri(dvmdii(x,y, &z));
  printf("rem:"); _voiri(z);
  printf("pow:\n");
  z = mulii(x,x); _voiri(z);
  z = mulii(z,z); _voiri(z);
  z = mulii(z,z); _voiri(z);
  z = mulii(z,z); _voiri(z);
  z = mulii(z,z); _voiri(z);
  printf("invmod:"); invmod(y,z,&r); _voiri(r);
  xr = itor(x, DEFAULTPREC);
  yr = itor(y, DEFAULTPREC);
  printf("\nREAL: %f\n", rtodbl(xr));
  printf("conv1:"); _voirr(xr);
  printf("conv2:"); _voirr(dbltor(rtodbl(xr)));
  printf("+:"); _voirr(addrr(xr,yr));
  printf("-:"); _voirr(subrr(xr,yr));
  printf("*:"); _voirr(mulrr(xr,yr));
  printf("/:"); _voirr(divrr(xr,yr));
  printf("gcc bug?:"); _voirr(divrs(dbltor(3.),2));
  return 0;
}
