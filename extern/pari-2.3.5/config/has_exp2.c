#include <math.h>
double (*f)(double) = exp2;
int main(){ return f != exp2; }
