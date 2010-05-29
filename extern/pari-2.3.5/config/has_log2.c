#include <math.h>
double (*f)(double) = log2;
int main(){ return f != log2; }
