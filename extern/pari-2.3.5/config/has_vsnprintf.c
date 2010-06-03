#include <stdio.h>
#include <stdarg.h>

main() {}
int f(int i,...) { char s[1]; va_list ap; va_start(ap,i);
  vsnprintf(s,1,"",ap); return 0;
}


