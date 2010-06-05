#include <stdio.h>
#include <sys/times.h>
#include <unistd.h>
main(){
  struct tms t;
  printf("%d%d", times(&t),
#ifdef _SC_CLK_TCK
                 sysconf(_SC_CLK_TCK)
#else
                 CLK_TCK
#endif
  );
}
