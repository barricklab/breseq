#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
main() {
  struct rlimit rip;
  getrlimit(RLIMIT_STACK, &rip);
  setrlimit(RLIMIT_STACK, &rip);
}

