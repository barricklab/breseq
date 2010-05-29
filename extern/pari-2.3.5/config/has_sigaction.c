#include <signal.h>
main()
{
  struct sigaction sa, oldsa;

  sa.sa_handler = SIG_DFL;
  sigemptyset(&sa.sa_mask);
  sa.sa_flags = SA_NODEFER;

  (void)(sigaction(SIGINT, &sa, &oldsa));
  return 0;
}
