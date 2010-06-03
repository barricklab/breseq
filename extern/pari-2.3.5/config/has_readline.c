#ifdef READLINE_LIBRARY
#  include <readline.h>
#else
#  include <readline/readline.h>
#endif
main() { char *s = readline("?"); }
