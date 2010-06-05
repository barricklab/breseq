#include <stdio.h>
#   ifdef READLINE_LIBRARY
#     include <readline.h>
#   else
#     include <readline/readline.h>
#   endif

#ifndef HAS_RL_LIBRARY_VERSION
extern char *rl_library_version; /* Might be mismatched header, try anyway! */
#endif

main(){ printf("%s", rl_library_version); }
