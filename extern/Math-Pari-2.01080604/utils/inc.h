#include <unistd.h>
#include <stdarg.h>
#include <setjmp.h>
#include <signal.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/time.h>
#include <memory.h>
#include <sys/resource.h>

int x;
int y = CLK_TCK * x;
