# include <sys/times.h>
# include <sys/time.h>
# include <time.h>
int main() {
  struct tms t; times(&t);
  return (int)t.tms_utime;
}
