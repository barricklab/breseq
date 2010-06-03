#include <stdlib.h>
#include <windows.h>

int
main(char **argv, int argc)
{
  HKEY handle;
  const char *key = "AcroExch.Document\\shell\\open\\command";
  const long SZ = 512;
  char str[SZ];
  long L = SZ, t, i;

  (void)RegOpenKeyEx(HKEY_CLASSES_ROOT, key, 0, KEY_READ, &handle);
  (void)RegQueryValueEx(handle, NULL, 0, &t, (LPBYTE)str, &L);
  RegCloseKey(handle);
  printf("%s\n", str); exit(0);
}

