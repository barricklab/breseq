#define RTLD_LAZY 1
#define RTLD_GLOBAL 1
void *dlopen(char *path, int mode);
void *dlsym(void *handle, char *symbol);
char *dlerror(void);
