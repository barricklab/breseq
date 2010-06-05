/* $Id: highlvl.c 7522 2005-12-09 18:14:24Z kb $

Copyright (C) 2000  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. */

/*******************************************************************/
/*                                                                 */
/*        SOME GP FUNCTION THAT MAY BE USEFUL OUTSIDE OF IT        */
/*                                                                 */
/*******************************************************************/
#include "pari.h"
#include "paripriv.h"
#include "../graph/rect.h"
#include "../language/anal.h"

#ifdef HAS_DLOPEN
#include <dlfcn.h>

void 
install0(char *name, char *code, char *gpname, char *lib)
{
  void *f, *handle;

  if (! *lib) lib = DL_DFLT_NAME;
  if (! *gpname) gpname = name;
  if (lib) lib = expand_tilde(lib);

#ifndef RTLD_GLOBAL /* OSF1 has dlopen but not RTLD_GLOBAL*/
#  define RTLD_GLOBAL 0
#endif
  handle = dlopen(lib, RTLD_LAZY|RTLD_GLOBAL);

  if (!handle)
  {
    const char *s = dlerror(); if (s) fprintferr("%s\n\n",s);
    if (lib) pari_err(talker,"couldn't open dynamic library '%s'",lib);
    pari_err(talker,"couldn't open dynamic symbol table of process");
  }
  f = dlsym(handle, name);
  if (!f)
  {
    if (lib) pari_err(talker,"can't find symbol '%s' in library '%s'",name,lib);
    pari_err(talker,"can't find symbol '%s' in dynamic symbol table of process",name);
  }
  if (lib) free(lib);
  install(f, gpname, code);
}
#else
#  ifdef _WIN32
#  include <windows.h>
void 
install0(char *name, char *code, char *gpname, char *lib)
{
  FARPROC f;
  HMODULE handle;
#ifdef WINCE
  short wlib[256], wname[256];

  MultiByteToWideChar(CP_ACP, 0, lib, strlen(lib)+1, wlib, 256);
  MultiByteToWideChar(CP_ACP, 0, name, strlen(name)+1, wname, 256);
  lib = wlib;
  name = wname;
#endif

#ifdef DL_DFLT_NAME
  if (! *lib) lib = DL_DFLT_NAME;
#endif
  if (! *gpname) gpname=name;
  if (lib) lib = expand_tilde(lib);
  
  handle = LoadLibrary(lib);
  if (!handle)
  {
    if (lib) pari_err(talker,"couldn't open dynamic library '%s'",lib);
    pari_err(talker,"couldn't open dynamic symbol table of process");
  }
  f = GetProcAddress(handle,name);
  if (!f)
  {
    if (lib) pari_err(talker,"can't find symbol '%s' in library '%s'",name,lib);
    pari_err(talker,"can't find symbol '%s' in dynamic symbol table of process",name);
  }
  if (lib) free(lib);
  install((void*)f,gpname,code);
}
#  else
void 
install0(char *name, char *code, char *gpname, char *lib) { pari_err(archer); }
#endif
#endif

void 
gpinstall(char *s, char *code, char *gpname, char *lib)
{
  if (GP_DATA->flags & SECURE)
  {
    fprintferr("[secure mode]: about to install '%s'. OK ? (^C if not)\n",s);
    hit_return();
  }
  install0(s, code, gpname, lib);
}

#include "highlvl.h"
