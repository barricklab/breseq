file="$objdir/paricfg.h"

echo Extracting $file
rm -f $file

case "$optimization" in
  full) ;;
  *) debuginfo=" -- $optimization";;
esac

if test -n "$perl"; then
  case "$osname" in
    os2) gphelp="perl -S gphelp -detex" ;;
    *) gphelp="$bindir/gphelp" ;;
  esac
fi

cat > $file << EOT
/*  This file was created by Configure. Any change made to it will be lost
 *  next time Configure is run. */
#ifndef __CONFIG_H__
#define __CONFIG_H__
EOT

# case "$osname" in
#   cygwin*) echo '#ifdef __CYGWIN32__' >> $file ;
# 	     echo '#include "libpari_dll.h"' >> $file ;
# 	     echo '#endif' >> $file ;;
# esac

cat >> $file << EOT
#define UNIX
#define GPHELP "$gphelp"
#define GPDATADIR "$datadir"
#define SHELL_Q '\\$shell_q'

#define PARIVERSION "GP/PARI CALCULATOR Version ${version}.${patch} (${status})"
#ifdef __cplusplus
# define PARIINFO "C++ ${pretty}${debuginfo}"
#else
# define PARIINFO "${pretty}${debuginfo}"
#endif
#define PARI_VERSION_CODE `expr $VersionMajor \\* 65536 + $VersionMinor \\* 256 + $patch`
#define PARI_VERSION(a,b,c) (((a) << 16) + ((b) << 8) + (c))
#define PARI_VERSION_SHIFT 8

#define PARI_DOUBLE_FORMAT ${doubleformat}
EOT

if test -n "$__gnuc__"; then
  echo "#define GCC_VERSION \"$__gnuc__\"" >> $file
fi

case $asmarch in
  none) echo '#define __HAS_NO_ASM__' >> $file;;
  hppa) echo '#define __HPPA__' >> $file;;
esac

if test -n "$ASMINLINE"; then
  echo '#define ASMINLINE' >> $file
fi

if test -n "$gzip"; then
  cat >> $file << EOT

/*  Location of GNU gzip program (enables reading of .Z and .gz files). */
#define GNUZCAT
#define ZCAT "$gzip -dc"

EOT
else if test -n "$zcat"; then
  cat >> $file << EOT

/*  Location of zcat program (enables reading of .Z files). */
#define ZCAT "$zcat"

EOT
  fi
fi

if test -n "$readline"; then
  cat >> $file <<EOT
/*  Use GNU readline library
 *  This is incompatible with SUN commandtools (not with shelltools).
 */
#define READLINE "$readline_version"
EOT
  if test -n "$rl_appendchar"; then
    echo '#define HAS_COMPLETION_APPEND_CHAR' >> $file
  fi
  if test -n "$rl_genericbind"; then
    echo '#define HAS_RL_GENERIC_BIND' >> $file
  fi
  if test -n "$rl_refresh_line_oldproto"; then
    echo '#define OLD_RL_REFRESH_LINE' >> $file
  fi
  if test -n "$rl_save_prompt"; then
    echo '#define HAS_RL_SAVE_PROMPT' >> $file
  else
    if test -n "$_rl_save_prompt"; then
      echo '#define HAS_UNDERSCORE_RL_SAVE_PROMPT' >> $file
    fi
  fi 
  if test -n "$rl_message"; then
    echo '#define HAS_RL_MESSAGE' >> $file
  fi
  if test -n "$rl_completion_matches"; then
    echo '#define HAS_RL_COMPLETION_MATCHES' >> $file
  fi
  if test -n "$rl_completion_func_t"; then
    echo '#define HAS_RL_COMPLETION_FUNC_T' >> $file
  fi
  if test -n "$rl_history"; then
    echo '#define HAS_HISTORY_H' >> $file
  fi
  if test -n "$rl_attempted_completion_over"; then
    echo '#define HAS_RL_ATTEMPTED_COMPLETION_OVER' >> $file
  fi
  if test -n "$rl_completion_query_items"; then
    echo '#define HAS_RL_COMPLETION_QUERY_ITEMS' >> $file
  fi
  if test -n "$rl_bind_key_in_map"; then
    echo '#define HAS_RL_BIND_KEY_IN_MAP' >> $file
  fi
  if test -n "$rl_library_version"; then
    echo '#define HAS_RL_LIBRARY_VERSION' >> $file
  fi
  if test -n "$rl_fake_library_version"; then
    echo '#define FAKE_RL_LIBRARY_VERSION' >> $file
  fi
  echo >> $file
fi

case "$CPPF_defined" in
yes) echo '#define CPPFunction_defined' >> $file;;
esac

case "$sizeof_long" in
8) echo '#define LONG_IS_64BIT' >> $file;;
esac

case "$has_exp2" in
yes) echo '#define HAS_EXP2' >> $file;;
esac
case "$has_log2" in
yes) echo '#define HAS_LOG2' >> $file;;
esac

case "$has_getrusage" in
yes) echo '#define USE_GETRUSAGE 1' >> $file;;
  *) case "$has_times" in
     yes) echo '#define USE_TIMES 1' >> $file;;
       *) case "$has_ftime" in
          yes) echo '#define USE_FTIME 1' >> $file;;
          esac;;
     esac;;
esac

case $has_sigaction in
yes) echo '#define HAS_SIGACTION' >> $file;;
esac

case $has_dlopen in
yes|builtin) cat >> $file << EOT
#define HAS_DLOPEN
EOT
;;
esac

case $has_getrlimit in
yes) echo '#define STACK_CHECK' >> $file;;
esac

case "$osname" in
os2)  echo '#define STACK_CHECK' >> $file;;
esac

case $has_vsnprintf in
yes) echo '#define HAS_VSNPRINTF' >> $file;;
esac

case $has_TIOCGWINSZ in
yes) echo '#define HAS_TIOCGWINSZ' >> $file;;
esac

case $has_strftime in
yes) echo '#define HAS_STRFTIME' >> $file;;
esac

case $has_opendir in
yes) echo '#define HAS_OPENDIR' >> $file;;
esac

case $has_stat in
yes) echo '#define HAS_STAT' >> $file;;
esac

echo '#endif' >> $file
