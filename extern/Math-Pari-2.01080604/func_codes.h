#define PERL_constant_NOTFOUND	1
#define PERL_constant_NOTDEF	2
#define PERL_constant_ISIV	3
#define PERL_constant_ISNO	4
#define PERL_constant_ISNV	5
#define PERL_constant_ISPV	6
#define PERL_constant_ISPVN	7
#define PERL_constant_ISSV	8
#define PERL_constant_ISUNDEF	9
#define PERL_constant_ISUV	10
#define PERL_constant_ISYES	11

#ifndef NVTYPE
typedef double NV; /* 5.6 and later define NVTYPE, and typedef NV to it.  */
#endif
#ifndef aTHX_
#define aTHX_ /* 5.6 or later define this for threading support.  */
#endif
#ifndef pTHX_
#define pTHX_ /* 5.6 or later define this for threading support.  */
#endif
static int
func_ord_by_type_1 (pTHX_ const char *name, IV *iv_return) {
  /* When generated this function returned values for the list of names given
     here.  However, subsequent manual editing may have added or removed some.
     G L p */
  /* Offset 0 gives the best switch position.  */
  switch (name[0]) {
  case 'G':
    {
      *iv_return = 18;
      return PERL_constant_ISIV;
    }
    break;
  case 'L':
    {
      *iv_return = 11;
      return PERL_constant_ISIV;
    }
    break;
  case 'p':
    {
      *iv_return = 0;
      return PERL_constant_ISIV;
    }
    break;
  }
  return PERL_constant_NOTFOUND;
}
static int
func_ord_by_type_2 (pTHX_ const char *name, IV *iv_return) {
  /* When generated this function returned values for the list of names given
     here.  However, subsequent manual editing may have added or removed some.
     GG GL Gp LG lG ls vS */
  /* Offset 1 gives the best switch position.  */
  switch (name[1]) {
  case 'G':
    if (name[0] == 'G') {
      *iv_return = 2;
      return PERL_constant_ISIV;
    }
    if (name[0] == 'L') {
      *iv_return = 24;
      return PERL_constant_ISIV;
    }
    if (name[0] == 'l') {
      *iv_return = 10;
      return PERL_constant_ISIV;
    }
    break;
  case 'L':
    if (name[0] == 'G') {
      *iv_return = 23;
      return PERL_constant_ISIV;
    }
    break;
  case 'S':
    if (name[0] == 'v') {
      *iv_return = 85;
      return PERL_constant_ISIV;
    }
    break;
  case 'p':
    if (name[0] == 'G') {
      *iv_return = 1;
      return PERL_constant_ISIV;
    }
    break;
  case 's':
    if (name[0] == 'l') {
      *iv_return = 16;
      return PERL_constant_ISIV;
    }
    break;
  }
  return PERL_constant_NOTFOUND;
}
static int
func_ord_by_type_3 (pTHX_ const char *name, IV *iv_return) {
  /* When generated this function returned values for the list of names given
     here.  However, subsequent manual editing may have added or removed some.
     GDn GGG GGL GGp GVE GVI GnG GnP lGG vLL vLs */
  /* Offset 2 gives the best switch position.  */
  switch (name[2]) {
  case 'E':
    if (name[0] == 'G' && name[1] == 'V') {
      *iv_return = 22;
      return PERL_constant_ISIV;
    }
    break;
  case 'G':
    if (name[0] == 'G' && name[1] == 'G') {
      *iv_return = 3;
      return PERL_constant_ISIV;
    }
    if (name[0] == 'G' && name[1] == 'n') {
      *iv_return = 26;
      return PERL_constant_ISIV;
    }
    if (name[0] == 'l' && name[1] == 'G') {
      *iv_return = 20;
      return PERL_constant_ISIV;
    }
    break;
  case 'I':
    if (name[0] == 'G' && name[1] == 'V') {
      *iv_return = 22;
      return PERL_constant_ISIV;
    }
    break;
  case 'L':
    if (name[0] == 'G' && name[1] == 'G') {
      *iv_return = 32;
      return PERL_constant_ISIV;
    }
    if (name[0] == 'v' && name[1] == 'L') {
      *iv_return = 19;
      return PERL_constant_ISIV;
    }
    break;
  case 'P':
    if (name[0] == 'G' && name[1] == 'n') {
      *iv_return = 12;
      return PERL_constant_ISIV;
    }
    break;
  case 'n':
    if (name[0] == 'G' && name[1] == 'D') {
      *iv_return = 14;
      return PERL_constant_ISIV;
    }
    break;
  case 'p':
    if (name[0] == 'G' && name[1] == 'G') {
      *iv_return = 29;
      return PERL_constant_ISIV;
    }
    break;
  case 's':
    if (name[0] == 'v' && name[1] == 'L') {
      *iv_return = 57;
      return PERL_constant_ISIV;
    }
    break;
  }
  return PERL_constant_NOTFOUND;
}
static int
func_ord_by_type_4 (pTHX_ const char *name, IV *iv_return) {
  /* When generated this function returned values for the list of names given
     here.  However, subsequent manual editing may have added or removed some.
     GGGG lGGG vGVE vGVI vLGG vLLL */
  /* Offset 3 gives the best switch position.  */
  switch (name[3]) {
  case 'E':
    if (memEQ(name, "vGV", 3)) {
    /*                  E     */
      *iv_return = 84;
      return PERL_constant_ISIV;
    }
    break;
  case 'G':
    if (memEQ(name, "GGG", 3)) {
    /*                  G     */
      *iv_return = 4;
      return PERL_constant_ISIV;
    }
    if (memEQ(name, "lGG", 3)) {
    /*                  G     */
      *iv_return = 30;
      return PERL_constant_ISIV;
    }
    if (memEQ(name, "vLG", 3)) {
    /*                  G     */
      *iv_return = 35;
      return PERL_constant_ISIV;
    }
    break;
  case 'I':
    if (memEQ(name, "vGV", 3)) {
    /*                  I     */
      *iv_return = 84;
      return PERL_constant_ISIV;
    }
    break;
  case 'L':
    if (memEQ(name, "vLL", 3)) {
    /*                  L     */
      *iv_return = 34;
      return PERL_constant_ISIV;
    }
    break;
  }
  return PERL_constant_NOTFOUND;
}
static int
func_ord_by_type_5 (pTHX_ const char *name, IV *iv_return) {
  /* When generated this function returned values for the list of names given
     here.  However, subsequent manual editing may have added or removed some.
     GDVDE GDVDI "V=GEp" "V=GIp" */
  /* Offset 3 gives the best switch position.  */
  switch (name[3]) {
  case 'D':
    if (memEQ(name, "GDVDE", 5)) {
    /*                  ^       */
      *iv_return = 28;
      return PERL_constant_ISIV;
    }
    if (memEQ(name, "GDVDI", 5)) {
    /*                  ^       */
      *iv_return = 28;
      return PERL_constant_ISIV;
    }
    break;
  case 'E':
    if (memEQ(name, "V=GEp", 5)) {
    /*                  ^       */
      *iv_return = 27;
      return PERL_constant_ISIV;
    }
    break;
  case 'I':
    if (memEQ(name, "V=GIp", 5)) {
    /*                  ^       */
      *iv_return = 27;
      return PERL_constant_ISIV;
    }
    break;
  }
  return PERL_constant_NOTFOUND;
}
static int
func_ord_by_type_6 (pTHX_ const char *name, IV *iv_return) {
  /* When generated this function returned values for the list of names given
     here.  However, subsequent manual editing may have added or removed some.
     "V=GGEp" "V=GGIp" vLGGGG "vV=GGE" "vV=GGI" */
  /* Offset 5 gives the best switch position.  */
  switch (name[5]) {
  case 'E':
    if (memEQ(name, "vV=GG", 5)) {
    /*                    E     */
      *iv_return = 83;
      return PERL_constant_ISIV;
    }
    break;
  case 'G':
    if (memEQ(name, "vLGGG", 5)) {
    /*                    G     */
      *iv_return = 59;
      return PERL_constant_ISIV;
    }
    break;
  case 'I':
    if (memEQ(name, "vV=GG", 5)) {
    /*                    I     */
      *iv_return = 83;
      return PERL_constant_ISIV;
    }
    break;
  case 'p':
    if (memEQ(name, "V=GGE", 5)) {
    /*                    p     */
      *iv_return = 37;
      return PERL_constant_ISIV;
    }
    if (memEQ(name, "V=GGI", 5)) {
    /*                    p     */
      *iv_return = 37;
      return PERL_constant_ISIV;
    }
    break;
  }
  return PERL_constant_NOTFOUND;
}
static int
func_ord_by_type_7 (pTHX_ const char *name, IV *iv_return) {
  /* When generated this function returned values for the list of names given
     here.  However, subsequent manual editing may have added or removed some.
     "GDGDGD&" "GGD0,L," "LGD0,L," "V=GGEDG" "V=GGIDG" "vV=GGGE" "vV=GGGI" */
  /* Offset 6 gives the best switch position.  */
  switch (name[6]) {
  case '&':
    if (memEQ(name, "GDGDGD", 6)) {
    /*                     &     */
      *iv_return = 31;
      return PERL_constant_ISIV;
    }
    break;
  case ',':
    if (memEQ(name, "GGD0,L", 6)) {
    /*                     ,     */
      *iv_return = 25;
      return PERL_constant_ISIV;
    }
    if (memEQ(name, "LGD0,L", 6)) {
    /*                     ,     */
      *iv_return = 45;
      return PERL_constant_ISIV;
    }
    break;
  case 'E':
    if (memEQ(name, "vV=GGG", 6)) {
    /*                     E     */
      *iv_return = 86;
      return PERL_constant_ISIV;
    }
    break;
  case 'G':
    if (memEQ(name, "V=GGED", 6)) {
    /*                     G     */
      *iv_return = 47;
      return PERL_constant_ISIV;
    }
    if (memEQ(name, "V=GGID", 6)) {
    /*                     G     */
      *iv_return = 47;
      return PERL_constant_ISIV;
    }
    break;
  case 'I':
    if (memEQ(name, "vV=GGG", 6)) {
    /*                     I     */
      *iv_return = 86;
      return PERL_constant_ISIV;
    }
    break;
  }
  return PERL_constant_NOTFOUND;
}
static int
func_ord_by_type_8 (pTHX_ const char *name, IV *iv_return) {
  /* When generated this function returned values for the list of names given
     here.  However, subsequent manual editing may have added or removed some.
     "GD0,L,DG" GGDVDVDE GGDVDVDI */
  /* Offset 7 gives the best switch position.  */
  switch (name[7]) {
  case 'E':
    if (memEQ(name, "GGDVDVD", 7)) {
    /*                      E     */
      *iv_return = 49;
      return PERL_constant_ISIV;
    }
    break;
  case 'G':
    if (memEQ(name, "GD0,L,D", 7)) {
    /*                      G     */
      *iv_return = 13;
      return PERL_constant_ISIV;
    }
    break;
  case 'I':
    if (memEQ(name, "GGDVDVD", 7)) {
    /*                      I     */
      *iv_return = 49;
      return PERL_constant_ISIV;
    }
    break;
  }
  return PERL_constant_NOTFOUND;
}
static int
func_ord_by_type_17 (pTHX_ const char *name, IV *iv_return) {
  /* When generated this function returned values for the list of names given
     here.  However, subsequent manual editing may have added or removed some.
     "GD0,G,D0,G,D0,L,p" "LV=GGEpD0,L,D0,L," "LV=GGIpD0,L,D0,L," */
  /* Offset 5 gives the best switch position.  */
  switch (name[5]) {
  case ',':
    if (memEQ(name, "GD0,G,D0,G,D0,L,p", 17)) {
    /*                    ^                  */
      *iv_return = 62;
      return PERL_constant_ISIV;
    }
    break;
  case 'E':
    if (memEQ(name, "LV=GGEpD0,L,D0,L,", 17)) {
    /*                    ^                  */
      *iv_return = 73;
      return PERL_constant_ISIV;
    }
    break;
  case 'I':
    if (memEQ(name, "LV=GGIpD0,L,D0,L,", 17)) {
    /*                    ^                  */
      *iv_return = 73;
      return PERL_constant_ISIV;
    }
    break;
  }
  return PERL_constant_NOTFOUND;
}
static int
func_ord_by_type (pTHX_ const char *name, STRLEN len, IV *iv_return) {
  /* Initially switch on the length of the name.  */
  /* When generated this function returned values for the list of names given
     in this section of perl code.  Rather than manually editing these functions
     to add or remove constants, which would result in this comment and section
     of code becoming inaccurate, we recommend that you edit this section of
     code, and use it to regenerate a new set of constant functions which you
     then use to replace the originals.

     Regenerate these constant functions by feeding this entire source file to
     perl -x

#!i:/emx.add/BIN/perl.exe -w
use ExtUtils::Constant qw (constant_types C_constant XS_constant);

my $types = {map {($_, 1)} qw(IV)};
my @names = (qw(),
            {name=>"G", type=>"IV", macro=>"1", value=>"18"},
            {name=>"GD0,G,D0,G,D0,L,p", type=>"IV", macro=>"1", value=>"62"},
            {name=>"GD0,L,D0,G,", type=>"IV", macro=>"1", value=>"13"},
            {name=>"GD0,L,DG", type=>"IV", macro=>"1", value=>"13"},
            {name=>"GD0,L,DGp", type=>"IV", macro=>"1", value=>"96"},
            {name=>"GDGDGD&", type=>"IV", macro=>"1", value=>"31"},
            {name=>"GDGDGD0,L,p", type=>"IV", macro=>"1", value=>"62"},
            {name=>"GDVDE", type=>"IV", macro=>"1", value=>"28"},
            {name=>"GDVDI", type=>"IV", macro=>"1", value=>"28"},
            {name=>"GDn", type=>"IV", macro=>"1", value=>"14"},
            {name=>"GG", type=>"IV", macro=>"1", value=>"2"},
            {name=>"GGD0,L,", type=>"IV", macro=>"1", value=>"25"},
            {name=>"GGDVDVDE", type=>"IV", macro=>"1", value=>"49"},
            {name=>"GGDVDVDI", type=>"IV", macro=>"1", value=>"49"},
            {name=>"GGG", type=>"IV", macro=>"1", value=>"3"},
            {name=>"GGGD0,L,p", type=>"IV", macro=>"1", value=>"33"},
            {name=>"GGGG", type=>"IV", macro=>"1", value=>"4"},
            {name=>"GGL", type=>"IV", macro=>"1", value=>"32"},
            {name=>"GGp", type=>"IV", macro=>"1", value=>"29"},
            {name=>"GL", type=>"IV", macro=>"1", value=>"23"},
            {name=>"GVE", type=>"IV", macro=>"1", value=>"22"},
            {name=>"GVI", type=>"IV", macro=>"1", value=>"22"},
            {name=>"GnG", type=>"IV", macro=>"1", value=>"26"},
            {name=>"GnP", type=>"IV", macro=>"1", value=>"12"},
            {name=>"Gp", type=>"IV", macro=>"1", value=>"1"},
            {name=>"L", type=>"IV", macro=>"1", value=>"11"},
            {name=>"LG", type=>"IV", macro=>"1", value=>"24"},
            {name=>"LGD0,L,", type=>"IV", macro=>"1", value=>"45"},
            {name=>"LV=GGEpD0,L,D0,L,", type=>"IV", macro=>"1", value=>"73"},
            {name=>"LV=GGIpD0,L,D0,L,", type=>"IV", macro=>"1", value=>"73"},
            {name=>"V=GEp", type=>"IV", macro=>"1", value=>"27"},
            {name=>"V=GGEDG", type=>"IV", macro=>"1", value=>"47"},
            {name=>"V=GGEp", type=>"IV", macro=>"1", value=>"37"},
            {name=>"V=GGIDG", type=>"IV", macro=>"1", value=>"47"},
            {name=>"V=GGIp", type=>"IV", macro=>"1", value=>"37"},
            {name=>"V=GIp", type=>"IV", macro=>"1", value=>"27"},
            {name=>"lG", type=>"IV", macro=>"1", value=>"10"},
            {name=>"lGG", type=>"IV", macro=>"1", value=>"20"},
            {name=>"lGGG", type=>"IV", macro=>"1", value=>"30"},
            {name=>"ls", type=>"IV", macro=>"1", value=>"16"},
            {name=>"p", type=>"IV", macro=>"1", value=>"0"},
            {name=>"vGVE", type=>"IV", macro=>"1", value=>"84"},
            {name=>"vGVI", type=>"IV", macro=>"1", value=>"84"},
            {name=>"vLGG", type=>"IV", macro=>"1", value=>"35"},
            {name=>"vLGGGG", type=>"IV", macro=>"1", value=>"59"},
            {name=>"vLL", type=>"IV", macro=>"1", value=>"19"},
            {name=>"vLLL", type=>"IV", macro=>"1", value=>"34"},
            {name=>"vLs", type=>"IV", macro=>"1", value=>"57"},
            {name=>"vS", type=>"IV", macro=>"1", value=>"85"},
            {name=>"vV=GED0,L,", type=>"IV", macro=>"1", value=>"87"},
            {name=>"vV=GGE", type=>"IV", macro=>"1", value=>"83"},
            {name=>"vV=GGGE", type=>"IV", macro=>"1", value=>"86"},
            {name=>"vV=GGGI", type=>"IV", macro=>"1", value=>"86"},
            {name=>"vV=GGI", type=>"IV", macro=>"1", value=>"83"},
            {name=>"vV=GID0,L,", type=>"IV", macro=>"1", value=>"87"});

print constant_types(); # macro defs
foreach (C_constant ("Math::Pari::func_type", 'func_ord_by_type', 'IV', $types, undef, 3, @names) ) {
    print $_, "\n"; # C constant subs
}
print "#### XS Section:\n";
print XS_constant ("Math::Pari::func_type", $types);
__END__
   */

  switch (len) {
  case 1:
    return func_ord_by_type_1 (aTHX_ name, iv_return);
    break;
  case 2:
    return func_ord_by_type_2 (aTHX_ name, iv_return);
    break;
  case 3:
    return func_ord_by_type_3 (aTHX_ name, iv_return);
    break;
  case 4:
    return func_ord_by_type_4 (aTHX_ name, iv_return);
    break;
  case 5:
    return func_ord_by_type_5 (aTHX_ name, iv_return);
    break;
  case 6:
    return func_ord_by_type_6 (aTHX_ name, iv_return);
    break;
  case 7:
    return func_ord_by_type_7 (aTHX_ name, iv_return);
    break;
  case 8:
    return func_ord_by_type_8 (aTHX_ name, iv_return);
    break;
  case 9:
    /* Names all of length 9.  */
    /* "GD0,L,DGp" "GGGD0,L,p" */
    /* Offset 1 gives the best switch position.  */
    switch (name[1]) {
    case 'D':
      if (memEQ(name, "GD0,L,DGp", 9)) {
      /*                ^             */
        *iv_return = 96;
        return PERL_constant_ISIV;
      }
      break;
    case 'G':
      if (memEQ(name, "GGGD0,L,p", 9)) {
      /*                ^             */
        *iv_return = 33;
        return PERL_constant_ISIV;
      }
      break;
    }
    break;
  case 10:
    /* Names all of length 10.  */
    /* "vV=GED0,L," "vV=GID0,L," */
    /* Offset 4 gives the best switch position.  */
    switch (name[4]) {
    case 'E':
      if (memEQ(name, "vV=GED0,L,", 10)) {
      /*                   ^            */
        *iv_return = 87;
        return PERL_constant_ISIV;
      }
      break;
    case 'I':
      if (memEQ(name, "vV=GID0,L,", 10)) {
      /*                   ^            */
        *iv_return = 87;
        return PERL_constant_ISIV;
      }
      break;
    }
    break;
  case 11:
    /* Names all of length 11.  */
    /* "GD0,L,D0,G," "GDGDGD0,L,p" */
    /* Offset 7 gives the best switch position.  */
    switch (name[7]) {
    case ',':
      if (memEQ(name, "GDGDGD0,L,p", 11)) {
      /*                      ^          */
        *iv_return = 62;
        return PERL_constant_ISIV;
      }
      break;
    case '0':
      if (memEQ(name, "GD0,L,D0,G,", 11)) {
      /*                      ^          */
        *iv_return = 13;
        return PERL_constant_ISIV;
      }
      break;
    }
    break;
  case 17:
    return func_ord_by_type_17 (aTHX_ name, iv_return);
    break;
  }
  return PERL_constant_NOTFOUND;
}

