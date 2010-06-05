cat >> $file << EOT
L1OBJS=$kern1/tune.h $kern1/int.h $knone/level1.h
parilvl1.h: \$(L1OBJS) $src/headers/paritune.h
	cat \$(L1OBJS) > parilvl1.h
MP_C=$kern1/mp.c $knone/cmp.c $knone/gcdll.c $knone/ratlift.c\
  $knone/invmod.c $kern1/gcd.c $knone/mp_indep.c $knone/add.c
mp.c: \$(MP_C)
	cat \$(MP_C) > mp.c
mp\$(_O): .headers mp.c
	\$(CC) -c \$(CFLAGS) \$(DLCFLAGS) \$(CPPFLAGS) -o mp\$(_O) mp.c
EOT
