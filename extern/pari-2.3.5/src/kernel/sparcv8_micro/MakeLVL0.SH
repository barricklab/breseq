cat >> $file << EOT
parilvl0.h: $kern0/asm0-common.h \$(L0MODS)
	$cfg/genkernel $kern0/asm0-common.h $kern0/asm0.h > parilvl0.h
EOT
