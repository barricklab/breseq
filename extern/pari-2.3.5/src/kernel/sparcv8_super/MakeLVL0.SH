cat >> $file << EOT
parilvl0.h: $src/kernel/sparcv8_micro/asm0-common.h \$(L0MODS)
	$cfg/genkernel $src/kernel/sparcv8_micro/asm0-common.h $kern0/asm0.h > parilvl0.h
EOT
