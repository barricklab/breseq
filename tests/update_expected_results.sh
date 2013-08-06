#!/usr/bin/env bash

for i in $(ls); do
  expected=$i/expected.gd
  new_expected=$i/output/evidence/annotated.gd
  if [ -e $expected ] && [ -e $new_expected ]; then
    echo "Updating $i"
    hg rm $expected
    cp $new_expected $expected
    hg add $expected
  fi
done
