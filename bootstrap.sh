#!/bin/sh -e

## make directories
if [ ! -d aux_build ] ; then
mkdir aux_build
fi

if [ ! -d aux_build/m4 ] ; then
mkdir aux_build/m4
fi


test -n "$srcdir" || srcdir=`dirname "$0"`
test -n "$srcdir" || srcdir=.
autoreconf --force --install --verbose "$srcdir"
#test -n "$NOCONFIGURE" || "$srcdir/configure" "$@"

