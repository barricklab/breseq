#!/bin/bash
#   This script is meant to run on TACC's Ranger cluster and assumes that the 
# environmental variables $SCRATCH, $WORK and $HOME are set. If it is found 
# that additional files are needed then add their url to the $FILE_URLS array. 
# For any new files that aren't in the compressed formats of *.tar.gz or *.tgz 
# the function untar_files() will need to be updated. If the files to be 
# downloaded are already in an executable format and aren't available in the 
# immediate parent directory then install_files() will need to be updated.

#   PLEASE, update this script if changes are needed. Troubleshooting 
# installation is not a pleasant use of time.

#-----------------------------------------------------------------------------
# GLOBAL VARIABLES  
#-----------------------------------------------------------------------------
FILE_URLS=(
    http://ftpmirror.gnu.org/libtool/libtool-2.4.2.tar.gz
    http://ftp.gnu.org/gnu/m4/m4-1.4.16.tar.gz
    http://ftp.gnu.org/gnu/autoconf/autoconf-latest.tar.gz
    http://ftp.gnu.org/gnu/automake/automake-1.11.1.tar.gz
    http://downloads.ghostscript.com/public/ghostscript-9.04.tar.gz
    http://zlib.net/zlib-1.2.5.tar.gz
    ftp://ftp.sanger.ac.uk/pub4/resources/software/ssaha2/ssaha2_v2.5.5_x86_64.tgz
)
BUILD_DIR="$SCRATCH/build"
LOCAL_DIR="$HOME/local"
PROFILE_USER_FILE="$HOME/.profile_user"
BRESEQ_BUILD_DIR="$WORK/src/breseq"
#-----------------------------------------------------------------------------

write_profile_user_file() {
cat <<"EOF" > $PROFILE_USER_FILE
#Load any default modules. 
module load R 
module swap pgi gcc 
#Breseq ./configure settings and flags. 
LOCAL_DIR="$HOME/local" 
PATH="$LOCAL_DIR/bin:$PATH" 
LD_LIBRARY_PATH="$LOCAL_DIR/lib:$LD_LIBRARY_PATH" 
CPPFLAGS="$CPPFLAGS -I$LOCAL_DIR/include" 
CFLAGS="$CFLAGS -I$LOCAL_DIR/include -L$LOCAL_DIR/lib" 
CPPFLAGS="$CPPFLAGS -I$LOCAL_DIR/include -L$LOCAL_DIR/lib" 
LDFLAGS="$LDFLAGS -L$LOCAL_DIR/lib"
EOF
}

download_files() {
    for i in "${FILE_URLS[@]}"; do
      if ls $(basename $i)
        then
            echo $i already exists!
        else
            wget $i
        fi
    done
}

untar_files() {
    for i in $(ls *.tar.gz); do 
        tar xvzf $i
    done
    
    for i in $(ls *.tgz); do
        tar xvzf $i
    done
}

install_files() {
    for i in $(ls -d */); do 
        pushd $i #CD in
        # Has configure file and needs to be built
        if ls configure
        then
            ./configure --prefix $LOCAL_DIR
            make install
        else
        # Has executables to be copied to a bin directory
            if find . -perm +111 -type f
            then
                for i in $(find . -perm +111 -type f); do 
                    cp $i $LOCAL_DIR/bin
                done
            fi
        fi
        popd #CD out
    done
}

install_breseq() {
    pushd $BRESEQ_BUILD_DIR
    hg clone http://breseq.googlecode.com/hg/ breseq
    pushd breseq
    ./bootstrap.sh
    ./configure --prefix $LOCAL_DIR
    make install
    make test
}

main() {
    mkdir -p $BUILD_DIR
    mkdir -p $LOCAL_DIR
    mkdir -p $BRESEQ_BUILD_DIR

    pushd $BUILD_DIR
        download_files
        untar_files
        install_files
    popd
    rm -rf $BUILD_DIR
    
    write_profile_user_file
    source $PROFILE_USER_FILE

    install_breseq
}

main
