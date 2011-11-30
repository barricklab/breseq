#!/bin/bash

_usage_download_sequences() {
cat << EOF
    DOWNLOAD -u <user name> -p <password> -o <output dir> < file1.gd file2.gd file3.gd ... >
EOF
}

#-----------------------------------------------------------------------------
#   This script is used to download the necessary reference and read sequences 
# given a genome diff with the proper #=REFSEQ and #=READSEQ tags in the header. 
# The script then parses a $KEY:$VALUE from the header to determine where that 
# sequence should be downloaded from.
# 
# Usage: ./gd_download_sequence.sh -u <user name> -p <password> -o <output dir> < file1.gd file2.gd file3.gd ... >"
#   Example: ./gd_download_sequence.sh -o 02_Downloads $(ls 01_Data/*.gd)"
#   Example: ./gd_download_sequence.sh -o downloads -u john -p "12)34" RJW1129.gd
# 
# *** You may need to quote the pass word if non common character are used.
# 
#   The user name and password options are used for downloading sequences from 
# ftp://backup.barricklab.org, you may store your user name and password in 
# the respective $DCAMP_BACKUP_USER and $DCAMP_BACKUP_PASSWORD environmental 
# variables. If no output directory is set then the default will be the 
# current working directory.
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# GLOBAL VARIABLES  
#-----------------------------------------------------------------------------
declare VERBOSE=false
declare OUTPUT_DIR=$(pwd)
declare USER_NAME=$DCAMP_BACKUP_USER
declare PASS_WORD=$DCAMP_BACKUP_PASSWORD

KEY_GENBANK="Genbank:"
KEY_SRA="SRA:"
KEY_BARRICK_PUBLIC="BarrickLab-Public:"
KEY_BARRICK_PRIVATE="BarrickLab-Private:"
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#::parse_genome_diff()
#
#   Reads through a given genome diff file to locate lines that begin 
# with #=REFSEQ or #=READSEQ and obtain the KEY:VALUE string within.
# 
# An example genome diff header is as follows:
#       #=GENOME_DIFF 1.0
#       #=AUTHOR    Barrick JE 
#       #=REFSEQ    BarrickLab-Public:release/reference/REL606.5.gbk
#       #=READSEQ   SRA:SRR066595
# 
# From reading in the above lines we would parse the following:
#       BarrickLab-Public:release/reference/REL606.5.gbk
#       SRA:SRR066595
#-----------------------------------------------------------------------------
parse_genome_diff() {
  local gd_file_name=$1
  #     Awk is given the two tags as the conditions to match, awk defaults 
  # to any whitespace delimeter such that column 2 is the KEY:VALUE we are 
  # looking for.
  for i in $(awk '/#=REFSEQ/{ print $2;} /#=READSEQ/{ print $2}' $gd_file_name); do
    if $VERBOSE ;then
      echo "parse_genome_diff(): $i in file: $gd_file_name"
    fi
    _parse_seq_key_value $i
    echo ""
    echo ""
  done
}
#-----------------------------------------------------------------------------
#::_parse_seq_key_value()
#
#   Receives KEY:VALUE pairs and attempts to find a matching key 
# from the know KEY_*** constants declared above. After which, the VALUE is 
# sent to the appropriate download function where the output file name and url 
# are built and the sequence is ultimately downloaded.
#
# ***Will require additions if new KEY:VALUES are introduced.
#-----------------------------------------------------------------------------
_parse_seq_key_value() {
  key_value=$1
  if $VERBOSE
  then
    echo "_parse_seq_key_value(): $key_value"
  fi

  #     If the KEY matches any of our known keys strip the KEY:VALUE of KEY:
  # and send to the proper function.
  if [[ "$key_value" =~ "$KEY_GENBANK" ]]; then
    _download_key_genbank ${i##$KEY_GENBANK}
  elif [[ "$key_value" =~ "$KEY_SRA" ]]; then
    _download_key_sra ${i##$KEY_SRA}
  elif [[ "$key_value" =~ "$KEY_BARRICK_PUBLIC" ]]; then
    _download_barrick_public ${i##$KEY_BARRICK_PUBLIC}
  elif [[ "$key_value" =~ "$KEY_BARRICK_PRIVATE" ]]; then
    _download_barrick_private ${i##$KEY_BARRICK_PRIVATE}
  fi
}
#-----------------------------------------------------------------------------
#::_download_key_***()
# 
#   Each is resposible for building the file name and url for the respective KEY
# and then downloading it. Currently only BarrickLab-Private needs a 
# user:password to loging to an ftp server.
#
# ***Will require additions if new KEY:VALUES are introduced.
#-----------------------------------------------------------------------------
_download_key_genbank() {
  local ref_seq=$1
  local new_file_path=$OUTPUT_DIR/$1.gbk
 
  if   _file_exists $new_file_path; then
    return 0
  fi

  local url="http://www.ncbi.nlm.nih.gov/sviewer/?db=nuccore&val=$ref_seq&report=gbwithparts&retmode=text"
  if $VERBOSE ; then
    echo Genbank URL: $url
  else
    $(wget -O $new_file_path $url)
  fi
}

_download_key_sra() {
  local read_seq=$1
  local new_file_path=$OUTPUT_DIR/$1.fastq.gz

  if  _file_exists $new_file_path; then
    return 0
  fi

  #String voodoo to find correct directories.
  local first=${read_seq:0:6}
  local second=${read_seq:0:9}
  local url="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$first/$second/$read_seq.fastq.gz"
  if $VERBOSE ; then
    echo SRA URL: $url
  else
    $(wget -O $new_file_path $url)
  fi
}

_download_barrick_public() {
  local url_file_path=$1
  local new_file_path=$OUTPUT_DIR/$(basename $1)
  
  if   _file_exists $new_file_path; then
    return 0
  fi                                        
  
  local url="http://barricklab.org/$url_file_path"
  if $VERBOSE ; then
    echo BarrickLab-Public URL:$url
  else
    $(wget -O $new_file_path $url)
  fi
}

_download_barrick_private() {
  local url_file_path=$1
  local new_file_path=$OUTPUT_DIR/$(basename $1)

  if _file_exists $new_file_path; then
    return $FILE_ALREADY_EXISTS
  fi

  if [ -z "$USER_NAME" ] || [ -z "$PASS_WORD" ]; then
    echo "WARNING! Command line options -u <user name> and/or -p<pass_word> not set."
    echo "Access to backup.barricklab.org denied."
    return 
  fi

  local url="ftp://$USER_NAME:$PASS_WORD@backup.barricklab.org/$url_file_path"
  if $VERBOSE ; then
    echo BarrickLab-Private URL: $url
  else
    $(wget -O $new_file_path $url)
  fi
}

#-----------------------------------------------------------------------------
#::_file_exists()
#   Check that the file already hasnt been downloaded. Read files are commonly
# compressed in *.gz format so this function also checks for uncompressed
# files with the same name.
#
# ***May want to add checks for *.tar.gz
#-----------------------------------------------------------------------------
_file_exists() {
  local new_file_path=$1
  local unzipped_file_path=${new_file_path%%.gz}

  if $VERBOSE ; then
    echo "_file_exists()::new_file_path= $new_file_path"
    echo "_file_exists()::unzipped_file_path= $unzipped_file_path"
  fi
  
  if [[ -e "$new_file_path" ]] && [[ -s "$new_file_path" ]]; then
    echo WARNING! The following file aready exists: $new_file_path 
    return 0
  fi

  if [[ -e "$unzipped_file_path" ]] && [[ -s "$unzipped_file_path" ]]; then
    echo WARNING! The following file aready exists: $unzipped_file_path 
    return 0
  fi

  return 1
}


_do_download_sequences() {
  local OUTPUT_DIR=$(pwd)
  local USER_NAME=$DCAMP_BACKUP_USER
  local PASS_WORD=$DCAMP_BACKUP_PASSWORD
  local -a args=("$@")

  #Check for user input
  while getopts "ho:u:p:" option; do
    case $option in 
      o)
        OUTPUT_DIR=$OPTARG
        ;;
      u)
        USER_NAME=$OPTARG
        ;;
      p)
        PASS_WORD=$OPTARG
        ;;
      h)
        exit -1
        ;;
    esac
  done
  shift $(($OPTIND - 1))

  if $VERBOSE
  then
    echo OUTPUT_DIR:$OUTPUT_DIR
    echo USER_NAME:$USER_NAME
    echo PASS_WORD:$PASS_WORD
    echo OPTIND:$OPTIND
    echo "size(ARGS):"${#args[*]}
    echo ""
  fi

  #Create directory
  if [[ "$OUTPUT_DIR" != "$(pwd)" ]]; then
    mkdir -p $OUTPUT_DIR
  fi

  #The other arguments are genome_diffs
  for i in $*; do  
    parse_genome_diff ${args[$i]}
  done                                                  
  

  exit 0
}

_do_download_sequences $@
