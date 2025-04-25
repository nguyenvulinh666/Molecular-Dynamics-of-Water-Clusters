#!/bin/bash
# PATCH SCRIPT FOR GROMACS
#

# this script needs to be launched from the root directory of the host code
destination="$PWD"

# definitions specific to this code
CODE="qespresso-4.1.2"
LINKED_FILES="$plumedir/common_files/*.h $plumedir/common_files/*.c"
WHERE_LINKS="./clib/"

function to_do_before_patch () {
  echo > /dev/null
}

function to_do_after_patch () {
  {
    echo -n "PLUMED_OBJECTS="
    for file in $plumedir/common_files/*.c
      do f=${file##*/}
      echo " \\"
      echo -n "	${f%.c}.o"
    done
    echo
    echo -n "PLUMED_SRC="
    for file in $plumedir/common_files/*.c
      do f=${file##*/}
      echo " \\"
      echo -n "	${f%.c}.c"
    done
    echo
    echo
  } > ./clib/plumed.inc
  mv make.sys make.sys.original
  awk '{if($1=="DFLAGS")print $0,"-DPLUMED_QESPRESSO"; else print}' make.sys.original > make.sys
}

function to_do_before_revert () {
  rm ./clib/plumed.inc
}

function to_do_after_revert () {
  echo > /dev/null
  mv make.sys.original make.sys
}

#########

NAME="$0"
source $plumedir/patches/patch_tool.sh

