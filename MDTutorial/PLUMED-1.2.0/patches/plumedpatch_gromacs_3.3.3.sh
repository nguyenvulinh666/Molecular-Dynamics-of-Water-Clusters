#!/bin/bash
# PATCH SCRIPT FOR GROMACS
#

# this script needs to be launched from the root directory of the host code
destination="$PWD"

# definitions specific to this code
CODE="gromacs-3.3.3"
LINKED_FILES="$plumedir/common_files/*.h $plumedir/common_files/*.c"
WHERE_LINKS="./src/kernel/"

PO_FILES="$(
  for file in $plumedir/common_files/*.c
  do f=${file##*/}
  echo $WHERE_LINKS/.deps/${f%.c}.Po
  done
)"

function to_do_before_patch () {
  echo > /dev/null
}

function to_do_after_patch () {
  touch $PO_FILES
  {
    echo -n "PLUMED_OBJECTS="
    for file in $plumedir/common_files/*.c
      do f=${file##*/}
      echo " \\"
      echo -n "	${f%.c}.\$(OBJEXT)"
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
  } > ./src/kernel/plumed.inc
  for file in $plumedir/common_files/*.c
    do f=${file##*/}
    echo "include ./\$(DEPDIR)/${f%.c}.Po"
  done > ./src/kernel/plumed.Po.inc

}

function to_do_before_revert () {
  rm $PO_FILES
  rm ./src/kernel/plumed.inc
  rm ./src/kernel/plumed.Po.inc
}

function to_do_after_revert () {
  echo > /dev/null
}

#########

NAME="$0"
source $plumedir/patches/patch_tool.sh

