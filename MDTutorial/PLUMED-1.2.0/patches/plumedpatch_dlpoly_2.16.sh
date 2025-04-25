#!/bin/bash
# PATCH SCRIPT FOR DLPOLY
#

# this script needs to be launched from the root directory of the host code
destination="$PWD"

# definitions specific to this code
CODE="dlpoly_2.16"
LINKED_FILES="$plumedir/common_files/*.h $plumedir/common_files/*.c"
WHERE_LINKS="srcf90/Plumed/"

function to_do_before_patch () {
  mkdir srcf90/Plumed/
}

function to_do_after_patch () {
  {
    echo -n "PLUMED_OBJECTS="
    for file in $plumedir/common_files/*.c
      do f=${file##*/}
      echo " \\"
      echo -n "	Plumed/${f%.c}.o"
    done
    echo
  } > ./srcf90/Plumed/plumed.inc
}

function to_do_before_revert () {
  echo > /dev/null
}

function to_do_after_revert () {
  rm -fr srcf90/Plumed/
}


#########

NAME="$0"
source $plumedir/patches/patch_tool.sh

