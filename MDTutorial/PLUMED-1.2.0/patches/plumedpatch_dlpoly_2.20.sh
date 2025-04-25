#!/bin/bash
# PATCH SCRIPT FOR DLPOLY
#

# this script needs to be launched from the root directory of the host code
destination="$PWD"

# definitions specific to this code
CODE="dlpoly_2.20"
LINKED_FILES="$plumedir/common_files/*.h $plumedir/common_files/*.c"
WHERE_LINKS="srcmod/Plumed/"

function to_do_before_patch () {
  mkdir srcmod/Plumed/
}

function to_do_after_patch () {
  # to be removed whem a pbc routine common to all the codes will be written
  awk '/subroutine images/{print;for(i=0;i<230;i++){getline;print;if($2=="subroutine")exit}}' \
       srcmod/utility_module.f >  srcmod/images.f
  awk '/subroutine invert/{print;for(i=0;i<230;i++){getline;print;if($2=="subroutine")exit}}' \
       srcmod/utility_module.f >> srcmod/images.f
 {
    echo -n "PLUMED_OBJECTS="
    for file in $plumedir/common_files/*.c
      do f=${file##*/}
      echo " \\"
      echo -n "	Plumed/${f%.c}.o"
    done
    echo
  } > ./srcmod/Plumed/plumed.inc

}

function to_do_before_revert () {
  rm srcmod/images.f
}

function to_do_after_revert () {
  rm -fr srcmod/Plumed/
}


#########

NAME="$0"
source $plumedir/patches/patch_tool.sh

