#!/bin/bash
# PATCH SCRIPT FOR ACEMD 
#

# this script needs to be launched from the root directory of the host code
plumedir="/local/acellera/md_meta/"
destination="$PWD"
# definitions specific to this code
CODE="ACEMD"
LINKED_FILES="$plumedir/common_files/metadyn.h $plumedir/common_files/*.c"
WHERE_LINKS="$plumedir/ACEMD"


function to_do_before_patch () {
  echo > /dev/null
cp Makefile Makefile.preplumed
}

function to_do_after_patch () {
  echo > /dev/null 
}

function to_do_before_revert () {
  echo > /dev/null 
}

function to_do_after_revert () {
  echo > /dev/null
}

#########

NAME="$0"
source $plumedir/patches/patch_tool.sh

