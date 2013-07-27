#!/bin/bash
#
# 
#
#  (c) Copyright 2013 Jean-Olivier Irisson
#      GNU General Public License v3
#
#--------------------------------------------------------------------------


sourcedir="/Volumes/ISIIShydro"
bckpdir="/Users/jiho/Work/projects/PUF/VISUFRONT/data/"
mkdir -p $bckpdir
rsync -avz "$sourcedir" "$bckpdir"

echo "Done"

# wait for the user to press a key
read

exit 0
