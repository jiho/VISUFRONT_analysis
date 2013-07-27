#!/bin/bash
#
# 
#
#  (c) Copyright 2013 Jean-Olivier Irisson
#      GNU General Public License v3
#
#--------------------------------------------------------------------------


workdir="/Users/jiho/Work/projects/PUF/VISUFRONT/data"
bckpdir="/Volumes/archive/VISUFRONT_data_backup/"
mkdir -p $bckpdir
rsync -avz "$workdir" "$bckpdir"

echo "Done"

# wait for the user to press a key
read

exit 0
