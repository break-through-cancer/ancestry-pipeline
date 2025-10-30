#!/bin/bash -ue
if [ ! -f genetic_map_hg38_withX.txt.gz ]; then
    wget -q https://alkesgroup.broadinstitute.org/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz
fi
