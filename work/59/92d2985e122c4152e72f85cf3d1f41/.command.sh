#!/bin/bash -ue
if [ ! -f genetic_map_hg38_withX.txt.gz ]; then
    echo "Downloading genetic map with curl..."
    curl -s -L -o genetic_map_hg38_withX.txt.gz         https://alkesgroup.broadinstitute.org/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz
else
    echo "Genetic map already exists, skipping download."
fi
