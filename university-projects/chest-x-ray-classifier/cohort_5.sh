#!/bin/bash

# iterate over the folder /dmc/ml_storage/machine_learning/Final_ML_CBS/data/age_csvs and read the file contents
for file in /content/cohorts_5/*.csv
do
        # get the filename without the extension
        filename=$(basename "$file" .csv)
        mkdir /content/data/$filename
        cat $file | while read line
        do
        # copy the image with the filename to the folder with the filename
                echo $line
                cp /content/train/$line /content/data/$filename/
    done
done
