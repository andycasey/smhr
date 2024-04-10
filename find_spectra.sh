#!/bin/bash

#targets=$(find /Volumes/ext/ -type f -iname "*J00401252+2729247_rv*" 2>/dev/null)
#IFS=',' read -ra line <<<

while read line
do
  IFS=',' read -ra data <<< "$line"
  star=${data[0]}
  run=${data[1]}
  if [ "$run" == "Hol20" ]; then
    echo $star
    file="$(find /Volumes/ext/ -type f -iname "*${star}_rv.fits*" 2>/dev/null | tail -n 1)"
    cp $file /Volumes/My\ Passport/Holmbeck2020_stars/
    file="$(find /Volumes/ext/ -type f -iname "*${star}_rv.txt*" 2>/dev/null | tail -n 1)"
    cp $file /Volumes/My\ Passport/Holmbeck2020_stars/
  fi
done < /Users/holmbeck/Research/RPA_all_classes_new.csv 
