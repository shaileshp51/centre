#!/usr/bin/env sh

for d in $( ls -d example-? ); do 
  echo "removing all *.log *.info *.txt files from directory: $d"; 
  find $d -maxdepth 1 -type f \
    \( -name "?*.log" -o -name "?*.info" -o -name "?*.txt" \) -exec rm {} +; 
done

for d in $( ls -d example-?/inpfiles ); do 
  echo "removing all but *.cfg files from directory: $d"; 
  find $d -type f -not -name "*.cfg" -exec rm {} +; 
done

for d in $( ls -d example-?/outfiles-h ); do 
  echo "removing directory: $d"; rm -r $d; 
done

for d in $( ls -d example-?/outfiles-v ); do 
  echo "removing directory: $d"; rm -r $d; 
done
