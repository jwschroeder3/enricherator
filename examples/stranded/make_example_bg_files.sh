#!/usr/bin/bash

echo "genotype,sample,rep,strand,file,norm_factor" > example_info.csv
while IFS= read -r line; do
  IFS=',' read -r -a array <<< "$line"
  geno="${array[0]}"
  samp_type="${array[1]}"
  rep="${array[2]}"
  strand="${array[3]}"
  infile="${array[4]}"
  base=$(basename $infile)
  outfile="count_files/sample_${base}"
  echo "Grabbing region from $infile, placing into $outfile"
  grep CP006881.1 $infile \
    | grep -A 1400 -P "CP006881.1\t128000\t" \
    | cut -f 1-4 \
    > $outfile
  echo "Writing row to info file"
  echo "${geno},${samp_type},${rep},${strand},${outfile},1" >> example_info.csv
done < $1
