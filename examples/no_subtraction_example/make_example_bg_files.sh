#!/usr/bin/bash

echo "genotype,sample,rep,strand,file,input_group,control_group,extract_group" > example_info.csv
while IFS= read -r line; do
  IFS=',' read -r -a array <<< "$line"
  geno="${array[0]}"
  samp_type="${array[1]}"
  rep="${array[2]}"
  strand="${array[3]}"
  infile="${array[4]}"
  inp_group="${array[5]}"
  control_group="${array[6]}"
  ext_group="${array[7]}"
  base=$(basename $infile)
  outfile="count_files/sample_${base}"
  echo "Grabbing region from $infile, placing into $outfile"
  grep -A 2500 -B 2500 -P "2597750\t2597755" $infile \
    | cut -f 1-4 \
    > $outfile
  echo "Writing row to info file"
  echo "${geno},${samp_type},${rep},${strand},${outfile},${inp_group},${control_group},${ext_group}" >> example_info.csv
done < /corexfs/schroedj/TEST_IPODERATOR/nosub_IPODerator_info.csv
