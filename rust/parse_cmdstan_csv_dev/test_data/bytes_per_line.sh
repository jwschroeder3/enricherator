#!/bin/bash
file=$1
while read line
do
    a=$( echo "$line" | wc -c) # echo -n to prevent counting new line
    echo $a
done <"$file"
