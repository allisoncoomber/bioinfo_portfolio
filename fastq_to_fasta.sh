#!/bin/bash

for file in path/*.fq*
do
	gunzip $file
	echo "$file" unzipped.
	no_ext=`echo $file | cut -f1 -d.`
	sed -n '1~4s/^@/>/p;2~4p' "$no_ext".fq > "$no_ext".fasta
	echo "$no_ext".fasta generated.
	gzip "$no_ext".fasta
	echo "$no_ext".fasta zipped.
	rm "$no_ext".fq
done
