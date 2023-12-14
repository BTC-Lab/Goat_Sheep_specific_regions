################################### Sheep ###################################

## Calculate coverage for 1000 regions
for region in `cat sheep_specific_regions.bed | head -1000` 
do 
	total=0
	echo -n $region >> coverage_matrix_sheep.tsv
	echo -n -e '\t' >> coverage_matrix_sheep.tsv
	
	tmp=$(echo $region | cut -d ':' -f 2 )

	length=$(echo $((tmp*-1)))
	
	for bedfile in `ls samples_bed_files`
	do 
		sample=$(echo $bedfile | sed -e 's/.bed//g')

		average=$(samtools depth -r $region /path/to/sample/alignment/folder/${sample}/${sample}_deduped.bam | awk -v length2="$length" '{sum += $3} END {print sum/length2}' )

		total=$(awk "BEGIN {print $total + $average}")
	done

	result=$(awk "BEGIN {print $total / 15}")

	echo "$result" >> cov_matrix_sheep.tsv
done

## further step to decide the top 10 regions to view in IGV
cat cov_matrix_sheep.tsv | sort -k2 -n | head -10 | cut -f 1

################################### Goat ###################################

## Calculate coverage for 1000 regions
for region in `cat goat_specific_regions.bed | head -1000` 
do 
	total=0
	echo -n $region >> coverage_matrix_goat.tsv
	echo -n -e '\t' >> coverage_matrix_goat.tsv
	
	tmp=$(echo $region | cut -d ':' -f 2 )

	length=$(echo $((tmp*-1)))
	
	for bedfile in `ls samples_bed_files`
	do 
		sample=$(echo $bedfile | sed -e 's/.bed//g')

		average=$(samtools depth -r $region /path/to/sample/alignment/folder/${sample}/${sample}_deduped.bam | awk -v length2="$length" '{sum += $3} END {print sum/length2}' )

		total=$(awk "BEGIN {print $total + $average}")
	done

	result=$(awk "BEGIN {print $total / 15}")

	echo "$result" >> cov_matrix_goat.tsv
done

## further step to decide the top 10 regions to view in IGV
cat cov_matrix_goat.tsv | sort -k2 -n | head -10 | cut -f 1
