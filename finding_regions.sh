################################### Sheep ###################################
## get regions that has a coverage of < 10 for real_sheep samples 
bedtools genomecov -bga -ibam $bamfile | awk '$4<10' > ${sample}.bed

## merge bed files for all samples (training data = 15 samples)
cat ${sample}.bed | cut -f 1,2,3 >> cov.bed
sort -k 1,1 cov.bed > cov.sorted.bed
bedtools genomecov -bga -i cov.sorted.bed -g sheep.genome > cov.merged.bed
cat cov.merged.bed | awk '$4 > 0' > cov.merged.notzero.bed

## get regions that has a coverage of < 2 for goats aligned to sheep genome
bedtools genomecov -bga -ibam $bamfile | awk '$4<2' > ${sample}.bed

## merge bed files for all samples (training data = 15 samples)
cat ${sample}.bed | cut -f 1,2,3 >> cov.bed
sort -k 1,1 cov.bed > cov.sorted.bed
bedtools genomecov -bga -i cov.sorted.bed -g sheep.genome > cov.merged.bed
cat cov.merged.bed | awk '$4 == 15' > cov.merged15.bed

cat cov.merged.notzero.bed | cut -f 1,2,3 > real_sheep.bed
cat cov.merged15.bed | cut -f 1,2,3 > goat2sheep.bed

## get uniq regions from non sheep and get top regions based on length
bedtools intersect -a goat2sheep.bed -b real_sheep.bed -v | bedtools merge > sheep_uniq.bed
sort -k 1,1 sheep_uniq.bed > sheep_uniq.sorted.bed
bedtools genomecov -bga -i sheep_uniq.sorted.bed -g sheep.genome | awk '$4 >0' > sheep_uniq.merged.bed
cat sheep_uniq.merged.bed | awk '{print $1, $2, $3, $3 - $2}' | sed -e 's/ /\t/g' | sort -k4 -n -r | grep -v 'NW' | awk '{print $1":"$2"-"$3}' > sheep_specific_regions.bed

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
## After further manual checking, the top 10 regions are saved to "sheep_specific_regions_selected.bed" 

################################### Goat ###################################
## get regions that has a coverage of < 10 for real_goat samples
bedtools genomecov -bga -ibam $bamfile | awk '$4<10' > ${sample}.bed

## merge bed files for all samples (training data = 15 samples)
cat ${sample}.bed | cut -f 1,2,3 >> cov.bed
sort -k 1,1 cov.bed > cov.sorted.bed
bedtools genomecov -bga -i cov.sorted.bed -g goat.genome > cov.merged.bed
cat cov.merged.bed | awk '$4 > 0' > cov.merged.notzero.bed

## get regions that has a coverage of < 2 for sheep aligned to goat genome
bedtools genomecov -bga -ibam $bamfile | awk '$4<2' > ${sample}.bed

## merge bed files for all samples (training data = 15 samples)
cat ${sample}.bed | cut -f 1,2,3 >> cov.bed
sort -k 1,1 cov.bed > cov.sorted.bed
bedtools genomecov -bga -i cov.sorted.bed -g goat.genome > cov.merged.bed
cat cov.merged.bed | awk '$4 == 15' > cov.merged15.bed

cat cov.merged.notzero.bed | cut -f 1,2,3 > real_goat.bed
cat cov.merged15.bed | cut -f 1,2,3 > sheep2goat.bed

## get uniq regions from non goat and get top regions based on length
bedtools intersect -a sheep2goat.bed -b real_goat.bed -v | bedtools merge > goat_uniq.bed
sort -k 1,1 goat_uniq.bed > goat_uniq.sorted.bed
bedtools genomecov -bga -i goat_uniq.sorted.bed -g goat.genome | awk '$4 >0' > goat_uniq.merged.bed
cat goat_uniq.merged.bed | awk '{print $1, $2, $3, $3 - $2}' | sed -e 's/ /\t/g' | sort -k4 -n -r | grep -v 'NW' | awk '{print $1":"$2"-"$3}' > goat_specific_regions.bed

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
## After further manual checking, the top 10 regions are saved to "goat_specific_regions_selected.bed" 
