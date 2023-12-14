################################### Sheep ###################################
## get regions that has a coverage of < 10 for real_sheep 
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

################################### Goat ###################################
## get regions that has a coverage of < 10 for real_goat 
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
