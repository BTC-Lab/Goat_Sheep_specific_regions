regions_bed="$1"
bamFile="$2"

for region in `awk '{print $1":"$2"-"$3}' $regions_bed`
do
	tmp=$(echo $region | cut -d ':' -f 2 )
	length=$(echo $((tmp*-1)))
	average=$(samtools depth -r $region $bamFile | awk -v length2="$length" '{sum += $3} END {print sum/length2}' )

	echo -n $(echo $region | sed 's/-/:/g' | sed 's/:/\t/g') 
	echo -n -e "\t" 
	echo $average 
done


