regions_bed="$1"
bamFile="$2"

echo "## command= sh calculate_coverage.sh ${regions_bed} ${bamFile} "
echo -e "##chrom\tstart\tend\tcoverage" 

lines=$(cat $regions_bed | wc -l )
covavg=0

for region in `awk '{print $1":"$2"-"$3}' $regions_bed`
do
	tmp=$(echo $region | cut -d ':' -f 2 )
	length=$(echo $((tmp*-1)))
	average=$(samtools depth -r $region $bamFile | awk -v length2="$length" '{sum += $3} END {print sum/length2}' )

	echo -n $(echo $region | sed 's/-/:/g' | sed 's/:/\t/g') 
	echo -n -e "\t" 
	echo $average 
	
	covavg=$(awk "BEGIN {print $covavg + $average}")
done


result=$(awk "BEGIN {print $covavg / $lines}")

echo ""
echo "## average coverage for the given regions is: ${result}"

i=$(printf '%.0f' "$result")

if [ $i -le 1 ]; then
echo "## warning: average coverage is < 1, the sample does not seem to match with current reference genome, alignment to another reference genome is recommended" 
else
echo "## if average coverage is < 1, the sample does not seem to match with current reference genome, alignment to another reference genome is recommended" 
fi



