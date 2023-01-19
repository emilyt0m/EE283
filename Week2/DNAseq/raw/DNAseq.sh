SourceDir="/data/class/ecoevo283/public/RAWDATA/DNAseq"
DestDir="/pub/etom2/EE283/EE283/Week2/DNAseq/raw"
FILES=$(basename -a $SourceDir/*)
for f in $FILES
do
   echo "Processing $f file..."
   ln -s $SourceDir/$f $DestDir/$f
done


