SourceDir="/data/class/ecoevo283/public/RAWDATA/RNAseq/RNAseq384plex_flowcell01"
DestDir="/pub/etom2/EE283/EE283/Week2/RNAseq/raw"
Files=$(find ${SourceDir}/ -type f -iname "*.fastq.gz" | grep "Project_plex")
for F in $Files
do
	ln -s $F $DestDir/$(basename $F)
done
