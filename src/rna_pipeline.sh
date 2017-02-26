#!/bin/bash

# Copyright (C) 2016-2017 Bohdan Khomtchouk and Derek Van Booven
# This file is part of NGStoolkit.


# This file is a single-click automated RNA-seq pipeline starting with bz2 files through differential expression.
# USAGE : $ sh rna_pipeline.sh <pipeline_type> <genome> <paired/single> <config> <strandedness>

# There are 6 input parameters:
# 1.	Pipeline type { simple, complete }.  Simple is for individual flow cells where all people want are the stats for a given run, and complete is through from bz2 to differential expression.
# 2.	Genome type { hg19, hg38, mm10, rn6 }.  Specify the genome for the samples.
# 3.	Read type { paired, single }.  Paired/single ended read flag.
# 4.	Configuration file { <filename> }.  This file is a list of sampleIDs and the groups for differential expression.  First column is the sample ID and the second column is the group.
# 5.	Stranded information { u, s, r }.  Type of strandedness to take the counts from STAR output.  U = unstranded, s = stranded, r = reverse stranded as defined by STAR.

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Step 0 - Initialization
mkdir logs
mkdir logs/unzip
mkdir logs/align
mkdir logs/stats
mkdir logs/de

runflag=1

case $2 in
	hg19)
		reference="/hihg/ref/genomes/star_2.5.0a_hg19"
		gtf="/hihg/ref/gtf/gencode.v19.annotation.gtf"
		R4="/hihg/ref/gtf/forgene/gencode.forR.txt"
		;;
	hg38)
		reference="/hihg/ref/genomes/star_2.5.0a_hg38"
		gtf="/hihg/ref/gtf/gencode.v25.annotation.gtf"
		R4="/hihg/ref/gtf/forgene/gencode.v25.forR.txt"
		;;
	mm10)
		reference="/hihg/ref/genomes/star_2.5.0a_mm10"
		gtf="/hihg/ref/gtf/Mus_musculus.GRCm38.77.gtf"
		R4="/hihg/ref/gtf/forgene/mm10.forR.txt"
		;;
	rn6)
		reference="/hihg/ref/genomes/star_2.5.0a_rn6"
		gtf="/hihg/ref/gtf/Rattus_norvegicus.Rnor_6.0.87.gtf"
		R4="/hihg/ref/gtf/forgene/rn6.forR.txt"
		;;
esac

#Step 1 - Submit a job that will unzip and remove the final character in all bz2 files in the directory
# Get # of files that we are going to unzip for later
filenum=$(ls *bz2 | wc -l)
echo "number of files is $filenum"
i=0
for f in *bz2
do
	#get the job id that was returned from LSF
	#stringZ="Job is submitted to <aut> project. Job <8033051> is submitted to queue <general>"
	echo "bunzip2 $f" >> unzip$i.job
	f2=${f/.bz2/}
	echo "awk '{if(NR%2==0) print substr(\$0, 1, length(\$0)-1); else print \$0}' $f2 > $f2.tmp" >> unzip$i.job
	echo "mv $f2.tmp $f2" >> unzip$i.job
	VAR=$((bsub -o logs/unzip/unzip$i.out -e logs/unzip/unzip$i.err -q general -P hihgbioinf sh unzip$i.job; ) 2>&1)
	#cut the jobid from the stderr return from LSF call
	jobid=$(echo $VAR | awk '{ print $8 }')
	#strip off the first and last character
	jobid2=${jobid#"<"}
	jobid2=${jobid2%">"}
	echo $jobid2
	#Place the jobid into the array for steady checking later
	jobarray[$i]=$jobid2
	((i++))
done
#check job status continuously and if they are still running wait and check again.
echo "Printing the job array as finalized for the unzipping stage"
echo ${jobarray[@]}
sleep 5
echo "right before while"
while :
do
	runflag=1
        for ((j=0; j<$i; j++))
        do
                #check the bjobs and see if they exist
                running=$(bjobs ${jobarray[$j]} | grep 'RUN\|PEND' | wc -l)
                if [ $running == 1 ]; then
                        echo "job ${jobarray[$j]} is still running/pending"
                        runflag=0
                fi
        done
        if [ $runflag == 1 ]; then
                echo "we have hit the end of the line for the unzipping block of the pipeline now move the files around and check for EXIT status"
                break
        else
                sleep 30
        fi
done
#Now that everything has finished check the error logs for a failed job.  If failed job exit the pipeline and flag error.
mv *job logs/unzip
exit_status=$(grep Exit logs/unzip/*out | wc -l)
if [ $exit_status == 0 ]; then
	echo "We have had a successful run here through unzipping"
else
	echo "We have had an Exit somewhere in unzipping"
fi

#Step 2 - Submit alignment jobs to the cluster
# First make a file1 list and a file2 list

i=0
for g in *_1.txt
do

	g2=$(echo $g | sed -e 's/1.txt/2.txt/g')
	echo $g
	echo $g2
	# now that we have both files of a paired end run we need to extract sample ID from the 
	sampleID=$(echo $g | cut -f3 -d '_')
	echo $sampleID
	# Build the job submission and submit!
	# First test for the single/paired flag and change job accordingly
	if [ $3 == "single" ]; then
		echo "/hihg/smoke/applications/STAR/STAR-STAR_2.5.2a/bin/Linux_x86_64/STAR --limitIObufferSize 50000000 --outFileNamePrefix $sampleID --readFilesIn $g --runThreadN 8 --genomeDir $reference --quantMode TranscriptomeSAM GeneCounts --sjdbGTFfile $gtf --outSAMtype BAM SortedByCoordinate" >> align$i.job
	else
		echo "/hihg/smoke/applications/STAR/STAR-STAR_2.5.2a/bin/Linux_x86_64/STAR --limitIObufferSize 50000000 --outFileNamePrefix $sampleID --readFilesIn $g $g2 --runThreadN 8 --genomeDir $reference --quantMode TranscriptomeSAM GeneCounts --sjdbGTFfile $gtf --outSAMtype BAM SortedByCoordinate" >> align$i.job
	fi
	VAR=$((bsub -o logs/align/align$i.out -e logs/align/align$i.err -q hihg -x -P hihgbioinf sh align$i.job; ) 2>&1)
	#cut the jobid from the stderr return from LSF call
        jobid=$(echo $VAR | awk '{ print $8 }')
        #strip off the first and last character
        jobid2=${jobid#"<"}
        jobid2=${jobid2%">"}
        echo $jobid2
        #Place the jobid into the array for steady checking later
        jobarray[$i]=$jobid2
        ((i++))
done

#check job status continuously and if they are still running wait and check again.
echo "Printing the job array as finalized for the alignment stage"
echo ${jobarray[@]}
sleep 5
while :
do
        runflag=1
        for ((j=0; j<$i; j++))
        do
                #check the bjobs and see if they exist
                running=$(bjobs ${jobarray[$j]} | grep 'RUN\|PEND' | wc -l)
                if [ $running == 1 ]; then
                        echo "job ${jobarray[$j]} is still running/pending"
                        runflag=0
                fi
        done
        if [ $runflag == 1 ]; then
                echo "we have hit the end of the line for the alignment block of the pipeline now move the files around and check for EXIT status"
                break
        else
                sleep 30
        fi
done
#Now that everything has finished check the error logs for a failed job.  If failed job exit the pipeline and flag error.
mv *job logs/align
exit_status=$(grep Exit logs/align/*out | wc -l)
if [ $exit_status == 0 ]; then
        echo "We have had a successful run here through alignment"
else
        echo "We have had an Exit somewhere in alignment"
fi

# Step 3 - Duplicate calculation and stats summarization
i=0
for g in *Aligned.sortedByCoord.out.bam
do
	# extract the SampleID from the filename
	echo $g
	sampleID=${g/Aligned.sortedByCoord.out.bam/}
	echo $sampleID

	# Build the job submission and submit!
	echo "java -jar /hihg/smoke/applications/picard/picard-tools-1.42/MarkDuplicates.jar REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT AS=true I=$g O=$sampleID.rmdup.bam M=$sampleID.rmdup.picard.metrics.txt TMP_DIR=./tmp" >> stats$i.job
	echo "/hihg/smoke/applications/java/jdk1.8.0_77/bin/java -jar /hihg/smoke/applications/picard/picard-tools-2.1.1/picard.jar CollectInsertSizeMetrics I=$g O=$sampleID.IS.metrics.txt H=$sampleID.rmdup.picard.metrics.txt.histogram.pdf M=0.5" >> stats$i.job

	VAR=$((bsub -o logs/stats/stats$i.out -e logs/stats/stats$i.err -q general -R rusage[mem=7000] -P hihgbioinf sh stats$i.job; ) 2>&1)
	#cut the jobid from the stderr return from LSF call
        jobid=$(echo $VAR | awk '{ print $8 }')
        #strip off the first and last character
        jobid2=${jobid#"<"}
        jobid2=${jobid2%">"}
        echo $jobid2
        #Place the jobid into the array for steady checking later
        jobarray[$i]=$jobid2
        ((i++))
done

#check job status continuously and if they are still running wait and check again.
echo "Printing the job array as finalized for the alignment stage"
echo ${jobarray[@]}
sleep 5
while :
do
        runflag=1
        for ((j=0; j<$i; j++))
        do
                #check the bjobs and see if they exist
                running=$(bjobs ${jobarray[$j]} | grep 'RUN\|PEND' | wc -l)
                if [ $running == 1 ]; then
                        echo "job ${jobarray[$j]} is still running/pending"
                        runflag=0
                fi
        done
        if [ $runflag == 1 ]; then
                echo "we have hit the end of the line for the stats generation block of the pipeline now move the files around and check for EXIT status"
                break
        else
                sleep 30
        fi
done
#Now that everything has finished check the error logs for a failed job.  If failed job exit the pipeline and flag error.
mv *job logs/stats
exit_status=$(grep Exit logs/stats/*out | wc -l)
if [ $exit_status == 0 ]; then
        echo "We have had a successful run here up to stats"
else
        echo "We have had an Exit somewhere in rmdup and insert size stats collection"
fi

ls *Log.final.out > tmp1
echo -e "Sample ID\tNumber of reads\tUniquely mapped reads\tUnique Alignment %\tMulti-mapped Reads\tMulti Alignment %\tUnmapped (too short)%\tMean Insert Size\tDuplicate %" > tmp0
sed -i 's/Log.final.out//g' tmp1
grep "Number of input reads" *Log.final.out | cut -f2 > tmp2
grep "Uniquely mapped reads number" *Log.final.out | cut -f2 > tmp3
grep "Number of reads mapped to multiple loci" *Log.final.out | cut -f2 > tmp4
grep "% of reads unmapped: too short" *Log.final.out | cut -f2 > tmp5
paste tmp1 tmp2 tmp3 tmp4 tmp5 > star_stats.out
awk '{ print $1 "\t" $2 "\t" $3 "\t" $3/$2*100 "%\t" $4 "\t" $4/$2*100 "%\t" $5 }' star_stats.out > star1_stats.out
rm tmp1
rm tmp2
rm tmp3
rm tmp4
rm tmp5
rm star_stats.out
for f in *rmdup*txt
do
        head -8 $f | tail -n 1 | cut -f8 >> tmp1
done

for i in *IS*txt
do
        head -8 $i | tail -n 1 | cut -f5 >> tmp2
done
paste star1_stats.out tmp1 tmp2 > final_pipeline_stats.out
rm tmp1
rm tmp2
rm star1_stats.out
cat tmp0 final_pipeline_stats.out > final_stats.out
rm final_pipeline_stats.out
rm tmp0

if [[ $3 == "simple" ]]; then
	exit 0
fi

input=$4
i=0
if [[ $5 == "u" ]]; then
	j=2
fi
if [[ $5 == "s" ]]; then
	j=3
fi
if [[ $5 == "r" ]]; then
	j=4
fi
echo "cut -f$j"

while IFS= read -r line
do
	sampleID=$(echo $line | cut -d ' ' -f1)
	group=$(echo $line | cut -d ' ' -f2)
	#echo $sampleID
	echo $sampleID >> $sampleID.tmp0
	echo $group >> $sampleID.tmp0
	cut -f1 $sampleID*Gene.out.tab > genelist.tmp
	cut -f$j $sampleID*Gene.out.tab > $sampleID.tmp1
	geneWC=$(wc -l genelist.tmp)
	geneN=$(echo $geneWC | cut -d' ' -f1)
	(( geneN2 = $geneN - 4 ))
	tail -n $geneN2 $sampleID.tmp1 > $sampleID.tmp2
	cat $sampleID.tmp0 $sampleID.tmp2 > $sampleID.tmp3
	rm $sampleID.tmp0
	rm $sampleID.tmp1
	rm $sampleID.tmp2
	echo $geneN
	samplearray[$i]=$sampleID
	(( i++ ))
done < $input

echo "Symbol" >> header1
echo "Group" >> header1
echo $geneWC
geneN=$(echo $geneWC | cut -d' ' -f1)
(( geneN2 = $geneN - 4 ))
tail -n $geneN2 genelist.tmp > genelist2.tmp
cat header1 genelist2.tmp > firstcol
rm header1
rm genelist.tmp
cp firstcol running.tmp
rm genelist2.tmp
for ((j=0; j<$i; j++))
do
	echo "concatenating ${samplearray[$j]}.tmp"
	paste running.tmp ${samplearray[$j]}.tmp3 > running.tmp2
	mv running.tmp2 running.tmp
	rm ${samplearray[$j]}.tmp3
done
rm firstcol
mv running.tmp overall.txt
echo "R CMD BATCH '--args overall.txt $2 ALL' /hihg/smoke/applications/pipelines/rna/DE_pipeline.R" > de.job
bsub -q general -P hihgbioinf -o de.out -e de.err -R span[ptile=8] -n 8 -R rusage[mem=8000] sh de.job

echo "END OF PIPELINE"
