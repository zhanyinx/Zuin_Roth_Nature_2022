#!/bin/bash

echo "Total number of reads"
/tungstenfs/groups/gbioinfo/Appz/SAM_tools/samtools-current/samtools view   $1 | awk '{ print $1}' | sort | uniq | wc -l 


echo ""
echo "Total number of reads >5kb"
/tungstenfs/groups/gbioinfo/Appz/SAM_tools/samtools-current/samtools view   $1 | awk '{if(length($10)>5000) print $1}' | sort | uniq | wc -l


echo ""
echo "Total number of mapped reads"
/tungstenfs/groups/gbioinfo/Appz/SAM_tools/samtools-current/samtools view -F 260  $1 | awk '{ print $1}' | sort | uniq | wc -l



echo ""
echo "Total number of mapped reads (>5kb)"
/tungstenfs/groups/gbioinfo/Appz/SAM_tools/samtools-current/samtools view -F 260  $1 | awk '{ if(length($10)>5000)  print $1}' | sort | uniq | wc -l


echo ""
echo "Total number of mapped reads (>5kb) to cassette $2"
/tungstenfs/groups/gbioinfo/Appz/SAM_tools/samtools-current/samtools view -F 260  $1 | awk '{ if(length($10)>5000 && $3=="Southern_"'"$2"'".dna")  print $1}' | sort | uniq | wc -l


echo ""
echo "number of reads spanning cassette"
left=`awk '{if($1 ~ /'"$2"'/ && $4 ~ /left/ && $4 ~/homology/) print $2 }' /tungstenfs/scratch/ggiorget/zhan/2020/Nanopore//annotations/cassette_annotation.BED`
right=`awk '{if($1 ~ /'"$2"'/ && $4 ~ /right/ && $4 ~/homology/) print $3 }' /tungstenfs/scratch/ggiorget/zhan/2020/Nanopore//annotations/cassette_annotation.BED`
/tungstenfs/groups/gbioinfo/Appz/SAM_tools/samtools-current/samtools view -F 260  $1 | awk 'BEGIN{left="'"$left"'"+0.; right = "'"$right"'"+0.}{if($3=="Southern_"'"$2"'".dna" && $4<=left && $4+length($10)>=right)  print $0}' > appo
awk '{
	for(i=1;i<=length($6);i++) {
		if(!(substr($6,i,1) ~ /^[0-9]/)){
			a= substr($6,i,1); break
		}
	}
	left=1
	right=1
	if(a=="S") left=-1
	if(substr($6,length($6),1)=="S") right=-1 
	
	print $6,left,right
	
	}' appo | sed 's/S/ /g' | sed 's/M/ /g' | sed 's/I/ /g' | sed 's/D/ /g' > appo1

awk 'BEGIN{left="'"$left"'"+0.; right = "'"$right"'"+0.; n=0; fn=0; m=0}{
	if(FNR==1) fn++
	if(fn==1){
		len[n]=0
		if($(NF-1)==-1){
			lsoft_clipped[n]=$1
		}
		else{
			lsoft_clipped[n]=0
			len[n]+=$1
		}
		if($NF == -1) rsoft_clipped[n]=$(NF-2)
		else{
			rsoft_clipped[n]=0
			len[n]+=$(NF-2)
		}

		for(i=2;i<NF-2;i++) len[n]+=$i
		n++
	}
	if(fn==2){
		if($3=="Southern_"'"$2"'".dna" && $4<=left && ($4+len[m])>=right)  c++
		m++
	}
	}END{print c}' appo1 appo 



echo ""
echo "number of reads spanning WT allele"
left=`awk '{if($4 ~ /L_HA/ ) print $2 }' /tungstenfs/scratch/ggiorget/zhan/2020/Nanopore//annotations/guide_rna_HA_annotations_mm9.BED`
right=`awk '{if($4 ~ /R_HA/ ) print $3 }' /tungstenfs/scratch/ggiorget/zhan/2020/Nanopore//annotations/guide_rna_HA_annotations_mm9.BED`
chr=`awk '{if($4 ~ /R_HA/ ) print $1 }' /tungstenfs/scratch/ggiorget/zhan/2020/Nanopore//annotations/guide_rna_HA_annotations_mm9.BED`
/tungstenfs/groups/gbioinfo/Appz/SAM_tools/samtools-current/samtools view -F 260  $3 | awk 'BEGIN{left="'"$left"'"+0.; right = "'"$right"'"+0.}{if($3=="'"$chr"'" && $4<=left && $4+length($10)>=right)  print $0}'  > gappo
awk '{
	for(i=1;i<=length($6);i++) {
		if(!(substr($6,i,1) ~ /^[0-9]/)){
			a= substr($6,i,1); break
		}
	}
	left=1
	right=1
	if(a=="S") left=-1
	if(substr($6,length($6),1)=="S") right=-1 
	
	print $6,left,right
	
	}' gappo | sed 's/S/ /g' | sed 's/M/ /g' | sed 's/I/ /g' | sed 's/D/ /g' > gappo1

awk 'BEGIN{left="'"$left"'"+0.; right = "'"$right"'"+0.; n=0; fn=0; m=0}{
	if(FNR==1) fn++
	if(fn==1){
		len[n]=0
		if($(NF-1)==-1){
			lsoft_clipped[n]=$1
		}
		else{
			lsoft_clipped[n]=0
			len[n]+=$1
		}
		if($NF == -1) rsoft_clipped[n]=$(NF-2)
		else{
			rsoft_clipped[n]=0
			len[n]+=$(NF-2)
		}

		for(i=2;i<NF-2;i++) len[n]+=$i
		n++
	}
	if(fn==2){
		if($3=="'"$chr"'" && $4<=left && ($4+len[m])>=right)  c++
		m++
	}
	}END{print c}' gappo1 gappo 


rm *appo*
