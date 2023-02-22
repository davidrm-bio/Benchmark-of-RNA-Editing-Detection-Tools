#!/bin/bash

# Taken from https://github.com/sarah-ku/hyperTRIBER 

export PATH=/binf-isilon/rennie/gsn480/scratch/bin/:$PATH

for f in `cat samples.txt`
do
myname=${f}_filt_sortrmdup

name_fwd1=${myname}_fwd1.bam
name_fwd2=${myname}_fwd2.bam
name_fwd=${myname}_fwd.bam

name_rev1=${myname}_rev1.bam
name_rev2=${myname}_rev2.bam
name_rev=${myname}_rev.bam

samtools view -bh -f 99 ${myname}.bam > $name_fwd1
samtools index $name_fwd1

samtools view -bh -f 147 ${myname}.bam > $name_fwd2
samtools index $name_fwd2

samtools merge -f $name_fwd $name_fwd1 $name_fwd2
samtools index $name_fwd

samtools view -bh -f 83 ${myname}.bam > $name_rev1
samtools index $name_rev1

samtools view -bh -f 163 ${myname}.bam > $name_rev2
samtools index $name_rev2

samtools merge -f $name_rev $name_rev1 $name_rev2
samtools index $name_rev

rm $name_fwd1
rm $name_fwd2
rm $name_fwd1.bai
rm $name_fwd2.bai

rm $name_rev1
rm $name_rev2
rm $name_rev1.bai
rm $name_rev2.bai
done