
# Step 1: Trimming - fastx_trimmer
nohup \time -f "Program: %C\nReal time: \t%E\nUser time (s): %U\nSystem time (s): %S\nCPU: %P\nAverage total memory use: %K"  -o Resources.txt bash trimming.sh  &

# Step 2: Mapping - BWA, HISAT2 & STAR
nohup \time -f "Program: %C\nReal time: \t%E\nUser time (s): %U\nSystem time (s): %S\nCPU: %P\nAverage total memory use: %K"  -o Resources.txt bash map_bwa.sh  &
nohup \time -f "Program: %C\nReal time: \t%E\nUser time (s): %U\nSystem time (s): %S\nCPU: %P\nAverage total memory use: %K"  -o Resources.txt bash map_hisat2.sh &
nohup \time -f "Program: %C\nReal time: \t%E\nUser time (s): %U\nSystem time (s): %S\nCPU: %P\nAverage total memory use: %K"  -o Resources.txt bash map_STAR.sh &

# Step 3: QC on BAM - Samtools & Picard-Tools
nohup \time -f "Program: %C\nReal time: \t%E\nUser time (s): %U\nSystem time (s): %S\nCPU: %P\nAverage total memory use: %K"  -o Resources.txt bash BAM_QC.sh &

# Step 4: Variant Calling - GATK
nohup \time -f "Program: %C\nReal time: \t%E\nUser time (s): %U\nSystem time (s): %S\nCPU: %P\nAverage total memory use: %K"  -o Resources.txt bash Pre-gatk.sh &
nohup \time -f "Program: %C\nReal time: \t%E\nUser time (s): %U\nSystem time (s): %S\nCPU: %P\nAverage total memory use: %K"  -o Resources.txt bash gatk.sh &
