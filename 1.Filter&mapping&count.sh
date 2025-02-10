# Cheak MD5
MD5sum *.fq.gz

# Fastqc
fastqc -o ./report -t 16 ./*.gz

# Remove connectors and low quality reads
for ((i=1; i<=n; i++))
do
    fastp -i /public/home/fuguiling/work/NAU/0.NAU/G"$i"_1.fq.gz \
          -I /public/home/fuguiling/work/NAU/0.NAU/G"$i"_2.fq.gz \
          -o ./fp_G"$i"_1.fq.gz \
          -O ./fp_G"$i"_2.fq.gz \
          --thread=16 -h ./G"$i".html
done

# Fastqc
fastqc -o ./report -t 16 ./*.gz

# Biuld ref_index
hisat2-build ./AD1_ref.fa TM-1index

# Mapping
for ((i=1; i<=n; i++))
do
    hisat2 -x ./TM-1index \
          -S G"$i".sam \
		  -p 16 -q \
          -1 ./fp_G"$i"_1.fq.gz \
		  -2 ./fp_G"$i"_2.fq.gz
done

# Sam2bam
for ((i=1; i<=n; i++))
do
    samtools sort \
          -o G"$i".bam \
		  ./G"$i".sam
done

# Calculate count
for ((i=1; i<=n; i++))
do
    Rscript ./run-featurecounts.R \
          -b ./G"$i".bam \
		  -g ./Ghirsutum_527_v2.1_geneid.gtf \
		  -o G"$i"
done
