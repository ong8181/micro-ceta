####
#### Quality filtering for Illumina data
####

# Mu31F - Dc320R (RevComp)
#^GACACTGAAAATGTCTAGATGG...CACCGCGGTCATACGATTRA
#^TYAATCGTATGACCGCGGTG...CCATCTAGACATTTTCAGTGTC
# Dc671F - Dc1015R (RevComp)
#^GCTACTYCAGTCTATATACC...GGTAAGCRTACYGGAARGTGTG
#^CACACYTTCCRGTAYGCTTACC...GGTATATAGACTGRAGTAGC
# Mu2084F - Dc2438R (RevComp)
#^ATGAAYGGCCACACGAGGGTTTTA...CCTCGATGTTGGATCAGGACA
#^TGTCCTGATCCAACATCGAGG...TAAAACCCTCGTGTGGCCRTTCAT
# Mu9459F - Mu9822R (RevComp)
#^CTGACTTCCAATCAGTTRGTTTCGG...CYCAARARGGYYTAGAATG
#^CATTCTARRCCYTYTTGRG...CCGAAACYAACTGATTGGAAGTCAG

# -------------------------------------------- #
# Step 0. Demultiplex and preparations
# -------------------------------------------- #
# Specify primer sequence (forward...reverse_revcomp, reverse...forward_revcomp)
PRIMER_FOR_READ=^CTGACTTCCAATCAGTTRGTTTCGG...CYCAARARGGYYTAGAATG
PRIMER_REV_READ=^CATTCTARRCCYTYTTGRG...CCGAAACYAACTGATTGGAAGTCAG
# Set parameters
WORKING_DIR=/home/rizabics/Desktop/NOV-001/Cet4_Mu9459F

# Set other parameters
INPUT_DIR=${WORKING_DIR}/seqdata_demultiplexed_Run3
OUTPUT_DIR=01_QualityFiltering_Run3Out

cd ${INPUT_DIR}
mkdir ../${OUTPUT_DIR}
mkdir ../${OUTPUT_DIR}/00_temp

# Check and save MD5
md5sum *.gz > 0_md5sum.txt


# -------------------------------------------- #
# Step 1. Global quality filtering
# -------------------------------------------- #
## Length > 100 bp
## No TruSeq adaptor
for file in *_L2_1.fq.gz; do
fastp \
--thread=16 \
--length_required=100 \
--adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
--adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
--in1 ${file} \
--in2 ${file%_L2_1.fq.gz}_L2_2.fq.gz \
--out1 ../${OUTPUT_DIR}/00_temp/${file%_L2_1.fq.gz}_flt_R1.fastq.gz \
--out2 ../${OUTPUT_DIR}/00_temp/${file%_L2_1.fq.gz}_flt_R2.fastq.gz \
-j ../${OUTPUT_DIR}/00_temp/report_${file%_L2_1.fq.gz}_fastp.json \
-h ../${OUTPUT_DIR}/00_temp/report_${file%_L2_1.fq.gz}_fastp.html
done

# Delete fastp summary
rm ../${OUTPUT_DIR}/00_temp/*fastp.html
rm ../${OUTPUT_DIR}/00_temp/*fastp.json


# -------------------------------------------- #
# Step 2.Trim primers
# -------------------------------------------- #
# Prepare directories
cd ../${OUTPUT_DIR}/00_temp
mkdir ../01_temp

# ----------------- cutadapt version
# Single Primer removal
for file in *_R1.fastq.gz; do
cutadapt -j 36 \
-a $PRIMER_FOR_READ \
-A $PRIMER_REV_READ \
-n 2 \
--discard-untrimmed \
-o ../01_temp/${file%_flt_R1.fastq.gz}_trimmed_R1.fastq.gz \
-p ../01_temp/${file%_flt_R1.fastq.gz}_trimmed_R2.fastq.gz \
${file} \
${file%_R1.fastq.gz}_R2.fastq.gz
done


# -------------------------------------------- #
# Step 3. Get summary table 
# -------------------------------------------- #
cd ${WORKING_DIR}/seqdata
cd ${INPUT_DIR}
seqkit stats -a *.gz > ${WORKING_DIR}/${OUTPUT_DIR}/01_1_demultiplexed.txt
seqkit stats -a *L2_1.fq.gz > ${WORKING_DIR}/${OUTPUT_DIR}/01_1_demultiplexedR1.txt
cd ../${OUTPUT_DIR}/00_temp
seqkit stats -a *.gz > ${WORKING_DIR}/${OUTPUT_DIR}/01_2_after_flt.txt
cd ../01_temp
seqkit stats -a *.gz > ${WORKING_DIR}/${OUTPUT_DIR}/01_3_after_cutadapt.txt
seqkit stats -a *R1.fastq.gz > ${WORKING_DIR}/${OUTPUT_DIR}/01_3_after_cutadaptR1.txt


# -------------------------------------------- #
# Step 4. Clean up 
# -------------------------------------------- #
cd ..
mv 01_temp/*.gz ./
rm -r 00_temp
rm -r 01_temp
