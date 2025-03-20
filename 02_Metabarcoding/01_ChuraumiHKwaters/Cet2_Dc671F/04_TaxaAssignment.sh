####
#### Taxa assignment using claident
####

# --------------------------------------------------------------------------------------- #
# OTU analysis
# --------------------------------------------------------------------------------------- #

# --------------------------------------------------------------------------------------- #
# Set file and folder names
# --------------------------------------------------------------------------------------- #
WORKING_FOLDER="XXXXXX/Cet2_Dc671F"

FASTA_FILE="../03_OTUClusteringOut/OTU_seqs.fa"
OUTPUT_FOLDER="04_TaxaAssignmentOut"

cd ${WORKING_FOLDER}
mkdir ${OUTPUT_FOLDER}
cd ${OUTPUT_FOLDER}

# overall_genus
clmakecachedb --blastdb=overall_genus --numthreads=36 ${FASTA_FILE} overall_genus_cache
clidentseq --blastdb=overall_genus_cache --numthreads=36 ${FASTA_FILE} overall_genus_clidentseq
classigntax --taxdb=overall_genus --maxpopposer=0.10 --minsoratio=9 overall_genus_clidentseq overall_genus_classigntax

# overall_class
clmakecachedb --blastdb=overall_class --numthreads=36 ${FASTA_FILE} overall_class_cache
clidentseq --blastdb=overall_class_cache --numthreads=36 ${FASTA_FILE} overall_class_clidentseq
classigntax --taxdb=overall_class --maxpopposer=0.10 --minsoratio=9 overall_class_clidentseq overall_class_classigntax

# Merge identification results (overall_class + overall_genus)
clmergeassign --priority=descend overall_genus_classigntax overall_class_classigntax final_classigntax_otu

# Delete large files
rm -r overall_class_cache
rm -r overall_genus_cache



# --------------------------------------------------------------------------------------- #
# ASV analysis
# --------------------------------------------------------------------------------------- #

# --------------------------------------------------------------------------------------- #
# Set file and folder names
# --------------------------------------------------------------------------------------- #
FASTA_FILE_ASV="../../02_DADA2Out/ASV_seqs.fa"
OUTPUT_FOLDER_ASV="TaxaAssignmentOut_ASV"

cd ${WORKING_FOLDER}/${OUTPUT_FOLDER}
mkdir ${OUTPUT_FOLDER_ASV}
cd ${OUTPUT_FOLDER_ASV}

# Overall_genus
clmakecachedb --blastdb=overall_genus --numthreads=36 ${FASTA_FILE_ASV} overall_genus_cache_asv
clidentseq --blastdb=overall_genus_cache_asv --numthreads=36 ${FASTA_FILE_ASV} overall_genus_clidentseq_asv
classigntax --taxdb=overall_genus --maxpopposer=0.10 --minsoratio=9 overall_genus_clidentseq_asv overall_genus_classigntax_asv

# overall_class
clmakecachedb --blastdb=overall_class --numthreads=36 ${FASTA_FILE_ASV} overall_class_cache_asv
clidentseq --blastdb=overall_class_cache_asv --numthreads=36 ${FASTA_FILE_ASV} overall_class_clidentseq_asv
classigntax --taxdb=overall_class --maxpopposer=0.10 --minsoratio=9 overall_class_clidentseq_asv overall_class_classigntax_asv

# Merge identification results (overall_class + overall_genus)
clmergeassign --priority=descend overall_genus_classigntax_asv overall_class_classigntax_asv final_classigntax_asv

# Delete large files
rm -r overall_class_cache_asv
rm -r overall_genus_cache_asv


