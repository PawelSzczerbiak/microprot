#!/bin/bash -l

#SBATCH -J test_micro
#SBATCH -N 1
#SBATCH --ntasks-per-node=24
#SBATCH --mem=120GB
#SBATCH --time=00:30:00
#SBATCH -A tomaszlabmcb1
#SBATCH -p plgrid-testing
#SBATCH -C memfs

#=============================================
#           UTILITY FUNCTIONS
#=============================================

timestamp() {
  date --rfc-3339=seconds
}

showtime() {
  DIFF=$(echo "$(date +%s.%N) - $START" | bc)
  echo $1": " $DIFF " s"
}

printparams(){
  echo "Number of seqs:  " $1
  echo "Parallelization: " $2
  echo "Number of cores: " $3
  echo "OMP_NUM_THREADS: " $4
  echo "E-value:         " $5
}
#=============================================

echo "Job ID: " $SLURM_JOB_ID

# Arguments
BATCH=$1_$2; SEQ_ids=$2; UHGP_NAME=$3
echo "Batch name:   " $BATCH
echo "Sequence ids: " $SEQ_ids
echo "UHGP DB name: " $UHGP_NAME

#=========== MODIFIABLE SECTION ==============

echo -e "\nSetup environment..."
timestamp; START=$(date +%s.%N)

# hh-suite module import
module load plgrid/tools/hh-suite/3.1.0
# conda activation
#source /net/archive/groups/plggtomlab/tools/miniconda3/etc/profile.d/conda.sh
#conda activate microprot_old # TODO
export PATH="/net/archive/groups/plggtomlab/tools/miniconda3/bin:$PATH"
source activate microprot_old # TODO
showtime "Setup time"

# Directories
export MICRO_DIR=/net/archive/groups/plggtomlab/tools/microprot/microprot
INP_DIR=$PLG_GROUPS_STORAGE/plggtomlab/microprot/input/$BATCH
export OUT_DIR=$SCRATCH/plggtomlab/microprot/output/$BATCH
export LOG_DIR=$SCRATCH/plggtomlab/microprot/log/$BATCH

export CM_DIR_NAME="01-CM"
export PFAM_DIR_NAME="02-Pfam"
export MSA_DIR_NAME="03-MSA_hhblits"
export NOTRIPE_DIR_NAME="04-not_ripe"
export ROSETTA_DIR_NAME="05-Rosetta"

# Databases
UHGP_DIR=$PLG_GROUPS_STORAGE/plggtomlab/dbs/uhgp/
PDB_DIR=$PLG_GROUPS_STORAGE/plggtomlab/dbs/pdb/
PFAMDB_DIR=$PLG_GROUPS_STORAGE/plggtomlab/dbs/pfam/
UNIDB_DIR=$PLG_GROUPS_STORAGE/plggtomlab/dbs/uniclust/uniclust30_2018_08/

#======== END OF MODIFIABLE SECTION ==========

echo -e "\nCreating input files / output & log folders..."
timestamp; START=$(date +%s.%N)

rm -rf $LOG_DIR $INP_DIR $OUT_DIR
mkdir -p $LOG_DIR $INP_DIR $OUT_DIR

# Create input files with fasta sequences
python $MICRO_DIR/scripts/parse_inputs.py -f \
$UHGP_DIR/$UHGP_NAME -d $SEQ_ids -i $INP_DIR -o $LOG_DIR
echo "Number of extracted sequences: " $(find $INP_DIR/*.fasta | wc -l)

# Create separate output folder for every sequence
cd $MICRO_DIR/slurm/
find $INP_DIR/*.fasta | xargs -t -d "\n" -P 24 -n 1 ./create_outputs.sh

showtime "Creation time"

#=============================================
#                    CM
#=============================================

export STEP="CM"
export DIR_NAME=$CM_DIR_NAME
export NEXT_DIR_NAME=$PFAM_DIR_NAME

#=========================
echo -e "\nLoading PDB..."
#=========================
timestamp; START=$(date +%s.%N)

mkdir -p $MEMFS/pdb
cp $PDB_DIR/*.ffindex $PDB_DIR/*.tsv $PDB_DIR/*dat* $MEMFS/pdb
echo "DB content: " 
find $MEMFS/pdb -type f
export MEMDB=$MEMFS/pdb/pdb70

showtime "PDB loading time" 

#======================
echo -e "\nHHSearch..."
#======================
timestamp; START=$(date +%s.%N)

SEQS=$(find $OUT_DIR/*/$DIR_NAME/*.fasta | wc -l)
P=$SEQS
if [ $P -gt $OMP_NUM_THREADS ] ; then  P=$OMP_NUM_THREADS ; fi
export CORES=$(($OMP_NUM_THREADS/$P))
export EVALUE=0.1
printparams $SEQS $P $CORES $OMP_NUM_THREADS $EVALUE

cd $MICRO_DIR/slurm/
find $OUT_DIR/*/$DIR_NAME/*.fasta | xargs -t -d "\n" -P $P -n 1 ./run_hhsearch.sh

showtime "Computation time"

#==========================
echo -e "\nRemoving PDB..."
#==========================
timestamp; START=$(date +%s.%N)

rm -rf $MEMFS/pdb

showtime "PDB removing time"

#==========================
echo -e "\nCM splitting..."
#==========================
timestamp; START=$(date +%s.%N)

export MATCH_EXP="non_match"  # create .fasta from non_match files
export MAX_EVALUE=0.1
export MIN_FRAG_LEN=40
echo "Max e-value:         " $MAX_EVALUE
echo "Min fragment length: " $MIN_FRAG_LEN

cd $MICRO_DIR/slurm/
find $OUT_DIR/*/$DIR_NAME/*.out | xargs -t -d "\n" -P 24 -n 1 ./split_search.sh

showtime "CM splitting time"

#=============================================
#                   Pfam
#=============================================

export STEP="Pfam"
export DIR_NAME=$PFAM_DIR_NAME
export NEXT_DIR_NAME=$MSA_DIR_NAME

#==========================
echo -e "\nLoading Pfam..."
#==========================
timestamp; START=$(date +%s.%N)

mkdir -p $MEMFS/pfam
cp $PFAMDB_DIR/*.ffindex $PFAMDB_DIR/*ffdata $MEMFS/pfam
echo "DB content: "
find $MEMFS/pfam -type f
export MEMDB=$MEMFS/pfam/pfam

showtime "Pfam loading time"

#======================
echo -e "\nHHSearch..."
#======================
timestamp; START=$(date +%s.%N)

SEQS=$(find $OUT_DIR/*/$DIR_NAME/*.fasta | wc -l)
P=$SEQS
if [ $P -gt $OMP_NUM_THREADS ] ; then P=$OMP_NUM_THREADS ; fi
export CORES=$(($OMP_NUM_THREADS/$P))
export EVALUE=0.1
printparams $SEQS $P $CORES $OMP_NUM_THREADS $EVALUE

cd $MICRO_DIR/slurm/
find $OUT_DIR/*/$DIR_NAME/*.fasta | xargs -t -d "\n" -P $P -n 1 ./run_hhsearch.sh

showtime "Computation time"

#===========================
echo -e "\nRemoving Pfam..."
#===========================
timestamp; START=$(date +%s.%N)

rm -rf $MEMFS/pfam

showtime "Pfam removing time"

#============================
echo -e "\nPfam splitting..."
#============================
timestamp; START=$(date +%s.%N)

export MATCH_EXP="match"  # create .fasta from both match and non_match files
export MAX_EVALUE=0.01
export MIN_PROB=90.0
export MIN_FRAG_LEN=40
echo "Max e-value:         " $MAX_EVALUE
echo "Min probability:     " $MIN_PROB
echo "Min fragment length: " $MIN_FRAG_LEN

cd $MICRO_DIR/slurm/
find $OUT_DIR/*/$DIR_NAME/*.out | xargs -t -d "\n" -P 24 -n 1 ./split_search.sh

showtime "Pfam splitting time"

#=============================================
#                MSA HHblits
#=============================================

export STEP="MSA"
export DIR_NAME=$MSA_DIR_NAME

#==============================
echo -e "\nLoading Uniclust..."
#==============================
timestamp; START=$(date +%s.%N)

mkdir -p $MEMFS/uniclust
cp $UNIDB_DIR/*index $UNIDB_DIR/*cs219 $UNIDB_DIR/*sizes $UNIDB_DIR/*dat* $MEMFS/uniclust
echo "DB content: "
find $MEMFS/uniclust -type f
export MEMDB=$MEMFS/uniclust/uniclust30_2018_08

showtime "Uniclust loading time"

#=====================
echo -e "\nHHblits..."
#=====================
timestamp; START=$(date +%s.%N)

SEQS=$(find $OUT_DIR/*/$DIR_NAME/*.fasta | wc -l)
P=$SEQS
if [ $P -gt 6 ] ; then P=6 ; fi
export CORES=$(($OMP_NUM_THREADS/$P))
export EVALUE=0.001
printparams $SEQS $P $CORES $OMP_NUM_THREADS $EVALUE

cd $MICRO_DIR/slurm/
find $OUT_DIR/*/$DIR_NAME/*.fasta | xargs -t -d "\n" -P $P -n 1 ./run_hhblits.sh

showtime "Computation time"

#===============================
echo -e "\nRemoving Uniclust..."
#===============================
timestamp; START=$(date +%s.%N)

rm -rf $MEMFS/uniclust

showtime "Uniclust removing time"
