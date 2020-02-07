#!/bin/bash -l

#SBATCH -J test_micro
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=120GB
#SBATCH --time=00:05:00
#SBATCH -A tomaszlabmcb1
#SBATCH -p plgrid-testing
#SBATCH -C memfs

#=============================================
#       UTILITY FUNCTIONS & VARIABLES
#=============================================

readonly start_total=$(date +%s.%N)
readonly format=" %-25s %s\n"
readonly format_long=" %-50s %s\n"

timestamp() {
  date --rfc-3339=seconds
}

showtime() {
  if [ $# = 2 ]; then start="${2}"; fi
  local diff=$(echo "$(date +%s.%N) - ${start}" | bc)
  printf "%s: %s s\n" "${1}" "${diff}"
}

printparams() {
  printf "${format}" \
    "Number of seqs :" "${1}" \
    "Parallelization:" "${2}" \
    "Number of cores:" "${3}" \
    "OMP_NUM_THREADS:" "${4}" \
    "hh-suite params:" "${5}"
}

header() {
  local chrlen=${#1}
  printf "\n"
  printf '=%.0s' $(seq 1 "${chrlen}")
  printf "\n%s\n" "${1}"
  printf '=%.0s' $(seq 1 "${chrlen}")
  printf "\n\n"
}
#=========== MODIFIABLE SECTION ==============

# 1) Arguments

# number of expected arguments
readonly no_of_args=12

if [ $# -ne ${no_of_args} ]; then
  printf "%s %s %s %s\n" "ERROR: expected" "${no_of_args}" "arguments got:" $#
  args=("$@")
  for i in "${!args[@]}"; do
      printf "%s:\t%s\n" "$((i+1))" "${args[$i]}"
  done
  exit 1
fi

readonly BATCH="${1}"
readonly SEQ_IDs="${2}"
readonly UHGP_NAME="${3}"
readonly PDB_NAME="${4}"
readonly PFAM_NAME="${5}"
readonly UNICLUST_NAME="${6}"
readonly MAX_CORES_CM="${7}"
readonly MAX_CORES_PFAM="${8}"
readonly MAX_CORES_MSA="${9}"
readonly MAX_CORES_NEFF="${10}"
readonly CONFIG="${11}"
readonly DEST_DIR="${12}"

export BATCH # used in organize_output script

printf "$format" \
  "Slurm job ID:" "${SLURM_JOB_ID}" \
  "Batch name:" "${BATCH}" \
  "Sequence ids:" "${SEQ_IDs}" \
  "UHGP DB name:" "${UHGP_NAME}" \
  "PDB DB name:" "${PDB_NAME}" \
  "Pfam DB name:" "${PFAM_NAME}" \
  "Uniclust DB name:" "${UNICLUST_NAME}" \
  "Max cores for CM:" "${MAX_CORES_CM}" \
  "Max cores for Pfam:" "${MAX_CORES_PFAM}" \
  "Max cores for MSA:" "${MAX_CORES_MSA}" \
  "Max cores for NEFF:" "${MAX_CORES_NEFF}" \
  "Config file:" "${CONFIG}" \
  "Destination directory:" "${DEST_DIR}"

printf "\n%s\n" "Setup environment..."
timestamp; start=$(date +%s.%N)

# 2) Directories

# main directories
MAIN_DIR="${PLG_GROUPS_STORAGE}/plggtomlab/"
WORK_DIR="${SCRATCH}/plggtomlab/"
# microprot directory
export MICRO_DIR="${MAIN_DIR}/tools/microprot/microprot"
# input directory (where the input data will be stored)
INP_DIR="${MAIN_DIR}/microprot/input/${BATCH}"
# Result directory for files staged for sending to DEST_DIR
export RES_DIR="${MAIN_DIR}/microprot"
# main working & log directories
export OUT_DIR="${WORK_DIR}/microprot/output/${BATCH}"
export LOG_DIR="${MAIN_DIR}/microprot/log/${BATCH}"
# folder names for the subsequent microprot steps
export CM_DIR_NAME="01-CM"
export PFAM_DIR_NAME="02-Pfam"
export MSA_DIR_NAME="03-MSA_hhblits"
export NOTRIPE_DIR_NAME="04-not_ripe"
export ROSETTA_DIR_NAME="05-Rosetta"
# step names
CM_STEP="CM"
PFAM_STEP="Pfam"
MSA_STEP="MSA"
RIPE_STEP="Ripe"
# databases
UHGP_DIR="${MAIN_DIR}/dbs/uhgp/"
PDB_DIR="${MAIN_DIR}/dbs/pdb/${PDB_NAME}"
PFAMDB_DIR="${MAIN_DIR}/dbs/pfam/${PFAM_NAME}"
UNIDB_DIR="${MAIN_DIR}/dbs/uniclust/${UNICLUST_NAME}"

# 3) Modules

# import config variables
source "${CONFIG}"
# hh-suite module import
module load plgrid/tools/hh-suite/3.1.0
# conda activation (choose prefered way)
#source "${MAIN_DIR}/tools/miniconda3/etc/profile.d/conda.sh"
#conda activate microprot_old # TODO
export PATH="${MAIN_DIR}/tools/miniconda3/bin:$PATH"
source activate microprot_old # TODO

# 4) Parameters

CM_MATCH_EXP="non_match" # create .fasta from .non_match files
PFAM_MATCH_EXP="match" # create .fasta from both .match and .non_match files

showtime "Setup time"

#======== END OF MODIFIABLE SECTION ==========

printf "\n%s\n" "Creating input files / output & log folders..."
timestamp; start=$(date +%s.%N)

rm -rf "${LOG_DIR}" "${INP_DIR}" "${OUT_DIR}"
mkdir -p "${LOG_DIR}" "${INP_DIR}" "${OUT_DIR}"

# Create input files with fasta sequences
python "${MICRO_DIR}"/scripts/parse_inputs.py -f \
  "${UHGP_DIR}/${UHGP_NAME}" -d "${SEQ_IDs}" -i "${INP_DIR}" -o "${LOG_DIR}"
printf "${format_long}" "Number of extracted sequences:" "$(find "${INP_DIR}"/*.fasta | wc -l)"

# Create separate output folder for every sequence
# and move input sequences there
(
  cd "${MICRO_DIR}"/slurm/
  find "${INP_DIR}"/*.fasta | xargs -d "\n" -P ${OMP_NUM_THREADS} -n 1 ./create_folders
)
printf "${format_long}" "Number of sequences distributed for processing:" \
 "$(find "${OUT_DIR}"/*/"${CM_DIR_NAME}"/*.fasta | wc -l)"
showtime "Creation time"

#=============================================
#                    CM
#=============================================

export STEP="${CM_STEP}"
DIR_NAME="${CM_DIR_NAME}"
export NEXT_DIR_NAME="${PFAM_DIR_NAME}"

#======================
header "Loading PDB..."
#======================
timestamp; start=$(date +%s.%N)

mkdir -p "${MEMFS}"/pdb
cp "${PDB_DIR}"/*.ffindex "${PDB_DIR}"/*.tsv "${PDB_DIR}"/*dat* "${MEMFS}"/pdb
printf "%s\n" "DB content:"
find "${MEMFS}"/pdb -type f
export MEMDB="${MEMFS}/pdb/${PDB_NAME/_*/}"  # e.g. pdb70_2019_12 -> pdb70
export DB_NAME="${PDB_NAME}"

showtime "PDB loading time"

#===================
header "HHSearch..."
#===================
timestamp; start=$(date +%s.%N)

SEQS=$(find "${OUT_DIR}"/*/"${DIR_NAME}"/*.fasta | wc -l)
P="${SEQS}"
if [ ${P} -gt ${MAX_CORES_CM} ]; then P=${MAX_CORES_CM}; fi
export CORES=$((${OMP_NUM_THREADS} / ${P}))
export PARAMS="${CM_PARAMS}"
printparams "${SEQS}" "${P}" "${CORES}" "${OMP_NUM_THREADS}" "${PARAMS}"
(
  cd "${MICRO_DIR}"/slurm/
  find "${OUT_DIR}"/*/"${DIR_NAME}"/*.fasta | xargs -d "\n" -P "${P}" -n 1 ./run_hhsearch
)
showtime "Computation time"

#=======================
header "Removing PDB..."
#=======================
timestamp; start=$(date +%s.%N)

rm -rf "${MEMFS}"/pdb

showtime "PDB removing time"

#=======================
header "CM splitting..."
#=======================
timestamp; start=$(date +%s.%N)

export MATCH_EXP="${CM_MATCH_EXP}"
export MAX_EVALUE="${CM_MAX_EVALUE}"
export MIN_FRAG_LEN="${CM_MIN_FRAG_LEN}"
printf "${format}" "Max e-value:" "${MAX_EVALUE}" "Min fragment length:" "${MIN_FRAG_LEN}"
(
  cd "${MICRO_DIR}"/slurm/
  find "${OUT_DIR}"/*/"${DIR_NAME}"/*.out | xargs -d "\n" -P "${OMP_NUM_THREADS}" -n 1 ./split_search
)
showtime "CM splitting time"

#=============================================
#                   Pfam
#=============================================

export STEP="${PFAM_STEP}"
DIR_NAME="${PFAM_DIR_NAME}"
SEQS=$(find "${OUT_DIR}"/*/"${DIR_NAME}"/*.fasta | wc -l)

# Pfam and further steps are performed only
# when the CM step didnt exclude all sequences
if [ ${SEQS} -gt 0 ]; then

  export NEXT_DIR_NAME="${MSA_DIR_NAME}"

  #=======================
  header "Loading Pfam..."
  #=======================
  timestamp; start=$(date +%s.%N)

  mkdir -p "${MEMFS}"/pfam
  cp "${PFAMDB_DIR}"/*.ffindex "${PFAMDB_DIR}"/*ffdata "${MEMFS}"/pfam
  printf "%s\n" "DB content:"
  find "${MEMFS}"/pfam -type f
  export MEMDB="${MEMFS}"/pfam/pfam
  export DB_NAME="${PFAM_NAME}"

  showtime "Pfam loading time"

  #===================
  header "HHSearch..."
  #===================
  timestamp; start=$(date +%s.%N)

  P="${SEQS}"
  if [ ${P} -gt ${MAX_CORES_PFAM} ]; then P=${MAX_CORES_PFAM}; fi
  export CORES=$((${OMP_NUM_THREADS} / ${P}))
  export PARAMS="${PFAM_PARAMS}"
  printparams "${SEQS}" "${P}" "${CORES}" "${OMP_NUM_THREADS}" "${PARAMS}"
  (
    cd "${MICRO_DIR}"/slurm/
    find "${OUT_DIR}"/*/"${DIR_NAME}"/*.fasta | xargs -d "\n" -P "${P}" -n 1 ./run_hhsearch
  )
  showtime "Computation time"

  #========================
  header "Removing Pfam..."
  #========================
  timestamp; start=$(date +%s.%N)

  rm -rf "${MEMFS}"/pfam

  showtime "Pfam removing time"

  #=========================
  header "Pfam splitting..."
  #=========================
  timestamp; start=$(date +%s.%N)

  export MATCH_EXP="${PFAM_MATCH_EXP}" # create .fasta from both match and non_match files
  export MAX_EVALUE="${PFAM_MAX_EVALUE}"
  export MIN_PROB="${PFAM_MIN_PROB}"
  export MIN_FRAG_LEN="${PFAM_MIN_FRAG_LEN}"
  printf "${format}" "Max e-value:" "${MAX_EVALUE}" \
    "Min probability:" "${MIN_PROB}" "Min fragment length:" "${MIN_FRAG_LEN}"
  (
    cd "${MICRO_DIR}"/slurm/
    find "${OUT_DIR}"/*/"${DIR_NAME}"/*.out | xargs -d "\n" -P "${OMP_NUM_THREADS}" -n 1 ./split_search
  )
  showtime "Pfam splitting time"

  #=============================================
  #                MSA HHblits
  #=============================================

  export STEP="${MSA_STEP}"
  DIR_NAME="${MSA_DIR_NAME}"

  #===========================
  header "Loading Uniclust..."
  #===========================
  timestamp; start=$(date +%s.%N)

  mkdir -p "${MEMFS}"/uniclust
  cp "${UNIDB_DIR}"/*index "${UNIDB_DIR}"/*cs219 "${UNIDB_DIR}"/*sizes "${UNIDB_DIR}"/*dat* "${MEMFS}"/uniclust
  printf "%s\n" "DB content:"
  find "${MEMFS}"/uniclust -type f
  export MEMDB="${MEMFS}/uniclust/${UNICLUST_NAME}"
  export DB_NAME="${UNICLUST_NAME}"

  showtime "Uniclust loading time"

  #==================
  header "HHblits..."
  #==================
  timestamp; start=$(date +%s.%N)

  SEQS=$(find "${OUT_DIR}"/*/"${DIR_NAME}"/*.fasta | wc -l)
  P="${SEQS}"
  if [ ${P} -gt ${MAX_CORES_MSA} ]; then P=${MAX_CORES_MSA}; fi
  export CORES=$((${OMP_NUM_THREADS} / ${P}))
  # NOTE: HHblits parameters are imported as MSA_PARAMS (see config)
  printparams "${SEQS}" "${P}" "${CORES}" "${OMP_NUM_THREADS}" "${MSA_PARAMS}"
  (
    cd "${MICRO_DIR}"/slurm/
    find "${OUT_DIR}"/*/"${DIR_NAME}"/*.fasta | xargs -d "\n" -P "${P}" -n 1 ./run_hhblits
  )
  showtime "Computation time"

  #============================
  header "Removing Uniclust..."
  #============================
  timestamp; start=$(date +%s.%N)

  rm -rf "${MEMFS}"/uniclust

  showtime "Uniclust removing time"

  #=============================================
  #                 MSA ripe
  #=============================================

  export STEP="${RIPE_STEP}"
  DIR_NAME="${MSA_DIR_NAME}"

  #============================
  header "Calculating N_eff..."
  #============================
  timestamp; start=$(date +%s.%N)

  # NOTE: MSA cutoff and min neff are imported as MSA_CUTOFF and MSA_MIN_NEFF (see config)
  printf "${format}" "Cutoff:" "${MSA_CUTOFF}" "Min N_eff:" "${MSA_MIN_NEFF}"
  (
    cd "${MICRO_DIR}"/slurm/
    find "${OUT_DIR}"/*/"${DIR_NAME}"/*.a3m | xargs -d "\n" -P "${MAX_CORES_NEFF}" -n 1 ./calculate_Neff
  )
  showtime "Calculation time"
else
  printf "\n%s\n" "Pfam and MSA steps omitted - no sequences after the CM step."
fi

#=============================================
#               Rules check
#=============================================

# CM and Pfam rules are the same except the regexp for (non_)match files

# CM rules:
# 1)  no. of .fasta files = no. of .a3m/.out files in CM dir (look for errors in _CM_hhsearch.err)
# 2)  no. of lines in .non_match = no. of .fasta files in Pfam dir (look for errors in _CM_split.err)
# 3)  _CM_split.err filies empty

# Pfam rules:
# 1)  no. of .fasta files = no of .a3m/.out files in Pfam dir (look for errors in _Pfam_hhsearch.err)
# 2)  no. of lines in .non_match + .match = no of .fasta files in MSA dir (look for errors in _Pfam_split.err)
# 3)  _Pfam_split.err filies empty

# MSA rules:
# 1)  no. of .fasta files = no of .a3m/.out files in MSA dir (look for errors in _MSA_hhblits.err)
# 2)  no. of .fasta files = no of .neff files in MSA dir (look for errors in _Ripe.err)
# 3)  no. of .fasta files = no of .a3m files in NOTRIPE + ROSETTA dirs (look for errors in _Ripe.err)
# 4)  _Ripe.err filies empty

#=========================
header "Checking rules..."
#=========================
timestamp; start=$(date +%s.%N)

export ERR_DIR="${LOG_DIR}/problematic_sequences"
rm -rf "${ERR_DIR}"
mkdir -p "${ERR_DIR}"
(
  cd "${MICRO_DIR}"/slurm/

  export STEP="${CM_STEP}"
  export MATCH_EXP="${CM_MATCH_EXP}"
  export NEXT_DIR_NAME="${PFAM_DIR_NAME}"
  find "${OUT_DIR}"/*/"${CM_DIR_NAME}"/*.fasta | xargs -d "\n" -P "${OMP_NUM_THREADS}" -n 1 ./check_standard

  export STEP="${PFAM_STEP}";
  export MATCH_EXP="${PFAM_MATCH_EXP}"
  export NEXT_DIR_NAME="${MSA_DIR_NAME}"
  find "${OUT_DIR}"/*/"${PFAM_DIR_NAME}"/*.fasta | xargs -d "\n" -P "${OMP_NUM_THREADS}" -n 1 ./check_standard

  export STEP="${MSA_STEP}"
  find "${OUT_DIR}"/*/"${MSA_DIR_NAME}"/*.fasta | xargs -d "\n" -P "${OMP_NUM_THREADS}" -n 1 ./check_MSA
)

showtime "Checking time"

# Do not copy files and update logs when testing
if [[ ${BATCH,,} != *"test"* ]]; then

  #============================
  header "Organizing output..."
  #============================
  timestamp; start=$(date +%s.%N)

  # 1) Create lists:
  #  - problematic_sequences_A.out with sequences that failed to be checked
  #  - problematic_sequences_B.out with sequences that failed when being checked
  # 2) Create archives for the rest of sequences and move them to result directory
  #    along with Rosetta .a3m files
  # 3) Update *db.index file(s)
  # 4) Update *submited_sequences.out file(s)
  (
    cd "${MICRO_DIR}"/slurm/
    find "${OUT_DIR}"/*/"${CM_DIR_NAME}"/*.fasta | xargs -d "\n" -P "${OMP_NUM_THREADS}" -n 1 ./organize_output
  )
  # TODO One more list:
  # problematic_sequences_C.out with unprocessed sequences

  # Update destination directory
  rsync -rzq "${RES_DIR}/tmp_output/${BATCH}/"* "${DEST_DIR}"

  # Remove output folder if empty
  if [ -z "$(ls -A "${OUT_DIR}")" ]; then
    rm -rf "${OUT_DIR}"
    test $? -eq 0 && printf "\n%s %s\n\n" "Removed:" "${OUT_DIR}"
  else
    printf "%s %s. %s\n" "Cannot remove" "${OUT_DIR}" "Folder not empty!"
  fi

  # Save SLURM_JOB_ID
  echo "${BATCH} # ${SLURM_JOB_ID} # $(date +%Y-%m-%d_%H:%M:%S)" >> "${RES_DIR}/job_mapping.out"
  showtime "Organizing output time"
fi

showtime "Total time" "${start_total}"
