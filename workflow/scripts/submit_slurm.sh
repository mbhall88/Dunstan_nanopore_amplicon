#!/usr/bin/env bash
set -eu

module load singularity/3.8.5

JOB_NAME="snakemake_master_process."$(date "+%s")
LOG_DIR="logs"

if [[ ! -d "$LOG_DIR" ]]; then
    echo "Error: Log directory $LOG_DIR does not exist"
    exit 1
fi

MEMORY="1G"
TIME="3d"
THREADS=2
PROFILE="slurm.punim1637"
BINDS="/data/scratch/projects/punim1637/"
SINGULARITY_ARGS="--nv -B $BINDS"
CMD="snakemake --profile $PROFILE --rerun-incomplete --local-cores $THREADS $* --singularity-args '$SINGULARITY_ARGS'; chmod +r .snakemake/"

ssubmit -t "$TIME" -m "$MEMORY" -o "$LOG_DIR"/"$JOB_NAME".o \
    -e "$LOG_DIR"/"$JOB_NAME".e "$JOB_NAME" "$CMD" -- -c "$THREADS"
