#!/bin/bash
#SBATCH --job-name=banff
#SBATCH --partition=gpu-a100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --array=0-9
#SBATCH --output=slurm-%A_%a.out
set -euo pipefail

# Submit from Feedforward_network. TASK is the same name accepted locally.
# A comma-separated partition list requests ONE allocation, not duplicate jobs.
export BANFF_TASK="${BANFF_TASK:-}"
if [[ -z "$BANFF_TASK" && "${BANFF_MODE:-train}" != publication ]]; then
    echo 'Set BANFF_TASK (for example toyota or Lorenz), or BANFF_MODE=publication.' >&2
    exit 2
fi
package_root="${BANFF_ROOT:-${SLURM_SUBMIT_DIR}}"
test -f "$package_root/setup.m"
module load "${MATLAB_MODULE:-matlab/r2023a}"
export BANFF_ROOT="$package_root"
export BANFF_SEED="${BANFF_SEED:-${SLURM_ARRAY_TASK_ID:-0}}"
export BANFF_MODE="${BANFF_MODE:-train}"
matlab -batch "addpath(fullfile(getenv('BANFF_ROOT'),'slurm')); run_job(getenv('BANFF_TASK'),str2double(getenv('BANFF_SEED')),getenv('BANFF_MODE'));"
