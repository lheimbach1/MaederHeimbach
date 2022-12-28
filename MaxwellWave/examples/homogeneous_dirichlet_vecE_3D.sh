#!/bin/bash -l
#SBATCH --job-name="homogeneous_dirichlet_vecE_3D"
#SBATCH --output=homogeneous_dirichlet_vecE_3D.o
#SBATCH --error=homogeneous_dirichlet_vecE_3D.e
#SBATCH --time=12:00:00
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=1
#SBATCH --partition=normal
#SBATCH --constraint=gpu
#SBATCH --account class04

module load daint-gpu
module load Julia/1.7.2-CrayGNU-21.09-cuda

export MPICH_RDMA_ENABLED_CUDA=1
export IGG_CUDAAWARE_MPI=1

srun -n8 bash -c 'LD_PRELOAD="/usr/lib64/libcuda.so:/usr/local/cuda/lib64/libcudart.so" julia --project --check-bounds=no -O3 examples/homogeneous_dirichlet_vecE_3D.jl'