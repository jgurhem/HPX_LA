#@ class            = clallmds+
#@ job_name         = hpx-run
#@ total_tasks      = NNODES
#@ node             = NNODES
#@ wall_clock_limit = TIMEWALL
#@ output           = $(job_name).$(jobid).log
#@ error            = $(job_name).$(jobid).err
#@ job_type         = mpich
#@ environment      = COPY_ALL
#@ node_usage       = not_shared
#@ queue
#

module purge
module load gnu/7.3.0 cmake/3.14.1-gnu54 openmpi/2.1.2_intel15.0.0_tm
rm -rf core.*
bash scripts/run_mpi.sh Poincare blockLU build/lu_tiled_dist NNODES NCORES NBLOCKS DATASIZE BLOCKSIZE test_DATASIZE.json
