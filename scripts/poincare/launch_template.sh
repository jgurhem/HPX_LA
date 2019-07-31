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

. scripts/poincare/load_env.sh
rm -rf core.*
bash scripts/run_mpi.sh Poincare blockLU build/lu_tiled_dist NNODES NCORES NBLOCKS DATASIZE BLOCKSIZE test_DATASIZE.json
