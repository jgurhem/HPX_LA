#!/bin/bash

timewall="00:30:00"
nodes=$1
cores=$(($nodes * 16))
datasize=$2
blocks=$3
blocksize=$(($datasize / $blocks))

sed "s/TIMEWALL/${timewall}/;
     s/NNODES/${nodes}/;
     s/DATASIZE/${datasize}/g;
     s/NBLOCKS/${blocks}/;
     s/BLOCKSIZE/${blocksize}/;
     s/NCORES/${cores}/" scripts/poincare/launch_template.sh > submit_s${datasize}_n${nodes}_b${blocks}.sh

. scripts/poincare/load_env.sh
llsubmit submit_s${datasize}_n${nodes}_b${blocks}.sh
