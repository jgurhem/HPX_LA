#!/bin/bash

if [ $# -ne 9 ]
then
        echo "Wrong number of parameters"
        exit 1
fi

echo bash $0 $*

machine=${1}
app=${2}
exe=${3}
nb_nodes=${4}
nb_cores=${5}
nb_blocks=${6}
datasize=${7}
blocksize=${8}
res_file=${9}

d1=$(date +%s.%N)
t=$(mpirun -n $nb_nodes build/lu_tiled_dist --T $nb_blocks --N $blocksize -l $nb_nodes)
res_=$?
d2=$(date +%s.%N)
t_app=$(echo "$d2 $d1" | awk '{printf "%f", $1 - $2}')
echo ${t[@]}

success="false"
if [ $res_ -eq 0 ]
then
  success="true"
else
  t=0
fi

{
cat << EOF
{"machine":"$machine",\
"nb_cores":"$nb_cores",\
"nb_nodes":"$nb_nodes",\
"test":"$app",\
"lang":"HPX",\
"nb_blocks":"$nb_blocks",\
"blocksize":"$blocksize",\
"datasize":"$datasize",\
"date":"$(date +%Y%m%d-%H%M%S)",\
"time_app":"$t_app",\
"time_calc":"$t",\
"success":"$success"}
EOF
} | tee -a $res_file

