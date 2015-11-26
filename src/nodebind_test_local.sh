#!/bin/bash
ROOT=$HOME/scratch/cwilks
MAX_THREADS=`grep 'processor\s*:' /proc/cpuinfo | wc -l`
#MAX_THREADS=10
DR=`dirname $0`
READS="/home/cwilks3/bowtie2/seqs_by_100.fq"

echo "Max threads = $MAX_THREADS"
#MAX_THREADS=120

NODES=`numactl --hardware | grep available: | awk '{print $2}'`

echo "Max threads = $MAX_THREADS"
echo "NUMA nodes = $NODES"

# This will spawn a bowtie process for each NUMA node and will continue to
# increase the number of threads for each iteration. Therefore this function
# will run n*t threads at any step, where n is the total number of NUMA nodes
# and t is the thread step from 1 to MAX_THREADS.
run_th () {
    local TIMING_DIR="$HOME/results/elephant6/raw"
    declare -a INPUT_READS=()
    declare -a OUTPUT_SAMFILE=()
    local timing_file
    local cmd
    #mkdir -p $TIMING_DIR
    #h=0
    #for mode in "-v 1" "-n" "-n -e50000" "-v 2" ; do
    #for mode in "-k 30" ; do
      #h=$((h + 1))
      #for ((t=$NODES; t<=$MAX_THREADS; t+=$NODES)); do
      for ((t=$MAX_THREADS; t<=$MAX_THREADS; t+=$NODES)); do
        ((nthread=$t/$NODES))
        #cmd="./${1} $cmd_tmpl --$mode -U "
        #cmd="./${1} -p $t $mode --best -t $HG19_INDEX "
        #cmd="./classify -M -d minikraken/database.kdb -i minikraken/database.idx -n minikraken/taxonomy/nodes.dmp "
        cmd="./${1} -t $t -M -d minikraken/database.kdb -i minikraken/database.idx -n minikraken/taxonomy/nodes.dmp -f "
        data_file="$HOME/results/elephant6/raw/${2}${t}.out"
        
        echo "mode: $mode, threads: $t"
        echo "Concatenating input reads"
        for ((i=0; i<$NODES; i++)); do
            INPUT_READS[$i]=$(mktemp -p /tmp bowtie_test_XXXX_${i}.fq)
            #OUTPUT_SAMFILE[$i]=$(mktemp -p /tmp bowtie_test_XXXX_${i}.sam)
            for ((j=0; j<$nthread; j++)); do
                cat $READS >> ${INPUT_READS[$i]}
            done
        done
        # make sure input and output are on a local filesystem, not NFS
        echo "Running bowtie"
        pids=""
        for ((i=0; i<$NODES; i++)); do
            timing_file="${TIMING_DIR}/${2}${t}_${i}.out"
            #mkdir -p "${TIMING_DIR}/$mode"
            #(numactl -N $i $cmd ${INPUT_READS[$i]} -p $nthread -S ${OUTPUT_SAMFILE[$i]} | grep "thread:" > $timing_file) &
            echo "(numactl -m $i -N $i $cmd -f ${INPUT_READS[$i]} -t $nthread | grep "thread:" > $timing_file)"
            (numactl -m $i -N $i $cmd -f ${INPUT_READS[$i]} -t $nthread | grep "thread:" > $timing_file) &
            echo "  spawned node $i process with pid $!"
            pids="$! $pids"
        done
        for pid in $pids ; do
            echo "  waiting for PID $pid"
            wait $pid
        done
        # cleanup
        for fs in "${INPUT_READS[@]}"; do
            echo " removing $fs"
            rm $fs
        done
        for fs in "${OUTPUT_SAMFILE[@]}"; do
            echo "removing $fs"
            rm $fs
        done
      done
}
run_th classify kraken_no_io_tbbreal_mem_per_node_no_kraken_output
exit 0

# Normal (all synchronization enabled), no TBB
if [ ! -f "classify" ] ; then
  git checkout master
  rm -f classify
  make -j $MAX_THREADS
fi
run_th classify kraken_master
exit 0
