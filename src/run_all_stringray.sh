#!/bin/bash

ROOT=$HOME/scratch/cwilks
MAX_THREADS=`grep 'processor\s*:' /proc/cpuinfo | wc -l`
#MAX_THREADS=10
DR=`dirname $0`
READS="/home/cwilks3/bowtie2/seqs_by_100.fq"

echo "Max threads = $MAX_THREADS"
MAX_THREADS=120
run_th () {
#for mode in very-fast fast sensitive very-sensitive ; do
  for ((t=$MAX_THREADS; t<=$MAX_THREADS; t=${t}+10)); do
    cmd="./${1} -t $t -M -d minikraken/database.kdb -i minikraken/database.idx -n minikraken/taxonomy/nodes.dmp -f "
    data_file="$HOME/results/elephant6/raw/${2}${t}.out"
    #data_file="/dev/shm/raw/$mode/${2}${t}.out"
    echo "mode: $mode, threads: $t"
    echo "Concatenating input reads"
    cp $READS /tmp/.run_all_reads.fq
    #cp $READS /dev/shm/.run_all_reads.fq
    for ((i=1;i<${t};i++)); do cat $READS >> /tmp/.run_all_reads.fq; done
    #for ((i=1;i<${t};i++)); do cat $READS >> /dev/shm/.run_all_reads.fq; done
    # make sure input and output are on a local filesystem, not NFS
    echo "Running kraken $cmd /tmp/.run_all_reads.fq"
    $cmd /tmp/.run_all_reads.fq | grep "thread:" > $data_file
  done
}

# Normal (all synchronization enabled), no TBB
#if [ ! -f "classify" ] ; then
#  git checkout tbb
#  make -j $MAX_THREADS
#fi
#run_th classify krakentbb120
#exit 0

run_th classify kraken_no_io_tbbreal_no_kraken_output2_no_classifyloop
exit 0

make clean    
if [ ! -f "classify" ] ; then
  git checkout no_io_tbb
  make -j $MAX_THREADS
  #mv lighter lighter-no_io
fi
run_th classify kraken_no_io_tbb2
exit 0

