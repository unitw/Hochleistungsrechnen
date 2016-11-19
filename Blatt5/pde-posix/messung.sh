#!/bin/bash

#SBATCH -N 1
#SBATCH -c 12
#SBATCH --ntasks-per-node 1
#SBATCH --output=results/messung.out

run_output="results/messung.log"
TIMEFORMAT='%E'

NR_OF_RUNS=${1:-3}
START_THREADS=${2:-1}

# for (( i = 99; i >= 0; i-- )); do
# 	src_out="results/messung.out"
# 	src_log="results/messung.log"
	
# 	if (($i > 0)); then
# 		src_out="results/messung"$i".out"
# 		src_log="results/messung"$i".log"
# 	fi

# 	if [[ -f $src_out ]]; then
# 		mv $src_out "results/messung$((i+1)).out"
# 		mv $src_log "results/messung$((i+1)).log"
# 	fi
# done

for (( i = $START_THREADS; i <= 12; i++ )); do
	echo "$i Thread(s):" >> $run_output
	echo "$i Thread(s): "

	for (( j = 1; j <= $NR_OF_RUNS; j++ )); do
		time_for_run=$(time (./partdiff-posix $i 2 512 2 2 1024 >> $run_output 2>&1) 2>&1)
		echo Lauf $j: $time_for_run seconds
	done
done
