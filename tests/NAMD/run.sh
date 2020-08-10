for i in `seq 1 16`; do
  echo $i
  time julia -t $i test_parallel.jl >& /dev/null
done
wait_for.tcl leandro julia 1
