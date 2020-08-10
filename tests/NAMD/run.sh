for i in `seq 1 4`; do
  echo $i
  julia error2.jl &
done
wait_for.tcl leandro julia 1 
