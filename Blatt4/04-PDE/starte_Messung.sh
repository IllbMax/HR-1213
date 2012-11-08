#!/bin/bash

messzahl=2     # Anzahl der Messungen pro Thread und Interline

minthreads=1   # unter Grenze der Threaditeration 
maxthreads=2   # obere Grenze der Threaditeration

mininterlines=512   # unter Grenze der Interlinesiteration
maxinterlines=512   # obere Grenze der Interlinesiteration
deltainterlines=10  # Schrittweite der Interlinesiteration  


threads=$minthreads

method=2
methodName=""
if [ $method -eq 1 ] 
  then methodName="GaussSeidel"
  else methodName="Jacobi"
fi

interlines=$mininterlines

inf_func=2
termination=2
termparam=500



while [ $threads -le $maxthreads ]
do
  interlines=$mininterlines
  while [ $interlines -le $maxinterlines ]
  do
    i=0
    while [ $i -lt $messzahl ]
    do
      s=_ # sepatator
      threadsname=$(printf "%02d" $threads)
      iname=$(printf "%02d" $i)
      # threads_method_interlines_infFunc_Termination_Termparam#messung.out
      name="output/$threadsname$s$methodName$s$interlines$s$inf_func$s$termination$s$termparam#$iname.out"   
      sbatch --output $name job_script $threads $method $interlines $inf_func $termination $termparam 
      #echo --output $name job_script $threads $method $interlines $inf_func $termination $termparam
      i=`expr $i + 1`
    done
    interlines=`expr $interlines + $deltainterlines`
  done
  threads=`expr $threads + 1`
done

