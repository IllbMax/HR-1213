#!/bin/bash

########################
#####  Messrahmen  #####
########################
# Ausgabeverzeichnis und Wiederholung der Messungen bestimmen  

outputdir="output"      # Verzeichnis, indem die Ausgabe gespeichert werden
job_script="job_script" # Name des Jobscripts
coresPerNode=12 # Anzahl der verfügbaren Cores auf einem Node
messzahl=5     # Anzahl der Messungen pro Knoten und pro Thread und Interlines

########################
##### Messparamter #####
########################
# Aufrufparameter (und deren Increment) vom PDE-Programm bestimmen

minthreads=1  # unter Grenze der Threaditeration 
maxthreads=12   # obere Grenze der Threaditeration

mininterlines=512    # unter Grenze der Interlinesiteration
maxinterlines=512    # obere Grenze der Interlinesiteration
deltainterlines=10   # Schrittweite der Interlinesiteration  
interlineIncrMeth=0  # Schrittincrementart: 
                     # 0: additiv, 1: multiplikativ

threads=$minthreads

method=2             # Berechnungsmethode 
methodName=""        # 1: Gauss Seidel,  2: Jacobi
if [ $method -eq 1 ] 
  then methodName="GaussSeidel"
  else methodName="Jacobi"
fi

interlines=$mininterlines

inf_func=2           # Stroerfunktion:  1: f_0,  2: sin
termination=2        # Art des Abbruchs: 1: Genauigkeit,  2: Iterationen
termparam=1024       # Parameter des Abbruchs  


########################################
#####   Jobberechnung / vergabe    #####  
########################################

jobcount=0   # zählt die Anzahl der gestarteten Jobs
while [ $threads -le $maxthreads ]
do
  interlines=$mininterlines
  while [ $interlines -le $maxinterlines ]
  do
    i=0
    mleft=$messzahl  # Anzal der noch fehlenden Messungen
    ppn=`expr $coresPerNode / $threads ` # Anzahl der pde Anwendungen, die auf einem Node passen
    if [ $ppn -le 0 ]
      then 
        echo "Es können nicht mehr Threads gestartet werden als Cores auf dem Node vorhanden sind."
        exit 2
    fi
    while [ $mleft -gt 0 ]
    do
      nopt=$ppn   # pde Anwedenungen, die gestartet werden sollen: -n nopt 
      if [ $mleft -lt $ppn ] # es passen alle übrigen Messungen auf einen Node
        then nopt=$mleft      
      fi
      mleft=`expr $mleft - $nopt`
      s=_ # sepatator
      threadsname=$(printf "%02d" $threads)
      iname=$(printf "%02d" $i)
      # Nomenklatur:  threads_method_interlines_infFunc_Termination_Termparam#messung.out
      name="$outputdir/$threadsname$s$methodName$s$interlines$s$inf_func$s$termination$s$termparam#$iname.out"   
      sbatch --output $name -n $nopt $job_script $threads $method $interlines $inf_func $termination $termparam 
      #echo --output $name -n $nopt $job_script $threads $method $interlines $inf_func $termination $termparam
                      
      jobcount=`expr $jobcount + 1`
      i=`expr $i + 1`      
    done
    
    # Increment Interlines
    if [ $interlineIncrMeth -eq 0 ]
      then interlines=`expr $interlines + $deltainterlines`   
      else interlines=`expr $interlines \* $deltainterlines`      
    fi
  done
  # Increment Threads
  threads=`expr $threads + 1`
done

echo "Es wurden $jobcount Jobs gestartet.\n"
