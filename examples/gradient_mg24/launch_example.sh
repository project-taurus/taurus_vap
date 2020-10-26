#!/bin/bash

# For debugging purposes
#set -x

heredir=$(pwd)
outdir=$heredir/out

if [ ! -d $outdir ]; then mkdir $outdir; fi 

cp $heredir/auxiliaries/launch.sh .

echo "* Eta = factor for the gradient " > results_hfb
echo "* Mu = factor for the momentum " >> results_hfb
echo "* Iter = Number of iterations "  >> results_hfb
echo "* Ener = Energy at the last iteration " >> results_hfb
echo " " > results_hfb

echo "* Eta = factor for the gradient " > results_vap
echo "* Mu = factor for the momentum " >> results_vap
echo "* Iter = Number of iterations "  >> results_vap
echo "* Ener = Energy at the last iteration " >> results_vap
echo " " > results_vap


# Loop for HFB/VAPNP
for ggg in 1 5
do

  if [ $ggg -eq 1 ]; then 
    result=results_hfb
  else
    result=results_vap
  fi

  # Loop for mu (-2 and -1 are the automatic recipes)
  for imu in {-2..3}
  do

    if [ $imu -lt 0 ]; then
      mu=000 
      feta=0
      let xo=$imu+3
    elif [ $imu -eq 0 ]; then
      mu=000
      feta=1
      xo=0
    else
      let mu=300*$imu 
      feta=1
      xo=0
    fi 
 
    # Loop for eta
    for ieta in {0..30}      
    do

      cp $heredir/auxiliaries/template_input.txt temp
      
      # Defining eta in the correct format
      let eta=5*$ieta*$feta
    
      if [ $eta -eq 0 ] && [ $xo -eq 0 ]; then
        eta=1
      fi 

      if [ $eta -lt 10 ]; then
        eta=00$eta
      elif [ $eta -lt 100 ]; then
        eta=0$eta
      fi
      
      # launching the script 
      cat temp | sed s/ETA/$eta/ | sed s/MUU/$mu/ \
               | sed s/XGGG/$ggg/ | sed s/XO/$xo/ > input.txt
     
      bash launch.sh > output.$xo.$eta.$mu.$ggg.txt
              
      # Determine the last iteration 
      line0=$(grep -n "converged" output.$xo.$eta.$mu.$ggg.txt)
      if [ -z "$line0" ]; then
        line0=$(grep -n "reached" output.$xo.$eta.$mu.$ggg.txt)
      fi
      line0=$(echo $line0 | sed "s/:.*//" )
      let line0=$line0-2
      
      iter=$(awk 'NR=='$line0'{print $1}' output.$xo.$eta.$mu.$ggg.txt)
  
      # Consistent formatting
      if [ $iter -lt 10 ]; then
        iter=$(echo "   $iter")
      elif [ $iter -lt 100 ]; then
        iter=$(echo "  $iter")
      elif [ $iter -lt 1000 ]; then
        iter=$(echo " $iter")
      fi
      
      # Determine the energy
      line1=$(grep -m 1 -n "Full H   " output.$xo.$eta.$mu.$ggg.txt)
      line1=$(echo $line1 | sed "s/:.*//" )
  
      ener=$(awk 'NR=='$line1'{print $6}' output.$xo.$eta.$mu.$ggg.txt)

      if [ $imu -lt -0 ]; then       
        echo Recipe $xo: Iter = $iter  Ener = $ener >> $result
      elif [ $imu -eq 0 ]; then       
        if [ $ieta -eq 0 ]; then       
          echo " " >> $result
          echo "          Mu=0.0            Mu=0.3            Mu=0.6            Mu=0.9 "   >> $result
          echo " Eta   Iter     Ener     Iter     Ener     Iter     Ener     Iter     Ener" >> $result
        fi
        echo 0.$eta  " $iter"  " $ener"  I300$eta  E300$eta  \
                                       I600$eta  E600$eta  \
                                       I900$eta  E900$eta  >> $result
      else
        sed "s/I${mu}${eta}/ $iter/g" $result > temp2
        mv temp2 $result
        sed "s/E${mu}${eta}/ $ener/g" $result > temp2
        mv temp2 $result
      fi
  
      mv output.$xo.$eta.$mu.$ggg.txt $outdir/
 
      if [ $imu -lt 0 ]; then break; fi
 
    done
    echo xo=$xo mu=$mu Mphi=$ggg finished 
  done
done
    
#Find the minimum
   
echo The results can be found in the files: results_hfb  results_vap

# Clean up
rm -f launch.sh temp input.txt
