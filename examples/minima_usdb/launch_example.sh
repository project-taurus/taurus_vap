#!/bin/bash

# For debugging purposes
#set -x

heredir=$(pwd)
outdir=$heredir/out

if [ ! -d $outdir ]; then mkdir $outdir; fi 

cp $heredir/auxiliaries/launch.sh .

# Get the line numbers and then the values
echo "* Z = number of protons " > results
echo "* N = number of neutrons" >> results
echo "* M = number of discretization points in the PNP " >> results
echo "  if M=1 : HFB calculations                      " >> results
echo "  if M>1 : VAPNP calculations (the values are the expectation" >> results
echo "  values of the projected state) " >> results
echo "* Energy = total energy" >> results
echo "* Epair = particle-particle part of the energy" >> results
echo "* Beta = triaxial beta deformation"   >> results
echo "* Gamma = triaxial gamma deformation" >> results
echo "  (the orientation of the nucleus is not fixed such that you"    >> results
echo "  might obtain same results for Gamma=X<60 and Gamma=n*60 +- X)" >> results
echo " " >> results
echo "  Z      N    M    Energy      Epair      Beta      Gamma  " >> results

# Loop for nuclei
for nucl in 1 2 3
do

  cp $heredir/auxiliaries/template_input.txt temp

  if [ $nucl -eq 1 ]; then
    zzz=4.00
    nnn=4.00
    bbb=0
    bb1=0 
    bb2=0 
    sed '19,20d' temp > temp2
    mv temp2 temp
  elif [ $nucl -eq 2 ]; then 
    zzz=4.00
    nnn=5.00
    bbb=1
    bb1=13
    bb2=0 
    sed '20d' temp > temp2
    mv temp2 temp
  else
    zzz=5.00
    nnn=5.00
    bbb=2
    bb1=1 
    bb2=13
  fi

  # Loop for HFB/VAPNP
  for ggg in 1 7
  do

    enermin=666.0
    epaimin=666.0
    betamin=666.0
    gammmin=666.0

    cat temp | sed s/XZZZ/$zzz/ | sed s/XNNN/$nnn/ \
             | sed s/XGGG/$ggg/ | sed s/XBBB/$bbb/ \
             | sed s/XBB1/$bb1/ | sed s/XBB2/$bb2/ > input.txt
    

    # Loop for the number of runs
    for i in {1..10}
    do

      bash launch.sh > output.$zzz.$nnn.$ggg.$i.txt

      # Get the line numbers and then the values
      line0=$(grep -m 1 -n "Full H   " output.$zzz.$nnn.$ggg.$i.txt)
      line0=$(echo $line0 | sed "s/:.*//" )
      let line1=$line0-2
      let line2=$line0+24
      let line3=$line0+25

      enertmp=$(awk 'NR=='$line0'{print $6}' output.$zzz.$nnn.$ggg.$i.txt)
      epaitmp=$(awk 'NR=='$line1'{print $6}' output.$zzz.$nnn.$ggg.$i.txt)
      betatmp=$(awk 'NR=='$line2'{print $2}' output.$zzz.$nnn.$ggg.$i.txt)
      gammtmp=$(awk 'NR=='$line3'{print $2}' output.$zzz.$nnn.$ggg.$i.txt)

      # Store if it is the minimum
      if (( $(echo "$enertmp < $enermin" | bc -l) )); then
        enermin=$enertmp
        epaimin=$epaitmp
        betamin=$betatmp
        gammmin=$gammtmp
      fi 
      
      mv output.$zzz.$nnn.$ggg.$i.txt $outdir/
    done
 
    echo $zzz ' '$nnn ' '$ggg  ' '$enermin ' '$epaimin ' '$betamin ' '$gammmin >> results

    echo Z=$zzz N=$nnn M=$ggg finished 
  done
done

echo The results can be found in the file: results

# Clean up
rm -f launch.sh temp input.txt
