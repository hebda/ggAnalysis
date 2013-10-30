#!/bin/bash

masses=(300 500 700 1000)

#for i in {0..215}; do
for i in {0..21}; do
    filename=submit/submit_$i
    echo "#!/bin/sh" > $filename

    echo "cd $PWD" >> $filename
    echo "export SCRAM_ARCH=slc5_amd64_gcc462" >> $filename
    echo "eval \`scramv1 runtime -sh\`" >> $filename

    for j in {0..9}; do
	optimIndex=$((i*10+j))
	if [ $optimIndex -gt 215 ]; then
	    continue 
	fi

	echo "root -l -b -q \"Regression/TMVARegression.C($optimIndex,0)\"" >> $filename
	echo "root -l -b -q \"Regression/TMVARegression.C($optimIndex,1)\"" >> $filename

	for mass in "${masses[@]}"; do
	    echo "./src/AddRegVars $optimIndex $mass 1" >> $filename
	    echo "./src/RegVal2 $optimIndex $mass 1" >> $filename
	done
	
	echo "rm temp_${optimIndex}.root" >> $filename
    done
    chmod 755 $filename
    bsub -q 1nd -J job$i $PWD/$filename
    #source $PWD/$filename
done
