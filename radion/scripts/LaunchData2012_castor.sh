#!/bin/bash


usage() {
    echo "`basename $0` -f filelist -c configfile -n njobs(default=1) [-i]"
    echo "  -f filelist: filelist of input (e.g. scripts/ExampleDataDiphoton8TeVSkimEOS.list)"
    echo "  -c configfile: config file (e.g. configFilesDir/mvaAnalysisTest.config)"
    echo "  -i: interactive mode (by default job sent to batch system)"
}

mainConfigFile=notSpecified
inFileList=notSpecified
nJobs=1
interactive=0

if ! options=$( getopt -o hc:f:n:i -l help -- "$@" )
then
    # something went wrong, getopt will put out an error message for us
    exit 1
fi

eval set -- $options
while [ $# -gt 0 ]; do
    case "$1" in
	-h | --help) usage; exit 0;;
	-c) mainConfigFile=$2; shift;;
	-f) inFileList=$2; shift;;
	-n) nJobs=$2; shift;;
	-i) interactive=1;;
	(--) shift; break;;
	(-*) echo "$0: error - unrecognized option $1" 1>&2; usage >> /dev/stderr; exit 1;;
	(*)  break;;
    esac
    shift
done

###### print and check configuration
echo "================================="
echo " - config file: " $mainConfigFile
echo " - file list  : " $inFileList
echo " - njobs      : " $nJobs
echo " FIX ME: user has to verify that nJobs x nEvtPerJobs (in config) > nTot"

if [ $mainConfigFile == "notSpecified" ]; then
    echo no good: $mainConfigFile
    usage; 
    exit 1;
elif [ ! -f $mainConfigFile ]; then
    echo " config file does not exist... bailout"
    exit 1;
fi

if [ $inFileList == "notSpecified" ]; then
    usage; 
    exit 1;
elif [ ! -f $inFileList ]; then
    echo " input filelist does not exist... bailout"
    exit 1;
fi
echo "================================="

#exit 0;

optionsub="-q 1nd "


dirafs=`pwd`
dirscript=tmpBatchOut/
eosPref=root://castorcms/


config() {
    file=$1
    script=${2}
    ijob=$3
    config=$4
    
    filelist="unknown"
    if [ $# -ge 5 ]; then
	filelist=${5}
    fi

    exe="runAna.C+(${ijob},\"$file\",\"${config}\")"
    echo "$exe"

    cat > $script<<EOF 
#!/bin/bash
cd $dirafs
source scripts/env532.sh
echo "castor: $dircastor" 
echo "ROOTSYS: \$ROOTSYS"  
gcccom="\`which gcc\`"
echo "gcc:" \$gcccom
echo "where am I:\`pwd\`"
echo  root -b -q rootlogon.C '$exe'   
root -b -q '$exe'         
EOF
    chmod +x $script
}


mkdir -p $dirscript

configfile=$mainConfigFile
let "nJobs=nJobs-1"
for f in $( cat $inFileList ); do
    for ijob in  $( seq 0 1 $nJobs ); do
#	echo " test file $f exist"
#	cmsLs $f
#	retval=$?
#	if [ $retval -eq 0 ]; then
	filename=`basename $f`
	script=script$$
	script=${script}_${filename/root/sh}_${ijob}
	cd $dirscript/
	config ${eosPref}${f} $script $ijob  $configfile
	echo "-> created script: "$script
	if [ $interactive -eq 1 ]; then
	    source $script
	else
	    bsub $optionsub $script 
	fi
	cd -
#	fi
    done
done


