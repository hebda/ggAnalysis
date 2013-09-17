if [ ! -d VertexAnalysis ]; then
    ln -s 53x_VertexAnalysis VertexAnalysis
fi

ROOTSYS=/afs/cern.ch/sw/lcg/app/releases/ROOT/5.32.00/x86_64-slc5-gcc46-opt/root/
XRDCP=/afs/cern.ch/sw/lcg/external/xrootd/3.1.0p2/x86_64-slc5-gcc46-opt/
. /afs/cern.ch/sw/lcg/external/gcc/4.6.2/x86_64-slc5/setup.sh
. $ROOTSYS/bin/thisroot.sh
LD_LIBRARY_PATH=`pwd`/lib:$XRDCP/lib64:$LD_LIBRARY_PATH:.
LD_LIBRARY_PATH=/afs/cern.ch/cms/slc5_amd64_gcc462/external/boost/1.47.0/lib/:$LD_LIBRARY_PATH
LD_LIBRARY_PATH=/afs/cern.ch/cms/slc5_amd64_gcc462/cms/cmssw/CMSSW_5_2_5/lib/slc5_amd64_gcc462/:$LD_LIBRARY_PATH
LD_LIBRARY_PATH=/afs/cern.ch/cms/slc5_amd64_gcc462/cms/cmssw-patch/CMSSW_5_2_5_patch1/external/slc5_amd64_gcc462/lib/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH
export PATH=`pwd`/bin/:$XRDCP/bin/:$PATH

### environment with gcc 4.6.2 