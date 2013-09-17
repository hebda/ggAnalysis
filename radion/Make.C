void Make(int arg1, char* arg2, char* arg3){
  gROOT->ProcessLine(".L VertexAnalysis/lib/libh2gglobeVertexAnalysis.so");
  gROOT->ProcessLine(".L JetMETObjects/lib/libCondFormatsJetMETObjects.so");
  char command[400];
  sprintf(command,".x runAna.C+(%i,\"%s\",\"%s\")",arg1,arg2,arg3);
  gROOT->ProcessLine(command);
}

