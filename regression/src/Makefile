CPPFLAGS=`root-config --cflags`
LDFLAGS=-L${ROOTSYS}/lib -L${ROOFITSYS}/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic -lHtml -lMinuit -lRooFitCore -lRooFit -lTMVA -lRooStats -lFoam -lMathCore -lMathMore -g

all: RegVal2 AddRegVars ROC

RegVal2: RegVal2.C
	g++ $(CPPFLAGS) $(LDFLAGS) -I/${ROOFITSYS}/include -o $@ $^

AddRegVars: AddRegVars.C
	g++ $(CPPFLAGS) $(LDFLAGS) -o $@ $^

ROC: ROC.C
	g++ $(CPPFLAGS) $(LDFLAGS) -o $@ $^

CSplots2012: CSplots2012.C
	g++ $(CPPFLAGS) $(LDFLAGS) -o $@ $^

calcSys: calcSys.C
	g++ $(CPPFLAGS) $(LDFLAGS) -o $@ $^

combineSys: combineSys.C
	g++ $(CPPFLAGS) $(LDFLAGS) -o $@ $^

calcSigEff: calcSigEff.C
	g++ $(CPPFLAGS) $(LDFLAGS) -o $@ $^

compareSigShapes: compareSigShapes.C
	g++ $(CPPFLAGS) $(LDFLAGS) -o $@ $^

eventScan: eventScan.C
	g++ $(CPPFLAGS) $(LDFLAGS) -o $@ $^

makeSignifPlot2012: makeSignifPlot2012.C
	g++ $(CPPFLAGS) $(LDFLAGS) -o $@ $^

AddVariables: AddVariables.C
	g++ $(CPPFLAGS) $(LDFLAGS) -o $@ $^
