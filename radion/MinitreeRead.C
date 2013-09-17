#include "MinitreeOut.h"
#include <TH1.h>
#include <TString.h>
#include <TCanvas.h>
#include <TCut.h>

#include <vector>
#include <fstream>
#include <iomanip>
#include <string>
#include <iostream>
using namespace std;


//------------------------------ synchronization code ------------------------------//
void Synchronize(char* minitreeFileName, string fileout = "ggAna_eventlist.dump") {
  MiniTree m(minitreeFileName,1);
  m.setSynchVariables();
  vector<int> nEvtPerCat; nEvtPerCat.resize(9,0);
  vector<TH1F*> hmass;
  ofstream outf( fileout.c_str() );
  bool unblind = false;

  for( unsigned ih = 0 ; ih < nEvtPerCat.size(); ih++ ) {
    TString hname = TString::Format("hmass_cat_%d",ih);
    int nbin = 80;
    if( ih >= 4 ) nbin = 20;
    hmass.push_back( new TH1F( hname, hname, nbin,100,180) );
  }
  
  for( int ientry = 0; ientry < m.mtree_mainTree->GetEntries(); ientry++ ) {
    m.mtree_mainTree->GetEntry(ientry);
    if( m.mtree_catMva < 0 &&  ( m.mtree_pit1 < 1./3. || m.mtree_pit2 < 1./4.) ) continue;
    if( m.mtree_mass >= 100 && m.mtree_mass <= 180 &&m.mtree_catMva >= 0 ) {
      m.mtree_massResoRightVtx *= m.mtree_mass;
      m.mtree_massResoRandVtx  *= m.mtree_mass;
      m.mtree_massResoEng      *= m.mtree_mass;
      m.mtree_massResoAng      *= m.mtree_mass;
      for( unsigned i = 0 ; i < 2; i++ ) {
	 m.mtree_relResOverE[i] *=  m.mtree_eg[i];
	 m.mtree_relSmearing[i] *=  m.mtree_eg[i];
	 m.mtree_relSmearing[i] = sqrt( m.mtree_relSmearing[i] * m.mtree_relSmearing[i] +
					m.mtree_relResOverE[i] * m.mtree_relResOverE[i] );
      }      
      
      for( unsigned i = 0 ; i < m.sync_iName.size(); i++ )
	outf << m.sync_iName[i] << ":" << *m.sync_iVal[i]  << "  ";
      for( unsigned i = 0 ; i < m.sync_lName.size(); i++ )
	outf << m.sync_lName[i] << ":" << *m.sync_lVal[i]  << "  ";
      for( unsigned i = 0 ; i < m.sync_fName.size(); i++ )
	outf << m.sync_fName[i] << ":" << *m.sync_fVal[i]  << "  ";
      outf << endl;	
    }
    
    if( m.mtree_mass >= 100 && m.mtree_mass <= 180 ) {
      bool massblindedok =  m.mtree_mass < 115 && m.mtree_mass > 150;
      if( m.mtree_diphoMva < -0.05  ) continue;
      for ( unsigned ic = 0 ; ic < 4; ++ic )   
	if( m.mtree_catMva == int(ic) && m.mtree_pit1>1./3.&&m.mtree_pit2>0.25) {
	  nEvtPerCat[ic]++; 
	  if( unblind || massblindedok ) hmass[ic]->Fill(m.mtree_mass);
	}
      
      if( m.mtree_vbfTag == 1 && m.mtree_vbfCat == 0 && m.mtree_lepTag == 0 ) {
	if( unblind || massblindedok) hmass[4]->Fill(m.mtree_mass);
	nEvtPerCat[4]++;
      }
      if( m.mtree_vbfTag == 1 && m.mtree_vbfCat == 1 && m.mtree_lepTag == 0 ) {
	if( unblind || massblindedok ) hmass[5]->Fill(m.mtree_mass);
	nEvtPerCat[5]++;
      }
      if( m.mtree_lepTag == 1 && m.mtree_lepCat == 0 ) {
	if( unblind || massblindedok ) hmass[7]->Fill(m.mtree_mass);
	nEvtPerCat[7]++;
      }
      if( m.mtree_lepTag == 1 && m.mtree_lepCat == 1  ) {
	if( unblind || massblindedok ) hmass[6]->Fill(m.mtree_mass);
	nEvtPerCat[6]++;
      }
      if( m.mtree_metTag == 1  ) {
	if( unblind || massblindedok ) hmass[8]->Fill(m.mtree_mass);
	nEvtPerCat[8]++;
      }
    }  
  }
  
  
  vector<string> catName; catName.resize(nEvtPerCat.size(),"unknown");
  catName[0] = " incl (mva) category 1: ";
  catName[1] = " incl (mva) category 2: ";
  catName[2] = " incl (mva) category 3: ";
  catName[3] = " incl (mva) category 4: ";
  catName[4] = "  vbf (mva) category 5: ";
  catName[5] = "  vbf (mva) category 6: ";
  catName[6] = "            muon tag  : ";
  catName[7] = "            elec tag  : ";
  catName[8] = "             met tag  : ";
  
  TCanvas *can = new TCanvas("canmass","canmass",1000,1000);
  can->Divide(3,3);
  for( unsigned ic = 0 ; ic < catName.size(); ++ic )  {
    cout << catName[ic] << nEvtPerCat[ic] << endl;
    can->cd(ic+1); hmass[ic]->DrawCopy("e");
  }
  outf.close();
}
//------------------------------ end synchronization code -----------------------------//



void Spin2_dist(void) {
  vector<TString> sFiles;
  vector<TString> sNames;
  vector<int> sColors;
  vector<int> sStyles;

  vector<MiniTree*> sM;
  vector<TString> hNames, hTitles;
  vector<int> hNbins;
  vector<float> hXmin,hXmax;
  vector<vector<TH1F*> > h;
  TString dirMC   = "diphoton2012_mvaSel_cms53x_v2/mc/";
  TString dirData = "diphoton2012_mvaSel_cms53x_v2/data/";
  sFiles.push_back( dirMC + "job_hgg_gg2evenm.root_0.root"); sNames.push_back("Spin 2^{+}m" ); sColors.push_back(kBlue);  sStyles.push_back(2);
  sFiles.push_back( dirMC + "job_hgg_gg0even.root_0.root" ); sNames.push_back("Spin 0"      ); sColors.push_back(kRed);   sStyles.push_back(1);
  sFiles.push_back( dirData + "MiniTreeData.root"         ); sNames.push_back("Data"        ); sColors.push_back(kBlack); sStyles.push_back(1);

  hNames.push_back("pt"); hTitles.push_back("p_{T}   [GeV/c ]"); hNbins.push_back(50); hXmin.push_back(0); hXmax.push_back(100);
  hNames.push_back("abs(cThetaStar_CS)"); hTitles.push_back("|cos #theta^{*}| CS"); hNbins.push_back(20); hXmin.push_back(0); hXmax.push_back(1);
  hNames.push_back("max(abs(sceta[0]),abs(sceta[1]))"); hTitles.push_back("Most Forward #gamma #eta_{SC}"); hNbins.push_back(10); hXmin.push_back(0); hXmax.push_back(2.5);
  hNames.push_back("max(abs(sceta[0]),abs(sceta[1]))"); hTitles.push_back("Most Forward #gamma #eta_{SC}"); hNbins.push_back(10); hXmin.push_back(0); hXmax.push_back(2.5);
  hNames.push_back("max(abs(sceta[0]),abs(sceta[1]))"); hTitles.push_back("Most Forward #gamma #eta_{SC}"); hNbins.push_back(10); hXmin.push_back(0); hXmax.push_back(2.5);
  hNames.push_back("max(abs(sceta[0]),abs(sceta[1]))"); hTitles.push_back("Most Forward #gamma #eta_{SC}"); hNbins.push_back(10); hXmin.push_back(0); hXmax.push_back(2.5);
  hNames.push_back("max(abs(sceta[0]),abs(sceta[1]))"); hTitles.push_back("Most Forward #gamma #eta_{SC}"); hNbins.push_back(10); hXmin.push_back(0); hXmax.push_back(2.5);
  hNames.push_back("diphoMva"); hTitles.push_back("di-#gamma MVA"); hNbins.push_back(20); hXmin.push_back(-1.0); hXmax.push_back(1.0);
  hNames.push_back("diphoMva"); hTitles.push_back("di-#gamma MVA"); hNbins.push_back(20); hXmin.push_back(-1.0); hXmax.push_back(1.0);
  hNames.push_back("diphoMva"); hTitles.push_back("di-#gamma MVA"); hNbins.push_back(20); hXmin.push_back(-1.0); hXmax.push_back(1.0);
  hNames.push_back("diphoMva"); hTitles.push_back("di-#gamma MVA"); hNbins.push_back(20); hXmin.push_back(-1.0); hXmax.push_back(1.0);
  hNames.push_back("diphoMva"); hTitles.push_back("di-#gamma MVA"); hNbins.push_back(20); hXmin.push_back(-1.0); hXmax.push_back(1.0);


  vector<TCut> cut; TCut cut0 = "mass>120&&mass<130&&vbfTag==0&&lepTag==0&&metTag==0&&ptg[0]/mass>0.333&&ptg[1]/mass>0.25&&diphoMva>-0.8";
  cut.push_back( cut0 );
  cut.push_back( cut0 );
  cut.push_back( cut0 && "abs(cThetaStar_CS)<0.20" );
  cut.push_back( cut0 && "abs(cThetaStar_CS)>0.2&&abs(cThetaStar_CS)<0.4" );
  cut.push_back( cut0 && "abs(cThetaStar_CS)>0.4&&abs(cThetaStar_CS)<0.6" );
  cut.push_back( cut0 && "abs(cThetaStar_CS)>0.5&&abs(cThetaStar_CS)<0.8" );
  cut.push_back( cut0 && "abs(cThetaStar_CS)>0.8&&abs(cThetaStar_CS)<1.0" );
  cut.push_back( cut0 && "abs(cThetaStar_CS)<0.20" );
  cut.push_back( cut0 && "abs(cThetaStar_CS)>0.20&&abs(cThetaStar_CS)<0.4" );
  cut.push_back( cut0 && "abs(cThetaStar_CS)>0.40&&abs(cThetaStar_CS)<0.6" );
  cut.push_back( cut0 && "abs(cThetaStar_CS)>0.50&&abs(cThetaStar_CS)<0.8" );
  cut.push_back( cut0 && "abs(cThetaStar_CS)>0.8&&abs(cThetaStar_CS)<1.0" );

  for( unsigned im = 0 ; im < sFiles.size(); im++ ) {
    sM.push_back( new MiniTree(sFiles[im],1) );
    vector<TH1F*> htmp;
 
    for( unsigned iv = 0 ; iv < hNames.size(); iv++ ) {
      TString htmpName =  TString::Format( "htoto_%d_%d",iv,im);
      htmp.push_back( new TH1F( htmpName, hNames[iv],hNbins[iv],hXmin[iv],hXmax[iv]));
      htmp[iv]->SetLineColor( sColors[im] );
      htmp[iv]->SetLineStyle( sStyles[im] );
      htmp[iv]->SetLineWidth( 2 );
      htmp[iv]->GetXaxis()->SetTitle( hTitles[iv] );
      htmp[iv]->GetXaxis()->CenterTitle();
      htmp[iv]->Sumw2();
      sM[im]->mtree_mainTree->Draw( hNames[iv] + ">>" + htmpName, "wei"*cut[iv], "goff" );
 
    }
    h.push_back( htmp);
  }
  TString pdfout = "Spin2Output.pdf";
  vector<TCanvas*> can;
  for( unsigned iv = 0 ; iv < hNames.size(); iv++ ) {
    TString ctmpName =  TString::Format( "ctoto_%d",iv);
    can.push_back( new TCanvas( ctmpName, ctmpName, 700,600 ) );
    float ymin = 0;
    float ymax = -999;
    for( unsigned im = 0 ; im < sFiles.size(); im++ ) {
      h[im][iv]->Scale( 1./h[im][iv]->Integral() );
      if( h[im][iv]->GetMaximum() > ymax ) ymax = h[im][iv]->GetMaximum(); 
    }
    for( unsigned im = 0 ; im < sFiles.size(); im++ ) {
      TString opt = "hist same";
      if( im == 0 ) opt = "hist";
      if( sNames[im] == "Data" || sNames[im] == "data" ) opt = "e same";
      can[iv]->cd();
      h[im][iv]->SetMinimum(ymin);
      h[im][iv]->SetMaximum(ymax*1.20);
      h[im][iv]->DrawCopy(opt);
    }
    if( iv == 0 ) can[iv]->Print( pdfout + "[" );
    can[iv]->Print( pdfout );
  }
  can[0]->Print( pdfout + "]" );

}
