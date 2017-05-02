/********************************************
 * NT
 * 05/02/2017
 * Input: root file with phase space run, radiative OFF & wo any cut
 * Output: real solid angle and the 2D shape of theta vs phi at target
 * The cut applied is from sieve slit runs, the sieve hole center
 * from reconstructed quantities
 *********************************************/
#include <iomanip>
#include <vector>
#include "TMath.h"
#include <cassert>
#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TROOT.h"
#include "TSystem.h"

void rec_phase_space()
{
  /* Limit for the table */
  /* Define histograms */

  TH2F * hthph     = new TH2F("hthph","",200,-20,20,200,20,60);
  TH2F * hthphCut  = new TH2F("hthphCut","",200,-20,20,200,20,60);  



  /* start each root file */

  int bcsa=0;
  int scsa=0; 
  const int DIM = 4;

  double c_phiTg[DIM]={0.006951,0.000486,0.004929,0.012954};
    //{0.006959,0.001251,0.005041,0.011853};// for 0dp from run 1866
  double c_thetaTg[DIM]={0.046409,0.046302,0.033442,0.036281};
    //{0.045869,0.045380,0.032330,0.035429};// for 0dp from run 1866

  double edge[DIM];
  
  char command[120];
  sprintf(command,"ls *.root |cat > rootlist");
  gSystem->Exec(command);
  char rootfile[100];
  ifstream rootlist("rootlist",ios_base::in);


  while(true){
    rootlist>>rootfile;
    if(rootlist.eof()) break;
  
    Float_t phfoc,yfoc,phitg,thetatg,deltatg,xfoc,thfoc,wmm;
    TFile *rootlist1  = new TFile(rootfile);//TFile::Open(rootfile);
    // TFile *fSim  = new TFile("./rootfiles/center_foil_phase_space_wo_cut.root");

    TTree *ntup = (TTree*)gROOT->FindObject("h1");
    //TTree *ntup = (TTree*)fSim->Get("h1");
    Int_t n_entries = (Int_t)ntup->GetEntries();
    ntup->SetBranchAddress("phfoc",&(phfoc));
    ntup->SetBranchAddress("yfoc",&(yfoc));
    ntup->SetBranchAddress("xfoc",&(xfoc));
    ntup->SetBranchAddress("thfoc",&(thfoc));
    ntup->SetBranchAddress("deltat",&(deltatg));
    ntup->SetBranchAddress("thetat",&(thetatg));
    ntup->SetBranchAddress("phit",&(phitg));
    ntup->SetBranchAddress("wmm",&(wmm));
    cout<<"Number of events= "<<n_entries<<endl;
    for(int kk=0;kk<n_entries-1;kk++){
      ntup->GetEntry(kk);	  
      if(yfoc<=0.03&&yfoc>=-0.02
	 &&phfoc<=0.05&&phfoc>=-0.05
	 &&thfoc<=0.03&&thfoc>=-0.03&&wmm<=0.004&&wmm>=-0.001)
	{
	  double ph_rec =0.001692 -0.012293*xfoc -0.007802*yfoc-0.40254*thfoc +0.627065*phfoc
	    -5.927*yfoc*phfoc -6.172*yfoc*thfoc
	    -11.326*yfoc*yfoc +6.111*phfoc*phfoc;
	  double th_rec =0.037990 -0.004055*xfoc +0.849883*yfoc +0.07456*thfoc -0.163464*phfoc;
	  double slope, intercept;
	  for(int i=0; i<DIM; i++)
	    {
	      if(i==(DIM-1))
		{
		  slope =(c_thetaTg[0]-c_thetaTg[i])/(c_phiTg[0]-c_phiTg[i]);
		} else 
		{
		  slope= (c_thetaTg[i+1]-c_thetaTg[i])/(c_phiTg[i+1]-c_phiTg[i]);
		}
	      intercept = c_thetaTg[i] - slope*c_phiTg[i];
	      edge[i] = slope*ph_rec + intercept;
	    }
	  if(th_rec<=edge[0]&&th_rec>=edge[1]&&th_rec>=edge[2]&&th_rec<=edge[3])
	    {
	      hthph->Fill(phitg,thetatg);
	      bcsa++;
	      if(thetatg>=35&&thetatg<=45&&phitg>=4&&phitg<=6){
		scsa++;
		hthphCut->Fill(phitg,thetatg);
	      }
	    }// cut at reconstructed variables
	}// cut on focal plane variable limits
    }// loop events
    delete rootlist1;
  }
  TCanvas * c = new TCanvas("c","",800,600);
  c->Clear();
  hthph->Draw("colz");
  
  cout<<"**************"<<endl;
  cout<<"Real solid angle= "<<bcsa*20./scsa<<endl;
    
}
