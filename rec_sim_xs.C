/* Nguyen Ton
 * 05/02/2017
 * Calculate carbon cross section using a boundary define by apply a 2D cut on 
 * reconstructed quantities theta vs phi target.
 * This cut is applied on the original quantities, not on the reconstructed quantites to compare
 * with experimental data
 * and a W cut
 */
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
const int CORNERS = 7;
const int Nbin=100;
void rec_sim_xs()
{
  int dpp;
  int runnum = 1899;
  double p0;
  double dp_off ;
  double dpp_min,dpp_max,hdp_min,hdp_max;
  double w_min(-0.005),w_max(0.01);
  double x_min,x_max,y_min,y_max,ph_min,ph_max,th_min,th_max;
  ifstream dataLimit("/home/ton/Documents/GDH_analysis/gd/MC_1st/Lookup_table_w_collimator/4cm_target_center_cd/input_limit.dat",ios_base::in);
  while(dataLimit>>x_min>>x_max>>y_min>>y_max>>th_min>>th_max>>ph_min>>ph_max){
    cout<<"xfoc limit= "<<x_min<<"\t"<<x_max<<endl;
    cout<<"yfoc limit= "<<y_min<<"\t"<<y_max<<endl;
    cout<<"theta foc limit= "<<th_min<<"\t"<<th_max<<endl;
    cout<<"phi foc limit= "<<ph_min<<"\t"<<ph_max<<endl;
  }
  if(dataLimit.bad())
    {
      cerr<<"Error: file not found!"<<endl;
      exit(8);
    }


  if(runnum==1866||runnum==1805||runnum==1899||runnum==1900){
    dpp=0;
    p0=1148.85;
    dp_off = 0.0007;
    dpp_min = -0.004-dp_off;
    dpp_max = 0.001 -dp_off;
    hdp_min =-0.02; hdp_max =0.005;
  }

  if(runnum==1909){
    dpp=-2;
    p0=1125.88;
    dp_off = 0.00021;
    dpp_min = 0.016 - dp_off;
    hdp_min=0; hdp_max =0.025;
  }


  TFile * fout = new TFile("sim_1899_reconstructed_target.root","RECREATE");
  TH1F * hx      = new TH1F("hx","",Nbin,x_min,x_max);
  TH1F * hphi    = new TH1F("hphi","",Nbin,ph_min,ph_max);
  TH1F * hth     = new TH1F("hth","",Nbin,th_min,th_max);
  TH1F * hdp     = new TH1F("hdp","",Nbin,hdp_min,hdp_max);
  TH1F * hy     = new TH1F("hy","",Nbin,y_min,y_max);
  TH1F * hxs     = new TH1F("hxs","",Nbin,0,25000);
  TH2F * hthph     = new TH2F("hthph","",Nbin,-20,20,Nbin,0,70);


  Float_t deltatg,xs,phfoc,yfoc,phitg,thetatg,xfoc,thfoc,dpor,thor,phor,wor,wmm;
  TFile * f = new TFile("1899_elastic_wo_cut.root");
  TTree *ntup = (TTree*)gROOT->FindObject("h1");
  
  double xsc=0.0;
  int tacc=0;
  Int_t n_entries = (Int_t)ntup->GetEntries();
  ntup->SetBranchAddress("phfoc",&(phfoc));
  ntup->SetBranchAddress("yfoc",&(yfoc));
  ntup->SetBranchAddress("xfoc",&(xfoc));
  ntup->SetBranchAddress("thfoc",&(thfoc));
  ntup->SetBranchAddress("deltat",&(deltatg));
  ntup->SetBranchAddress("thetat",&(thetatg));
  ntup->SetBranchAddress("phit",&(phitg));
  ntup->SetBranchAddress("wmm",&(wmm));
  ntup->SetBranchAddress("wor",&(wor));
  ntup->SetBranchAddress("xs",&(xs));
  ntup->SetBranchAddress("thor",&(thor));
  ntup->SetBranchAddress("phor",&(phor));
  ntup->SetBranchAddress("dpor",&(dpor));
  cout<<"Number of events= "<<n_entries<<endl;
  double corner_phi[CORNERS],corner_theta[CORNERS];
  double edge[CORNERS];
  corner_phi[0] = 8.766;//8.554;//8.616;//8.591;
  corner_phi[1] = 8.766;//8.554;//9.052;//8.591;
  corner_phi[2] = 8.319;//7.229;//7.720;//7.219;
  corner_phi[3] = 3.694;//2.357;//2.049;//1.494;
  corner_phi[4] = 2.004;//1.474;//3.312;//3.251;
  corner_phi[5] = 3.081;//3.255;//2.623;//3.251;
  corner_phi[6] = 2.634;//3.255;//

  corner_theta[0] = 50.074;//49.510;//49.839;//49.218;
  corner_theta[1] = 40.331;//43.863;//40.78;//43.606;
  corner_theta[2] = 34.096;//33.867;//33.025;//33.452;
  corner_theta[3] = 32.529;//31.167;//32.391;//31.180;
  corner_theta[4] = 32.529;//31.167;//41.697;//43.606;
  corner_theta[5] = 41.693;//41.689;//50.756;//49.218;
  corner_theta[6] = 51.675;//49.510;// 
  const int DIM = 4;
  double c_phiTg[DIM]={0.006959,0.001251,0.005041,0.011853};
  double c_thetaTg[DIM]={0.045869,0.045380,0.032330,0.035429};
  double edgeCacc[DIM];
  for(int kk=0;kk<n_entries-1;kk++){
    ntup->GetEntry(kk);	
    double slope, intercept;
    for(int ll=0; ll<CORNERS; ll++)
      {
	if(ll==(CORNERS-1)){
	  slope = (corner_theta[ll]-corner_theta[0])/(corner_phi[ll]-corner_phi[0]);
	} else {
	  slope = (corner_theta[ll+1]-corner_theta[ll])/(corner_phi[ll+1]-corner_phi[ll]);
	}
	intercept = corner_theta[ll] - slope*corner_phi[ll];
	edge[ll] = slope*phitg + intercept;
      }  
    if(phor<=corner_phi[0]&&thor>=edge[1]&&thor>=edge[2]&&thor>=corner_theta[3]&&(thor<=edge[4]||thor>=edge[5])&&thor<=edge[6]
       &&wor<=w_max&&wor>=w_min)      
      { //cut on theta vs phi target
      tacc++;
      }

    if(yfoc<=0.03&&yfoc>=-0.02&&phfoc<=0.05&&phfoc>=-0.05&&thfoc<=0.03&&thfoc>=-0.03
       &&wmm<=w_max&&wmm>=w_min)
      {
	double ph_rec =0.001692 -0.012293*xfoc -0.007802*yfoc-0.40254*thfoc +0.627065*phfoc
	  -5.927*yfoc*phfoc -6.172*yfoc*thfoc
	  -11.326*yfoc*yfoc +6.111*phfoc*phfoc;
	double th_rec =0.037990 -0.004055*xfoc +0.849883*yfoc +0.07456*thfoc -0.163464*phfoc;
	double slopeCacc, interceptCacc;
	for(int i=0; i<DIM; i++)
	  {
	    if(i==(DIM-1))
	      {
		slopeCacc =(c_thetaTg[0]-c_thetaTg[i])/(c_phiTg[0]-c_phiTg[i]);
	      } else 
	      {
		slopeCacc= (c_thetaTg[i+1]-c_thetaTg[i])/(c_phiTg[i+1]-c_phiTg[i]);
	      }
	    interceptCacc = c_thetaTg[i] - slopeCacc*c_phiTg[i];
	    edgeCacc[i] = slopeCacc*ph_rec + interceptCacc;
	  }
	if(phor<=corner_phi[0]&&thor>=edge[1]&&thor>=edge[2]&&thor>=corner_theta[3]&&(thor<=edge[4]||phor>=corner_phi[5])&&thor<=corner_theta[6])
	  {

	  hx->Fill(xfoc,xs);
	  hphi->Fill(phfoc,xs);
	  hth->Fill(thfoc,xs);
	  hy->Fill(yfoc,xs);
	  hdp->Fill(deltatg/100.,xs);
	  hxs->Fill(xs,xs);

	  hthph->Fill(phor,thor);
	}// choose event
      }//4D cut
  }// loop events
  hthph->Draw("colz");
  double cacc = hx->GetEntries();//number of event survive in my 4D table cut
  double cross_section = (hxs->Integral())/tacc;

  cout<<"Number of event inside sieve cut = "<<tacc<<endl;
  cout<<"Number of event inside 4D cut = "<<cacc<<endl;
  cout<<"Cross section= "<<cross_section<<endl;
  cout<<"**************"<<endl;

  fout->Write();
  
}
