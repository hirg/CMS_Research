// NOTE: I am not the original editor of this macro. I only edited it to 
// generate events as needed and plan on cleaning it up

cout << "Hello world!00" << endl;

#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TRandom3.h"
#include "TTree.h"
#include <iostream.h>
#include <fstream.h>
#include <stdlib.h>
cout << "Hello world!0" << endl;
using namespace std;
#define PI2 2.0*3.141592653589793

const Int_t maxnh = 20000;
Int_t b_npg, b_n, n;
Float_t b_phirg;
Float_t b_ptg[maxnh], b_etag[maxnh], b_phig[maxnh];

void proSTEGvn(int seed, int events)
{

  int tot_num=events;  //60000 events
  double MeanMult=200;
  double RMSMult=10;
  
//  int ifile = atoi(getenv("IFILE"));
  cout << "Hello world!1" << endl;

  //simple toy event generator
//  TFile f(Form(";C:\root/macros/vndata_50k_%d.root", ifile), "RECREATE", "ROOT file with histograms & tree");
  TFile f(Form("testCluster%d.root", seed), "RECREATE", "ROOT file with histograms & tree");

  cout << "Hello world!2" << endl;
  TTree *tree = new TTree("tree","Event tree with a few branches");
  cout << "Hello world!3" << endl;
  tree->Branch("npg", &b_npg, "npg/I");   // # of particles;
  tree->Branch("phirg", &b_phirg, "phirg/F");  // RP angle;
  tree->Branch("n", &b_n, "n/I");          // same as npg;
  tree->Branch("ptg", &b_ptg, "ptg[n]/F");  // ;
  tree->Branch("etag", &b_etag, "etag[n]/F");
  tree->Branch("phig", &b_phig, "phig[n]/F");
  
  TF1 *EtaDistr = new TF1("EtaDistr","exp(-(x-2.1)^2/6.3)+exp(-(x+2.1)^2/6.3)",-2.4,2.4);
  //TF1 *PhiDistr = new TF1("PhiDistr","1+2*[0]*cos(x)+2*[1]*cos(2*x)+2*[2]*cos(3*x)+2*[3]*cos(4*x)+2*[4]*cos(5*x)+2*[5]*cos(6*x)",0,PI2);
  TF1 *PhiDistr = new TF1("PhiDistr","1+2*[0]*cos(2*x)+2*[1]*cos(3*x)+2*[2]*cos(4*x)+2*[3]*cos(5*x)+2*[4]*cos(6*x)",0,PI2);
  //TF1 *PtDistr  = new TF1("PtDistr","exp (-(x/.40))+0.0015*exp (-(x/1.5))", 0.2,10);	//V~0.12
  TF1 *PtDistr  = new TF1("PtDistr","exp (-(x/0.90))+0.15*exp (-(x/15))", 0.1,10);	//V~=0.06
  //  TF1 *PtDistr = new TF1("PtDistr","[0]*x*TMath::Power(1+(sqrt(x*x+[1]*[1])-[1]/[2]),-[3])",0.2,10);
	//PtDistr->SetParameters(118.836,-0.335972,0.759243,118.836);	//Real data fit with Tsallis
  //TF1 *V1vsEta = new TF1("V1vsEta","-exp(-(x+1)^2)/20-x/30+exp(-(x-1)^2)/20",-2.4,2.4); 
  //TF1 *V2vsPt   = new TF1("V2vsPt","((x/3)^1.8/(1+(x/3)^1.8))*(.00005+(1/x)^0.8)",0.2,10);
TF1 *V2vsPt = new TF1("V2vsPt","((x/[0])^[1]/(1+(x/[2])^[3]))*(.00005+(1/x)^[4])",0.1,10);
V2vsPt->SetParameters(4.81159,1.80783,3.69272,3.11889,0.931485);	//Real data V~0.05

/*  TF1 *V2vsPt = new TF1("V2vsPt", "[0]", 0.1, 10);
  TF1 *V3vsPt = new TF1("V3vsPt", "[0]", 0.1, 10);
  TF1 *V4vsPt = new TF1("V4vsPt", "[0]", 0.1, 10);
  TF1 *V5vsPt = new TF1("V5vsPt", "[0]", 0.1, 10);
  TF1 *V6vsPt = new TF1("V6vsPt", "[0]", 0.1, 10);
	//V2vsPt->SetParameters(5,1.8,3,1.8,0.8); //V~0.06

	V2vsPt->SetParameter(0, 0.3);
	V3vsPt->SetParameter(0, 0);
	V4vsPt->SetParameter(0, 0);
	V5vsPt->SetParameter(0, 0.0);
	V6vsPt->SetParameter(0, 0.0);*/
TF1 *V3vsPt = new TF1("V3vsPt","((x/3.2)^2.3/(1+(x/3.4)^2.1))*(.00005+(1/x)^1.4)",0.1,10);
TF1 *V4vsPt = new TF1("V4vsPt","((x/4.8)^2.1/(1+(x/3.4)^2.1))*(.00005+(1/x)^1.4)",0.1,10);
TF1 *V5vsPt = new TF1("V5vsPt","((x/6.0)^3.2/(1+(x/11.4)^2.1))*(.00005+(1/x)^1.4)",0.1,10);
TF1 *V6vsPt = new TF1("V6vsPt","((x/5.6)^2.4/(1+(x/4.7)^2.1))*(.00005+(1/x)^1.4)",0.1,10);
 
  rnd = new TRandom3(seed);
  
  double v1, v2, v3, v4, v5, v6, ph, myphi, mypt, myeta, phi0, Psi;
  int nnf,k;
  double neta, nphi, npt;
  int slicept;
  
  for(int i=0; i<tot_num; i++){ 
    
    Psi = rnd->Uniform(0.0,PI2);
    b_phirg = Psi; 
    b_npg = rnd->Gaus(MeanMult,RMSMult); 
    n = 0;
  
    for(int j=0; j<b_npg;j++ ){
     
      mypt = PtDistr->GetRandom();
      myeta =  EtaDistr->GetRandom();
     //v1=V1vsEta->Eval(myeta);
      v2=V2vsPt->Eval(mypt);
      v3=V3vsPt->Eval(mypt);
      v4=V4vsPt->Eval(mypt);
      v5=V5vsPt->Eval(mypt);
      v6=V6vsPt->Eval(mypt);

	 /* v1 = 0.3;
	  v2 = 0;
	  v3 = 0;
	  v4 = 0;
	  v5 = 0;
	  v6 = 0;*/

      b_etag[n] = myeta;
      b_ptg[n]  = mypt;
      //PhiDistr->SetParameters(v1,v2,v3,v4,v5,v6);
      PhiDistr->SetParameters(v2,v3,v4,v5,v6);
      
      myphi = PhiDistr->GetRandom(); // randon selection dn/dphi
      myphi = myphi+Psi; // angle in lab frame -- needed for correct cumulant v2
      if (myphi>PI2) myphi=myphi-PI2; // 0 - 2*Pi
      
      b_phig[n] = myphi; // save angle in lab frame
         
      n++;
      
    } // End of loop over particles
 
    if (i%50 == 0) cout << i << " " << "events processed" << endl;

    b_n = n;
    tree->Fill();
  } // End of loop over events
  
  cout << "writing tree" << endl;
  tree->Write();
  cout << "writing to file" << endl;
  f.Write();
  cout << "closing file" << endl;
  f.Close();
  cout << "THE END" << endl;
}
