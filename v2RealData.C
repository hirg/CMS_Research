// 0-199 centrality; should limit events to 60-79; 10000 minutes; let job run 24 hours

//#include <math.h>
//#include <iostream>
//#include "TH1.h"
//#define cent_min 10
//#define cent_max 20
#define max_particles 50000
#define MAX_EVENTS 100000

void v2RealData(int cent_min, int cent_max)
{ 

TFile fInput("/scratch/appelte1/STEG_run182798_hiGeneralAndPixelTracks_v0.root");
int c_min = cent_min * 2;
int c_max = cent_max * 2 - 1;

TFile * fOutput = new TFile(Form("cent_%d-%d_realData.root",cent_min,cent_max), "RECREATE", "description");
// Get the tree from the file
TTree * steg = fInput.Get("/ana/tree");

// Setup branches to read events from the tree
//int max_particles = 1000;cent_min
Int_t np; // number of particles in the event
Int_t cent;
float phi[max_particles]; // arrays for the kinematics
float eta[max_particles]; // of each particle
float pt[max_particles];  // in the event

// total number of events
Int_t nevent = (Int_t) steg->GetEntries();
cout << nevent << endl;

// Here the TTree object is given the memory addresses
// of our variables above as places to put the contents
// of each branch when the "GetEntry" method is called
steg->SetBranchAddress("npg",&np);
steg->SetBranchAddress("cent",&cent);
steg->SetBranchAddress("phig",&phi);
steg->SetBranchAddress("etag",&eta);
steg->SetBranchAddress("ptg",&pt);

// arrays with sizes equal to number of pT partitions
Float_t x[49];
Float_t xE[49];

Float_t v2_observed[49];
Float_t v2_observed2[49];
Float_t v2_observedE[49];

Float_t v2_corrected[49];
Float_t v2_corrected2[49];
Float_t v2_correctedE[49];

Int_t nParticlesObs[49];

for (int i = 0; i < 49; i++)
{
    x[i] = (Float_t)i*0.2 + 0.3; // 0.3, 0.5,.....
    xE[i] = 0.1;

    nParticlesObs[i] = 0;
    v2_observed[i] = 0.0;
    v2_observed2[i] = 0.0;
    v2_observedE[i] = 0.0;
}

TH1F * testCorrection = new TH1F("TestCorrection", "TestCorrection", 100, -.1, 1.1);

Int_t nEvents = 0;
//nevent = MAX_EVENTS;
float correction = 0;
// Loop over the events in the tree
for( Int_t i=0; i<nevent; i++)
{
   
    // print a message every 100 events
    if( i%100 == 0) cout << "Processing event " << i << "...\n";
  
    // Load the information from event i
    steg->GetEntry(i);

    if (cent >= c_min && cent <= c_max) {
        
        nEvents++;
    
        float sinPhiA = 0.0;
        float cosPhiA = 0.0;
        float sinPhiB = 0.0;
        float cosPhiB = 0.0;
 
        //Loop over all particles in the event to find event planes
        for( Int_t j=0; j < np; j++)
        {
            if (eta[j] > -2.4 && eta[j] < -1.6) {
                sinPhiA += sin(2*phi[j]); // weight = 1
                cosPhiA += cos(2*phi[j]); // weight = 1
            }
            if (eta[j] > 1.6 && eta[j] < 2.4) {
                sinPhiB += sin(2*phi[j]); // weight = 1
                cosPhiB += cos(2*phi[j]); // weight = 1
            }
        }
        float eventPlaneA = atan2(sinPhiA, cosPhiA) / 2;
        float eventPlaneB = atan2(sinPhiB, cosPhiB) / 2;

        float curCorrection = cos(2*(eventPlaneA - eventPlaneB));
        testCorrection->Fill(curCorrection);

        correction += curCorrection;

        //Loop over all particles in the event to find v2_obs and v2_theor
        for( Int_t j=0; j < np; j++)
        {
            // find appropriate pt range
            bool pTRangeFound = false;
            Float_t midPt = 0.3;
            Int_t index = 0;
            while (!pTRangeFound && index < 49)
            {
                if (pt[j] >= midPt - 0.1 && pt[j] < midPt + 0.1) {
                    pTRangeFound = true;
                } else {
                    midPt += 0.2;
                    index++;
                }
            }
            if (pTRangeFound) {
                if (fabs(eta[j]) < 1.0) {
                    v2_observed[index] += cos(2*(phi[j] - eventPlaneA));
                    v2_observed2[index] += pow(cos(2*(phi[j] - eventPlaneA)),2.0);
                    nParticlesObs[index]++;
                }
            }
        }
    }
}

correction /= nEvents;
correction = sqrt(correction);

for (int i = 0; i < 49; i++) {
    
   // cout << i << ": " << nParticlesObs[i] << endl;
   // cout << i << ": " << v2_observed[i] << endl;
    
//cout << i << ": " << nParticlesObs[i] << endl;
    if (nParticlesObs[i] != 0) {
       v2_observed[i] /= nParticlesObs[i];
       v2_observed2[i] /= nParticlesObs[i];

       v2_corrected[i] = v2_observed[i] / correction;
       v2_corrected2[i] = v2_observed2[i] / pow(correction, 2.0);


       v2_observedE[i] = sqrt((v2_observed2[i] - pow(v2_observed[i], 2.0)) / nParticlesObs[i]);

       v2_correctedE[i] = sqrt((v2_corrected2[i] - pow(v2_corrected[i], 2.0)) / nParticlesObs[i]);

    }
}

// cout << v2_observed[9] << endl << v2_corrected[9] << endl;

/*TCanvas * c1 = new TCanvas("c1","c1",600,600);
c1->cd();
TH1F * hDum = new TH1F("hDum",";p_{T} [GeV/c];v_{2}",10,0,10);
hDum->SetMinimum(-0.05);
hDum->SetMaximum(0.15);
hDum->GetYaxis()->CenterTitle();
hDum->GetXaxis()->CenterTitle();
hDum->Draw();*/


TGraphErrors * g2 = new TGraphErrors(49, x, v2_observed, xE, v2_observedE);
g2->SetName("v2_Obs");
g2->Write();

/*g2->SetMarkerStyle(20);
g2->SetLineColor(kGreen);
g2->SetMarkerColor(kGreen);
g2->Draw("p");*/

TGraphErrors * g3 = new TGraphErrors(49, x, v2_corrected, xE, v2_correctedE);
g3->SetName("v2_Corrected");
g3->Write();

fOutput->Write();
fOutput->Close();

/*g3->SetMarkerStyle(20);
g3->SetLineColor(kRed);
g3->SetMarkerColor(kRed);
g3->Draw("p");

//testCorrection->Draw("ap");
TLegend * leg = new TLegend(.7,.7,.9,.9);
leg->SetFillColor(0);
leg->AddEntry(g2,"v_{2}^{obs}\{EP\}","lp");
leg->AddEntry(g3,"v_{2}\{EP\}","lp");
leg->Draw();*/
}

