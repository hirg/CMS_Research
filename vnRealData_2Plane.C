// 0-199 centrality; should limit events to 60-79; 10000 minutes; let job run 24 hours

//#include <math.h>
//#include <iostream>
//#include "TH1.h"
//#define cent_min 10
//#define cent_max 20
#define max_particles 50000
#define MAX_EVENTS 10000

void vnRealData(int n, int cent_min, int cent_max)
{ 
// n = order of flow; 2 = elliptic, 3 = trianglular,...
TFile fInput("/scratch/appelte1/STEG_run182798_hiGeneralAndPixelTracks_v0.root");
int c_min = cent_min * 2;
int c_max = cent_max * 2 - 1;

TFile * fOutput = new TFile(Form("cent_%d-%d_v%d_STEG_run182798_hiGeneralAndPixelTracks_v0.root",cent_min,cent_max,n), "RECREATE", "description");
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

Float_t vn_observedA[49];
Float_t vn_observedA2[49];
Float_t vn_observedAE[49];

Float_t vn_correctedA[49];
Float_t vn_correctedA2[49];
Float_t vn_correctedAE[49];

Float_t vn_observedB[49];
Float_t vn_observedB2[49];
Float_t vn_observedBE[49];

Float_t vn_correctedB[49];
Float_t vn_correctedB2[49];
Float_t vn_correctedBE[49];

Int_t nParticlesObsA[49];
Int_t nParticlesObsB[49];

for (int i = 0; i < 49; i++)
{
    x[i] = (Float_t)i*0.2 + 0.3; // 0.3, 0.5,.....
    xE[i] = 0.1;

    nParticlesObsA[i] = 0;
    vn_observedA[i] = 0.0;
    vn_observedA2[i] = 0.0;
    vn_observedAE[i] = 0.0;

    nParticlesObsB[i] = 0;
    vn_observedB[i] = 0.0;
    vn_observedB2[i] = 0.0;
    vn_observedBE[i] = 0.0;

}

TH1F * testCorrection = new TH1F("TestCorrection", "TestCorrection", 100, -.1, 1.1);
TH1F * phiA = new TH1F("phiA","phiA",200,-6.28,6.28);
TH1F * phiB = new TH1F("phiB","phiB",200,-6.28,6.28);
TH1F * phiVn = new TH1F("phiVn","phiVn",200,-6.28,6.28);

TH1F * planeA = new TH1F("planeA","planeA",200,-6.28,6.28);
TH1F * planeB = new TH1F("planeB","planeB",200,-6.28,6.28);


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
                sinPhiA += sin(n*phi[j]); // weight = 1
                cosPhiA += cos(n*phi[j]); // weight = 1
                phiA->Fill(phi[j]);
            }
            if (eta[j] > 1.6 && eta[j] < 2.4) {
                sinPhiB += sin(n*phi[j]); // weight = 1
                cosPhiB += cos(n*phi[j]); // weight = 1
                phiB->Fill(phi[j]);
            }
        }
        float eventPlaneA = atan2(sinPhiA, cosPhiA) / n;
        float eventPlaneB = atan2(sinPhiB, cosPhiB) / n;

        planeA->Fill(eventPlaneA);
        planeB->Fill(eventPlaneB);

        float curCorrection = cos(n*(eventPlaneA - eventPlaneB));
        testCorrection->Fill(curCorrection);

        correction += curCorrection;

        //Loop over all particles in the event to find vn_obs and vn_theor
        for( Int_t j=0; j < np; j++)
        {
            if (fabs(eta[j]) < 1.0) {
                
                phiVn->Fill(phi[j]);
                // find appropriate pt range
                int index = (int) (pt[j]/.2) - 1;               

                if (index >= 0 && index <= 48) {
                    vn_observedA[index] += cos(n*(phi[j] - eventPlaneA));
                    vn_observedA2[index] += pow(cos(n*(phi[j] - eventPlaneA)),2.0);
                    nParticlesObsA[index]++;
                   
                    vn_observedB[index] += cos(n*(phi[j] - eventPlaneB));
                    vn_observedB2[index] += pow(cos(n*(phi[j] - eventPlaneB)),2.0);
                    nParticlesObsB[index]++;

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
    if (nParticlesObsA[i] != 0) {
       vn_observedA[i] /= nParticlesObsA[i];
       vn_observedA2[i] /= nParticlesObsA[i];

       vn_correctedA[i] = vn_observedA[i] / correction;
       vn_correctedA2[i] = vn_observedA2[i] / pow(correction, 2.0);


       vn_observedAE[i] = sqrt((vn_observedA2[i] - pow(vn_observedA[i], 2.0)) / nParticlesObsA[i]);

       vn_correctedAE[i] = sqrt((vn_correctedA2[i] - pow(vn_correctedA[i], 2.0)) / nParticlesObsA[i]);

    }
    if (nParticlesObsB[i] != 0) {
       vn_observedB[i] /= nParticlesObsB[i];
       vn_observedB2[i] /= nParticlesObsB[i];

       vn_correctedB[i] = vn_observedB[i] / correction;
       vn_correctedB2[i] = vn_observedB2[i] / pow(correction, 2.0);


       vn_observedBE[i] = sqrt((vn_observedB2[i] - pow(vn_observedB[i], 2.0)) / nParticlesObsB[i]);

       vn_correctedBE[i] = sqrt((vn_correctedB2[i] - pow(vn_correctedB[i], 2.0)) / nParticlesObsB[i]);

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


TGraphErrors * obsA = new TGraphErrors(49, x, vn_observedA, xE, vn_observedAE);
obsA->SetName(Form("v%d_ObsA",n));
obsA->Write();

TGraphErrors * obsB = new TGraphErrors(49, x, vn_observedB, xE, vn_observedBE);
obsB->SetName(Form("v%d_ObsB",n));
obsB->Write();


/*g2->SetMarkerStyle(20);
g2->SetLineColor(kGreen);
g2->SetMarkerColor(kGreen);
g2->Draw("p");*/

TGraphErrors * corA = new TGraphErrors(49, x, vn_correctedA, xE, vn_correctedAE);
corA->SetName(Form("v%d_CorrectedA",n));
corA->Write();

TGraphErrors * corB = new TGraphErrors(49, x, vn_correctedB, xE, vn_correctedBE);
corB->SetName(Form("v%d_CorrectedB",n));
corB->Write();

phiA->Write();
phiB->Write();
phiVn->Write();
planeA->Write();
planeB->Write();

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

