// 0-199 centrality; should limit events to 60-79; 10000 minutes; let job run 24 hours

//#include <math.h>
//#include <iostream>
//#include "TH1.h"
//#define cent_min 10
//#define cent_max 20
#define max_particles 50000
#define MAX_EVENTS 10000

#include <queue>
using namespace std;

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

Float_t vn_observed[49];
Float_t vn_observed2[49];
Float_t vn_observedE[49];

Float_t vn_corrected[49];
Float_t vn_corrected2[49];
Float_t vn_correctedE[49];

Int_t nParticlesObs[49];

for (int i = 0; i < 49; i++)
{
    x[i] = (Float_t)i*0.2 + 0.3; // 0.3, 0.5,.....
    xE[i] = 0.1;

    nParticlesObs[i] = 0;
    vn_observed[i] = 0.0;
    vn_observed2[i] = 0.0;
    vn_observedE[i] = 0.0;

    vn_corrected[i] = 0.0;
    vn_corrected2[i] = 0.0;
    vn_correctedE[i] = 0.0;
}

TH1F * testCorrection = new TH1F("TestCorrection", "TestCorrection", 100, -.1, 1.1);
TH1F * phiA = new TH1F("phiA","phiA",200,-6.28,6.28);
TH1F * phiB = new TH1F("phiB","phiB",200,-6.28,6.28);
TH1F * phiC = new TH1F("phiC","phiC",200,-6.28,6.28);

TH1F * planeA = new TH1F("planeA","planeA",200,-6.28,6.28);
TH1F * planeB = new TH1F("planeB","planeB",200,-6.28,6.28);
TH1F * planeC = new TH1F("planeC","planeC",200,-6.28,6.28);


Int_t nEvents = 0;
nevent = MAX_EVENTS;
float cosAB = 0;	// cos[n(planeA - planeB)] 
float cosBC = 0;	// cos[n(planeB - planeC)]
float cosAC = 0; 	// cos[n(planeA - planeC)]

queue<float> * queueA;
queue<float> * queueB;

// Loop over the events in the tree to calculate event planes and resolutions
for( Int_t i=0; i<nevent; i++)
{
   
    // print a message every 100 events
    if( i%100 == 0) cout << "Calculating event planes and resolutions. Processing event " << i << "...\n";
  
    // Load the information from event i
    steg->GetEntry(i);

    if (cent >= c_min && cent <= c_max) {
        
        nEvents++;
    
        float sinPhiA = 0.0;
        float cosPhiA = 0.0;
        float sinPhiB = 0.0;
        float cosPhiB = 0.0;
        float sinPhiC = 0.0;
        float cosPhiC = 0.0;
 
        //Loop over all particles in the event to determine event planes
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
            if (fabs(eta[j]) < 1) {
                sinPhiC += sin(n*phi[j]); // weight = 1
                cosPhiC += cos(n*phi[j]); // weight = 1
                phiC->Fill(phi[j]);
            }
        }
        float eventPlaneA = atan2(sinPhiA, cosPhiA) / n;
        float eventPlaneB = atan2(sinPhiB, cosPhiB) / n;
        float eventPlaneC = atan2(sinPhiC, cosPhiC) / n;

        queueA->push(eventPlaneA);
        queueB->push(eventPlaneB);
        planeC->Fill(eventPlaneC);

        cosAB += cos(n*(eventPlaneA - eventPlaneB));
        cosBC += cos(n*(eventPlaneB - eventPlaneC));
        cosAC += cos(n*(eventPlaneA - eventPlaneC));
    }
}
cosAB /= nEvents;
cosBC /= nEvents;
cosAC /= nEvents;

float resWrtA = cosAB * cosAC / cosBC; // resWrtA = resolution w.r.t. planeA
float resWrtB = cosAB * cosBC / cosAC; // resWrtB = resolution w.r.t. planeB

resWrtA = sqrt(resWrtA);
resWrtB = sqrt(resWrtB);

// Loop over the events in the tree to calculate event planes and resolutions
for( Int_t i=0; i<nevent; i++)
{
   
    // print a message every 100 events
    if( i%100 == 0) cout << "Calculating vn. Processing event " << i << "...\n";
  
    // Load the information from event i
    steg->GetEntry(i);

    if (cent >= c_min && cent <= c_max) {
        float curVn;
        float plane;
        float resolution;

        float eventPlaneA = queueA->front();
        float eventPlaneB = queueB->front();

        queueA->pop();
        queueB->pop();

        planeA->Fill(eventPlaneA);
        planeB->Fill(eventPlaneB);

        //Loop over all particles in the event to find vn_obs and vn_theor
        for( Int_t j=0; j < np; j++)
        {
            if (fabs(eta[j]) < 1.0) {
                
                phiC->Fill(phi[j]);
                // find appropriate pt range
                int index = (int) (pt[j]/.2) - 1;               

                if (index >= 0 && index <= 48) {
                    if (eta[j] > 0) {
                        plane = eventPlaneA;
                        resolution = resWrtA;
                    } else {
                        plane = eventPlaneB;
                        resolution = resWrtB;
                    }
                    
                    curVn = cos(n*(phi[j] - plane));

                    vn_observed[index] += curVn;
                    vn_observed2[index] += pow(curVn,2.0);
                    vn_corrected[index] += curVn / resolution;
                    vn_corrected2[index] += pow(curVn / resolution,2.0);

                    nParticlesObs[index]++;                                      
                }
            }
        }
    }
}


for (int i = 0; i < 49; i++) {
    
   // cout << i << ": " << nParticlesObs[i] << endl;
   // cout << i << ": " << v2_observed[i] << endl;
    
//cout << i << ": " << nParticlesObs[i] << endl;
    if (nParticlesObs[i] != 0) {
       vn_observed[i] /= nParticlesObs[i];
       vn_observed2[i] /= nParticlesObs[i];

       vn_corrected[i] /= nParticlesObs[i];
       vn_corrected2[i] /= nParticlesObs[i];


       vn_observedE[i] = sqrt((vn_observed2[i] - pow(vn_observed[i], 2.0)) / nParticlesObs[i]);

       vn_correctedE[i] = sqrt((vn_corrected2[i] - pow(vn_corrected[i], 2.0)) / nParticlesObs[i]);

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


TGraphErrors * obs = new TGraphErrors(49, x, vn_observed, xE, vn_observedE);
obs->SetName(Form("v%d_ObsA",n));
obs->Write();


/*g2->SetMarkerStyle(20);
g2->SetLineColor(kGreen);
g2->SetMarkerColor(kGreen);
g2->Draw("p");*/

TGraphErrors * cor = new TGraphErrors(49, x, vn_corrected, xE, vn_correctedE);
cor->SetName(Form("v%d_corrected",n));
cor->Write();

phiA->Write();
phiB->Write();
phiC->Write();
planeA->Write();
planeB->Write();
planeC->Write();

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

