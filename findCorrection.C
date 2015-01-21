{
//#include <math.h>
//#include <iostream>
//#include "TH1.h"
#define max_particles 1000
#define MAX_EVENTS 100000
//void findPsi_2() 

TFile f("testCluster.root");

// Get the tree from the file
TTree * steg = f.Get("tree");

// Setup branches to read events from the tree
//int max_particles = 1000;
Int_t np; // number of particles in the event
float phi[max_particles]; // arrays for the kinematics
float eta[max_particles]; // of each particle
float pt[max_particles];  // in the event

float phir; // phi of event plane

// total number of events
Int_t nevent = (Int_t) steg->GetEntries();
cout << nevent << endl;

// Here the TTree object is given the memory addresses
// of our variables above as places to put the contents
// of each branch when the "GetEntry" method is called
steg->SetBranchAddress("npg",&np);
steg->SetBranchAddress("phig",&phi);
steg->SetBranchAddress("etag",&eta);
steg->SetBranchAddress("ptg",&pt);
steg->SetBranchAddress("phirg",&phir);

// arrays with sizes equal to number of pT partitions
Float_t x[49];
Float_t xE[49];

Float_t v2_theoretical[49];
Float_t v2_theoretical2[49];
Float_t v2_theoreticalE[49];

Float_t v2_observed[49];
Float_t v2_observed2[49];
Float_t v2_observedE[49];

Float_t v2_corrected[49];
Float_t v2_corrected2[49];
Float_t v2_correctedE[49];

Int_t nParticlesTheor[49];
Int_t nParticlesObs[49];

for (int i = 0; i < 49; i++)
{
    x[i] = i*0.2 + 0.3; // 0.1, 0.3, 0.5,.....
    xE[i] = 0.1;
}

TH1F * testCorrection = new TH1F("TestCorrection", "TestCorrection", 100, -.1, 1.1);

float correction = 0;
// Loop over the events in the tree
for( Int_t i=0; i<nevent; i++)
{
    float sinPhiA = 0.0;
    float cosPhiA = 0.0;
    float sinPhiB = 0.0;
    float cosPhiB = 0.0;
    // print a message every 100 events
    if( i%100 == 0) cout << "Processing event " << i << "...\n";
  
    // Load the information from event i
    steg->GetEntry(i);

 
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
           v2_theoretical[index] += cos(2*(phi[j] - phir));
           v2_theoretical2[index] += pow(cos(2*(phi[j] - phir)),2.0);
           nParticlesTheor[index]++;

           if (abs(eta[j]) < 1) {
               v2_observed[index] += cos(2*(phi[j] - eventPlaneA));
               v2_observed2[index] += pow(cos(2*(phi[j] - eventPlaneA)),2.0);
               nParticlesObs[index]++;
           }
       }
    }
}

correction /= nevent;
correction = sqrt(correction);

for (int i = 0; i < 49; i++) {
    if (nParticlesTheor[i] != 0) {

       v2_theoretical[i] /= nParticlesTheor[i];
       v2_theoretical2[i] /= nParticlesTheor[i];
       v2_theoreticalE[i] = sqrt((v2_theoretical2[i] - pow(v2_theoretical[i], 2.0)) / nParticlesTheor[i]);
    }

    if (nParticlesObs[i] != 0) {
       v2_observed[i] /= nParticlesObs[i];
       v2_observed2[i] /= nParticlesObs[i];

       v2_corrected[i] = v2_observed[i] / correction;
       v2_corrected2[i] = v2_observed2[i] / pow(correction, 2.0);


       v2_observedE[i] = sqrt((v2_observed2[i] - pow(v2_observed[i], 2.0)) / nParticlesObs[i]);

       v2_correctedE[i] = sqrt((v2_corrected2[i] - pow(v2_corrected[i], 2.0)) / nParticlesObs[i]);

    }
}

TCanvas * c1 = new TCanvas("c1","c1",600,600);
c1->cd();
TH1F * hDum = new TH1F("hDum",";p_{T} [GeV/c];v_{2}",10,0,10);
hDum->SetMinimum(-0.05);
hDum->SetMaximum(0.15);
hDum->GetYaxis()->CenterTitle();
hDum->GetXaxis()->CenterTitle();
hDum->Draw();


TF1 *V2vsPt   = new TF1("V2vsPt","((x/4.81159)^1.80783/(1+(x/3.69272)^3.11889))*(.00005+(1/x)^0.931485)",0.2,10);
V2vsPt->Draw("same");

TGraphErrors * g1 = new TGraphErrors(49, x, v2_theoretical, xE, v2_theoreticalE);
g1->SetMarkerStyle(21);
g1->SetLineColor(kBlue);
g1->SetMarkerColor(kBlue);
g1->Draw("p");

TGraphErrors * g2 = new TGraphErrors(49, x, v2_observed, xE, v2_observedE);
g2->SetMarkerStyle(20);
g2->SetLineColor(kGreen);
g2->SetMarkerColor(kGreen);
g2->Draw("p");

TGraphErrors * g3 = new TGraphErrors(49, x, v2_corrected, xE, v2_correctedE);
g3->SetMarkerStyle(20);
g3->SetLineColor(kRed);
g3->SetMarkerColor(kRed);
g3->Draw("p");

//testCorrection->Draw("ap");
TLegend * leg = new TLegend(.7,.7,.9,.9);
leg->SetFillColor(0);
leg->AddEntry(V2vsPt,"Input Model","l");
leg->AddEntry(g1,"v_{2}\{RP\}","lp");
leg->AddEntry(g2,"v_{2}^{obs}\{EP\}","lp");
leg->AddEntry(g3,"v_{2}\{EP\}","lp");
leg->Draw();
}

