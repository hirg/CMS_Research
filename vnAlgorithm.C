// 0-199 centrality; should limit events to 60-79; 10000 minutes; let job run 24 hours

//#include <math.h>
//#include <iostream>
//#include "TH1.h"
//#define cent_min 10
//#define cent_max 20
#define max_particles 50000
#define MAX_EVENTS 40000
#define jMAX 21 // to be used in flattening algorithm

//#include <queue>
//using namespace std;

void vnAlgorithm(int n, bool isRealData, int cent_min, int cent_max)
{ 
// n = order of flow; 2 = elliptic, 3 = trianglular,...
int c_min = cent_min * 2;
int c_max = cent_max * 2 - 1;

TFile * fInput;
TFile * fOutput;
TTree * steg;
Int_t cent = 0;

if (isRealData) {
    fInput = new TFile("/scratch/appelte1/STEG_run182798_hiGeneralAndPixelTracks_v0.root");

    fOutput = new TFile(Form("cent_%d-%d_v%d_STEG_run182798_hiGeneralAndPixelTracks_v0.root",cent_min,cent_max,n), "RECREATE", "description");
    // Get the tree from the file
    steg = (TTree*) fInput->Get("/ana/tree");    
    steg->SetBranchAddress("cent",&cent);
} else {

    fInput = new TFile("dataFromSTEG.root");

    fOutput = new TFile(Form("v%d_dataFromSTEG.root",n), "RECREATE", "description");
    // Get the tree from the file
    steg = (TTree*) fInput->Get("tree");
}

// Setup branches to read events from the tree
//int max_particles = 1000;cent_min
Int_t np; // number of particles in the event
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

TH1F * planeAFlat = new TH1F("planeAFlat","planeAFlat",200,-6.28,6.28);
TH1F * planeBFlat = new TH1F("planeBFlat","planeBFlat",200,-6.28,6.28);
TH1F * planeCFlat = new TH1F("planeCFlat","planeCFlat",200,-6.28,6.28);

TH1F * planeAUnflat = new TH1F("planeAUnflat","planeAUnflat",200,-6.28,6.28);
TH1F * planeBUnflat = new TH1F("planeBUnflat","planeBUnflat",200,-6.28,6.28);
TH1F * planeCUnflat = new TH1F("planeCUnflat","planeCUnflat",200,-6.28,6.28);


Int_t nEvents = 0;
nevent = MAX_EVENTS;
float cosAB = 0;	// cos[n(planeAFlat - planeBFlat)]; to be used for calculating resolution
float cosBC = 0;	// cos[n(planeBFlat - planeCFlat)]; to be used for calculating resolution
float cosAC = 0; 	// cos[n(planeAFlat - planeCFlat)]; to be used for calculating resolution

vector<float> queueA;
vector<float> queueB;
vector<float> queueC;

vector<float> sinPhiAVector;
vector<float> cosPhiAVector;
vector<float> sinPhiBVector;
vector<float> cosPhiBVector;
vector<float> sinPhiCVector;
vector<float> cosPhiCVector;

float AVGsinPhiA = 0.0;
float AVGcosPhiA = 0.0;
float AVGsinPhiB = 0.0;
float AVGcosPhiB = 0.0;
float AVGsinPhiC = 0.0;
float AVGcosPhiC = 0.0;


//float aPlanes[nevent];
//float bPlanes[nevent];

// Loop over the events in the tree to calculate event planes and resolutions
for( Int_t i=0; i<nevent; i++)
{
   
    // print a message every 1000 events
    if( i%1000 == 0) cout << "Step 1: Calculating unflattened, recentered event planes (part 1). Processing event " << i << " out of " << nevent << "...\n";
  
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

        AVGsinPhiA += sinPhiA;
        AVGcosPhiA += cosPhiA;
        AVGsinPhiB += sinPhiB;
        AVGcosPhiB += cosPhiB;
        AVGsinPhiC += sinPhiC;
        AVGcosPhiC += cosPhiC;

        sinPhiAVector.push_back(sinPhiA);
        cosPhiAVector.push_back(cosPhiA);
        sinPhiBVector.push_back(sinPhiB);
        cosPhiBVector.push_back(cosPhiB);
        sinPhiCVector.push_back(sinPhiC);
        cosPhiCVector.push_back(cosPhiC);
    }
}

AVGsinPhiA /= nEvents;
AVGcosPhiA /= nEvents;
AVGsinPhiB /= nEvents;
AVGcosPhiB /= nEvents;
AVGsinPhiC /= nEvents;
AVGcosPhiC /= nEvents;

int counter = 0;

/*TH1F * aDEBUG = new TH1F("aDEBUG","aDEBUG",200,-6.28,6.28);
TH1F * bDEBUG = new TH1F("bDEBUG","bDEBUG",200,-6.28,6.28);
TH1F * cDEBUG = new TH1F("cDEBUG","cDEBUG",200,-6.28,6.28);*/

// loop over events to get unflattened, recentered event planes
for( Int_t i=0; i<nevent; i++)
{
    // print a message every 1000 events
    if( i%1000 == 0) cout << "Step 2: Calculating unflattened, recentered event planes (part 2). Processing event " << i << " out of " << nevent << "...\n";
  
    // Load the information from event i
    steg->GetEntry(i);
    if (cent >= c_min && cent <= c_max) {

        float eventPlaneA = atan2(sinPhiAVector[counter] - AVGsinPhiA, cosPhiAVector[counter] - AVGcosPhiA) / n;
        float eventPlaneB = atan2(sinPhiBVector[counter] - AVGsinPhiB, cosPhiBVector[counter] - AVGcosPhiB) / n;
        float eventPlaneC = atan2(sinPhiCVector[counter] - AVGsinPhiC, cosPhiCVector[counter] - AVGcosPhiC) / n;

        /*aDEBUG->Fill(eventPlaneA);
        bDEBUG->Fill(eventPlaneB);
        cDEBUG->Fill(eventPlaneC);*/

        queueA.push_back(eventPlaneA);
        queueB.push_back(eventPlaneB);
        queueC.push_back(eventPlaneC);

        counter++;
    }
}

/*aDEBUG->Write();
bDEBUG->Write();
cDEBUG->Write();*/

vector<float> AcosFlatteningAVG; // each element = avg. of cos(j*n*eventPlaneA) for j = 1, 2, ... , jMAX
vector<float> AsinFlatteningAVG; // each element = avg. of sin(j*n*eventPlaneA) for j = 1, 2, ... , jMAX
vector<float> BcosFlatteningAVG; // each element = avg. of cos(j*n*eventPlaneB) for j = 1, 2, ... , jMAX
vector<float> BsinFlatteningAVG; // each element = avg. of sin(j*n*eventPlaneB) for j = 1, 2, ... , jMAX
vector<float> CcosFlatteningAVG; // each element = avg. of cos(j*n*eventPlaneC) for j = 1, 2, ... , jMAX
vector<float> CsinFlatteningAVG; // each element = avg. of sin(j*n*eventPlaneC) for j = 1, 2, ... , jMAX

// initialize all of these averages
for (int j = 0; j < jMAX; j++)
{
    AcosFlatteningAVG.push_back(0);
    AsinFlatteningAVG.push_back(0);
    BcosFlatteningAVG.push_back(0);
    BsinFlatteningAVG.push_back(0);
    CcosFlatteningAVG.push_back(0);
    CsinFlatteningAVG.push_back(0);
}

counter = 0;
// loop over events to flatten event planes
for( Int_t i=0; i<nevent; i++)
{
    // print a message every 1000 events
    if( i%1000 == 0) cout << "Step 3: Calculating flattened event planes (part 1). Processing event " << i << " out of " << nevent << "...\n";
  
    // Load the information from event i
    steg->GetEntry(i);
    if (cent >= c_min && cent <= c_max) {

        for (int j = 1; j <= jMAX; j++) {
            AcosFlatteningAVG[j - 1] += cos(j*n*queueA[counter]);
            AsinFlatteningAVG[j - 1] += sin(j*n*queueA[counter]);
            BcosFlatteningAVG[j - 1] += cos(j*n*queueB[counter]);
            BsinFlatteningAVG[j - 1] += sin(j*n*queueB[counter]);
            CcosFlatteningAVG[j - 1] += cos(j*n*queueC[counter]);
            CsinFlatteningAVG[j - 1] += sin(j*n*queueC[counter]);
        }

        counter++;
    }
}


 cout << "Flattening Averages:\n";
// finds averages
for (int j = 0; j < jMAX; j++)
{
    AcosFlatteningAVG[j] /= nEvents;
    //cout << "\t" << j+1 << ": (Acos) " << AcosFlatteningAVG[j] << "\n";
    AsinFlatteningAVG[j] /= nEvents;
    //cout << "\t" << j+1 << ": (Asin) " << AsinFlatteningAVG[j] << "\n";
    BcosFlatteningAVG[j] /= nEvents;
    //cout << "\t" << j+1 << ": (Bcos) " << BcosFlatteningAVG[j] << "\n";
    BsinFlatteningAVG[j] /= nEvents;
    //cout << "\t" << j+1 << ": (Bsin) " << BsinFlatteningAVG[j] << "\n";
    CcosFlatteningAVG[j] /= nEvents;
    //cout << "\t" << j+1 << ": (Ccos) " << CcosFlatteningAVG[j] << "\n";
    CsinFlatteningAVG[j] /= nEvents;
    //cout << "\t" << j+1 << ": (Csin) " << CsinFlatteningAVG[j] << "\n";
}

counter = 0;
float flatteningConstantA;
float flatteningConstantB;
float flatteningConstantC;

// loop over events to flatten event planes
for( Int_t i=0; i<nevent; i++)
{
    // print a message every 1000 events
    if( i%1000 == 0) cout << "Step 4: Calculating flattened event planes (part 2) and resolution. Processing event " << i << " out of " << nevent << "...\n";
  
    // Load the information from event i
    steg->GetEntry(i);
    if (cent >= c_min && cent <= c_max) {
        flatteningConstantA = 1;
        flatteningConstantB = 1;
        flatteningConstantC = 1;
        for (int j = 1; j <= jMAX; j++) {
            flatteningConstantA += 2 / ((float)(j*n)) * (-AsinFlatteningAVG[j - 1] * cos(j*n*queueA[counter]) + 
                    AcosFlatteningAVG[j - 1] * sin(j*n*queueA[counter]));
            flatteningConstantB += 2 / ((float)(j*n)) * (-BsinFlatteningAVG[j - 1] * cos(j*n*queueB[counter]) + 
                    BcosFlatteningAVG[j - 1] * sin(j*n*queueB[counter]));
            flatteningConstantC += 2 / ((float)(j*n)) * (-CsinFlatteningAVG[j - 1] * cos(j*n*queueC[counter]) + 
                    CcosFlatteningAVG[j - 1] * sin(j*n*queueC[counter]));
        }

        cout << "event " << i << " or " << counter << " flattening constant A: " << flatteningConstantA << "\n";
        cout << "event " << i << " or " << counter << " flattening constant B: " << flatteningConstantB << "\n";
        cout << "event " << i << " or " << counter << " flattening constant C: " << flatteningConstantC << "\n";

        planeAUnflat->Fill(queueA[counter]);
        planeBUnflat->Fill(queueB[counter]);
        planeCUnflat->Fill(queueC[counter]);

        queueA[counter] *= flatteningConstantA;
        queueB[counter] *= flatteningConstantB;
        queueC[counter] *= flatteningConstantC;

        planeAFlat->Fill(queueA[counter]);
        planeBFlat->Fill(queueB[counter]);
        planeCFlat->Fill(queueC[counter]);

        cosAB += cos(n*(queueA[counter] - queueB[counter]));
        cosBC += cos(n*(queueB[counter] - queueC[counter]));
        cosAC += cos(n*(queueA[counter] - queueC[counter]));

        counter++;
    }
}
cosAB /= nEvents;
cosBC /= nEvents;
cosAC /= nEvents;

float resWrtA = cosAB * cosAC / cosBC; // resWrtA = resolution w.r.t. planeAFlat
float resWrtB = cosAB * cosBC / cosAC; // resWrtB = resolution w.r.t. planeBFlat

resWrtA = sqrt(resWrtA);
resWrtB = sqrt(resWrtB);

counter = 0;
// Loop over the events to calculate v_n
for( Int_t i=0; i<nevent; i++)
{
   
    // print a message every 1000 events
    if( i%1000 == 0) cout << "Step 5: Calculating vn. Processing event " << i << " out of " << nevent << "...\n";
  
    // Load the information from event i
    steg->GetEntry(i);

    if (cent >= c_min && cent <= c_max) {
        float curVn;
        float plane;
        float resolution;
        //Loop over all particles in the event to find vn_obs and vn_theor
        for( Int_t j=0; j < np; j++)
        {
           if (fabs(eta[j]) < 1.0) {
                
                // find appropriate pt range
                int index = (int) (pt[j]/.2) - 1;               

                if (index >= 0 && index <= 48) {
                    if (eta[j] > 0) {
                        plane = queueA[counter];
                        resolution = resWrtA;
                    } else {
                        plane = queueB[counter];
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

        counter++;
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
obs->SetName(Form("v%d_Observed",n));
obs->Write();


/*g2->SetMarkerStyle(20);
g2->SetLineColor(kGreen);
g2->SetMarkerColor(kGreen);
g2->Draw("p");*/

TGraphErrors * cor = new TGraphErrors(49, x, vn_corrected, xE, vn_correctedE);
cor->SetName(Form("v%d_Corrected",n));
cor->Write();

phiA->Write();
phiB->Write();
phiC->Write();

planeAFlat->Write();
planeBFlat->Write();
planeCFlat->Write();

planeAUnflat->Write();
planeBUnflat->Write();
planeCUnflat->Write();


//fOutput->Write();
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

