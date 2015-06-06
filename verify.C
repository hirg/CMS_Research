#include <math.h>
#include <iostream>
#include <vector>

using namespace std;

#include "TH1.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TFile.h"

#define PI 3.141592653589793
#define max_particles 20000

struct PhiDist
{
	float eta_MIN;
	float eta_MAX;

	TH1F * phi_GEN;
	TH1F * phi_RECO;

	PhiDist(float minEta, float maxEta)
	: eta_MIN(minEta), eta_MAX(maxEta)
	{
		phi_GEN = new TH1F(Form("phi_GEN_eta=[%.1f,%.1f]",minEta,maxEta),Form("phi_GEN_eta=[%.1f,%.1f]",minEta,maxEta),200,-PI,PI);
		phi_RECO = new TH1F(Form("phi_RECO_eta=[%.1f,%.1f]",minEta,maxEta),Form("phi_RECO_eta=[%.1f,%.1f]",minEta,maxEta),200,-PI,PI);
	}
};

struct UniversalData
{
    TTree * steg;
    Int_t nEvents;

    Int_t np;       // # of particles in event
    float phi[max_particles];
    float eta[max_particles];
    //float pt[max_particles];
    bool reco[max_particles]; // whether or not particle has been reconstructed
};

void fillEta(UniversalData * data, TH1F * eta_GEN, TH1F * eta_RECO);
void fillPhi(UniversalData * data, PhiDist * phiDist, Int_t cur, Int_t tot);
bool withinEtaRange(float eta, float min, float max);
//void plot(TH1F * eta_GEN, TH1F * eta_RECO, vector<PhiDist*> v);

void verify()
{
	UniversalData * data = new UniversalData;

	//TFile * fInput = new TFile("/home/palmerjh/CMSSW_5_3_20/STEGEff/STEGData.root");
	TFile * fInput = new TFile("/home/palmerjh/CMSSW_5_3_20/STEGEff/testCluster0.root");
	TFile * fOutput = new TFile("verify.root", "RECREATE", "description");

	TH1F * eta_GEN = new TH1F("eta_GEN","eta_GEN",200,-2.5,2.5);
	TH1F * eta_RECO = new TH1F("eta_RECO","eta_RECO",200,-2.5,2.5);

	data->steg = (TTree*) fInput->Get("tree");
    (data->steg)->SetBranchAddress("reconstructed",&(data->reco));


    // Here the TTree object is given the memory addresses
    // of our variables above as places to put the contents
    // of each branch when the "GetEntry" method is called
    (data->steg)->SetBranchAddress("npg",&(data->np));
    (data->steg)->SetBranchAddress("phig",&(data->phi));
    (data->steg)->SetBranchAddress("etag",&(data->eta));
    //(data->steg)->SetBranchAddress("ptg",&(data->pt));

    // total number of events
    data->nEvents = (Int_t) (data->steg)->GetEntries();

    PhiDist * A = new PhiDist(-2.4,-1.6);
    PhiDist * B = new PhiDist(1.6,2.4);
    PhiDist * C = new PhiDist(-1,1);

    vector<PhiDist*> v; v.push_back(A); v.push_back(B); v.push_back(C);

    fillEta(data, eta_GEN, eta_RECO);

    for (Int_t i = 0; i < v.size(); i++)
    {
    	fillPhi(data, v[i], i, v.size());
    	v[i]->phi_GEN->Write();
    	v[i]->phi_RECO->Write();
    }

    eta_GEN->Write();
    eta_RECO->Write();

    fInput->Close();
    fOutput->Close();

   	cout << "DONE" << endl;

    //plot(eta_GEN, eta_RECO, v);
}

void fillEta(UniversalData * data, TH1F * eta_GEN, TH1F * eta_RECO)
{
	for( Int_t i=0; i < data->nEvents; i++)
    {
        // print a message every 1000 events
        if( i%10000 == 0) cout << "Filling eta histograms. Processing event " << i << " out of " << data->nEvents << "...\n";

        // Load the information from event i
        (data->steg)->GetEntry(i);

        //Loop over all particles in the event to determine event planes
        for( Int_t j=0; j < data->np; j++)
        {
        	eta_GEN->Fill((data->eta)[j]);

        	if((data->reco)[j]) {
        		eta_RECO->Fill((data->eta)[j]);
        	}
        }
    }
}

void fillPhi(UniversalData * data, PhiDist * phiDist, Int_t cur, Int_t tot)
{
	for( Int_t i=0; i < data->nEvents; i++)
    {
        // print a message every 1000 events
        if( i%10000 == 0) cout << Form("Filling phi histograms. (Step %d of %d) Processing event %d out of %d...\n",cur,tot,i,data->nEvents);

        // Load the information from event i
        (data->steg)->GetEntry(i);

        //Loop over all particles in the event to determine event planes
        for( Int_t j=0; j < data->np; j++)
        {
        	if (withinEtaRange((data->eta)[j], phiDist->eta_MIN, phiDist->eta_MAX)) {
        		phiDist->phi_GEN->Fill((data->phi)[j]);

        		if((data->reco)[j]) {
        			phiDist->phi_RECO->Fill((data->phi)[j]);
        		}
        	}	
        }
    }
}

bool withinEtaRange(float eta, float min, float max)
{
	return eta >= min && eta <= max;
}
//void plot(TH1F * eta_GEN, TH1F * eta_RECO, vector<PhiDist*> v);