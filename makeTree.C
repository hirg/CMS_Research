#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include "stdlib.h"
#include "TFile.h"
#include "TTree.h"
using namespace std;

const Int_t maxnh = 2000;
Int_t b_npg;
Int_t b_pid[maxnh];
Float_t b_ptg[maxnh], b_etag[maxnh], b_phig[maxnh];

/*Int_t b_npg;
vector<Int_t> b_pid;
vector<Float_t> b_ptg, b_etag, b_phig;*/

void printVectors(vector<Int_t> b_pid, vector<Float_t> b_ptg, vector<Float_t> b_etag, vector<Float_t> b_phig);

void makeTree(Int_t job)
{
	//TFile * in = new TFile(Form("/store/user/palmerjh/Results/job-%d/particle_list.dat"))
	TFile out(Form("makeTree_%d.root", job), "RECREATE", "ROOT file with tree"); // output file 
	TTree *tree = new TTree("tree","Event tree with a few branches");
	tree->Branch("npg", &b_npg, "npg/I");   // # of particles in event
  	tree->Branch("ptg", &b_ptg, "ptg[npg]/F");  // transverse momentum of each particle;
  	tree->Branch("etag", &b_etag, "etag[npg]/F"); // eta for each particle
  	tree->Branch("phig", &b_phig, "phig[npg]/F"); // phi for each particle
  	tree->Branch("pid", &b_pid, "pid[npg]/I"); // pid for each particle; pid_lookup.dat in same directory as particle_list.dat

	ifstream infile(Form("/store/user/palmerjh/Results/centrality_30_35/vis_08/job-%d/particle_list.dat", job));
	//ifstream infile(Form("job-%d/particle_list.dat", job));
	string line;

	getline(infile, line); // ignore first line of file
	b_npg = 0;
	Int_t curEvent = 1;
	Int_t event, trashUrQMD, pid;
	Float_t trashTau, trashX, trashY, trashEta, pT, phi, trashRap, eta; 

	while (getline(infile, line))
	{
		istringstream iss(line);

    	if (!(iss >> event >> trashUrQMD >> pid >> trashTau >> trashX >> trashY >> trashEta >> pT >> phi >> trashRap >> eta)) { exit(1); } // error

    	//cout << Form("%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", event, trashUrQMD, pid,
		//	trashTau, trashX, trashY, trashEta, pT, phi, trashRap, eta);
    	//cout << "# of particles: " << b_npg << endl;

    	if (curEvent != event) {
    		if (curEvent % 10 == 0) cout << Form("Processed %d of 50 events", curEvent) << endl;
    		tree->Fill();
    		//printVectors(b_pid,b_ptg,b_etag,b_phig);
    		//b_ptg.clear(); b_etag.clear(); b_phig.clear(); b_pid.clear();
    		//printVectors(b_pid,b_ptg,b_etag,b_phig);
    		b_npg = 0;
    		curEvent = event;
		}
		b_ptg[b_npg] = pT;
		b_etag[b_npg] = eta;
		b_phig[b_npg] = phi;
		b_pid[b_npg] = pid;

		b_npg++;
		//b_ptg.push_back(pT); b_etag.push_back(eta); b_phig.push_back(phi); b_pid.push_back(pid);
	}

	//cout << "# of particles: " << b_npg << endl;

	//printVectors(b_pid,b_ptg,b_etag,b_phig);

	tree->Fill(); // fills data from last event

	tree->Write();
	out.Write();
	out.Close();

	//delete tree;
	//tree = 0;

}

void printVectors(vector<Int_t> b_pid, vector<Float_t> b_ptg, vector<Float_t> b_etag, vector<Float_t> b_phig)
{
	cout << "pid content: ";
	for(Int_t i = 0; i < b_pid.size(); i++) {
		cout << Form("%d\t", b_pid[i]);
	}
	cout << endl;

	cout << "pT content: ";
	for(Int_t i = 0; i < b_ptg.size(); i++) {
		cout << Form("%.1f\t", b_ptg[i]);
	}
	cout << endl;

	cout << "b_etag content: ";
	for(Int_t i = 0; i < b_etag.size(); i++) {
		cout << Form("%.1f\t", b_etag[i]);
	}
	cout << endl;

	cout << "b_phig content: ";
	for(Int_t i = 0; i < b_phig.size(); i++) {
		cout << Form("%.1f\t", b_phig[i]);
	}
	cout << endl;
}
