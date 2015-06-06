#define max_particles 50000

void getEta()
{
	TFile * fInput = new TFile("/scratch/appelte1/STEG_run182798_hiGeneralAndPixelTracks_v0.root");
	TFile * fOutput = new TFile("etaDist.root","RECREATE","description");

	TTree * steg = (TTree*) fInput->Get("/ana/tree");
	Int_t np;
	float eta[max_particles]; // of each particle

	// total number of events
	Int_t nevent = (Int_t) steg->GetEntries();
	
	steg->SetBranchAddress("npg",&np);
	steg->SetBranchAddress("etag",&eta);
	
	TH1F * etaDist = new TH1F("eta","eta",200,-2.4,2.4);
	
	for( Int_t i=0; i<nevent; i++)
	{	  
    	// Load the information from event i
    	steg->GetEntry(i);

        for( Int_t j=0; j < np; j++)
        {
        	etaDist->Fill(eta[j]);
        }

        if (i%10000 == 0) cout << i << " events processed out of " << nevent << endl;
    }

    etaDist->Write();
    fInput->Close();
    fOutput->Close();

    cout << "DONE" << endl;

}