void functionWriter()
{
	TFile * f = new TFile("/home/palmerjh/CMSSW_5_3_20/STEGEff/STEGFunctions.root", "RECREATE", "contains STEG functions");

	TF1 *V2vsPt = new TF1("V2vsPt","((x/[0])^[1]/(1+(x/[2])^[3]))*(.00005+(1/x)^[4])",0.2,10);
  	V2vsPt->SetParameters(4.81159,1.80783,3.69272,3.11889,0.931485);  //Real data V~0.05

  	TF1 *V3vsPt = new TF1("V3vsPt","((x/3.2)^2.3/(1+(x/3.4)^2.1))*(.00005+(1/x)^1.4)",0.2,10);
	TF1 *V4vsPt = new TF1("V4vsPt","((x/4.8)^2.1/(1+(x/3.4)^2.1))*(.00005+(1/x)^1.4)",0.2,10);
  	TF1 *V5vsPt = new TF1("V5vsPt","((x/6.0)^3.2/(1+(x/11.4)^2.1))*(.00005+(1/x)^1.4)",0.2,10);
  	TF1 *V6vsPt = new TF1("V6vsPt","((x/5.6)^2.4/(1+(x/4.7)^2.1))*(.00005+(1/x)^1.4)",0.2,10);

  	V2vsPt->Write();
  	V3vsPt->Write();
	V4vsPt->Write();
  	V5vsPt->Write();
  	V6vsPt->Write();

  	f->Close();
}