make_complete_histo(){

	
	TChain* genie = new TChain("MasterTree");
	genie->Add("Level7_genie_ic.1460.0000XX_LowE.root");

	//change these 3 numbers for the amount of bins: Warning, there is a n^2 efficiency so be careful 
	const int nbins_en = 10; 
	const int nbins_cos_theta = 10; 
	const int nbins_phi =4; 

	const int nbins_total = nbins_en*nbins_cos_theta*nbins_phi;
	double binedges_en[nbins_en+1];
	double binedges_cos_theta[nbins_cos_theta+1];
	double binedges_phi[nbins_phi+1];


	//The final matrix histogram
	TH2D* matrix_histo = new TH2D("matrix_histo","matrix_histo",nbins_total,0,nbins_total,nbins_total,0,nbins_total);


	//Initialise the bin edges
	binedges_en[0]=0.0;
	for (int i = 1; i < nbins_en+1; ++i){ binedges_en[i]= binedges_en[i-1] + 200/(double)nbins_en; }
	binedges_cos_theta[0]=-1.0;
	for (int i = 1; i < nbins_cos_theta+1; ++i){ binedges_cos_theta[i]= binedges_cos_theta[i-1] + 2/(double)nbins_cos_theta; }
	binedges_phi[0]=0.0;
	for (int i = 1; i < nbins_phi+1; ++i){ binedges_phi[i]= binedges_phi[i-1] + (2*TMath::Pi())/(double)nbins_phi; }	


	//Fill each bin with the # recon(ij) in true(mn)
	for (int true_e = 0; true_e < nbins_en; ++true_e){ //true energy bin loop
		for (int recon_e = 0; recon_e < nbins_en; ++recon_e){ //recon energy bin loop
			cout << endl<< true_e << " " << recon_e << endl;
			for (int true_theta = 0; true_theta < nbins_cos_theta; ++true_theta){
				for (int recon_theta = 0; recon_theta < nbins_cos_theta; ++recon_theta)
				{
					//create the histogram to read values from, setting energy cuts and theta bins
					TH2D* firstrun = new TH2D("firstrun","firstrun",nbins_phi,binedges_phi,nbins_phi,binedges_phi);
					//create the stupidly long cut 
					char* cut_e = Form("(MCMuon.energy >%8.7f && MultiNest8D_Track.energy >%8.7f && MCMuon.energy <=%8.7f && MultiNest8D_Track.energy <=%8.7f", binedges_en[true_e],binedges_en[recon_e], binedges_en[true_e+1], binedges_en[recon_e+1]);
					char* cut_theta = Form(" && TMath::Cos(MCMuon.zenith) >%8.7f && TMath::Cos(MultiNest8D_Track.zenith) >%8.7f && TMath::Cos(MCMuon.zenith) <=%8.7f && TMath::Cos(MultiNest8D_Track.zenith) <=%8.7f)", binedges_cos_theta[true_theta], binedges_cos_theta[recon_theta], binedges_cos_theta[true_theta+1], binedges_cos_theta[recon_theta+1] );
					char* cut_total[9999];	
					strcpy(cut_total, cut_e);


					//project all azimuth values for this cut
					genie->Project("firstrun","MCMuon.azimuth:MultiNest8D_Track.azimuth",strcat(cut_total, cut_theta));

					//move values into matrix_histo
					for (int t_bin = 0; t_bin < nbins_phi; ++t_bin){	
						for (int r_bin = 0; r_bin < nbins_phi; ++r_bin){
							//read for all azimuth angles for specific energy and zenith
							double tempval = firstrun->GetBinContent(r_bin+1,t_bin+1);
							matrix_histo->SetBinContent((recon_e*nbins_cos_theta*nbins_phi)+(true_theta*nbins_phi)+r_bin+1,(true_e*nbins_cos_theta*nbins_phi)+(recon_theta*nbins_phi)+t_bin+1, tempval);
						}
					}

					//delete histogram to free space for next loop
					delete firstrun;
					cout << true_theta << " " << recon_theta << endl;
				}
			}
		}
	}
	//save the histogram for later use
	matrix_histo->SaveAs(Form("icecube_energy_%d_cos_%d_phi_%d.root",nbins_en,nbins_cos_theta,nbins_phi));

}