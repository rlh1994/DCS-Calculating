ice_genie_loop(){


  	//change these 3 numbers for number of bins, must match original make matrix
	const int nbins_en = 10; 
	const int nbins_cos_theta = 10; 
	const int nbins_phi = 4; 

	//change the name of the two files that contain neutrino and anti-neutrino simulations
	//get the genie run
	TChain* genie = new TChain("gst");
	genie->Add("genie_energy_cos_3_195_numu_50k.root");
	//add the second run (the anti-neutrino)
	TChain* genie_anti = new TChain("gst");
	genie_anti->Add("genie_energy_cos_3_195_numubar_50k.root");


	const int nbins_total = nbins_en*nbins_cos_theta*nbins_phi;
	double binedges_en[nbins_en+1];
	double binedges_cos_theta[nbins_cos_theta+1];
	double binedges_phi[nbins_phi+1];


	//get the cut icecube data to compare to
	TFile* m_histo_file = new TFile(Form("icecube_energy_%d_cos_%d_phi_%d.root",nbins_en,nbins_cos_theta,nbins_phi));
	TH2D* m_histo =  (TH2D*)m_histo_file->Get("matrix_histo");
	double m_matrix[nbins_total][nbins_total];
	for (int i = 0; i < nbins_total; ++i){
		for (int j = 0; j < nbins_total; ++j){
			m_matrix[j][i] = m_histo->GetBinContent(i+1,j+1);			
		}
	}

	//recreate the true data 
	double MC_true[nbins_total];
	for (int i = 0; i < nbins_total; ++i){
		MC_true[i] = 0.0;
		for (int j = 0; j < nbins_total; ++j){
			MC_true[i] += m_matrix[i][j];
		}
	}


	TH1D* simulation_histo = new TH1D("simulation_histo","simulation_histo",nbins_total,0,nbins_total);
	//Initialise the bin edges
	binedges_en[0]=0.0;
	for (int i = 1; i < nbins_en+1; ++i){ binedges_en[i]= binedges_en[i-1] + 200/(double)nbins_en; }
	binedges_cos_theta[0]=-1.0;
	for (int i = 1; i < nbins_cos_theta+1; ++i){ binedges_cos_theta[i]= binedges_cos_theta[i-1] + 2/(double)nbins_cos_theta; }
	binedges_phi[0]=-TMath::Pi();
	for (int i = 1; i < nbins_phi+1; ++i){ binedges_phi[i]= binedges_phi[i-1] + (2*TMath::Pi())/(double)nbins_phi; }

	for (int true_e = 0; true_e < nbins_en; ++true_e){ //true energy bin loop
			for (int true_theta = 0; true_theta < nbins_cos_theta; ++true_theta){
				//create the histogram to read values from, setting energy cuts and theta bins

				TH1D* firstrun = new TH1D("firstrun","firstrun",nbins_phi,binedges_phi);
				TH1D* secondrun = new TH1D("secondrun","secondrun",nbins_phi,binedges_phi);
				//create the stupidly long cut 
				char* cut_e = Form("El[0] >%8.7f && El[0] <=%8.7f ", binedges_en[true_e], binedges_en[true_e+1]);
				char* cut_theta = Form(" && pzl[0]/El[0] >%8.7f && pzl[0]/El[0] <=%8.7f ", binedges_cos_theta[true_theta], binedges_cos_theta[true_theta+1]);
				char* cut_total[9999];	
				strcpy(cut_total, cut_e);

				//project all azimuth values for this cut
				genie->Project("firstrun","(TMath::ATan2(pxl,pyl))",strcat(cut_total, cut_theta));
				genie_anti->Project("secondrun","(TMath::ATan2(pxl,pyl))",strcat(cut_total, cut_theta));

				//move values into simulation_histo
				for (int bin = 0; bin < nbins_phi; ++bin){	
					//read for all azimuth angles for specific energy and zenith
					double tempval = firstrun->GetBinContent(bin+1);
					double tempval2 = secondrun->GetBinContent(bin+1);
					simulation_histo->SetBinContent((true_e*nbins_cos_theta*nbins_phi)+(true_theta*nbins_phi)+bin+1, tempval); 
					simulation_histo->AddBinContent((true_e*nbins_cos_theta*nbins_phi)+(true_theta*nbins_phi)+bin+1, tempval2);
				}

				//delete histogram to free space for next loop
				delete firstrun;
				delete secondrun;
				cout << true_e << " " << true_theta << endl;
		}
	}

	//change to scale to 30 million events, will depends on how many genie events you ran
	simulation_histo->Scale(300);	
	
	TH1D* icecube = new TH1D("icecube","icecube",nbins_total,0,nbins_total);


	for (int i = 0; i < nbins_total; ++i){
		icecube->SetBinContent(i+1,MC_true[i]);
	}

	//effieciency file
	ofstream effic;
	effic.open(Form("efficiency_energy_%d_cos_%d_phi_%d.dat",nbins_en,nbins_cos_theta,nbins_phi));
	effic << "#energy bins: " << nbins_en << " 0:200 " << "#cos bins " << nbins_cos_theta << " -1:1 " << endl;
	double efficiency[nbins_total];
	for (int i = 0; i < nbins_total; ++i){
		if(simulation_histo->GetBinContent(i+1) !=0){
			efficiency[i] = icecube->GetBinContent(i+1)/simulation_histo->GetBinContent(i+1);
			effic << efficiency[i] << endl;}
		else {effic << efficiency[i] << endl;}
	}
	effic.close();


  return;
}
