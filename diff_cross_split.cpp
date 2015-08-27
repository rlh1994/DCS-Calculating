diff_cross_split(){

	//prelim1
	gStyle->SetOptStat(0000);
	gStyle->SetOptFit(1);
	gStyle->SetOptTitle(0);
	gStyle->SetTitle(0);
	gStyle->SetPalette(1);
	gStyle->SetFillColor(1);
	gStyle->SetPadColor(0);
	gStyle->SetCanvasColor(0);
	gStyle->SetStatColor(0);
	gStyle->SetTitleColor(1);
	gStyle->SetTitleSize(0.06);
	gStyle->SetPadBorderMode(0);
	gStyle->SetPadBottomMargin(0.15);
	gStyle->SetPadLeftMargin(0.12);
	gStyle->SetPadRightMargin(0.12); 
	gStyle->SetFrameBorderMode(0);
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetTextFont(22);

  	//change these 3 numbers for number of bins, must match original make matrix
	const int nbins_en = 10;
	const int nbins_theta = 10;
	const int nbins_phi = 4;
	TCanvas* c101 = new TCanvas("c101","MC Unfolded",10,10,800,600);
	TCanvas* c102 = new TCanvas("c102","MC Original",10,10,800,600);
	//Please split the cos canvases how you wish for the number of cos bins you have chosen.
	c101->Divide(5,2);
	c102->Divide(5,2);

	const int nbins_angles = nbins_theta*nbins_phi;
	const int nbins_total = nbins_en*nbins_theta*nbins_phi;
	const int years_in_sec = 365.25*24*60*60;
	double int_flux;

	//read in the integated flux per bin as calculated in honda_integrate.cpp
	ifstream flux;
	flux.open("integrated_flux.dat");
	flux >> int_flux;
	flux.close();

	//Read in the saved histogram
	TFile* m_histo_file = new TFile(Form("icecube_energy_%d_cos_%d_phi_%d.root",nbins_en,nbins_theta,nbins_phi));
	TH2D* m_histo =  (TH2D*)m_histo_file->Get("matrix_histo");

	//create tranformation matrix, true MC and recon MC vectors
	double m_matrix[nbins_total][nbins_total];
	double u_matrix[nbins_total][nbins_total];
	for (int i = 0; i < nbins_total; ++i){
		for (int j = 0; j < nbins_total; ++j){
			m_matrix[j][i] = m_histo->GetBinContent(i+1,j+1);
		}
	}

	//read in the efficiency as calculated in ice_genie_loop.cpp
	double efficiency[nbins_total];
	ifstream effic;
	string junk;
	effic.open(Form("efficiency_energy_%d_cos_%d_phi_%d.dat",nbins_en,nbins_theta,nbins_phi));
	getline(effic,junk);
	for (int i = 0; i < nbins_total; ++i){
		effic >> efficiency[i];
	}

	double MC_true[nbins_total];
	for (int i = 0; i < nbins_total; ++i){
		MC_true[i] = 0.0;
		for (int j = 0; j < nbins_total; ++j)
		{
			MC_true[i] += m_matrix[i][j];
		}
	}

	double MC_recon[nbins_total];
	double MC_recon_smear[nbins_total];
	for (int i = 0; i < nbins_total; ++i){
		MC_recon[i]=0;
		MC_recon_smear[i]=0;
		for (int j = 0; j < nbins_total; ++j)
		{
			MC_recon[i] += m_matrix[j][i];
		}
	}


	//smear the recon data a bit
	TRandom* rand = new TRandom();
	for (int i = 0; i < nbins_total; ++i){
		double bin_content = MC_recon[i];
		double ran_change = TMath::Gaus(0,TMath::Sqrt(bin_content)); //set to 0 if you don't want smearing
		MC_recon_smear[i] = (ran_change) + bin_content;
		if (MC_recon_smear[i] < 0)MC_recon_smear[i] = 0;
		MC_true[i] += ran_change;
		if (MC_true[i] < 0)MC_true[i] = 0;
	}


	//get norms for each column
	double norm[nbins_total];
	double unf[nbins_total];
	double temp[nbins_total];


	for (int iteration = 0; iteration < 1; ++iteration) //change the upper bound to get that many iterations 
	{
		//normalise the columns
		for (int i = 0; i < nbins_total; ++i){ 
			norm[i] = 0;
			for (int j = 0; j < nbins_total; ++j){
				norm[i] += m_matrix[j][i]; 
				if(norm[i]==0) norm[i] = 1;
			}
		}
		for (int i = 0; i < nbins_total; ++i){
			for (int j = 0; j < nbins_total; ++j)
			{
				u_matrix[i][j] = m_matrix[i][j]/(norm[j]);
			}	
		}

		//put the last reconstructed/unfolded run into a vector
		for (int i = 0; i < nbins_total; ++i)
		{
			if(iteration==0){temp[i] = MC_recon_smear[i];}
			else {temp[i] = unf[i];}

		}
		//do the actual unfolding
		for (int i = 0; i < nbins_total; ++i){
			unf[i] =0.0;
			for (int j = 0; j < nbins_total; ++j)
			{
				unf[i] += u_matrix[i][j] * temp[j];
				if(unf[i]<0.0) unf[i] =0.0;
			}
			//add in all the other considerations (eff, #targets, bin width, time, integrated flux per bin, convert to cm)
			if(efficiency[i]<1E-50)unf[i] = 0;
			else {unf[i] = unf[i]/(efficiency[i] * 8.366E+37 * (200/nbins_en) * (2.0/nbins_theta) * (2*TMath::Pi()/nbins_phi) * 3 * years_in_sec * (int_flux/10000));}
		}

		//reweight the tranformation matrix
		double weight[nbins_total];
		for (int i = 0; i < nbins_total; ++i){
			if(temp[i]!=0.0 )weight[i] = unf[i] / temp[i];
			else{ weight[i]=0.0;} //?
		}
		for (int i = 0; i < nbins_total; ++i){
			for (int j = 0; j < nbins_total; ++j){
				m_matrix[i][j] = m_matrix[i][j] * weight[i];
			}
		}
	}

	TH1D* coshistos[nbins_angles];
	for (int i = 0; i < nbins_angles; ++i){
		coshistos[i] =  new TH1D(Form("cos_%d",i),Form("cos_%d",i),nbins_en,0,200);
			for (int k = 0; k < nbins_en; ++k){
				coshistos[i]->SetBinContent(k+1,unf[(nbins_angles*k)+i]);
		}
	}
	TH1D* truehistos[nbins_angles];
	for (int i = 0; i < nbins_angles; ++i){
		truehistos[i] =  new TH1D(Form("true_%d",i),Form("true_%d",i),nbins_en,0,200);
			for (int k = 0; k < nbins_en; ++k){
				truehistos[i]->SetBinContent(k+1,MC_true[(nbins_angles*k)+i]);
		}
	}

	//Let's look at the histogram shall we?
	TLegend* leg1[nbins_theta];
	for (int i = 0, j=0; i < nbins_angles; i+=nbins_phi, j++)
	{

		c101->cd(j+1);
		coshistos[i]->SetLineColor(1);
		coshistos[i]->SetLineStyle(1);
		coshistos[i]->GetXaxis()->SetTitle(Form("%2.1f < cos(theta) < %2.1f         E (GeV)",-1.0+j*2/double(nbins_theta), -1.0+(j+1)*2/double(nbins_theta)));
		coshistos[i]->GetYaxis()->SetTitle("Diff xsec (cm^2 Gev sr s)^-1");
		coshistos[i]->GetYaxis()->SetTitleOffset(1.4);
		double maxheight = 0;

		for (int k = 0; k < nbins_phi; ++k)
		{

			maxheight = TMath::Max(coshistos[i+k]->GetMaximum(), maxheight);
		}
		coshistos[i]->SetMaximum(maxheight + maxheight/20.0);
		coshistos[i]->Draw("");
		

		for (int k = 1; k < nbins_phi; ++k){
			coshistos[i+k]->SetLineColor((k+1)*2);
			coshistos[i+k]->SetLineStyle(k+1);
			coshistos[i+k]->Draw("SAME");
		}

		leg1[j] = new TLegend(0.3,0.72,0.85,0.87);//left,bottom,right,up(NDC)
		for (int k = 0; k < nbins_phi; ++k){
			leg1[j]->AddEntry(coshistos[i+k],Form("%3.0f < Phi < %3.0f",k*360.0/double(nbins_phi),(k+1)*360.0/double(nbins_phi)));
		}
		leg1[j]->SetFillColor(0);
		leg1[j]->SetTextFont(22);
		leg1[j]->SetTextSize(0.03);
		leg1[j]->Draw("");

	}

	for (int i = 0, j=0; i < nbins_angles; i+=nbins_phi, j++)
	{
		c102->cd(j+1);
		truehistos[i]->SetLineColor(1);
		truehistos[i]->SetLineStyle(1);
		truehistos[i]->GetXaxis()->SetTitle(Form("%2.1f < cos(theta) < %2.1f         E (GeV)",-1.0+j*(2/double(nbins_theta)), -1.0+(j+1)*(2/double(nbins_theta))));
		truehistos[i]->GetYaxis()->SetTitle("Number of Events");
		truehistos[i]->GetYaxis()->SetTitleOffset(1.4);
		double maxheight = 0;

		for (int k = 0; k < nbins_phi; ++k){
			maxheight = TMath::Max(truehistos[i+k]->GetMaximum(), maxheight);
		}

		truehistos[i]->SetMaximum(maxheight + maxheight/20.0);
		truehistos[i]->Draw("");

		for (int k = 1; k < nbins_phi; ++k){
			truehistos[i+k]->SetLineColor((k+1)*2);
			truehistos[i+k]->SetLineStyle(k+1);
			truehistos[i+k]->Draw("SAME");
		}

		leg1[j] = new TLegend(0.3,0.72,0.85,0.87);//left,bottom,right,up(NDC)
		for (int k = 0; k < nbins_phi; ++k){
			leg1[j]->AddEntry(truehistos[i+k],Form("%3.0f < Phi < %3.0f",k*360.0/double(nbins_phi),(k+1)*360.0/double(nbins_phi)));
		}
		leg1[j]->SetFillColor(0);
		leg1[j]->SetTextFont(22);
		leg1[j]->SetTextSize(0.03);
		leg1[j]->Draw("");

	}



}