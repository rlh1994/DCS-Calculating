flux_create_cos(){
	
	vector<double> energyedge, energymiddle, costheta, flux;

	//you need to do the middle otherwise it won't read into GENIE correctly

	//get the energy edges with log spacing
	for (double i = -1; i < 3.0; i+=0.05){
		energyedge.push_back(TMath::Power(10,i));
		flux.push_back(TMath::Power(TMath::Power(10,i),-2));
	}
	//get the middle of those bins for GENIE to work correctly
	for (int k = 0; k < energyedge.size()-1; ++k)
	{
		energymiddle.push_back((energyedge[k+1]-energyedge[k])/2+ energyedge[k]);
	}
	//Fill in the middle cos values
	for (double i = 0.975; i > -1; i-=0.05)
	{
		costheta.push_back(i);
	}

	ofstream fake_fluka;
	fake_fluka.open("e-2_including_cos.dat");
	//output the data in fluka format, just not as nice alignment 
	for (int k = 0; k < costheta.size(); ++k){
		for (int j = 0; j < energymiddle.size(); ++j){
			fake_fluka << energymiddle[j] << " & " << costheta[k]
			<< " & " << flux[j] << endl;
		}
	}


}