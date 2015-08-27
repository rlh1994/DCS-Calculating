honda_integrate(){
	
	//some stuff we will, fixed from running the matrix before.
	const int nbins_en = 10; //steps of 10Gev
	const int nbins_theta = 1; //steps of 180 degrees, Pi rads
	const int nbins_phi = 1; //steps of 360 degrees, Pi/0.5
	const int nbins_total = nbins_en*nbins_theta*nbins_phi;

	int ibin, section, subsection, line, j;
  	double energy, costheta, flux1, flux2, phi;
  	double energy_edge[38];
  	double bin_flux[37];
  	string junk;
  	section = subsection = line = 1;
  	j=0; 
  	//have to hard code the bin edges because creating and reading a file was too much work
	{	energy_edge[0]= 2.981;
		energy_edge[1]= 3.3436;
		energy_edge[2]= 3.7526;
		energy_edge[3]= 4.2096;
		energy_edge[4]= 4.724;
		energy_edge[5]= 5.2998;
		energy_edge[6]= 5.947;
		energy_edge[7]= 6.6722;
		energy_edge[8]= 7.4868;
		energy_edge[9]= 8.3998;
		energy_edge[10]= 9.4252;
		energy_edge[11]= 10.5748;
		energy_edge[12]= 11.8652;
		energy_edge[13]= 13.3128;
		energy_edge[14]= 14.9372;
		energy_edge[15]= 16.7608;
		energy_edge[16]= 18.8052;
		energy_edge[17]= 21.1008;
		energy_edge[18]= 23.6732;
		energy_edge[19]= 26.5648;
		energy_edge[20]= 29.8032;
		energy_edge[21]= 33.4428;
		energy_edge[22]= 37.5192;
		energy_edge[23]= 42.1028;
		energy_edge[24]= 47.2332;
		energy_edge[25]= 53.0048;
		energy_edge[26]= 59.4632;
		energy_edge[27]= 66.7288;
		energy_edge[28]= 74.8612;
		energy_edge[29]= 84.0048;
		energy_edge[30]= 94.2452;
		energy_edge[31]= 105.755;
		energy_edge[32]= 118.645;
		energy_edge[33]= 133.135;
		energy_edge[34]= 149.365;
		energy_edge[35]= 167.615;
		energy_edge[36]= 188.045;
		energy_edge[37]= 211.015;} //initialising some values
 	
 	double int_flux=0;
 	for (int i = 0; i < 36; ++i){bin_flux[i] = 0;}

 	ifstream flux_stream;
 	//if you plan to use a different flux file, you will need to change this here
	flux_stream.open("kam-nu-20-12-n3650.d");

	//read in the honda data
	while ( flux_stream ) {
      flux1 = flux2 = 0.0;
      if (line == 1 || line == 2){
        getline(flux_stream, junk);
        line++; //ignore these lines
      } else {
        flux_stream >> energy >> flux1 >> flux2 >> junk >> junk;
        if(energy > 3.0 && energy < 196.0){
        	bin_flux[j] += flux1 + flux2;
        	++j;
        }
        line++;
        if( line == 104 ){ //new phi range
          ++subsection;
          line = 1;
          j=0;
          getline(flux_stream, junk);
          if (subsection == 13) //new costheta range
          {
            ++section;
            subsection = 1;
          }
        }
      }
    }


    for (int i = 0; i < 37; ++i){
    	cout << bin_flux[i] << " " << energy_edge[i] << " " << int_flux << endl;
    	if (i==0){
    		int_flux += bin_flux[i] * (energy_edge[i+1]-3.0);
  		}else if (i==36){
    		int_flux += bin_flux[i] * (195.0-energy_edge[i]);
    	} else{
    		int_flux += bin_flux[i] * (energy_edge[i+1]-energy_edge[i]);
    	}
	}
    

    ofstream final_answer;
    final_answer.open("integrated_flux.dat");
    	final_answer << int_flux << endl;


}