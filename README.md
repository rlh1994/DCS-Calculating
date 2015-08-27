# DCS-Calculating
ROOT macros and some pre-made data files used to calculate the DCS from and IceCube simulation file.

The user will need to provide an IceCube simulation ROOT file as well as a Honda Flux file.

make_complete_histo.cpp calcualtes and outputs the migration matrix, this should be run first.

honda_integrate.cpp calculates and outputs the integrated Honda flux, this should be run second.

ice_genie_loop.cpp calculates and outputs the efficiency and requires GENIE runs to calcualte it, examples of these files have been provided, this should be run third.

diff_cross_split calculates and outputs histograms for all cos and phi ranges of the differential cross section using all the previous files, this should be run last.

flux_create_cos.cpp can be used to create a FLUKA looking flux file that GENIE can use to run an E^(-2) flux in Gevgen_atmo with a higher max energy (using some modifications to GENIE code, see the PDF file Appendix A for more details).

e-2_including_cos.dat is an example output of flux_create_cos.cpp up to 1000GeV.

genie_energy_cos_3_195_numu/numubar_50k.root are GENIE simulations run using an E^(-2) flux with energy [3,195] each with 50k events of numu or numubar respectively.

For any more information please see the PDF file for a more in depth discussion.
