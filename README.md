# Hydro-mechanically coupled explicit Vumat for Abaqus using the hypoplastic model with intergranular strain extension and considering effective contact stresses

This set of fortran files contains the hydro-mechanically coupled explicit Vumat for Abaqus using the hypoplastic model with intergranular strain extension. The implementations are according to the following papers:
   -  Hydro-mechanically coupling according to the papers [*Vibratory pile driving in water-saturated sand: Back-analysis of model tests using a hydro-mechanically coupled CEL method*](https://doi.org/10.1016/j.sandf.2020.11.005) and [*Impact of the installation on the long-term cyclic behaviour of piles in sand: a numerical study*](https://doi.org/10.1016/j.soildyn.2020.106223) 
   -  Consideration of effective contact stresses for calculation of friction according to [*Hydro-mechanically coupled CEL analyses with effective contact stresses*](https://doi.org/10.1002/nag.3725) 
   -  Incorporation of effects of cavitation according to [*Monopile installation in clay and subsequent response to millions of lateral load cycles*](https://doi.org/10.1016/j.compgeo.2022.105221)
   -  Implementation of the hypoplastic model with intergranular strain extension according to [*Contributions to the numerical modelling of pile installation processes and high-cyclic loading of soils*](https://www.bgu.ruhr-uni-bochum.de/bgu/mam/images/dissertationen/staubach__2022__heft_73_contributions_to_the_numerical_modelling_of_pile_installation_processes_and_high-cyclic_loading_of_soils_mit_db.pdf)

**Please cite these papers if using the subroutines in your research!**

Please see the PDF file for detailed explanations of the routines. 

All simulations and benchmarks have been performed with Abaqus 2023 and Intel OneAPI 2023. Compatibility with other compilers may not be guaranteed!

The routines for the consideration of effective contact stresses only support simulations with up to 8 CPUs.

The auxiliary routines contained in "tools.f" were written by Prof. A. Niemunis.

Please do not send me emails with general questions about the implementations or why 
you are having problems - this is a research software and I cannot offer individual support without a joint research project. However, pointing out shortcomings is always welcome!

I do not warrant the suitability, accuracy, or operational 
characteristics of the implementations. The implementations are provided on an "as-is" basis. The user of the 
software agrees that the distributor of the implementations has no liability whatever for any damages to the user or to any other party.

**Excess pore water pressure during pile installation**
![excess_pore](https://github.com/user-attachments/assets/50d8afcd-ce7e-4121-ad93-3c3adf735aea)
