# abaqusCEL

This set of fortran files contains the hydro-mechanically coupled explicit Vumat for Abaqus using the hypoplastic model with intergranular strain extension. The implementations are according to the following works
   -  Hydro-mechanically coupling according to the papers "Vibratory pile driving in water-saturated sand: Back-analysis of model tests using a hydro-mechanically coupled CEL method" [See here](https://doi.org/10.1016/j.sandf.2020.11.005) and "Impact of the installation on the long-term cyclic behaviour of piles in sand: a numerical study" [See here](https://doi.org/10.1016/j.soildyn.2020.106223) 
   -  Consideration of effective contact stresses for calculation of friction according to "Hydro-mechanically coupled CEL analyses with effective contact stresses" [See here](https://doi.org/10.1002/nag.3725) 
   -  Incorporation of effects of cavitation according to "Monopile installation in clay and subsequent response to millions of lateral load cycles" [See here](https://doi.org/10.1016/j.compgeo.2022.105221)
   -  Implementation of the hypoplastic model with intergranular strain extension according to "Contributions to the numerical modelling of pile installation processes and high-cyclic loading of soils" [See here](https://www.bgu.ruhr-uni-bochum.de/bgu/mam/images/dissertationen/staubach__2022__heft_73_contributions_to_the_numerical_modelling_of_pile_installation_processes_and_high-cyclic_loading_of_soils_mit_db.pdf)

**Please cite these papers if using the subroutines in your research!**

!!! note "Abaqus version"
  All simulations and benchmarks have been performed with Abaqus 2023 and Intel OneAPI 2023
