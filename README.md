A DFT + DMFT interface code, 
developed based on the projected local orbital formalism as described in: //
J. Phys.: Condens. Matter 23 085601 (2011)
 
The code includes the CT-QMC solver routines (ctqmc_*.f90 files) from the iQIST package: 

Computer Physics Communications, pages 140-160, Volume 195 (2015)



The rest of the routines are created to carry out the following tasks for a complete DMFT iteration:  

(1) constructing the Green's function from the Bloch Hamiltonian and self energy; 
(2) constructing the hybridyzation energy for the CT-QMC input; 
(3) downfolding and upfolding the self energy between Wannier orbital basis and Bloch basis; 
(4) performing Fermi lever fixing to ensure particle number conservation;  
(5) performing a few choices of double counting corrections; 

The code takes as input:  

(1) the diagonal Bloch Hamiltonian (band energies) within the correlation window; 
(2) the projectors that downfold the Bloch Hamiltonian to Wannier basis; 
(3) a simple input file containing the solver control parameters;

The code outputs the converged self energy and the full Green's function. 

The code is MPI paralleled. 
