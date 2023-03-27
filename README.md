A DFT + DMFT interface code, <br />
developed based on the projected local orbital formalism as described in: <br />
J. Phys.: Condens. Matter 23 085601 (2011) <br />
   <br />
The code includes the CT-QMC solver routines (ctqmc_*.f90 files) from the iQIST package: <br />
Computer Physics Communications, pages 140-160, Volume 195 (2015) <br />
  <br />
  <br />
The rest of the routines are created to carry out the following tasks for a complete DMFT iteration:  <br />
  <br />
(1) constructing the Green's function from the Bloch Hamiltonian and self energy; <br />
(2) constructing the hybridyzation energy for the CT-QMC input; <br />
(3) downfolding/upfolding the self energy and Green's function between Wannier orbital basis and Bloch basis; <br />
(4) performing Fermi lever fixing to ensure particle number conservation;  <br />
(5) performing a few choices of double counting corrections; <br />
  <br />
The code takes as input:  <br />
  <br />
(1) the diagonal Bloch Hamiltonian (band energies) within the correlation window; <br />
(2) the projectors that downfold the Bloch Hamiltonian to Wannier basis; <br />
(3) a simple input file containing the solver control parameters;<br />
The input files can be based on any DFT calculation, as long as the user ensures the correctness of the projectors. <br />
  <br />
The code outputs the converged self energy and the full Green's function. <br />
  <br />
The code is MPI paralleled. <br />
