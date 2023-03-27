!
! L.Zhang: modified, added ED parameters, also U matrix parameters
!
  module control
     use constants, only : dp, int64

     implicit none

!=========================================================================
!>>> character variables                                               <<<
!=========================================================================

! solver control flag
     character(5), public, save :: solvctrl  = 'QMCED'

!=========================================================================
!>>> integer variables                                                 <<<
!=========================================================================

! control flag: symmetry of bands
! if issun == 1, the bands are not symmetrized
! if issun == 2, the bands are symmetrized according to symmetry matrix
     integer, public, save :: issun  = 1

! control flag: symmetry of spin orientation
! if isspn == 1, enforce spin up =  spin down
! if isspn == 2, let spin up and spin down states evolve independently
     integer, public, save :: isspn  = 1

! control flag: impurity green's function binning mode
! if isbin == 1, without binning mode
! if isbin == 2, with binning mode
     integer, public, save :: isbin  = 1

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! number of correlated bands
     integer, public, save :: nband  = 1

! number of spin projection
     integer, public, save :: nspin  = 2

! number of correlated orbitals (= nband * nspin)
     integer, public, save :: norbs  = 2

! number of atomic states (= 2**norbs)
     integer, public, save :: ncfgs  = 4

! maximum number of continuous time quantum Monte Carlo quantum impurity
! solver plus dynamical mean field theory self-consistent iterations
     integer, public, save :: niter  = 20

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! maximum perturbation expansion order
     integer, public, save :: mkink  = 1024

! maximum number of matsubara frequency point
     integer, public, save :: mfreq  = 8193

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! number of matsubara frequency sampling by continuous time quantum Monte
! Carlo quantum impurity solver
! note: the rest (mfreq - nfreq + 1 points) values are evaluated by using
! high frequency expansions or Hubbard-I approximation
     integer, public, save :: nfreq  = 128

! number of imaginary time slice sampling by continuous time quantum Monte
! Carlo quantum impurity solver
     integer, public, save :: ntime  = 1024

! flip period for spin up and spin down states
! note: care must be taken to prevent the system from being trapped in a
! state which breaks a symmetry of local hamiltonian when it should not
! be. to avoid unphysical trapping, we introduce "flip" moves, which
! exchange the operators corresponding, for example, to up and down spins
! in a given orbital.
! note: 0 or negative value means infinite long period
     integer(int64), public, save :: nflip  = 20000

! maximum number of thermalization steps
     integer(int64), public, save :: ntherm = 200000

! maximum number of quantum Monte Carlo sampling steps
     integer(int64), public, save :: nsweep = 20000000

! output period for quantum impurity solver
     integer(int64), public, save :: nwrite = 2000000

! clean update period for quantum impurity solver
     integer(int64), public, save :: nclean = 100000

! how often to sampling the gmat and nmat (and paux)
     integer(int64), public, save :: nmonte = 10

! how often to sampling the gtau and prob
     integer(int64), public, save :: ncarlo = 10

!=========================================================================
!>>> real variables                                                    <<<
!=========================================================================

! average Coulomb interaction
     real(dp), public, save, allocatable :: Uc(:)

! Hund's exchange interaction in z axis (Jz = Js = Jp = J)
    real(dp), public, save, allocatable :: Jz(:)

! Slater integrals
     real(dp), public, save, allocatable :: Fslater(:,:)

! Control parameter for reading U-matrix
     integer, public, save :: switch = 0

! We want cfb rotations in evergy iteration or not
     integer, public, save :: cfb_switch = 0

! Starting Double counting in MET approach
     real(dp), public, save, allocatable :: dc_start(:)

! The total number of particles
     real(dp), public, save :: ntotal
! Type of DCC to use
     integer, public, save, allocatable :: dc_type(:)

! inversion of temperature
     real(dp), public, save :: beta  = 8.00_dp

! controlflag for manual chemical potential
     integer, public, save :: muman  = 0

! manual chemical potential
     real(dp), public, save :: manfermi  = 0.0_dp

! mixing parameter for dynamical mean field theory self-consistent engine
     real(dp), public, save :: alpha = 0.70_dp

!=========================================================================
!>>> MPI related common variables                                      <<<
!=========================================================================

! number of processors: default value 1
     integer, public, save :: nprocs = 1

! the id of current process: default value 0
     integer, public, save :: myid   = 0

! denote as the controller process: default value 0
     integer, public, save :: master = 0

! the id of current process in cartesian topology (cid == myid)
     integer, public, save :: cid    = 0

! the x coordinates of current process in cartesian topology
     integer, public, save :: cx     = 0

! the y coordinates of current process in cartesian topology
     integer, public, save :: cy     = 0

!=========================================================================
!>>> ED Lanczos variables                                              <<<
!=========================================================================

! ED bath sites(per atom)
    integer, public, save, allocatable :: nbath(:)

! control integer for real or imaginary time

    integer, public, save :: ctrl = 0

! Real energy mesh variables
    real(dp) :: emin = -10.0            !  Minimal energy
    real(dp) :: emax =  10.0            !  Maximal energy
    real(dp) :: eta  = 0.05             !  Offset to complex plane

! bath fitting stuff
  integer, save :: fit = 3              ! Define choice of mathematical routine to fit
                                        ! Type of uncertainty function to fit Delta
  integer, save :: weight = 1           ! Constant

  real(dp), save :: ufp(3)              ! Fit's parameters to define uncertainty function
                                        ! (amplitude,position,spread)

! lanczos solver controls
  integer, save :: nev_def = 50       ! Number of lowest eigenvalues to find in Arnoldi
  integer, save :: ed_opt = 2         ! optimization = 0, 1, 2
  integer, save :: verbos = 1         ! verbosity = 0, 10, 20, 30
  integer, save :: ed_ncv0 = 2        ! ncv = 2*newnev + ncv0

  real(dp), save :: ed_tol = 0.001    ! Boltzmann factor cut
  real(dp), save :: ed_evp = 0.25     ! Max # of evalues in one sector
  real(dp), save :: ed_decomp = 0.03  ! Probability cut for eigenvect decomposition

  logical, save :: lunn = .false.      !   Density-density U-matrix in ED

! output not implemented...
  logical, save :: lsisj = .false.    !   Calculate < s_i,s_j > correlators in ED

  end module control
