!
! L.Zhang: modified, to account for multiple atoms
!
!=========================================================================
!>>> module ctqmc_core                                                 <<<
!=========================================================================
!>>> containing core (internal) variables used by continuous time quantum
! Monte Carlo quantum impurity solver

  module ctqmc_core

     use constants, only : dp, zero

     implicit none

! current perturbation expansion order
     integer, public, save  :: ckink = 0

! current status for current flavor channel, used to sync with stts
! cstat = 0 : null occupation
! cstat = 1 : partial occupation, segment scheme
! cstat = 2 : partial occupation, anti-segment scheme
! cstat = 3 : full occupation
     integer, public, save  :: cstat = 0

!-------------------------------------------------------------------------
!::: core variables: real, insert action counter                       :::
!-------------------------------------------------------------------------

! insert kink (operators pair) statistics: total insert count
     real(dp), public, save :: insert_tcount = zero

! insert kink (operators pair) statistics: total accepted insert count
     real(dp), public, save :: insert_accept = zero

! insert kink (operators pair) statistics: total rejected insert count
     real(dp), public, save :: insert_reject = zero

!-------------------------------------------------------------------------
!::: core variables: real, remove action counter                       :::
!-------------------------------------------------------------------------

! remove kink (operators pair) statistics: total remove count
     real(dp), public, save :: remove_tcount = zero

! remove kink (operators pair) statistics: total accepted remove count
     real(dp), public, save :: remove_accept = zero

! remove kink (operators pair) statistics: total rejected remove count
     real(dp), public, save :: remove_reject = zero

!-------------------------------------------------------------------------
!::: core variables: real, lshift action counter                       :::
!-------------------------------------------------------------------------

! lshift kink (operators pair) statistics: total lshift count
     real(dp), public, save :: lshift_tcount = zero

! lshift kink (operators pair) statistics: total accepted lshift count
     real(dp), public, save :: lshift_accept = zero

! lshift kink (operators pair) statistics: total rejected lshift count
     real(dp), public, save :: lshift_reject = zero

!-------------------------------------------------------------------------
!::: core variables: real, rshift action counter                       :::
!-------------------------------------------------------------------------

! rshift kink (operators pair) statistics: total rshift count
     real(dp), public, save :: rshift_tcount = zero

! rshift kink (operators pair) statistics: total accepted rshift count
     real(dp), public, save :: rshift_accept = zero

! rshift kink (operators pair) statistics: total rejected rshift count
     real(dp), public, save :: rshift_reject = zero

!-------------------------------------------------------------------------
!::: core variables: real, reswap action counter                       :::
!-------------------------------------------------------------------------

! reswap kink (operators pair) statistics: total reswap count
     real(dp), public, save :: reswap_tcount = zero

! reswap kink (operators pair) statistics: total accepted reswap count
     real(dp), public, save :: reswap_accept = zero

! reswap kink (operators pair) statistics: total rejected reswap count
     real(dp), public, save :: reswap_reject = zero

!-------------------------------------------------------------------------
!::: core variables: real, reflip action counter                       :::
!-------------------------------------------------------------------------

! reflip kink (operators pair) statistics: total reflip count
     real(dp), public, save :: reflip_tcount = zero

! reflip kink (operators pair) statistics: total accepted reflip count
     real(dp), public, save :: reflip_accept = zero

! reflip kink (operators pair) statistics: total rejected reflip count
     real(dp), public, save :: reflip_reject = zero

  end module ctqmc_core

!=========================================================================
!>>> module ctqmc_clur                                                 <<<
!=========================================================================
!>>> containing perturbation expansion series related arrays (colour part)
! used by continuous time quantum Monte Carlo quantum impurity solver

  module ctqmc_clur

     use constants, only : dp
     use stack

     implicit none

! memory address index for the imaginary time \tau_s
     integer, public, save, allocatable :: index_s(:,:)

! memory address index for the imaginary time \tau_e
     integer, public, save, allocatable :: index_e(:,:)

! imaginary time \tau_s of create  operators
     real(dp), public, save, allocatable :: time_s(:,:)

! imaginary time \tau_e of destroy operators
     real(dp), public, save, allocatable :: time_e(:,:)

! exp(i\omega t), s means create  operators
     complex(dp), public, save, allocatable :: exp_s(:,:,:)

! exp(i\omega t), e means destroy operators
     complex(dp), public, save, allocatable :: exp_e(:,:,:)

! container for the empty (unoccupied) memory address index
     type (istack), public, save, allocatable :: empty_s(:)

! container for the empty (unoccupied) memory address index
     type (istack), public, save, allocatable :: empty_e(:)

  end module ctqmc_clur

!=========================================================================
!>>> module ctqmc_umat                                                 <<<
!=========================================================================
!>>> containing util-matrix related arrays used by continuous time quantum
! Monte Carlo quantum impurity solver

  module ctqmc_umat

     use constants, only : dp, int64

     implicit none
! atoms

! histogram for perturbation expansion series
     integer(int64),  public, save, allocatable :: hist(:,:)

! current perturbation expansion order for different flavor channel
     integer,  public, save, allocatable :: rank(:,:)

! current occupation status for different flavor channel
     integer,  public, save, allocatable :: stts(:,:)

!-------------------------------------------------------------------------
!::: input data variables                                              :::
!-------------------------------------------------------------------------

! atoms. they will find you.

! symmetry properties for correlated orbitals
     integer,  public, save, allocatable :: symm(:,:)

! impurity level for correlated orbitals
     real(dp), public, save, allocatable :: eimp(:,:)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!-------------------------------------------------------------------------
!::: physical observables                                              :::
!-------------------------------------------------------------------------
! probability of eigenstates of local hamiltonian matrix
     real(dp), public, save, allocatable :: prob(:,:)

! auxiliary physical observables
! paux(1) : total energy, Etot
! paux(2) : potential engrgy, Epot
! paux(3) : kinetic energy, Ekin
! paux(4) : magnetic moment, < Sz >

!MK statistical factor
     real(dp), public, save, allocatable :: statfact(:)
     real(dp), public, save, allocatable :: energy_dc(:)
     real(dp), public, save :: total_energy
!/MK

     real(dp), public, save, allocatable :: paux(:,:)

! impurity occupation number, < n_i >
     real(dp), public, save, allocatable :: nmat(:,:)

! impurity double occupation number matrix, < n_i n_j >
     real(dp), public, save, allocatable :: nnmat(:,:,:)

! reduced Coulomb interaction matrix, two-index version and atom index
     real(dp), public, save, allocatable :: uumat(:,:,:)

! reduced Coulomb interaction matrix, two-index version and atom index
     real(dp), public, save, allocatable :: full_u_mat(:,:,:,:,:)

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!-------------------------------------------------------------------------
!::: mesh data variables                                               :::
!-------------------------------------------------------------------------
! imaginary time mesh
     real(dp), public, save, allocatable :: tmesh(:)

! real matsubara frequency mesh
     real(dp), public, save, allocatable :: rmesh(:)

! complex matsubara frequency mesh
     complex(dp), public, save, allocatable :: cmesh(:)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! identity matrix
     complex(dp), public, save, allocatable :: unity(:,:)

  end module ctqmc_umat

!=========================================================================
!>>> module ctqmc_mmat                                                 <<<
!=========================================================================
!>>> containing M-matrix and G-matrix related arrays used by continuous
! time quantum Monte Carlo quantum impurity solver

  module ctqmc_mmat

     use constants, only : dp

     implicit none

!atoms...

! helper matrix for evaluating M & G matrices
     real(dp), public, save, allocatable    :: lspace(:,:,:)

! helper matrix for evaluating M & G matrices
     real(dp), public, save, allocatable    :: rspace(:,:,:)

! M matrix, $ \mathscr{M} $
     real(dp), public, save, allocatable    :: mmat(:,:,:,:)

! helper matrix for evaluating G matrix
     complex(dp), public, save, allocatable :: lsum(:,:,:)

! helper matrix for evaluating G matrix
     complex(dp), public, save, allocatable :: rsum(:,:,:)

! G matrix, $ \mathscr{G} $
     complex(dp), public, save, allocatable :: gmat(:,:,:,:)

  end module ctqmc_mmat

!MK
!=========================================================================
!>>> module plo_quantities                                             <<<
!=========================================================================
!>>> contains quantities used in the projection formalism
!
  module plo_quantities

     use constants, only : dp,zero,czero !double precision

     implicit none

! Number of atoms in the unit cell, dimension of correlated subspace (passed to nband),
! number of k-points, number of bands in lda hamiltonian
    integer, public, save :: natom, nkp, lda_nband

! Fermi level read from dft, only used for projector and Hamiltonian setup
    real(dp), public, save :: efermi_dft = zero

! Fermi level
    real(dp), public, save :: efermi = zero

! Average interaction parameters. Don't ask.
    real(dp), public, save, allocatable :: Uavg(:)
    real(dp), public, save, allocatable :: Javg(:)

! Coefficients for the rotation of the basis
    real(dp), public, save, allocatable :: cfb_coeff(:,:,:)
    real(dp), public, save, allocatable :: cfb_det(:,:,:)
! Hopping matrix
    real(dp), public, save, allocatable :: hopmat(:,:)

! the double counting per atom
    real(dp), public, save, allocatable :: dcc(:,:)

! Total Number of particles in the system per spin,
! the total number is in ntotal
! the atoms. they are near. beware.
    real(dp), public, save, allocatable :: ntot(:,:)

! Hamiltonian in Bloch basis from VASP. In Bloch basis
     real(dp), public, save, allocatable :: hbloch(:,:,:)

! PLOs from PROCAR or somewhere else
     complex(dp), public, save, allocatable :: locpro(:,:,:,:,:)

! Number of particles in the correlated subspace per atom
     real(dp), public, save, allocatable :: ncorr(:,:,:)
     
! Density for charge selfconsistency
    complex(dp), public, save, allocatable :: dmft_dens(:,:)
! Occupancy matrix per atom
     real(dp), public, save, allocatable :: occ(:,:,:,:)
!
! Green functions and self energy per atom
     complex(dp), public, save, allocatable :: sigma(:,:,:,:,:)
! Not per atom since it is in the bloch basis
     complex(dp), public, save, allocatable :: Gbloc(:,:,:)

     complex(dp), public, save, allocatable :: G0loc(:,:,:,:,:)
     complex(dp), public, save, allocatable :: Gloc(:,:,:,:,:)

  end module plo_quantities

!
!=========================================================================
!>>> module lanc_quantities                                             <<<
!=========================================================================
!>>> contains quantities used in the Lanczos ED

!  module lanc_quantities
!
!     use constants, only : dp,zero,czero !double precision
!
!     implicit none
!
! energy grid
!
!     complex(dp), public, save, allocatable :: zenergy(:)
!     real(dp), public, save, allocatable :: real_mesh(:)
!
!
! matrix hybridization function for ED
!     complex(dp), public, save, allocatable :: lanc_Delta(:,:,:,:,:)
!
!
! Green function
!     complex(dp), public, save, allocatable :: lanc_Green(:,:,:,:,:)
!
!
! Green0 function
!     complex(dp), public, save, allocatable :: lanc_Green0(:,:,:,:,:)
!
! hybridization function, segment format for dumping
!     complex(dp), public, save, allocatable :: lanc_hybf(:,:,:,:)
!
! Green function, segment format for calculating and dumping
!     complex(dp), public, save, allocatable :: lanc_grnf(:,:,:,:)
!
! Green function, segment format for calculating and dumping
!     complex(dp), public, save, allocatable :: lanc_grnf0(:,:,:,:)
!
!  end module lanc_quantities

!=========================================================================
!>>> module ctqmc_gmat                                                 <<<
!=========================================================================
!>>> containing green's function matrix related arrays used by continuous
! time quantum Monte Carlo quantum impurity solver

  module ctqmc_gmat

     use constants, only : dp

     implicit none
!atoms

! impurity green's function, in imaginary time axis, matrix form
     real(dp), public, save, allocatable    :: gtau(:,:,:,:)

! impurity green's function, in matsubara frequency axis, matrix form
     complex(dp), public, save, allocatable :: grnf(:,:,:,:)

  end module ctqmc_gmat

!=========================================================================
!>>> module ctqmc_wmat                                                 <<<
!=========================================================================
!>>> containing weiss's function and hybridization function matrix related
! arrays used by continuous time quantum Monte Carlo quantum impurity solver

  module ctqmc_wmat

     use constants, only : dp
     use plo_quantities, only : natom

! Atoms, atoms everywhere

     implicit none

! bath weiss's function, in imaginary time axis, matrix form
     real(dp), public, save, allocatable    :: wtau(:,:,:,:)

! bath weiss's function, in matsubara frequency axis, matrix form
     complex(dp), public, save, allocatable :: wssf(:,:,:,:)

! hybridization function, in imaginary time axis, matrix form
     real(dp), public, save, allocatable    :: htau(:,:,:,:)

! hybridization function, in matsubara frequency axis, matrix form
     complex(dp), public, save, allocatable :: hybf(:,:,:,:)

! second order derivates for hybridization function, used to interpolate htau
     real(dp), public, save, allocatable    :: hsed(:,:,:,:)

  end module ctqmc_wmat

!=========================================================================
!>>> module ctqmc_smat                                                 <<<
!=========================================================================
!>>> containing self-energy function matrix related arrays used by
! continuous time quantum Monte Carlo quantum impurity solver

  module ctqmc_smat

     use constants, only : dp

     implicit none

! self-energy function, in matsubara frequency axis, matrix form atoms atoms atoms
     complex(dp), public, save, allocatable :: sig1(:,:,:,:)

! self-energy function, in matsubara frequency axis, matrix form
     complex(dp), public, save, allocatable :: sig2(:,:,:,:)

  end module ctqmc_smat 


!=========================================================================
!>>>                                                                   <<<
!=========================================================================

  subroutine plo_allocate()

  use plo_quantities
  use control
!
! Allocate
!
  allocate( G0loc(nband,nband,mfreq,nspin,natom), Gloc(nband,nband,mfreq,nspin,natom), &
	    Gbloc(lda_nband,mfreq,nspin), sigma(nband,nband,mfreq,nspin,natom),&
	    ncorr(2,nspin,natom),occ(nband,nband,nspin,natom), ntot(nspin,natom), dcc(nspin,natom),&
	    dc_start(natom),dc_type(natom),dmft_dens(lda_nband,lda_nband),Uavg(natom),Javg(natom),&
	    cfb_coeff(natom,nband,nband),hopmat(nband*natom,nband*natom),&
	    Uc(natom),Jz(natom),Fslater(0:6,natom),cfb_det(natom,nband,nband))
!
! And set to zero
!

    G0loc(:,:,:,:,:) = czero
    Gloc(:,:,:,:,:) = czero
    Gbloc(:,:,:) = czero
    sigma(:,:,:,:,:) = czero
    dcc(:,:) = zero

! That has to be here for reasons of the input files...
   allocate(nbath(natom))

  end subroutine plo_allocate

!=========================================================================
!>>>                                                                   <<<
!=========================================================================
  
  subroutine plo_deallocate()

  use plo_quantities
  use control
!
! Deallocate
!
  deallocate( G0loc, Gloc, Gbloc, sigma, ncorr, occ, ntot, dcc, dc_start, dc_type ,hbloch, locpro, &
  dmft_dens,Uavg,Javg,cfb_coeff,hopmat,Uc,Jz,Fslater,nbath)

  end subroutine plo_deallocate

!=========================================================================
!>>>                                                                   <<<
!=========================================================================

!/MK

  !subroutine lanc_allocate()
!
!   use control
!   use plo_quantities
!   use lanc_quantities
!  
!   allocate(zenergy(mfreq))
!   allocate(real_mesh(mfreq))
! 
!   allocate(lanc_Delta(nband,nband,mfreq,nspin,natom))
!   allocate(lanc_Green(nband,nband,mfreq,nspin,natom))
!   allocate(lanc_Green0(nband,nband,mfreq,nspin,natom))
! 
!   allocate(lanc_hybf(mfreq,norbs,norbs,natom),  stat=istat)
!   allocate(lanc_grnf(mfreq,norbs,norbs,natom),  stat=istat)
!   allocate(lanc_grnf0(mfreq,norbs,norbs,natom), stat=istat)
! 
!   lanc_Delta = czero
!   lanc_Green = czero
!   lanc_Green0 = czero
!   lanc_hybf = czero
!   lanc_grnf = czero
!   lanc_grnf0 = czero
! 
!  end subroutine lanc_allocate

!=========================================================================
!>>>                                                                   <<<
!=========================================================================

!  subroutine lanc_deallocate()
!
!   use control
!   use lanc_quantities
!
!   deallocate(lanc_Delta,zenergy,lanc_Green,lanc_Green0,lanc_hybf,lanc_grnf,lanc_grnf0,real_mesh)
!
!  end subroutine lanc_deallocate

!/MK

!=========================================================================
!>>> module context                                                    <<<
!=========================================================================
!>>> containing memory management subroutines and define global variables

  module context

     use constants
     use control

     use ctqmc_core
     use ctqmc_clur

     use ctqmc_umat
     use ctqmc_mmat

     use ctqmc_gmat
     use ctqmc_wmat
     use ctqmc_smat

     implicit none

! status flag
     integer, private :: istat

! declaration of module procedures: allocate memory
     public :: ctqmc_allocate_memory_clur
     public :: ctqmc_allocate_memory_umat
     public :: ctqmc_allocate_memory_mmat
     public :: ctqmc_allocate_memory_gmat
     public :: ctqmc_allocate_memory_wmat
     public :: ctqmc_allocate_memory_smat

! declaration of module procedures: deallocate memory
     public :: ctqmc_deallocate_memory_clur
     public :: ctqmc_deallocate_memory_umat
     public :: ctqmc_deallocate_memory_mmat
     public :: ctqmc_deallocate_memory_gmat
     public :: ctqmc_deallocate_memory_wmat
     public :: ctqmc_deallocate_memory_smat

     contains

!=========================================================================
!>>> allocate memory subroutines                                       <<<
!=========================================================================

!>>> allocate memory for clur-related variables

     subroutine ctqmc_allocate_memory_clur()

         implicit none

! loop index
         integer :: i

! allocate memory
         allocate(index_s(mkink,norbs),     stat=istat)
         allocate(index_e(mkink,norbs),     stat=istat)

         allocate(time_s(mkink,norbs),      stat=istat)
         allocate(time_e(mkink,norbs),      stat=istat)

         allocate(exp_s(nfreq,mkink,norbs), stat=istat)
         allocate(exp_e(nfreq,mkink,norbs), stat=istat)

         allocate(empty_s(norbs),           stat=istat)
         allocate(empty_e(norbs),           stat=istat)

! check the status
         if ( istat /= 0 ) then
             call ctqmc_print_error('ctqmc_allocate_memory_clur','can not allocate enough memory')
         endif

! initialize them
         index_s = 0
         index_e = 0

         time_s  = zero
         time_e  = zero

         exp_s   = czero
         exp_e   = czero

         do i=1,norbs
             empty_s(i) = istack_create(mkink)
             empty_e(i) = istack_create(mkink) 
         enddo ! over i={1,norbs} loop

         return
     end subroutine ctqmc_allocate_memory_clur


!=========================================================================
!>>> module ctqmc_allocate_memory_umat()                               <<<
!=========================================================================

!>>> allocate memory for umat-related variables
     subroutine ctqmc_allocate_memory_umat()
     
         use plo_quantities, only : natom
     
         implicit none

! allocate memory
         allocate(hist(mkink,natom),        stat=istat)
         allocate(rank(norbs,natom),        stat=istat)
         allocate(stts(norbs,natom),        stat=istat)

         allocate(symm(norbs,natom),        stat=istat)

         allocate(eimp(norbs,natom),        stat=istat)

         allocate(prob(ncfgs,natom),        stat=istat)
         allocate(paux(  4,natom  ),        stat=istat)
         allocate(nmat(norbs,natom),        stat=istat)

         allocate(statfact(natom),        stat=istat)
         allocate(energy_dc(natom),        stat=istat)

         allocate(nnmat(norbs,norbs,natom), stat=istat)
         allocate(uumat(norbs,norbs,natom), stat=istat)

         allocate(full_u_mat(nband,nband,nband,nband,natom), stat=istat)

         allocate(tmesh(ntime),       stat=istat)
         allocate(rmesh(mfreq),       stat=istat)

         allocate(cmesh(mfreq),       stat=istat)

         allocate(unity(norbs,norbs), stat=istat)

! check the status
         if ( istat /= 0 ) then
             call ctqmc_print_error('ctqmc_allocate_memory_umat','can not allocate enough memory')
         endif

! initialize them
         hist  = 0
         rank  = 0
         stts  = 0

         symm  = 1

         eimp  = zero

         prob  = zero
         paux  = zero
         nmat  = zero

         statfact = zero
         energy_dc = zero

         nnmat = zero
         uumat = zero

         full_u_mat = zero

         tmesh = zero
         rmesh = zero

         cmesh = czero

         unity = czero

         return
     end subroutine ctqmc_allocate_memory_umat

!=========================================================================
!>>>                                                                   <<<
!=========================================================================

!>>> allocate memory for mmat-related variables
     subroutine ctqmc_allocate_memory_mmat()

        use plo_quantities, only : natom

         implicit none

! allocate memory
         allocate(lspace(mkink,norbs,natom),     stat=istat)
         allocate(rspace(mkink,norbs,natom),     stat=istat)

         allocate(mmat(mkink,mkink,norbs,natom), stat=istat)

         allocate(lsum(nfreq,norbs,natom),       stat=istat)
         allocate(rsum(nfreq,norbs,natom),       stat=istat)

         allocate(gmat(nfreq,norbs,norbs,natom), stat=istat)

! check the status
         if ( istat /= 0 ) then
             call ctqmc_print_error('ctqmc_allocate_memory_mmat','can not allocate enough memory')
         endif

! initialize them
         lspace = zero
         rspace = zero

         mmat   = zero

         lsum   = czero
         rsum   = czero

         gmat   = czero

         return
     end subroutine ctqmc_allocate_memory_mmat

!=========================================================================
!>>>                                                                   <<<
!=========================================================================

!>>> allocate memory for gmat-related variables

     subroutine ctqmc_allocate_memory_gmat()

        use plo_quantities, only : natom

         implicit none

! allocate memory
         allocate(gtau(ntime,norbs,norbs,natom), stat=istat)

         allocate(grnf(mfreq,norbs,norbs,natom), stat=istat)

! check the status
         if ( istat /= 0 ) then
             call ctqmc_print_error('ctqmc_allocate_memory_gmat','can not allocate enough memory')
         endif

! initialize them
         gtau = zero

         grnf = czero

         return
     end subroutine ctqmc_allocate_memory_gmat

!=========================================================================
!>>>                                                                   <<<
!=========================================================================

!>>> allocate memory for wmat-related variables

     subroutine ctqmc_allocate_memory_wmat()

         use plo_quantities, only : natom
         implicit none

! allocate memory
         allocate(wtau(ntime,norbs,norbs,natom), stat=istat)
         allocate(htau(ntime,norbs,norbs,natom), stat=istat)
         allocate(hsed(ntime,norbs,norbs,natom), stat=istat)

         allocate(wssf(mfreq,norbs,norbs,natom), stat=istat)
         allocate(hybf(mfreq,norbs,norbs,natom), stat=istat)

! check the status
         if ( istat /= 0 ) then
             call ctqmc_print_error('ctqmc_allocate_memory_wmat','can not allocate enough memory')
         endif

! initialize them
         wtau = zero
         htau = zero
         hsed = zero

         wssf = czero
         hybf = czero

         return
     end subroutine ctqmc_allocate_memory_wmat

!=========================================================================
!>>>                                                                   <<<
!=========================================================================

!>>> allocate memory for smat-related variables

     subroutine ctqmc_allocate_memory_smat()
     
         use plo_quantities, only : natom
         implicit none

! allocate memory
         allocate(sig1(mfreq,norbs,norbs,natom), stat=istat)
         allocate(sig2(mfreq,norbs,norbs,natom), stat=istat)

! check the status
         if ( istat /= 0 ) then
             call ctqmc_print_error('ctqmc_allocate_memory_smat','can not allocate enough memory')
         endif

! initialize them
         sig1 = czero
         sig2 = czero

         return
     end subroutine ctqmc_allocate_memory_smat

!=========================================================================
!>>> deallocate memory subroutines                                     <<<
!=========================================================================

!>>> deallocate memory for clur-related variables
     subroutine ctqmc_deallocate_memory_clur()
         implicit none

! loop index
         integer :: i

         do i=1,norbs
             call istack_destroy(empty_s(i))
             call istack_destroy(empty_e(i))
         enddo ! over i={1,norbs} loop

         if ( allocated(index_s) ) deallocate(index_s)
         if ( allocated(index_e) ) deallocate(index_e)

         if ( allocated(time_s)  ) deallocate(time_s )
         if ( allocated(time_e)  ) deallocate(time_e )

         if ( allocated(exp_s)   ) deallocate(exp_s  )
         if ( allocated(exp_e)   ) deallocate(exp_e  )

         if ( allocated(empty_s) ) deallocate(empty_s)
         if ( allocated(empty_e) ) deallocate(empty_e)

         return
     end subroutine ctqmc_deallocate_memory_clur

!=========================================================================
!>>>                                                                   <<<
!=========================================================================

!>>> deallocate memory for umat-related variables
     subroutine ctqmc_deallocate_memory_umat()
         implicit none

         if ( allocated(hist)  )   deallocate(hist )
         if ( allocated(rank)  )   deallocate(rank )
         if ( allocated(stts)  )   deallocate(stts )

         if ( allocated(symm)  )   deallocate(symm )

         if ( allocated(eimp)  )   deallocate(eimp )

         if ( allocated(prob)  )   deallocate(prob )
         if ( allocated(paux)  )   deallocate(paux )
         if ( allocated(nmat)  )   deallocate(nmat )

         if ( allocated(statfact)  )   deallocate(statfact )
         if ( allocated(energy_dc)  )   deallocate(energy_dc )

         if ( allocated(nnmat) )   deallocate(nnmat)
         if ( allocated(uumat) )   deallocate(uumat)
         
         if ( allocated(full_u_mat) )   deallocate(full_u_mat)

         if ( allocated(tmesh) )   deallocate(tmesh)
         if ( allocated(rmesh) )   deallocate(rmesh)

         if ( allocated(cmesh) )   deallocate(cmesh)

         if ( allocated(unity) )   deallocate(unity)

         return
     end subroutine ctqmc_deallocate_memory_umat

!=========================================================================
!>>>                                                                   <<<
!=========================================================================

!>>> deallocate memory for mmat-related variables
     subroutine ctqmc_deallocate_memory_mmat()
         implicit none

         if ( allocated(lspace) )  deallocate(lspace)
         if ( allocated(rspace) )  deallocate(rspace)

         if ( allocated(mmat)   )  deallocate(mmat  )

         if ( allocated(lsum)   )  deallocate(lsum  )
         if ( allocated(rsum)   )  deallocate(rsum  )

         if ( allocated(gmat)   )  deallocate(gmat  )

         return
     end subroutine ctqmc_deallocate_memory_mmat

!=========================================================================
!>>>                                                                   <<<
!=========================================================================

!>>> deallocate memory for gmat-related variables
     subroutine ctqmc_deallocate_memory_gmat()
         implicit none

         if ( allocated(gtau) )    deallocate(gtau)

         if ( allocated(grnf) )    deallocate(grnf)

         return
     end subroutine ctqmc_deallocate_memory_gmat

!=========================================================================
!>>>                                                                   <<<
!=========================================================================

!>>> deallocate memory for wmat-related variables
     subroutine ctqmc_deallocate_memory_wmat()
         implicit none

         if ( allocated(wtau) )    deallocate(wtau)
         if ( allocated(htau) )    deallocate(htau)
         if ( allocated(hsed) )    deallocate(hsed)

         if ( allocated(wssf) )    deallocate(wssf)
         if ( allocated(hybf) )    deallocate(hybf)

         return
     end subroutine ctqmc_deallocate_memory_wmat

!=========================================================================
!>>>                                                                   <<<
!=========================================================================

!>>> deallocate memory for smat-related variables
     subroutine ctqmc_deallocate_memory_smat()
         implicit none

         if ( allocated(sig1) )    deallocate(sig1)
         if ( allocated(sig2) )    deallocate(sig2)

         return
     end subroutine ctqmc_deallocate_memory_smat

  end module context
