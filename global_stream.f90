!
! L.Zhang: modified, input file format changed to account for the solver,
!

! ============================================================================================ 
! ============================================================================================ 

!>>> setup key parameters for continuous time quantum Monte Carlo quantum
! impurity solver and dynamical mean field theory kernel
  subroutine ctqmc_config_global()
     use constants
     use control

     use mmpi

     implicit none

! local variables
! used to check whether the input files (solver.ctqmc.in, solver.hyb.in,solver.eimp.in) exists
     logical :: exists

!=========================================================================
! setup common variables
!=========================================================================
     issun  = 2            ! without symmetry    (1) or with symmetry   mode (2)
     isspn  = 1            ! spin projection, PM (1) or AFM             mode (2)
     isbin  = 2            ! without binning     (1) or with binning    mode (2)
!-------------------------------------------------------------------------
     nband  = 1            ! number of correlated bands
     nspin  = 2            ! number of spin projection
     niter  = 10           ! maximum number of DMFT + CTQMC self-consistent iterations
!-------------------------------------------------------------------------
     ntotal = 1.00_dp      ! Total number of particles
     beta   = 10.00_dp     ! inversion of temperature
     muman  = 0            ! manual chemical potential control
     manfermi  = 0.0       ! manual chemical potential setting
     alpha  = 0.70_dp      ! mixing parameter for self-consistent engine
     solvctrl = 'QMC  '    ! solver switch
     mfreq = 1024          ! number of frequencies
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!=========================================================================
! setup continuous time quantum Monte Carlo quantum impurity solver and some common variables
!=========================================================================
     mkink  = 1024         ! maximum perturbation expansions order
     nfreq  = 128          ! maximum number of matsubara frequency sampling by quantum impurity solver
     ntime  = 1024         ! number of time slice
     nflip  = 20000        ! flip period for spin up and spin down states
     ntherm = 200000       ! maximum number of thermalization steps
     nsweep = 20000000     ! maximum number of quantum Monte Carlo sampling steps
     nwrite = 2000000      ! output period
     nclean = 100000       ! clean update period
     nmonte = 10           ! how often to sampling the gmat and nmat
     ncarlo = 10           ! how often to sampling the gtau and prob
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!=========================================================================
! setup Lanczos ED solver variables
!=========================================================================
     eta = 0.01            ! imaginary offset for real axis
     fit = 3               ! fitting routine, 0 no fit, 1 2 3 some routines
     weight = 1            ! weight for the fit, 1: constant
     nev_def = 50          ! Number of lowest eigenvalues to find in Arnoldi
     ed_tol = 0.001        ! Boltzmann factor cutoff
     lunn = .false.        ! density-density interaction in ED

! read in input file if possible, only master node can do it
     if ( myid == master ) then
         exists = .false.

! inquire file status
         inquire (file = 'solver.global.in', exist = exists)

! read in parameters, default setting should be overrided
         if ( exists .eqv. .true. ) then
             open(mytmp, file='solver.global.in', form='formatted', status='unknown')

             read(mytmp,*)
             read(mytmp,*)
             read(mytmp,*)
!------------------------------------------------------------------------+
             read(mytmp,*) issun                                         !
             read(mytmp,*) isspn                                         !
             read(mytmp,*) isbin                                         !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+
             read(mytmp,*)
!------------------------------------------------------------------------+
             read(mytmp,*) nband                                         !
             read(mytmp,*) nspin                                         !
             read(mytmp,*) niter                                         !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+
             read(mytmp,*)
!------------------------------------------------------------------------+
             read(mytmp,*) ntotal                                        !
             read(mytmp,*) beta                                          !
             read(mytmp,*) muman                                         !
             read(mytmp,*) manfermi                                      !
             read(mytmp,*) alpha                                         !
             read(mytmp,*) solvctrl                                      !
             read(mytmp,*) mfreq                                         !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+
             read(mytmp,*)
!------------------------------------------------------------------------+
             read(mytmp,*) mkink                                         !
             read(mytmp,*) nfreq                                         !
             read(mytmp,*) ntime                                         !
             read(mytmp,*) nflip                                         !
             read(mytmp,*) ntherm                                        !
             read(mytmp,*) nsweep                                        !
             read(mytmp,*) nwrite                                        !
             read(mytmp,*) nclean                                        !
             read(mytmp,*) nmonte                                        !
             read(mytmp,*) ncarlo                                        !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+
             read(mytmp,*)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+
             read(mytmp,*) eta                                           !
             read(mytmp,*) fit                                           !
             read(mytmp,*) weight                                        !
             read(mytmp,*) nev_def                                       !
             read(mytmp,*) ed_tol                                        !
             read(mytmp,*) lunn                                          !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+
             read(mytmp,*)

             close(mytmp)
         endif ! back if ( exists .eqv. .true. ) block
! trim the empty spaces
        solvctrl = trim(adjustl(solvctrl))

! Adjust norbs and ncfgs automatically
    norbs = nspin*nband
    ncfgs = 2**norbs

     endif ! back if ( myid == master ) block

! since config parameters may be updated in master node, it is important
! to broadcast config parameters from root to all children processes
# if defined (MPI)

!------------------------------------------------------------------------+
     call mp_bcast( issun , master )                                     !
     call mp_bcast( isspn , master )                                     !
     call mp_bcast( isbin , master )                                     !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+
     call mp_barrier()

!------------------------------------------------------------------------+
     call mp_bcast( nband , master )                                     !
     call mp_bcast( nspin , master )                                     !
     call mp_bcast( norbs , master )                                     !
     call mp_bcast( ncfgs , master )                                     !
     call mp_bcast( niter , master )                                     !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+
     call mp_barrier()

!------------------------------------------------------------------------+

     call mp_barrier()

!------------------------------------------------------------------------+
     call mp_bcast( ntotal  , master )                                   !
     call mp_bcast( beta  , master )                                     !
     call mp_bcast( muman  , master )                                    !
     call mp_bcast( manfermi  , master )                                 !
     call mp_bcast( alpha , master )                                     !
     call mp_bcast( solvctrl , master )                                  !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+
     call mp_barrier()

!------------------------------------------------------------------------+
     call mp_bcast( mkink , master )                                     !
     call mp_bcast( mfreq , master )                                     !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+
     call mp_barrier()

!------------------------------------------------------------------------+
     call mp_bcast( nfreq , master )                                     !
     call mp_bcast( ntime , master )                                     !
     call mp_bcast( nflip , master )                                     !
     call mp_bcast( ntherm, master )                                     !
     call mp_bcast( nsweep, master )                                     !
     call mp_bcast( nwrite, master )                                     !
     call mp_bcast( nclean, master )                                     !
     call mp_bcast( nmonte, master )                                     !
     call mp_bcast( ncarlo, master )                                     !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+
     call mp_barrier()
!------------------------------------------------------------------------+
     call mp_bcast( eta, master )                                        !
     call mp_bcast( fit, master )                                        !
     call mp_bcast( weight, master )                                     !
     call mp_bcast( nev_def, master )                                    !
     call mp_bcast( ed_tol, master )                                     !
     call mp_bcast( lunn, master )                                       !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+
     call mp_barrier()


# endif  /* MPI */

     return
  end subroutine ctqmc_config_global


! ============================================================================================ 
! ============================================================================================ 

!>>> setup key parameters for continuous time quantum Monte Carlo quantum
! impurity solver that are impurity dependent

  subroutine ctqmc_config_impurity(iatom)

     use constants
     use control
     use mmpi

     implicit none

! local variables
     character(20) :: filename
     integer, intent(in   ) :: iatom

!-------------------------------------------------------------------------
     switch = 1            ! 0: we do the U-2J model, 1: we read from file, 2: we compute from F integrals
     Uc(iatom)  = 5.0      ! U aka F0
     Jz(iatom) = 1.0       ! J aka 1/14(F2+F4) in 5 band
     cfb_switch = 0        ! 0: no CFB rotation, 1: we rotate in every iteration
!-------------------------------------------------------------------------
     dc_type(iatom) = 0           ! Type of double counting to use
     dc_start(iatom) = 0.0        ! Starting double counting for MET approach
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
     nbath(iatom) = 1	  ! bath sites per atom
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! read in input file if possible, only master node can do it
     if ( myid == master ) then

    write(filename,*) iatom
    filename = 'solver.local.in_'//trim(adjustl(filename))

    open(mytmp, file=filename, form='formatted', status='unknown')

             read(mytmp,*)
!------------------------------------------------------------------------+
             read(mytmp,*) switch                                        !
             read(mytmp,*) Uc(iatom)                                     !
             read(mytmp,*) Jz(iatom)                                     !
             read(mytmp,*) cfb_switch                                    !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+
             read(mytmp,*)                                               !
             read(mytmp,*) dc_type(iatom)                                !
             read(mytmp,*) dc_start(iatom)                               !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+
             read(mytmp,*)                                               !
             read(mytmp,*) nbath(iatom)                                  !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+

             read(mytmp,*)
!------------------------------------------------------------------------+
             close(mytmp)
     endif ! back if ( myid == master ) block

! since config parameters may be updated in master node, it is important
! to broadcast config parameters from root to all children processes
# if defined (MPI)

!------------------------------------------------------------------------+
     call mp_bcast( switch    , master )                                 !
     call mp_bcast( Uc(iatom)    , master )                             !
     call mp_bcast( Jz(iatom)    , master )                             !
     call mp_bcast( cfb_switch    , master )                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+
     call mp_barrier()

!------------------------------------------------------------------------+
     call mp_bcast( dc_type  , master )                                  !
     call mp_bcast( dc_start  , master )                                 !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+
     call mp_bcast( nbath  , master )                                    !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+
     call mp_barrier()

# endif  /* MPI */

     return
  end subroutine ctqmc_config_impurity


! ============================================================================================ 
! ============================================================================================ 

!>>> allocate memory for global variables and then initialize them 

  subroutine ctqmc_setup_array() 

     use context

     implicit none

! allocate memory for context module
     call ctqmc_allocate_memory_clur()

     call ctqmc_allocate_memory_umat()
     call ctqmc_allocate_memory_mmat()

     call ctqmc_allocate_memory_gmat()
     call ctqmc_allocate_memory_wmat()
     call ctqmc_allocate_memory_smat()

     return
  end subroutine ctqmc_setup_array

! ============================================================================================ 
! ============================================================================================ 
!>>> initialize the meshes...

  subroutine ctqmc_meshes_init()

     use constants
     use control
     use context

     use mmpi

     implicit none

! local variables
! loop index
     integer  :: i
     integer  :: j
     integer  :: k

! build identity: unity
     unity = czero
     do i=1,norbs
         unity(i,i) = cone
     enddo ! over i={1,norbs} loop

! build imaginary time tau mesh: tmesh
     do i=1,ntime
         tmesh(i) = zero + ( beta - zero ) / real(ntime - 1) * real(i - 1)
     enddo ! over i={1,ntime} loop

! build matsubara frequency mesh: rmesh
     do j=1,mfreq
         rmesh(j) = ( two * real(j - 1) + one ) * ( pi / beta )
     enddo ! over j={1,mfreq} loop

! build matsubara frequency mesh: cmesh
     do k=1,mfreq
         cmesh(k) = czi * ( two * real(k - 1) + one ) * ( pi / beta )
     enddo ! over k={1,mfreq} loop

end subroutine ctqmc_meshes_init


!=========================================================================
!>>>                                                                   <<<
!=========================================================================

!  subroutine lanc_meshes_init(e_control)

!     use constants
!     use control
!     use context
!     use lanc_quantities
!     use plo_quantities

!     use mmpi

!     implicit none

!     integer :: i

! controlchar controlling if we want real axis or matsubara etc.
!     integer, intent(in) :: e_control

!     if (e_control==0)then

! matsubara mesh without mu
!      zenergy(:) = cmplx(zero, rmesh(:))

!     elseif (e_control==1)then
! real energy mesh without mu
!     forall( i = 1:mfreq )        &
!	real_mesh(i) = emin + ( emax - emin )*( i - 1 ) / ( mfreq - 1 )

!      zenergy(:) = cmplx(real_mesh(:),eta)

!     elseif(e_control==2)then
! real energy mesh WITH mu
!     forall( i = 1:mfreq )        &
!	real_mesh(i) = efermi + emin + ( emax - emin )*( i - 1 ) / ( mfreq - 1 )
!
!      zenergy(:) = cmplx(real_mesh(:),eta)

!     end if

!end subroutine lanc_meshes_init

!=========================================================================
!>>>                                                                   <<<
!=========================================================================

!>>> initialize the continuous time quantum Monte Carlo quantum impurity
! solver plus dynamical mean field theory self-consistent engine

  subroutine ctqmc_selfer_init()

     use constants
     use control
     use context
     use plo_quantities

     use mmpi

     implicit none

! local variables
! loop index
     integer :: iat

! used to check whether the self-energy exists
     logical  :: exists_4
!
! Check if a self-energy exists
!
!-------------------------------------------------------------------------
     if ( myid == master ) then ! only master node can do it

! Addition:Inquire about input self energy (useful for the charge selfoconsistent
! implementation)
         inquire (file = 'solver.sgm.in', exist = exists_4)

!==================================================================
!
! reset hybridization to zero
     hybf = czero

! setup initial symm
     symm = 1

! setup initial eimp
     eimp = zero

        if(exists_4 .eqv. .true.) then
         write(*,*) 'Self-energy input found!'
         write(*,*)''
         call plo_read_selfenergy(sigma)
        end if

     endif ! back if ( myid == master ) block

    call mp_barrier()
    call mp_bcast(exists_4, master)

! check for self-energy
    if(exists_4 .eqv. .true.) then

! broadcast selfenergy
         call mp_bcast(sigma, master)
         
         call plo_rewrite_array(sig2,sigma,'plo2solver')

         call plo_fermi()

         call hybfunc_and_levels_sigma(hybf)

         call find_symmetries()
    else

! calculate hybridization function etc.

        call plo_fermi()

        do iat=1,natom
            call hybfunc_and_levels(iat)
        enddo


      if(myid==master)then

        do iat=1,natom
            call find_symmetries(iat)
        enddo

      end if ! master

     end if

    call mp_barrier()

! write out the hybridization function
     if ( myid == master ) then ! only master node can do it
        do iat=1,natom
         call ctqmc_dump_hybf_init(rmesh,hybf,iat)
        enddo
     endif

! since the hybridization function may be updated in master node, it is
! important to broadcast it from root to all children processes

! broadcast hybridization data
     call mp_bcast(hybf, master)

! broadcast level data
     call mp_bcast(eimp, master)

! broadcast symmetry data
     call mp_bcast(symm, master)

! Produce Coulomb interaction matrix
!
     if ( myid == master ) then ! only master node can do it
	 call ctqmc_make_uumat(full_u_mat,uumat,switch,cfb_switch)
     end if

! broadcast umatrix data
     call mp_bcast(full_u_mat, master)
     
     call mp_bcast(uumat, master)

! block until all processes have reached here
     call mp_barrier()

  end subroutine ctqmc_selfer_init

! ============================================================================================ 
! ============================================================================================ 

!>>> initialize the continuous time quantum Monte Carlo quantum impurity solver

  subroutine ctqmc_solver_init(atom)
     use constants
     use control
     use context

     use stack
     use spring

     implicit none

! local variables
! loop index
     integer :: i
     integer :: j

     integer, intent(in) :: atom

! system time since 1970, Jan 1, used to generate the random number seed
     integer :: system_time

! random number seed for twist generator
     integer :: stream_seed

! init random number generator
     call system_clock(system_time)
     stream_seed = abs( system_time - ( myid * 1981 + 2008 ) * 951049 )
     call spring_sfmt_init(stream_seed)

! init empty_s and empty_e stack structure
     do i=1,norbs
         call istack_clean( empty_s(i) )
         call istack_clean( empty_e(i) )
     enddo ! over i={1,norbs} loop

     do i=1,norbs
         do j=mkink,1,-1
             call istack_push( empty_s(i), j )
             call istack_push( empty_e(i), j )
         enddo ! over j={mkink,1} loop
     enddo ! over i={1,norbs} loop

! init statistics variables
     insert_tcount = zero
     insert_accept = zero
     insert_reject = zero

     remove_tcount = zero
     remove_accept = zero
     remove_reject = zero

     lshift_tcount = zero
     lshift_accept = zero
     lshift_reject = zero

     rshift_tcount = zero
     rshift_accept = zero
     rshift_reject = zero

     reswap_tcount = zero
     reswap_accept = zero
     reswap_reject = zero

     reflip_tcount = zero
     reflip_accept = zero
     reflip_reject = zero

! init global variables
     ckink   = 0
     cstat   = 0

! init hist  array
     hist    = 0

! init rank  array
     rank    = 0

! init stts  array
! stts = 0 : null occupation
! stts = 1 : partial occupation, segment scheme
! stts = 2 : partial occupation, anti-segment scheme
! stts = 3 : full occupation
     stts    = 0

! init index array
     index_s = 0
     index_e = 0

! init time  array
     time_s  = zero
     time_e  = zero

! no no no, initialization not

! init probability for atomic states
!     prob    = zero

! init auxiliary physical observables
!     paux    = zero

! init occupation number array
!     nmat    = zero
!     nnmat   = zero

! init M-matrix related array
!     mmat    = zero
!     lspace  = zero
!     rspace  = zero

! init imaginary time impurity green's function array
!     gtau    = zero

! init imaginary time bath weiss's function array
!     wtau    = zero

! init exponent array exp_s and exp_e
     exp_s   = czero
     exp_e   = czero

! init G-matrix related array
!     gmat    = czero
!     lsum    = czero
!     rsum    = czero

! init impurity green's function array
!     grnf    = czero

! init bath weiss's function array
!     wssf    = czero

! init self-energy function array
! note: sig1 should not be reinitialized here, since it is used to keep
! the persistency of self-energy function
!<     sig1    = czero
!     sig2    = czero

! fourier transformation hybridization function from matsubara frequency
! space to imaginary time space

     call ctqmc_fourier_hybf(hybf,htau,atom)

! symmetrize the hybridization function on imaginary time axis if needed
     if ( issun == 2 .or. isspn == 1 ) then
         call ctqmc_symm_gtau(symm, htau,atom)
     endif

! calculate the 2nd-derivates of htau, which is used in spline subroutines
     call ctqmc_make_hsed(tmesh, htau, hsed)

! write out the hybridization function on imaginary time axis
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_htau(tmesh,htau,atom)
     endif

! write out the seed for random number stream, it is useful to reproduce
! the calculation process once fatal error occurs.
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a,i11)') 'seed:', stream_seed
     endif

     return
  end subroutine ctqmc_solver_init

!>>> garbage collection for this program, please refer to ctqmc_setup_array
  subroutine ctqmc_final_array()
     use context

     implicit none

! deallocate memory for context module
     call ctqmc_deallocate_memory_clur()

     call ctqmc_deallocate_memory_umat()
     call ctqmc_deallocate_memory_mmat()

     call ctqmc_deallocate_memory_gmat()
     call ctqmc_deallocate_memory_wmat()
     call ctqmc_deallocate_memory_smat()

! deallocate random other stuff

     return
  end subroutine ctqmc_final_array
