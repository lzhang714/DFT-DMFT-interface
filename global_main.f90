!
!
! DFT+DMFT code with following segment solver:
!
!=========+=========+=========+=========+=========+=========+=========+>>>
! A test program for dynamical mean field theory (DMFT) self-consistent  !
! engine plus hybridization expansion version continuous time quantum    !
! Monte Carlo (CTQMC) quantum impurity solver                            !
! author  : li huang                                                     !
! version : v2015.12.15T                                                 !
! status  : WARNING: IN TESTING STAGE, USE IT IN YOUR RISK               !
! comment : this impurity solver is based on segment picture formalism   !
!           any question, please contact with huangli712@yahoo.com.cn    !
!=========+=========+=========+=========+=========+=========+=========+>>>
!
!
! built starting from the solver codes by Long Zhang 
!

  program dft_dmft 

     use constants
     use control
     use plo_quantities
     use context, only : grnf,eimp,paux,statfact,energy_dc,nmat,norbs,uumat,total_energy,wssf

     use mmpi

     implicit none

! local variables
     integer :: iter,iatom
     logical :: convergence

! initialize mpi envirnoment
# if defined (MPI)
     call mp_init()
     call mp_comm_rank(myid)
     call mp_comm_size(nprocs)
# endif  /* MPI */

    ! print header 
    if ( myid == master ) then 
        call ctqmc_print_header()
    endif
    
    ! ---------------------------------------
    ! read Bloch Hamitonian and Projectors 
    ! ---------------------------------------
    call plo_read_projectors()

    ! setup the important parameters 
    call ctqmc_config_global()

    ! print out QMC runtime parameters in summary, only for checking
    select case (solvctrl)
      case('QMC  ')
      if ( myid == master ) then 
         call ctqmc_print_summary()
      endif
    end select

    ! ---------------------------------------
    ! allocate Arrays used in the PLO interface
    ! ---------------------------------------
    call plo_allocate()

    ! allocate memory and initialize ONCE for the time being
    call ctqmc_setup_array()

    ! initialize matsubara and imaginary time meshes ONCE
    call ctqmc_meshes_init()

    ! ---------------------------------------
    ! prepare initial hybridization function, init self-consistent iteration
    ! ONCE per atom. and read atom specific inputs
    ! ---------------------------------------
    do iatom=1,natom
     call ctqmc_config_impurity(iatom)
    enddo

    ! Initialize Hybridization Function
    call mp_barrier()
    call ctqmc_selfer_init()
    call mp_barrier()

!=========================================================================
!>>> DMFT ITERATION BEGIN                                              <<<
!=========================================================================

    DMFT_CTQMC_ITERATION: do iter=1,niter

         
      if ( myid == master ) then 
          write(mystd,'(2X,a,i3,a)') ' DFT+DMFT >>> DMFT iter:', iter
          write(*,*)' '
      endif
         
      ! 1/5   double counting correction
      call plo_dcc(iter)
         
      ! 2/5   call solver for each atom 
      do iatom=1,natom

          if ( myid == master ) then ! only master node can do it
             write(*,*)'================================================================'
             write(mystd,'(2X,a,i3)') 'segment QMC solver is applied to Atom', iatom
             write(*,*)'================================================================'
          endif

          call ctqmc_impurity_solver(iter,iatom)

      enddo 

      ! 3/5   find new fermi level etc.
      call plo_fermi()

      ! 4/5   build new hbdy 
      call ctqmc_dmft_selfer(iter)

      ! 5/5   check convergence 
      convergence = .false.
      call ctqmc_dmft_conver(convergence)
      if ( convergence .eqv. .true. ) then
          EXIT DMFT_CTQMC_ITERATION ! jump out the iteration
      endif

    enddo DMFT_CTQMC_ITERATION ! over iter={1,niter} loop

!=========================================================================
!>>> DMFT ITERATION END                                                <<<
!=========================================================================

# if defined (MPI)
      if ( myid == master ) then ! only master node can do it
        call plo_gbloc_gbtau(Gbloc)
      end if
      call mp_barrier()
      call ctqmc_final_array()
      call plo_deallocate()
# endif  /* MPI */

     if ( myid == master ) then ! only master node can do it
         call ctqmc_print_footer()
     endif

! finalize mpi envirnoment
# if defined (MPI)
     call mp_barrier()
     call mp_finalize()
# endif  /* MPI */

  end program dft_dmft 









