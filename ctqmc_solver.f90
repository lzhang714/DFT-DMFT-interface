!
! MK: modified, at least to account for multiple atoms, potentially also in other ways.
! Original header:
!
!-------------------------------------------------------------------------
! project : azalea
! program : ctqmc_impurity_solver
!           ctqmc_diagram_warmming
!           ctqmc_diagram_sampling
!           ctqmc_diagram_templing
!           ctqmc_diagram_checking
!           ctqmc_impurity_tester
! source  : ctqmc_solver.f90
! type    : subroutines
! author  : li huang (email:huangli712@yahoo.com.cn)
! history : 09/16/2009 by li huang
!           09/20/2009 by li huang
!           09/24/2009 by li huang
!           09/26/2009 by li huang
!           10/20/2009 by li huang
!           10/29/2009 by li huang
!           11/01/2009 by li huang
!           11/17/2009 by li huang
!           11/22/2009 by li huang
!           12/02/2009 by li huang
!           12/04/2009 by li huang
!           12/06/2009 by li huang
!           12/17/2009 by li huang
!           12/22/2009 by li huang
!           12/26/2009 by li huang
!           12/30/2009 by li huang
!           01/13/2010 by li huang
!           02/27/2010 by li huang
!           06/09/2010 by li huang
!           06/21/2010 by li huang
! purpose : the main subroutine for the hybridization expansion version
!           continuous time quantum Monte Carlo (CTQMC) quantum impurity
!           solver
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

!>>> core engine for hybridization expansion version continuous time
! quantum Monte Carlo quantum impurity solver
  subroutine ctqmc_impurity_solver(iter,atom)
     use constants
     use control
     use context
     use plo_quantities, only : natom

     implicit none

! external arguments
! current self-consistent iteration number and guess what the other is
     integer, intent(in) :: iter
     integer, intent(in) :: atom

! local variables
! loop index
     integer(int64)  :: i
     integer(int64)  :: j
     integer  :: m
     integer  :: n

! status flag
     integer  :: istat

! current QMC sweeping steps
     integer(int64)  :: cstep

! control flag, whether the solver is checked periodically
! cflag = 0  , do not check the quantum impurity solver
! cflag = 1  , check the quantum impurity solver periodically
! cflag = 99 , the quantum impurity solver is out of control
! cflag = 100, the quantum impurity solver has reached convergence
     integer  :: cflag

! starting time
     real(dp) :: time_begin

! ending time
     real(dp) :: time_end

! time consuming by current iteration 
     real(dp) :: time_iter

! time consuming by total iteration
     real(dp) :: time_niter

! histogram for perturbation expansion series, for mpi case
     integer(int64), allocatable  :: hist_mpi(:,:)

! impurity occupation number matrix, for mpi case
     real(dp), allocatable :: nmat_mpi(:,:)

! probability of atomic states, for mpi case
     real(dp), allocatable :: prob_mpi(:,:)

! impurity double occupation number matrix, for mpi case
     real(dp), allocatable :: nnmat_mpi(:,:,:)

! impurity green's function, imaginary time axis, for mpi case
     real(dp), allocatable :: gtau_mpi(:,:,:,:)

! impurity green's function, matsubara frequency axis, for mpi case
     complex(dp), allocatable :: grnf_mpi(:,:,:,:)

! allocate memory
     allocate(hist_mpi(mkink,natom),             stat=istat)
     if ( istat /= 0 ) then
         call ctqmc_print_error('ctqmc_impurity_solver','can not allocate enough memory')
     endif

     allocate(nmat_mpi(norbs,natom),             stat=istat)
     if ( istat /= 0 ) then
         call ctqmc_print_error('ctqmc_impurity_solver','can not allocate enough memory')
     endif

     allocate(prob_mpi(ncfgs,natom),             stat=istat)
     if ( istat /= 0 ) then
         call ctqmc_print_error('ctqmc_impurity_solver','can not allocate enough memory')
     endif

     allocate(nnmat_mpi(norbs,norbs,natom),      stat=istat)
     if ( istat /= 0 ) then
         call ctqmc_print_error('ctqmc_impurity_solver','can not allocate enough memory')
     endif

     allocate(gtau_mpi(ntime,norbs,norbs,natom), stat=istat)
     if ( istat /= 0 ) then
         call ctqmc_print_error('ctqmc_impurity_solver','can not allocate enough memory')
     endif

     allocate(grnf_mpi(mfreq,norbs,norbs,natom), stat=istat)
     if ( istat /= 0 ) then
         call ctqmc_print_error('ctqmc_impurity_solver','can not allocate enough memory')
     endif

! setup cstep
     cstep = 0

! setup cflag, check the status of quantum impurity solver periodically
     cflag = 1

! setup timer
     time_iter = zero
     time_niter = zero

! setup nsweep
! whether it is time to enter QMC data accumulating mode
     if ( iter == 999 ) then
         nsweep = nsweep * 10
         nwrite = nwrite * 10
     endif

!=========================================================================
!>>> starting quantum impurity solver                                  <<<
!=========================================================================

! print the header of continuous time quantum Monte Carlo quantum impurity solver
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(2X,a)') 'AZALEA >>> CTQMC quantum impurity solver running'
         write(mystd,*)''
         write(mystd,'(4X,a,i10)') 'nband :', nband
         write(mystd,'(4X,a,i10)') 'nspin :', nspin
         write(mystd,*)
     endif

!=========================================================================
!>>> initializing quantum impurity solver                              <<<
!=========================================================================

! init the continuous time quantum Monte Carlo quantum impurity solver
! setup the key variables
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a)') 'quantum impurity solver initializing'
     endif

     call cpu_time(time_begin) ! record starting time
! careful here, arrays are reset etc...
     call ctqmc_solver_init(atom)
     call cpu_time(time_end)   ! record ending   time

! print the time information
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
         write(mystd,*)
     endif

!=========================================================================
!>>> retrieving quantum impurity solver                                <<<
!=========================================================================

! init the continuous time quantum Monte Carlo quantum impurity solver further
! retrieving the time series information produced by previous running
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a)') 'quantum impurity solver retrieving'
     endif

     call cpu_time(time_begin) ! record starting time
     call ctqmc_retrieve_status(atom)
     call cpu_time(time_end)   ! record ending   time

! print the time information
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
         write(mystd,*)
     endif

!=========================================================================
!>>> warmming quantum impurity solver                                  <<<
!=========================================================================

! warmup the continuous time quantum Monte Carlo quantum impurity solver,
! in order to achieve equilibrium state
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a)') 'quantum impurity solver warmming'
     endif

     call cpu_time(time_begin) ! record starting time
     call ctqmc_diagram_warmming(atom)
     call cpu_time(time_end)   ! record ending   time

! print the time information
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
         write(mystd,*)
     endif

!=========================================================================
!>>> beginning main iteration                                          <<<
!=========================================================================

! start simulation
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a)') 'quantum impurity solver sampling'
         write(mystd,*)
     endif

     CTQMC_MAIN_ITERATION: do i=1, nsweep, nwrite

! record start time
         call cpu_time(time_begin)

         CTQMC_DUMP_ITERATION: do j=1, nwrite

!=========================================================================
!>>> sampling perturbation expansion series                            <<<
!=========================================================================

! increase cstep by 1
             cstep = cstep + 1

! sampling the perturbation expansion feynman diagrams randomly
             if ( beta > one ) then
                 call ctqmc_diagram_sampling(cstep,atom) ! at normal region
             else
                 call ctqmc_diagram_templing(cstep,atom) ! at high temperature region 
             endif ! back if ( beta > one ) block

!=========================================================================
!>>> sampling the physical observables                                 <<<
!=========================================================================

! record the histogram for perturbation expansion series
             call ctqmc_record_hist(atom)

! record the impurity (double) occupation number matrix and other
! auxiliary physical observables
             if ( mod(cstep, nmonte) == 0 ) then
                 call ctqmc_record_nmat(atom)
             endif

! record the impurity green's function in matsubara frequency space
             if ( mod(cstep, nmonte) == 0 ) then
                 call ctqmc_record_grnf(atom)
             endif

! record the probability of eigenstates
             if ( mod(cstep, ncarlo) == 0 ) then
                 call ctqmc_record_prob(atom)
             endif

! record the impurity green's function in imaginary time space
             if ( mod(cstep, ncarlo) == 0 ) then
                 call ctqmc_record_gtau(atom)
             endif

         enddo CTQMC_DUMP_ITERATION ! over j={1,nwrite} loop

!=========================================================================
!>>> reporting quantum impurity solver                                 <<<
!=========================================================================

! it is time to write out the statistics results
         if ( myid == master ) then ! only master node can do it

! about iteration number
             write(mystd,'(2X,a,i3,2(a,i15))') 'AZALEA >>> iter:',       &
                                   iter, ' sweep:', cstep, ' of ', nsweep

! about auxiliary physical observables
             istat = cstep / nmonte
             statfact(atom) = istat

             write(mystd,'(4X,a)')        'auxiliary system observables:'
             write(mystd,'(2(4X,a,f10.5))') 'etot :', paux(1,atom) / istat,   &
                                            'epot :', paux(2,atom) / istat
             write(mystd,'(2(4X,a,f10.5))') 'ekin :', paux(3,atom) / istat,   &
                                            '<Sz> :', paux(4,atom) / istat

! about insert action
             if ( insert_tcount <= half ) insert_tcount = -one ! if insert is disable
             write(mystd,'(4X,a)')        'insert kink statistics:'
             write(mystd,'(4X,a,3i10)')   'count:',                      &
                                           int( insert_tcount ),         &
                                           int( insert_accept ),         &
                                           int( insert_reject )
             write(mystd,'(4X,a,3f10.5)') 'ratio:', one,                 &
                                           insert_accept / insert_tcount,&
                                           insert_reject / insert_tcount

! about remove action
             if ( remove_tcount <= half ) remove_tcount = -one ! if remove is disable
             write(mystd,'(4X,a)')        'remove kink statistics:'
             write(mystd,'(4X,a,3i10)')   'count:',                      &
                                           int( remove_tcount ),         &
                                           int( remove_accept ),         &
                                           int( remove_reject )
             write(mystd,'(4X,a,3f10.5)') 'ratio:', one,                 &
                                           remove_accept / remove_tcount,&
                                           remove_reject / remove_tcount

! about lshift action
             if ( lshift_tcount <= half ) lshift_tcount = -one ! if lshift is disable
             write(mystd,'(4X,a)')        'lshift kink statistics:'
             write(mystd,'(4X,a,3i10)')   'count:',                      &
                                           int( lshift_tcount ),         &
                                           int( lshift_accept ),         &
                                           int( lshift_reject )
             write(mystd,'(4X,a,3f10.5)') 'ratio:', one,                 &
                                           lshift_accept / lshift_tcount,&
                                           lshift_reject / lshift_tcount

! about rshift action
             if ( rshift_tcount <= half ) rshift_tcount = -one ! if rshift is disable
             write(mystd,'(4X,a)')        'rshift kink statistics:'
             write(mystd,'(4X,a,3i10)')   'count:',                      &
                                           int( rshift_tcount ),         &
                                           int( rshift_accept ),         &
                                           int( rshift_reject )
             write(mystd,'(4X,a,3f10.5)') 'ratio:', one,                 &
                                           rshift_accept / rshift_tcount,&
                                           rshift_reject / rshift_tcount

! about reswap action
             if ( reswap_tcount <= half ) reswap_tcount = -one ! if reswap is disable
             write(mystd,'(4X,a)')        'global swap statistics:'
             write(mystd,'(4X,a,3i10)')   'count:',                      &
                                           int( reswap_tcount ),         &
                                           int( reswap_accept ),         &
                                           int( reswap_reject )
             write(mystd,'(4X,a,3f10.5)') 'ratio:', one,                 &
                                           reswap_accept / reswap_tcount,&
                                           reswap_reject / reswap_tcount

! about reflip action
             if ( reflip_tcount <= half ) reflip_tcount = -one ! if reflip is disable
             write(mystd,'(4X,a)')        'global flip statistics:'
             write(mystd,'(4X,a,3i10)')   'count:',                      &
                                           int( reflip_tcount ),         &
                                           int( reflip_accept ),         &
                                           int( reflip_reject )
             write(mystd,'(4X,a,3f10.5)') 'ratio:', one,                 &
                                           reflip_accept / reflip_tcount,&
                                           reflip_reject / reflip_tcount

         endif ! back if ( myid == master ) block

!=========================================================================
!>>> reducing immediate results                                        <<<
!=========================================================================

! collect the histogram data from hist to hist_mpi
         call ctqmc_reduce_hist(hist_mpi)

! collect the impurity green's function data from gtau to gtau_mpi
         call ctqmc_reduce_gtau(gtau_mpi)

!=========================================================================
!>>> symmetrizing immediate results                                    <<<
!=========================================================================

! symmetrize the impurity green's function over spin or over bands
         if ( issun == 2 .or. isspn == 1 ) then
             call ctqmc_symm_gtau(symm, gtau_mpi,atom)
         endif

!=========================================================================
!>>> writing immediate results                                         <<<
!=========================================================================

! write out the histogram data, hist_mpi
         if ( myid == master ) then ! only master node can do it
             call ctqmc_dump_hist(hist_mpi,atom)
         endif

! gtau_mpi need to be scaled properly before written
         do m=1,norbs
             do n=1,ntime
                 gtau_mpi(n,m,m,atom) = gtau_mpi(n,m,m,atom) * dble(ncarlo) / dble(cstep)
             enddo ! over n={1,ntime} loop
         enddo ! over m={1,norbs} loop

! write out the impurity green's function, gtau_mpi
         if ( myid == master ) then ! only master node can do it
             if ( iter /= 999 ) then
                 call ctqmc_dump_gtau(tmesh, gtau_mpi,atom)
             else
                 call ctqmc_dump_gbin(cstep / nwrite, tmesh, gtau_mpi)
                 write(mystd,'(4X,a)') '>>> quantum impurity solver status: binned'
             endif
         endif ! back if ( myid == master ) block

!=========================================================================
!>>> checking quantum impurity solver                                  <<<
!=========================================================================

         call ctqmc_diagram_checking(cflag,atom)

!=========================================================================
!>>> timing quantum impurity solver                                    <<<
!=========================================================================

! record ending time for this iteration
         call cpu_time(time_end)

! calculate timing information
         time_iter = time_end - time_begin
         time_niter = time_niter + time_iter
         time_begin = time_end

! print out the result
         if ( myid == master ) then ! only master node can do it
             call ctqmc_time_analyzer(time_iter, time_niter)
             write(mystd,*)
         endif

!=========================================================================
!>>> escaping quantum impurity solver                                  <<<
!=========================================================================

! if the quantum impurity solver is out of control or reaches convergence
         if ( cflag == 99 .or. cflag == 100 ) then
             EXIT CTQMC_MAIN_ITERATION ! jump out the iteration
         endif

     enddo CTQMC_MAIN_ITERATION ! over i={1,nsweep} loop

!=========================================================================
!>>> ending main iteration                                             <<<
!=========================================================================

!=========================================================================
!>>> reducing final results                                            <<<
!=========================================================================

! collect the histogram data from hist to hist_mpi
     call ctqmc_reduce_hist(hist_mpi)

! collect the (double) occupation matrix data from nmat to nmat_mpi, from
! nnmat to nnmat_mpi
     call ctqmc_reduce_nmat(nmat_mpi, nnmat_mpi)

! collect the probability data from prob to prob_mpi
     call ctqmc_reduce_prob(prob_mpi)

! collect the impurity green's function data from gtau to gtau_mpi
     call ctqmc_reduce_gtau(gtau_mpi)

! collect the impurity green's function data from grnf to grnf_mpi
     call ctqmc_reduce_grnf(grnf_mpi)

! update original data and calculate the averages simultaneously
     hist(:,:)  = hist_mpi(:,:)

     nmat(:,:)  = nmat_mpi(:,:)  * dble(nmonte) / dble(nsweep)
     prob(:,:)  = prob_mpi(:,:)  * dble(ncarlo) / dble(nsweep)

     do m=1,norbs
         do n=1,ntime
             gtau(n,m,m,atom) = gtau_mpi(n,m,m,atom) * dble(ncarlo) / dble(nsweep)
         enddo ! over n={1,ntime} loop
     enddo ! over m={1,norbs} loop

     do m=1,norbs
         do n=1,nfreq
             grnf(n,m,m,atom) = grnf_mpi(n,m,m,atom) * dble(nmonte) / dble(nsweep)
         enddo ! over n={1,nfreq} loop
     enddo ! over m={1,norbs} loop

     do m=1,norbs
         do n=1,norbs
             nnmat(n,m,atom) = nnmat_mpi(n,m,atom)   * dble(nmonte) / dble(nsweep)
         enddo ! over n={1,norbs} loop
     enddo ! over m={1,norbs} loop

!=========================================================================
!>>> symmetrizing final results                                        <<<
!=========================================================================

! symmetrize the occupation number matrix (nmat) over spin or over bands
     if ( issun == 2 .or. isspn == 1 ) then
         call ctqmc_symm_nmat(symm, nmat,atom)
     endif

! symmetrize the impurity green's function (gtau) over spin or over bands
     if ( issun == 2 .or. isspn == 1 ) then
         call ctqmc_symm_gtau(symm, gtau,atom)
     endif

! build high frequency self-energy function using moment expansion, and
! then adjust impurity green's function by using dyson's equation. it is
! suitable for weak correlation region
!     call ctqmc_build_selfenergy1(atom)

! build atomic green's function and self-energy function using improved
! Hubbard-I approximation, and then make interpolation for self-energy
! function between low frequency QMC data and high frequency Hubbard-I
! approximation data, the impurity green's function can be obtained by
! using dyson's equation finally
     call ctqmc_build_selfenergy2(atom)

! symmetrize the impurity green's function (grnf) over spin or over bands
     if ( issun == 2 .or. isspn == 1 ) then
         call ctqmc_symm_grnf(symm, grnf, atom)
     endif

! symmetrize the impurity self-energy function (sig2) over spin or over bands
     if ( issun == 2 .or. isspn == 1 ) then
         call ctqmc_symm_grnf(symm, sig2, atom)
     endif

!=========================================================================
!>>> writing final results                                             <<<
!=========================================================================

! write out the final histogram data, hist
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_hist(hist,atom)
     endif

! write out the final (double) occupation matrix data, nmat and nnmat
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_nmat(nmat, nnmat,atom)
     endif

! write out the final probability data, prob
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_prob(prob,atom)
     endif

! write out the final impurity green's function data, gtau
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_gtau(tmesh, gtau,atom)
     endif

! write out the final impurity green's function data, grnf
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_grnf(rmesh, grnf,atom)
     endif

! write out the final self-energy function data, sig2
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_sigf(rmesh, sig2,atom)
     endif

!=========================================================================
!>>> saving quantum impurity solver                                    <<<
!=========================================================================

! save the perturbation expansion series information to the disk file
     if ( myid == master ) then ! only master node can do it
         call ctqmc_save_status(atom)
     endif

!=========================================================================
!>>> finishing quantum impurity solver                                 <<<
!=========================================================================

! print the footer of continuous time quantum Monte Carlo quantum impurity solver
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(2X,a)') 'AZALEA >>> CTQMC quantum impurity solver shutdown'
         write(mystd,*)
     endif

! deallocate memory
     deallocate(hist_mpi)
     deallocate(nmat_mpi)
     deallocate(prob_mpi)
     deallocate(gtau_mpi)
     deallocate(grnf_mpi)
     deallocate(nnmat_mpi)

     return
  end subroutine ctqmc_impurity_solver

!>>> perform thermalization on perturbation expansion series to achieve
! thermodynamics equilibrium state
  subroutine ctqmc_diagram_warmming(atom)
     use constants, only : zero,int64
     use control, only : ntherm
     use context

     implicit none

! local variables
! loop index
     integer(int64) :: i

     integer, intent(in) :: atom

! warm up the diagram series
     do i=1,ntherm
         call ctqmc_diagram_sampling(i,atom)
     enddo ! over i={1,ntherm} loop

! reinit statistics variables
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

     return
  end subroutine ctqmc_diagram_warmming

!>>> visit the perturbation expansion diagrams randomly
  subroutine ctqmc_diagram_sampling(cstep,atom)
     use constants, only : dp, int64
     use control, only : nflip, nclean

     use spring

     implicit none

! external arguments
! current QMC sweep steps
     integer(int64), intent(in) :: cstep

! we need
     integer, intent(in) :: atom

! change the order of perturbation expansion series
     if ( spring_sfmt_stream() < 0.9_dp ) then
         if ( spring_sfmt_stream() > 0.5_dp ) then
             call ctqmc_insert_kink(atom)  ! insert one new kink
         else
             call ctqmc_remove_kink(atom)  ! remove one old kink
         endif
! do not change the order of perturbation expansion series
     else
         if ( spring_sfmt_stream() > 0.5_dp ) then
             call ctqmc_lshift_kink(atom)  ! shift the left  endpoints
         else
             call ctqmc_rshift_kink(atom)  ! shift the right endpoints
         endif
     endif ! back if ( spring_sfmt_stream() < 0.9_dp ) block

! numerical trick: perform global spin flip periodically
! note: ctqmc_reflip_kink(1), flip inter-orbital spins randomly
     if ( nflip > 0  .and. mod(cstep,  nflip) == 0 ) then
         if ( spring_sfmt_stream() > 0.5_dp ) then
             call ctqmc_reflip_kink(2,atom) ! flip intra-orbital spins one by one
         else
             call ctqmc_reflip_kink(3,atom) ! flip intra-orbital spins globally
         endif
     endif

! numerical trick: perform global update periodically
     if ( nclean > 0 .and. mod(cstep, nclean) == 0 ) then
         call ctqmc_reload_kink(atom)
     endif

     return
  end subroutine ctqmc_diagram_sampling

!>>> visit the perturbation expansion diagrams randomly at very high temperature
  subroutine ctqmc_diagram_templing(cstep,atom)
     use constants, only : dp, int64
     use control, only : nflip, nclean

     use spring

     implicit none

! external arguments
! current QMC sweep steps
     integer(int64), intent(in) :: cstep
     
     integer, intent(in) :: atom
     
! change the order of perturbation expansion series
     if ( spring_sfmt_stream() < 0.1_dp ) then
         if ( spring_sfmt_stream() > 0.5_dp ) then
             call ctqmc_insert_kink(atom)  ! insert one new kink
         else
             call ctqmc_remove_kink(atom)  ! remove one old kink
         endif
! do not change the order of perturbation expansion series
     else
         if ( spring_sfmt_stream() > 0.8_dp ) then
             call ctqmc_lshift_kink(atom)  ! shift the left  endpoints
         else
             call ctqmc_reswap_kink(atom)  ! swap creator and destroyer
         endif
     endif ! back if ( spring_sfmt_stream() < 0.9_dp ) block

! numerical trick: perform global spin flip periodically
! note: ctqmc_reflip_kink(1), flip inter-orbital spins randomly
     if ( nflip > 0  .and. mod(cstep,  nflip) == 0 ) then
         if ( spring_sfmt_stream() > 0.5_dp ) then
             call ctqmc_reflip_kink(2,atom) ! flip intra-orbital spins one by one
         else
             call ctqmc_reflip_kink(3,atom) ! flip intra-orbital spins globally
         endif
     endif

! numerical trick: perform global update periodically
     if ( nclean > 0 .and. mod(cstep, nclean) == 0 ) then
         call ctqmc_reload_kink(atom)
     endif

     return
  end subroutine ctqmc_diagram_templing

!>>> checking whether the quantum impurity solver is consistent internally
  subroutine ctqmc_diagram_checking(cflag,atom)
     use constants
     use control
     use context

     implicit none

! external arguments
! control flag
     integer, intent(inout) :: cflag

     integer, intent(in) :: atom

! local variables
! loop index
     integer :: i
     integer :: j

     if ( cflag == 1 ) then

! check perturbation expansion order
         do i=1,norbs
             if ( stts(i,atom) == 0 ) then
                 if ( rank(i,atom) /= 0 ) cflag = 99
             endif

             if ( stts(i,atom) == 1 ) then
                 if ( rank(i,atom) == 0 ) cflag = 99
             endif

             if ( stts(i,atom) == 2 ) then
                 if ( rank(i,atom) == 0 ) cflag = 99
             endif

             if ( stts(i,atom) == 3 ) then
                 if ( rank(i,atom) /= 0 ) cflag = 99
             endif
         enddo ! over i={1,norbs} loop

! check time order of operators
         do i=1,norbs
             do j=1,rank(i,atom)-1
                 if ( time_s( index_s(j, i), i ) > time_s( index_s(j+1, i), i ) ) then
                     cflag = 99
                 endif
                 if ( time_e( index_e(j, i), i ) > time_e( index_e(j+1, i), i ) ) then
                     cflag = 99
                 endif
             enddo ! over j={1,rank(i,atom)-1} loop
         enddo ! over i={1,norbs} loop

! check segment and anti-segment
         do i=1,norbs
             if ( stts(i,atom) == 1 ) then
                 if ( time_s( index_s(1, i), i ) > time_e( index_e(1, i), i ) ) then
                     cflag = 99
                 endif
             endif ! back if ( stts(i,atom) == 1 ) block

             if ( stts(i,atom) == 2 ) then
                 if ( time_s( index_s(1, i), i ) < time_e( index_e(1, i), i ) ) then
                     cflag = 99
                 endif
             endif ! back if ( stts(i,atom) == 2 ) block
         enddo ! over i={1,norbs} loop

! write the results, only master node can do it
         if ( myid == master ) then
             if ( cflag == 99 ) then
                 write(mystd,'(4X,a)') '>>> quantum impurity solver status: error?'
                 write(mystd,'(4X,a)') '>>> please check the status file: solver.status.dat'
                 call ctqmc_save_status(atom)
                 call ctqmc_print_error('ctqmc_diagram_checking','unknown fatal error occur')
             else
                 write(mystd,'(4X,a)') '>>> quantum impurity solver status: normal'
             endif
         endif ! back if ( myid == master ) block

     endif ! back if ( cflag == 1 ) block

     return
  end subroutine ctqmc_diagram_checking

!>>> testing subroutine, please active it on ctqmc_diagram_sampling()
  subroutine ctqmc_impurity_tester()
     use constants
     use control
     use context

     implicit none

!-------------------------------------------------------------------------
! insert your debug code here
!-------------------------------------------------------------------------

     call ctqmc_make_display(2)
     call ctqmc_print_error('ctqmc_impurity_tester','in debug mode')

     return
  end subroutine ctqmc_impurity_tester
