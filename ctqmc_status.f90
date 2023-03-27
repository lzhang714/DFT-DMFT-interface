!
! MK: modified, at least to account for multiple atoms, potentially also in other ways.
! Original header:
!
!-------------------------------------------------------------------------
! project : azalea
! program : ctqmc_save_status
!           ctqmc_retrieve_status
! source  : ctqmc_status.f90
! type    : subroutine
! author  : li huang (email:huangli712@yahoo.com.cn)
! history : 09/23/2009 by li huang
!           09/26/2009 by li huang
!           10/10/2009 by li huang
!           10/20/2009 by li huang
!           10/29/2009 by li huang
!           11/01/2009 by li huang
!           11/10/2009 by li huang
!           11/29/2009 by li huang
!           12/01/2009 by li huang
!           12/02/2009 by li huang
!           12/26/2009 by li huang
!           02/21/2010 by li huang
! purpose : save or retrieve the perturbation expansion series information
!           to or from the status file for hybridization expansion version
!           continuous time quantum Monte Carlo (CTQMC) quantum impurity
!           solver.
!           it can be used to save the computational time to achieve the
!           equilibrium state
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

!>>> save the current perturbation expansion series information for the
! continuous time quantum Monte Carlo quantum impurity solver
  subroutine ctqmc_save_status(atom)
     use constants
     use control
     use context

     implicit none
     
     integer,intent(in) :: atom

! local variables
! loop index over orbitals
     integer :: i

! loop index over segments
     integer :: j
     
     character(20) :: filename
     
        write(filename,*) atom
        filename = 'solver.status_'//trim(adjustl(filename))//'.dat'
        open(mytmp, file=filename, form='formatted', status='unknown')

! dump the segments
     do i=1,norbs
         write(mytmp,'(a,i4)') '# flavor     :', i
         write(mytmp,'(a,i4)') '# status     :', stts(i,atom)

! write out start point values for segments (create  operators)
         write(mytmp,'(a,i4)') '# time_s data:', rank(i,atom)
         do j=1,rank(i,atom)
             write(mytmp,'(2i4,f12.6)') i, j, time_s( index_s(j, i), i )
         enddo ! over j={1,rank(i,atom)} loop

! write out end   point values for segments (destroy operators)
         write(mytmp,'(a,i4)') '# time_e data:', rank(i,atom)
         do j=1,rank(i,atom)
             write(mytmp,'(2i4,f12.6)') i, j, time_e( index_e(j, i), i )
         enddo ! over j={1,rank(i,atom)} loop

         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,norbs} loop

! close the file handler
     close(mytmp)

     return
  end subroutine ctqmc_save_status

!>>> retrieve the perturbation expansion series information to initialize
! the continuous time quantum Monte Carlo quantum impurity solver
  subroutine ctqmc_retrieve_status(atom)
     use constants
     use control
     use context

     use mmpi

     implicit none

     integer, intent(in) :: atom

! local variables
! loop index
     integer  :: i
     integer  :: j

! dummy integer variables
     integer  :: i1
     integer  :: j1

! used to check whether the input file (solver.status.dat) exists
     logical  :: exists

! dummy character variables
     character(14) :: chr

! determinant ratio for insert segments
     real(dp) :: deter_ratio

! dummy variables, used to store imaginary time points
     real(dp) :: tau_s(mkink,norbs)
     real(dp) :: tau_e(mkink,norbs)
     
     character(20) :: filename

! initialize variables
     exists = .false.

     tau_s = zero
     tau_e = zero

        write(filename,*) atom
        filename = 'solver.status_'//trim(adjustl(filename))//'.dat'

! inquire file status: solver.status.dat, only master node can do it
     if ( myid == master ) then
         inquire (file = filename, exist = exists)
     endif

! broadcast exists from master node to all children nodes
# if defined (MPI)

! broadcast data
     call mp_bcast( exists, master )

! block until all processes have reached here
     call mp_barrier()

# endif  /* MPI */

! if solver.status.dat does not exist, return parent subroutine immediately
     if ( exists .eqv. .false. ) return

! read solver.status.dat, only master node can do it
     if ( myid == master ) then

! open the status file
        open(mytmp, file=filename, form='formatted', status='unknown')

! read in key data
         do i=1,norbs
             read(mytmp, '(a14,i4)') chr, i1
             read(mytmp, '(a14,i4)') chr, cstat

             read(mytmp, '(a14,i4)') chr, ckink
             do j=1,ckink
                 read(mytmp,*) i1, j1, tau_s(j, i)
             enddo ! over j={1,ckink} loop

             read(mytmp, '(a14,i4)') chr, ckink
             do j=1,ckink
                 read(mytmp,*) i1, j1, tau_e(j, i)
             enddo ! over j={1,ckink} loop

             read(mytmp,*) ! skip two lines
             read(mytmp,*)

             stts(i,atom) = cstat
             rank(i,atom) = ckink
         enddo ! over i={1,norbs} loop

! close the status file
         close(mytmp)

     endif ! back if ( myid == master ) block

! broadcast rank, stts, tau_s, and tau_e from master node to all children nodes
# if defined (MPI)

! broadcast data
     call mp_bcast( rank,  master )
     call mp_bcast( stts,  master )

! block until all processes have reached here
     call mp_barrier()

! broadcast data
     call mp_bcast( tau_s, master )
     call mp_bcast( tau_e, master )

! block until all processes have reached here
     call mp_barrier()

# endif  /* MPI */

! check the validity of tau_s
     if ( maxval(tau_s) > beta ) then
         call ctqmc_print_error('ctqmc_retrieve_status','the retrieved tau_s data are not correct')
     endif

! check the validity of tau_e
     if ( maxval(tau_e) > beta ) then
         call ctqmc_print_error('ctqmc_retrieve_status','the retrieved tau_e data are not correct')
     endif

! restore all the segments or anti-segments
     do i=1,norbs

! segment scheme
         if ( stts(i,atom) == 1 ) then
             do j=1,rank(i,atom)
                 ckink = j - 1 ! update ckink simultaneously
                 call cat_insert_detrat(i, tau_s(j, i), tau_e(j, i), deter_ratio,atom)
                 call cat_insert_matrix(i, j, j, tau_s(j, i), tau_e(j, i), deter_ratio,atom)
             enddo ! over j={1,rank(i,atom)} loop
         endif ! back if ( stts(i,atom) == 1 ) block

! anti-segment scheme
         if ( stts(i,atom) == 2 ) then
             do j=1,rank(i,atom)-1
                 ckink = j - 1 ! update ckink simultaneously
                 call cat_insert_detrat(i, tau_s(j, i), tau_e(j+1, i), deter_ratio,atom)
                 call cat_insert_matrix(i, j, j, tau_s(j, i), tau_e(j+1, i), deter_ratio,atom)
             enddo ! over j={1,rank(i,atom)-1} loop
             ckink = rank(i,atom) - 1
             call cat_insert_detrat(i, tau_s(ckink+1, i), tau_e(1, i), deter_ratio,atom)
             call cat_insert_matrix(i, ckink+1, 1, tau_s(ckink+1, i), tau_e(1, i), deter_ratio,atom)
         endif ! back if ( stts(i,atom) == 1 ) block

     enddo ! over i={1,norbs} loop

     return
  end subroutine ctqmc_retrieve_status
