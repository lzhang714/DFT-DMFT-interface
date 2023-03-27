!
! MK: modified, at least to account for multiple atoms, potentially also in other ways.
! Original header:
!
!-------------------------------------------------------------------------
! project : azalea
! program : ctqmc_print_header
!           ctqmc_print_footer
!           ctqmc_print_summary
!           ctqmc_print_error
!           ctqmc_print_exception
! source  : ctqmc_print.f90
! type    : subroutines
! author  : li huang (email:huangli712@yahoo.com.cn)
! history : 09/15/2009 by li huang
!           09/20/2009 by li huang
!           12/01/2009 by li huang
!           02/21/2010 by li huang
! purpose : provide printing infrastructure for hybridization expansion
!           version continuous time quantum Monte Carlo (CTQMC) quantum
!           impurity solver
! input   :
! output  :
! status  : very unstable
! comment :
!-------------------------------------------------------------------------

!>>> print the startup information for continuous time quantum Monte Carlo
! quantum impurity solver plus dynamical mean field theory self-consistent
! engine
  subroutine ctqmc_print_header()
     use constants
     use control, only : nprocs

     implicit none

     write(mystd,'(2X,a)') 'DFT + DMFT interface code.'    '
     write(mystd,'(2X,a)') 'Segment solver version: 17 July 2015  '    '

     write(mystd,'(2X,a)') 'Interface developed by Long Zhang @ UF'    '

     write(mystd,'(2X,a)') 'Segment Solver developed by Li Huang, CAEP & IOP'
     write(mystd,'(2X,a)') 'Segment Solver license: GPL2 and later versions' 

     write(mystd,*)

# if defined (MPI)

     write(mystd,'(2X,a,i3)') '>>> parallelism: Yes >>> processors:', nprocs

# else   /* MPI */

     write(mystd,'(2X,a,i3)') '>>> parallelism: No  >>> processors:', 1

# endif  /* MPI */

     write(mystd,*)

     return
  end subroutine ctqmc_print_header

!>>> print the ending information for continuous time quantum Monte Carlo
! quantum impurity solver plus dynamical mean field theory self-consistent
! engine
  subroutine ctqmc_print_footer()
     use constants

     implicit none

! used to record the time information
     real(dp) :: tot_time

! obtain time information
     call cpu_time(tot_time)

     write(mystd,*)
     write(mystd,'(2X,a,f10.2,a)') 'total time spent:', tot_time, 's'
     write(mystd,*)

     write(mystd,'(2X,a)') '>>> I am tired and want to go to bed. You should do the same ;-) Bye!'
     write(mystd,'(2X,a)') '>>> ending'

     return
  end subroutine ctqmc_print_footer

!>>> print the running parameters, only for reference
  subroutine ctqmc_print_summary()
     use constants
     use control

     implicit none

     write(mystd,'(2X,a)') 'QMC >>> parameters list:'

     write(mystd,'(1(4X,a,i15))')   'isbin :', isbin
     write(mystd,'(2(4X,a,i15))')   'issun :', issun  , 'isspn :', isspn

     write(mystd,'(2(4X,a,i15))')   'mkink :', mkink  , 'mfreq :', mfreq
     write(mystd,'(2(4X,a,i15))')   'nband :', nband  , 'nspin :', nspin
     write(mystd,'(2(4X,a,i15))')   'norbs :', norbs  , 'ncfgs :', ncfgs
     write(mystd,'(2(4X,a,i15))')   'niter :', niter  , 'nfreq :', nfreq
     write(mystd,'(2(4X,a,i15))')   'ntime :', ntime  , 'nflip :', nflip

     write(mystd,'(2(4X,a,i15))')   'ntherm:', ntherm , 'nsweep:', nsweep
     write(mystd,'(2(4X,a,i15))')   'nclean:', nclean , 'nwrite:', nwrite
     write(mystd,'(2(4X,a,i15))')   'nmonte:', nmonte , 'ncarlo:', ncarlo

     write(mystd,'(1(4X,a,f15.5))')  'temp  :', ev2k/beta

     write(mystd,*)

     return
  end subroutine ctqmc_print_summary

!>>> print the error information and STOP the program
  subroutine ctqmc_print_error(sub, msg)
     use constants

     implicit none

! external arguments
! subroutine name
     character(len=*), intent(in) :: sub

! error message
     character(len=*), intent(in) :: msg

! print error information
     write(mystd,'(2X,4a)') 'fatal error occurred in ', sub, ': ', msg

! TERMINATE THE PROGRAM
!-------------------------------------------------------------------------
     STOP
!-------------------------------------------------------------------------

     return
  end subroutine ctqmc_print_error

!>>> print normal runtime exceptional information, and continue
  subroutine ctqmc_print_exception(sub, msg)
     use constants

     implicit none

! external arguments
! subroutine name
     character(len=*), intent(in) :: sub

! exception message
     character(len=*), intent(in) :: msg

! print error information
     write(mystd,'(2X,4a)') 'runtime exception occurred in ', sub, ': ', msg

! CONTINUE/PAUSE THE PROGRAM
!-------------------------------------------------------------------------
     CONTINUE ! OR PAUSE
!-------------------------------------------------------------------------

     return
  end subroutine ctqmc_print_exception
