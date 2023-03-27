!
! MK: modified, multiple atoms, ED things, umatrix etc.etc.
! Original header:
!
!-------------------------------------------------------------------------
! project : azalea
! program : ctqmc_dump_gtau
!           ctqmc_dump_wtau
!           ctqmc_dump_htau
!           ctqmc_dump_gbin
!           ctqmc_dump_grnf
!           ctqmc_dump_wssf
!           ctqmc_dump_hybf
!           ctqmc_dump_sigf
!           ctqmc_dump_hist
!           ctqmc_dump_nmat
!           ctqmc_dump_prob
! source  : ctqmc_dump.f90
! type    : subroutine
! author  : li huang (email:huangli712@yahoo.com.cn)
! history : 09/16/2009 by li huang
!           09/17/2009 by li huang
!           09/18/2009 by li huang
!           09/20/2009 by li huang
!           09/22/2009 by li huang
!           10/25/2009 by li huang
!           11/01/2009 by li huang
!           11/30/2009 by li huang
!           12/01/2009 by li huang
!           12/04/2009 by li huang
!           12/09/2009 by li huang
!           12/26/2009 by li huang
!           12/30/2009 by li huang
!           02/28/2010 by li huang
!           03/04/2010 by li huang
!           08/23/2010 by li huang
! purpose : dump key observables produced by the hybridization expansion
!           version continuous time quantum Monte Carlo (CTQMC) quantum
!           impurity solver and dynamical mean field theory (DMFT) self
!           -consistent engine to disk files
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

!>>> write out impurity green's function in imaginary time space
  subroutine ctqmc_dump_gtau(tmesh,gtau,atom)
     use constants
     use control
     use plo_quantities, only : natom

     implicit none

! external arguments
! imaginary time mesh
     real(dp), intent(in) :: tmesh(ntime)

     integer, intent(in) :: atom

     character(20) :: filename

! impurity green's function
     real(dp), intent(in) :: gtau(ntime,norbs,norbs,natom)

! local variables
! loop index
     integer  :: i
     integer  :: j

! dummy variables
     real(dp) :: raux

! scaled impurity green's function
     real(dp) :: gaux(ntime,norbs,norbs)

! evaluate gaux first
     raux = real(ntime) / (beta * beta)

     do i=1,norbs
         do j=1,ntime
             gaux(j,i,i) = gtau(j,i,i,atom) * raux
         enddo ! over j={1,ntime} loop
     enddo ! over i={1,norbs} loop

        write(filename,*) atom
        filename = 'solver.green_'//trim(adjustl(filename))//'.dat'
        open(mytmp, file=filename, form='formatted', status='unknown')

! write it
     do i=1,nband
         do j=1,ntime
             write(mytmp,'(2i5,3f12.6)') i, j, tmesh(j), gaux(j,i,i), gaux(j,i+nband,i+nband)
         enddo ! over j={1,ntime} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,nband} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_gtau


!>>> write out Bloch green's function in imaginary time space
!  subroutine ctqmc_dump_gbtau(tmesh, gbtau)
!     use constants
!     use control

!     implicit none

! external arguments
! imaginary time mesh
!     real(dp), intent(in) :: tmesh(ntime)

! Bloch Greens function
!     real(dp), intent(in) :: gbtau(ntime,norbs,norbs)

! local variables
! loop index
!     integer  :: i
!     integer  :: j

! dummy variables
!     real(dp) :: raux

! scaled impurity green's function
!     real(dp) :: gaux(ntime,norbs,norbs)

! evaluate gaux first
!     raux = real(ntime) / (beta * beta)
!     do i=1,norbs
!         do j=1,ntime
!             gaux(j,i,i) = gtau(j,i,i) * raux
!         enddo ! over j={1,ntime} loop
!     enddo ! over i={1,norbs} loop

! open data file: solver.green.dat
!     open(mytmp, file='solver.gbloch.dat', form='formatted', status='unknown')

! write it
!     do i=1,nband
!         do j=1,ntime
!             write(mytmp,'(2i5,3f12.6)') i, j, tmesh(j), gaux(j,i,i), gaux(j,i+nband,i+nband)
!         enddo ! over j={1,ntime} loop
!         write(mytmp,*) ! write empty lines
!         write(mytmp,*)
!     enddo ! over i={1,nband} loop

! close data file
!     close(mytmp)

!     return
!  end subroutine ctqmc_dump_gbtau

!>>> write out bath weiss's function in imaginary time space
  subroutine ctqmc_dump_wtau(tmesh, wtau,atom)
     use constants
     use control
     use plo_quantities, only : natom

     implicit none

! external arguments
! imaginary time mesh
     real(dp), intent(in) :: tmesh(ntime)
     integer, intent(in) :: atom
     
     character(20) :: filename

! bath weiss's function
     real(dp), intent(in) :: wtau(ntime,norbs,norbs,natom)

! local variables
! loop index
     integer :: i
     integer :: j

        write(filename,*) atom
        filename = 'solver.weiss_'//trim(adjustl(filename))//'.dat'
        open(mytmp, file=filename, form='formatted', status='unknown')

! write it
     do i=1,nband
         do j=1,ntime
             write(mytmp,'(2i5,3f12.6)') i, j, tmesh(j), wtau(j,i,i,atom), wtau(j,i+nband,i+nband,atom)
         enddo ! over j={1,ntime} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,nband} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_wtau

!>>> write out hybridization function in imaginary time space
  subroutine ctqmc_dump_htau(tmesh,htau,atom)
     use constants
     use control
     use plo_quantities, only : natom

     implicit none

! external arguments
! imaginary time mesh
     real(dp), intent(in) :: tmesh(ntime)
     
     integer, intent(in) :: atom
     
     character(20) :: filename

! hybridization function
     real(dp), intent(in) :: htau(ntime,norbs,norbs,natom)

! local variables
! loop index
     integer :: i
     integer :: j

    write(filename,*) atom
    filename = 'solver.hybri_'//trim(adjustl(filename))//'.dat'
    open(mytmp, file=filename, form='formatted', status='unknown')

! write it
     do i=1,nband
         do j=1,ntime
             write(mytmp,'(2i5,3f12.6)') i, j, tmesh(j), htau(j,i,i,atom), htau(j,i+nband,i+nband,atom)
         enddo ! over j={1,ntime} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,nband} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_htau

!>>> write out impurity green's function in imaginary time space (binning mode)
  subroutine ctqmc_dump_gbin(ibin, tmesh, gtau)
     use constants
     use control

     implicit none

! external arguments
! current bin index, integer representation
     integer, intent(in)  :: ibin

! imaginary time mesh
     real(dp), intent(in) :: tmesh(ntime)

! impurity green's function
     real(dp), intent(in) :: gtau(ntime,norbs,norbs)

! local variables
! loop index
     integer  :: i
     integer  :: j

! dummy variables
     real(dp) :: raux

! scaled impurity green's function
     real(dp) :: gaux(ntime,norbs,norbs)

! current bin index, string representation
     character(len=10) :: sbin

! evaluate gaux first
     raux = real(ntime) / (beta * beta)
     do i=1,norbs
         do j=1,ntime
             gaux(j,i,i) = gtau(j,i,i) * raux
         enddo ! over j={1,ntime} loop
     enddo ! over i={1,norbs} loop

! open data file: solver.green.bin.x
     write(sbin,'(i10)') ibin ! convert ibin to sbin
     open(mytmp, file='solver.green.bin.'//trim(adjustl(sbin)), form='formatted', status='unknown')

! write it
     do i=1,nband
         do j=1,ntime
             write(mytmp,'(2i5,3f12.6)') i, j, tmesh(j), gaux(j,i,i), gaux(j,i+nband,i+nband)
         enddo ! over j={1,ntime} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,nband} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_gbin

!>>> write out impurity green's function in matsubara frequency space
  subroutine ctqmc_dump_grnf(rmesh, grnf,atom)
     use constants
     use control
     use plo_quantities, only : natom
     implicit none

! external arguments
! matsubara frequency mesh
     real(dp), intent(in) :: rmesh(mfreq)

     integer, intent(in) :: atom

    character(20) :: filename

! impurity green's function
     complex(dp), intent(in) :: grnf(mfreq,norbs,norbs,natom)

! local variables
! loop index
     integer :: i
     integer :: j

        write(filename,*) atom
        filename = 'solver.grn_'//trim(adjustl(filename))//'.dat'
        open(mytmp, file=filename, form='formatted', status='unknown')

! write it
     do i=1,nband
         do j=1,mfreq
             write(mytmp,'(i5,5f16.8)') i, rmesh(j), &
                                  real(grnf(j,i,i,atom)), &
                                 aimag(grnf(j,i,i,atom)), &
                      real(grnf(j,i+nband,i+nband,atom)), &
                     aimag(grnf(j,i+nband,i+nband,atom))
         enddo ! over j={1,mfreq} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,nband} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_grnf

!=========================================================================
!>>>                                                                   <<<
!=========================================================================

!  subroutine lanc_dump_grnf_real(rmesh, grnf,atom)
!  end subroutine lanc_dump_grnf_real

!=========================================================================
!>>>                                                                   <<<
!=========================================================================

!  subroutine lanc_dump_spectra(rmesh,grnf_loc,grnf_bloch,atom)
!  end subroutine lanc_dump_spectra

!=========================================================================
!>>>                                                                   <<<
!=========================================================================

!  subroutine lanc_dump_grnf_lanc(rmesh, grnf,atom)
!  end subroutine lanc_dump_grnf_lanc

!=========================================================================
!>>>                                                                   <<<
!=========================================================================

!
! code bloat at its best, maybe compactify sometime...
!
!  subroutine lanc_dump_grnf_lanc_real(rmesh, grnf,atom)
!  end subroutine lanc_dump_grnf_lanc_real

!=========================================================================
!>>>                                                                   <<<
!=========================================================================

!>>> write out bath weiss's function in matsubara frequency space

  subroutine ctqmc_dump_wssf(rmesh, wssf,atom)

     use constants
     use control
     use plo_quantities, only : natom

     implicit none

! external arguments
! matsubara frequency mesh
     real(dp), intent(in) :: rmesh(mfreq)

     integer, intent(in) :: atom
     
     character(20) :: filename

! bath weiss's function
     complex(dp), intent(in) :: wssf(mfreq,norbs,norbs,natom)

! local variables
! loop index
     integer :: i
     integer :: j

        write(filename,*) atom
        filename = 'solver.wss_'//trim(adjustl(filename))//'.dat'
        open(mytmp, file=filename, form='formatted', status='unknown')

! write it
     do i=1,nband
         do j=1,mfreq
             write(mytmp,'(i5,5f16.8)') i, rmesh(j), &
                                  real(wssf(j,i,i,atom)), &
                                 aimag(wssf(j,i,i,atom)), &
                      real(wssf(j,i+nband,i+nband,atom)), &
                     aimag(wssf(j,i+nband,i+nband,atom))
         enddo ! over j={1,mfreq} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,nband} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_wssf

!=========================================================================
!>>>                                                                   <<<
!=========================================================================

  subroutine ctqmc_dump_wssf_real(rmesh, wssf,atom)

     use constants
     use control
     use plo_quantities, only : natom

     implicit none

! external arguments
! matsubara frequency mesh
     real(dp), intent(in) :: rmesh(mfreq)

     integer, intent(in) :: atom
     
     character(20) :: filename

! bath weiss's function
     complex(dp), intent(in) :: wssf(mfreq,norbs,norbs,natom)

! local variables
! loop index
     integer :: i
     integer :: j

        write(filename,*) atom
        filename = 'solver.wss_re_'//trim(adjustl(filename))//'.dat'
        open(mytmp, file=filename, form='formatted', status='unknown')

! write it
     do i=1,nband
         do j=1,mfreq
             write(mytmp,'(i5,5f16.8)') i, rmesh(j), &
                                  real(wssf(j,i,i,atom)), &
                                 aimag(wssf(j,i,i,atom)), &
                      real(wssf(j,i+nband,i+nband,atom)), &
                     aimag(wssf(j,i+nband,i+nband,atom))
         enddo ! over j={1,mfreq} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,nband} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_wssf_real

!=========================================================================
!>>>                                                                   <<<
!=========================================================================


!>>> write out hybridization function in matsubara frequency space

  subroutine ctqmc_dump_hybf(rmesh,hybf,atom)

     use constants
     use control
     use plo_quantities, only : natom

     implicit none

! external arguments
! matsubara frequency mesh
     real(dp), intent(in) :: rmesh(mfreq)

     integer, intent(in) :: atom
     
     character(20) :: filename

! hybridization function
     complex(dp), intent(in) :: hybf(mfreq,norbs,norbs,natom)

! local variables
! loop index
     integer :: i
     integer :: j

        write(filename,*) atom
        filename = 'solver.hyb_'//trim(adjustl(filename))//'.dat'
        open(mytmp, file=filename, form='formatted', status='unknown')

! write it
     do i=1,nband
         do j=1,mfreq
             write(mytmp,'(i5,5f16.8)') i, rmesh(j), &
                                  real(hybf(j,i,i,atom)), &
                                 aimag(hybf(j,i,i,atom)), &
                      real(hybf(j,i+nband,i+nband,atom)), &
                     aimag(hybf(j,i+nband,i+nband,atom))
         enddo ! over j={1,mfreq} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,nband} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_hybf

!=========================================================================
!>>>                                                                   <<<
!=========================================================================

!>>> write out initial hybridization function in matsubara frequency space

  subroutine ctqmc_dump_hybf_init(rmesh,hybf,atom)

     use constants
     use control
     use plo_quantities, only : natom

     implicit none

! external arguments
! matsubara frequency mesh
     real(dp), intent(in) :: rmesh(mfreq)

     integer, intent(in) :: atom
     
     character(20) :: filename

! hybridization function
     complex(dp), intent(in) :: hybf(mfreq,norbs,norbs,natom)

! local variables
! loop index
     integer :: i
     integer :: j

        write(filename,*) atom
        filename = 'solver.hyb_ini_'//trim(adjustl(filename))//'.dat'
        open(mytmp, file=filename, form='formatted', status='unknown')

! write it
     do i=1,nband
         do j=1,mfreq
             write(mytmp,'(i5,5f16.8)') i, rmesh(j), &
                                  real(hybf(j,i,i,atom)), &
                                 aimag(hybf(j,i,i,atom)), &
                      real(hybf(j,i+nband,i+nband,atom)), &
                     aimag(hybf(j,i+nband,i+nband,atom))
         enddo ! over j={1,mfreq} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,nband} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_hybf_init

!=========================================================================
!>>>                                                                   <<<
!=========================================================================


!>>> write out fitted hybridization function in matsubara frequency space

!  subroutine lanc_dump_hyb_fit(rmesh,hybf,atom)
!  end subroutine lanc_dump_hyb_fit

!=========================================================================
!>>>                                                                   <<<
!=========================================================================

!>>> write out self-energy function in matsubara frequency space
  subroutine ctqmc_dump_sigf(rmesh, sigf,atom)
     use constants
     use control
     use plo_quantities

     implicit none

! external arguments
! matsubara frequency mesh
     real(dp), intent(in) :: rmesh(mfreq)
     integer, intent(in) :: atom
     
     character(20) :: filename

! self-energy function
     complex(dp), intent(in) :: sigf(mfreq,norbs,norbs,natom)

! local variables
! loop index
     integer :: i
     integer :: j
     integer :: is
     integer :: iom

        write(filename,*) atom
        filename = 'solver.sgm_'//trim(adjustl(filename))//'.dat'
        open(mytmp, file=filename, form='formatted', status='unknown')

! write it
     do i=1,nband
         do j=1,mfreq
             write(mytmp,'(i5,5f16.8)') i, rmesh(j), &
                                  real(sigf(j,i,i,atom)), &
                                 aimag(sigf(j,i,i,atom)), &
                      real(sigf(j,i+nband,i+nband,atom)), &
                     aimag(sigf(j,i+nband,i+nband,atom))
         enddo ! over j={1,mfreq} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,nband} loop

! close data file
     close(mytmp)

! open data file: solver.sgm.in, for reusing it later
!        write(filename,*) atom
!        filename = 'solver.sgm_'//trim(adjustl(filename))
!        write(*,*) filename
!        open(mytmp, file=filename, form='formatted', status='unknown')

        write(filename,*) atom
        filename = 'solver.sgm.in_'//trim(adjustl(filename))//'.dat'
        open(mytmp, file=filename, form='formatted', status='unknown')

! write the full Sigma

    call plo_rewrite_array(sigf,sigma,'solver2plo')

     do is = 1, nspin
      do i = 1, nband
        do j = i, nband
          write(mytmp,*) ( sigma(i,j,iom,is,atom), iom=1,mfreq )
        end do
      end do
    end do

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_sigf

!=========================================================================
!>>>                                                                   <<<
!=========================================================================

!  subroutine lanc_dump_sigf_real(rmesh, sigf,atom)
!  end subroutine lanc_dump_sigf_real

!=========================================================================
!>>>                                                                   <<<
!=========================================================================


!>>> write out the Monte Carlo sampling histogram for perturbation expansion series

  subroutine ctqmc_dump_hist(hist,atom)

     use constants
     use control
     use plo_quantities, only : natom

     implicit none

! external arguments
! histogram data
     integer(int64), intent(in) :: hist(mkink,natom)

     integer, intent(in) :: atom

! local variables
! loop index
     integer  :: i

     character(20) :: filename

! dummy variables
     real(dp) :: raux

! scaled histogram data
     real(dp) :: haux(mkink)

! evaluate haux at first
     raux = real( sum(hist(:,atom)) )
     do i=1,mkink
         haux(i) = real( hist(i,atom) ) / raux
     enddo ! over i={1,mkink} loop

    write(filename,*) atom
    filename = 'solver.hist_'//trim(adjustl(filename))//'.dat'
    open(mytmp, file=filename, form='formatted', status='unknown')

! write it
     do i=1,mkink
         write(mytmp,'(i5,i10,f10.5)') i, hist(i,atom), haux(i)
     enddo ! over i={1,mkink} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_hist

!>>> write out the occupation matrix and double occupation matrix
  subroutine ctqmc_dump_nmat(nmat, nnmat, atom)
     use constants
     use control
     use plo_quantities, only : natom

     implicit none

! external arguments
! occupation matrix data
     real(dp), intent(in) :: nmat(norbs,natom)

     integer, intent(in) :: atom

! double occupation matrix data
     real(dp), intent(in) :: nnmat(norbs,norbs,natom)

! local variables
! loop index
     integer :: i
     integer :: j

     character(20) :: filename

    write(filename,*) atom
    filename = 'solver.nmat_'//trim(adjustl(filename))//'.dat'
    open(mytmp, file=filename, form='formatted', status='unknown')

! write it
     write(mytmp,'(a)') '  < n_i >   data:'
     do i=1,norbs
         write(mytmp,'(i5,f12.6)') i, nmat(i,atom)
     enddo ! over i={1,norbs} loop
     write(mytmp,'(a5,f12.6)') 'sup', sum( nmat(1:nband,atom) )
     write(mytmp,'(a5,f12.6)') 'sdn', sum( nmat(nband+1:norbs,atom) )
     write(mytmp,'(a5,f12.6)') 'sum', sum( nmat(1:norbs,atom) )

     write(mytmp,'(a)') '< n_i n_j > data:'
     do i=1,norbs
         do j=1,norbs
             write(mytmp,'(2i5,f12.6)') i, j, nnmat(i,j,atom)
         enddo ! over j={1,norbs} loop
     enddo ! over i={1,norbs} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_nmat

!>>> write out the probability of eigenstates of local hamiltonian matrix
  subroutine ctqmc_dump_prob(prob,atom)
     use constants
     use control
     use plo_quantities, only : natom

     implicit none

! external arguments
! probability data of eigenstates
     real(dp), intent(in) :: prob(ncfgs,natom)

     integer, intent(in) :: atom

! local variables
! loop index
     integer  :: i
     integer  :: j

! occupation number of eigenstates
     integer  :: noccs(ncfgs)

! net spin of eigenstates
     integer  :: soccs(ncfgs)

! atomic basis sets
     integer  :: basis(ncfgs,norbs)

! probability of occupation number distribution
     real(dp) :: oprob(0:norbs)

! probability of net spin distribution
     real(dp) :: sprob(-nband:nband)
     
     character(20) :: filename

! build atomic basis set, we do not order them according to their
! occupation numbers
     do i=1,ncfgs
         do j=1,norbs
             if ( btest(i-1,j-1) .eqv. .true. ) then
                 basis(i,j) = 1
             else
                 basis(i,j) = 0
             endif
         enddo ! over j={1,norbs} loop
     enddo ! over i={1,ncfgs} loop

! build occupation numbers for atomic basis set
     do i=1,ncfgs
         noccs(i) = sum( basis(i,:) )
     enddo ! over i={1,ncfgs} loop

! build net spin for eigenstates
     do i=1,ncfgs
         soccs(i) = ( sum( basis(i,1:nband) ) - sum( basis(i,nband+1:norbs) ) )
     enddo ! over i={1,ncfgs} loop

! evaluate oprob
     oprob = zero
     do i=1,ncfgs
         j = noccs(i)
         oprob(j) = oprob(j) + prob(i,atom)
     enddo ! over i={1,ncfgs} loop

! evaluate sprob
     sprob = zero
     do i=1,ncfgs
         j = soccs(i)
         sprob(j) = sprob(j) + prob(i,atom)
     enddo ! over i={1,ncfgs} loop

        write(filename,*) atom
        filename = 'solver.prob_'//trim(adjustl(filename))//'.dat'
        open(mytmp, file=filename, form='formatted', status='unknown')

! write it
     write(mytmp,'(a)') '# state probability: index | prob | occupy | spin'
     do i=1,ncfgs
         write(mytmp,'(i5,3f12.6)') i, prob(i,atom), real(noccs(i)), real(soccs(i)) * half
     enddo ! over i={1,ncfgs} loop

     write(mytmp,'(a)') '# orbital probability: index | occupy | prob'
     do i=0,norbs
         write(mytmp,'(i5,2f12.6)') i+1, real(i), oprob(i)
     enddo ! over i={0,norbs} loop
     write(mytmp,'(a5,12X,f12.6)') 'sum', sum(oprob)

     write(mytmp,'(a)') '# spin probability: index | spin | prob'
     do i=-nband,nband
         write(mytmp,'(i5,2f12.6)') i+nband+1, i*half, sprob(i)
     enddo ! over i={-nband,nband} loop
     write(mytmp,'(a5,12X,f12.6)') 'sum', sum(sprob)

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_prob
