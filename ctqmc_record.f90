!
! MK: modified, at least to account for multiple atoms, potentially also in other ways.
! Original header:
!
!-------------------------------------------------------------------------
! project : azalea
! program : ctqmc_record_gtau
!           ctqmc_record_grnf
!           ctqmc_record_hist
!           ctqmc_record_nmat
!           ctqmc_record_prob <<<---
!           ctqmc_reduce_gtau
!           ctqmc_reduce_grnf
!           ctqmc_reduce_hist
!           ctqmc_reduce_nmat
!           ctqmc_reduce_prob <<<---
!           ctqmc_symm_nmat
!           ctqmc_symm_gtau
!           ctqmc_symm_grnf
!           ctqmc_smth_sigf   <<<---
!           ctqmc_build_selfenergy1
!           ctqmc_build_selfenergy2
! source  : ctqmc_record.f90
! type    : subroutine
! author  : li huang (email:huangli712@yahoo.com.cn)
! history : 09/16/2009 by li huang
!           09/18/2009 by li huang
!           09/20/2009 by li huang
!           09/25/2009 by li huang
!           09/27/2009 by li huang
!           10/29/2009 by li huang
!           11/01/2009 by li huang
!           11/03/2009 by li huang
!           11/10/2009 by li huang
!           11/19/2009 by li huang
!           11/30/2009 by li huang
!           12/06/2009 by li huang
!           12/09/2009 by li huang
!           12/18/2009 by li huang
!           12/22/2009 by li huang
!           12/26/2009 by li huang
!           12/29/2009 by li huang
!           01/14/2010 by li huang
!           02/01/2010 by li huang
!           02/24/2010 by li huang
!           02/27/2010 by li huang
!           09/29/2010 by li huang
! purpose : measure, record, and postprocess the key observables produced
!           by the hybridization expansion version continuous time quantum
!           Monte Carlo (CTQMC) quantum impurity solver
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

!>>> record the impurity green's function in imaginary time axis
  subroutine ctqmc_record_gtau(atom)
     use constants
     use control
     use context

     implicit none

! local variables
! loop indices for start and end points
     integer  :: is
     integer  :: ie

! loop index for flavor channel
     integer  :: flvr

     integer,intent(in)  :: atom

! index for imaginary time \tau
     integer  :: curr

! used to store the element of mmat matrix
     real(dp) :: maux

! imaginary time for start and end points
     real(dp) :: taus
     real(dp) :: taue

! length betweem taus and taue
     real(dp) :: dtau

! interval for imaginary time slice
     real(dp) :: step

! evaluate step at first
     step = real(ntime - 1) / beta

     CTQMC_FLAVOR_LOOP: do flvr=1,norbs

! get imaginary time value for segments
         do is=1,rank(flvr,atom)
             taus = time_s( index_s(is, flvr), flvr )

             do ie=1,rank(flvr,atom)
                 taue = time_e( index_e(ie, flvr), flvr )

! evaluate dtau
                 dtau = taue - taus

! get matrix element from mmat, pay special attention to the sign of dtau
                 maux = mmat(ie, is, flvr,atom) * sign(one, dtau)

! adjust dtau, keep it stay in (zero, beta)
                 if ( dtau < zero ) then
                     dtau = dtau + beta
                 endif

! determine index for imaginary time
                 curr = nint( dtau * step ) + 1

! special tricks for the first point and the last point
                 if ( curr == 1 .or. curr == ntime ) then
                     maux = two * maux
                 endif

! record gtau, we normalize gtau in ctqmc_dump_gtau() subroutine
                 gtau(curr, flvr, flvr,atom) = gtau(curr, flvr, flvr,atom) - maux

             enddo ! over ie={1,rank(flvr,atom)} loop
         enddo ! over is={1,rank(flvr,atom)} loop

     enddo CTQMC_FLAVOR_LOOP ! over flvr={1,norbs} loop

     return
  end subroutine ctqmc_record_gtau

!>>> record the impurity green's function in matsubara frequency space
  subroutine ctqmc_record_grnf(atom)
     use constants
     use control
     use context

     implicit none

! local variables
! loop index over matsubara frequencies
     integer :: ifrq

! loop index for flavor channel
     integer :: flvr

     integer,intent(in) :: atom

! note: only the first nfreq points of grnf are modified
     do flvr=1,norbs
         do ifrq=1,nfreq
             grnf(ifrq, flvr, flvr,atom) = grnf(ifrq, flvr, flvr,atom) + gmat(ifrq, flvr, flvr,atom)
         enddo ! over ifrq={1,nfreq} loop
     enddo ! over flvr={1,norbs} loop

     return
  end subroutine ctqmc_record_grnf

!>>> record the histogram of perturbation expansion series
  subroutine ctqmc_record_hist(atom)
     use context

     implicit none
     
     integer,intent(in) :: atom

! note: if ckink == 0, we record its count in hist(mkink)
     if ( ckink > 0 ) then
         hist(ckink,atom) = hist(ckink,atom) + 1
     else
         hist(mkink,atom) = hist(mkink,atom) + 1
     endif

     return
  end subroutine ctqmc_record_hist

!>>> record the occupation matrix, double occupation matrix, and auxiliary
! physical observables simulataneously
  subroutine ctqmc_record_nmat(atom)
     use constants
     use control
     use context

     implicit none

! local variables
! loop index over segments
     integer  :: i
     integer, intent(in) :: atom

! loop index for flavor channel
     integer  :: flvr

! imaginary time for start and end points
     real(dp) :: ts
     real(dp) :: te

! total length of segments
     real(dp) :: sgmt(norbs)

! used to record overlaps between two segments
     real(dp) :: oaux(norbs)
     real(dp) :: ovlp(norbs,norbs)

! evaluate occupation matrix: < n_i >
!-------------------------------------------------------------------------
     do flvr=1,norbs

! case 1: null occupation
         if ( stts(flvr,atom) == 0 ) then
             sgmt(flvr) = zero

! case 2: partial occupation, segment scheme
         else if ( stts(flvr,atom) == 1 ) then
             sgmt(flvr) = zero
             do i=1,rank(flvr,atom)
                 ts = time_s(index_s(i, flvr), flvr)
                 te = time_e(index_e(i, flvr), flvr)
                 sgmt(flvr) = sgmt(flvr) + abs( te - ts )
             enddo ! over i={1,rank(flvr,atom)} loop

! case 3: partial occupation, anti-segment scheme
         else if ( stts(flvr,atom) == 2 ) then
             sgmt(flvr) = beta
             do i=1,rank(flvr,atom)
                 ts = time_s(index_s(i, flvr), flvr)
                 te = time_e(index_e(i, flvr), flvr)
                 sgmt(flvr) = sgmt(flvr) - abs( ts - te )
             enddo ! over i={1,rank(flvr,atom)} loop

! case 4: full occupation
         else if ( stts(flvr,atom) == 3 ) then
             sgmt(flvr) = beta

         endif ! back if ( stts(flvr,atom) == 0 ) block

         nmat(flvr,atom) = nmat(flvr,atom) + sgmt(flvr) / beta
     enddo ! over flvr={1,norbs} loop
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! evaluate double occupation matrix: < n_i n_j >
!-------------------------------------------------------------------------
     do flvr=1,norbs

! case 1: null occupation
         if ( stts(flvr,atom) == 0 ) then
             ovlp(flvr,:) = zero

! case 2: partial occupation, segment scheme
         else if ( stts(flvr,atom) == 1 ) then
             ovlp(flvr,:) = zero
             do i=1,rank(flvr,atom)
                 ts = time_s(index_s(i, flvr), flvr)
                 te = time_e(index_e(i, flvr), flvr)
                 call ctqmc_make_overlap(flvr, ts, te, oaux,atom)
                 ovlp(flvr,:) = ovlp(flvr,:) + oaux
             enddo ! over i={1,rank(flvr,atom)} loop

! case 3: partial occupation, anti-segment scheme
! pay special attention to the head and tail parts
         else if ( stts(flvr,atom) == 2 ) then
             ovlp(flvr,:) = zero
             do i=1,rank(flvr,atom)-1
                 ts = time_s(index_s(i,   flvr), flvr)
                 te = time_e(index_e(i+1, flvr), flvr)
                 call ctqmc_make_overlap(flvr, ts, te, oaux,atom)
                 ovlp(flvr,:) = ovlp(flvr,:) + oaux
             enddo ! over i={1,rank(flvr,atom)-1} loop

             te = time_e(index_e(1, flvr), flvr)
             call ctqmc_make_overlap(flvr, zero, te, oaux,atom)
             ovlp(flvr,:) = ovlp(flvr,:) + oaux

             ts = time_s(index_s(rank(flvr,atom), flvr), flvr)
             call ctqmc_make_overlap(flvr, ts, beta, oaux,atom)
             ovlp(flvr,:) = ovlp(flvr,:) + oaux

! case 4: full occupation
         else if ( stts(flvr,atom) == 3 ) then
             call ctqmc_make_overlap(flvr, zero, beta, oaux,atom)
             ovlp(flvr,:) = oaux

         endif ! back if ( stts(flvr,atom) == 0 ) block

         nnmat(flvr,:,atom) = nnmat(flvr,:,atom) + ovlp(flvr,:) / beta
     enddo ! over flvr={1,norbs} loop
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! evaluate spin magnetization: < Sz >
!-------------------------------------------------------------------------
     do flvr=1,nband
         paux(4,atom) = paux(4,atom) + ( sgmt(flvr) - sgmt(flvr+nband) ) / beta
     enddo ! over flvr={1,nband} loop
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! evaluate kinetic energy: ekin
! equation : -T < k >
!-------------------------------------------------------------------------
     paux(3,atom) = paux(3,atom) - real(ckink * norbs) / beta
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! evaluate potential energy: epot
!-------------------------------------------------------------------------
     do flvr=1,norbs
         do i=1,flvr
             paux(2,atom) = paux(2,atom) + uumat(flvr,i,atom) * ovlp(flvr,i) / beta
         enddo ! over i={1,flvr} loop
     enddo ! over flvr={1,norbs} loop
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! evaluate total energy: etot
!-------------------------------------------------------------------------
     paux(1,atom) = paux(2,atom) + paux(3,atom)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

     return
  end subroutine ctqmc_record_nmat

!>>> record the probability of atomic states
  subroutine ctqmc_record_prob(atom)
     use constants
     use control
     use context

     implicit none

! local variables
! current flavor channel
     integer :: flvr
     integer, intent(in   ) :: atom

! atomic state index
     integer :: pstat

! current atomic state for segment representation
     integer :: state(norbs)

! generate current atomic state
     do flvr=1,norbs
         select case ( stts(flvr,atom) )

             case (0:1)
                 state(flvr) = 0

             case (2:3)
                 state(flvr) = 1

         end select
     enddo ! over flvr={1,norbs} loop

! convert atomic state array to index
     call ctqmc_make_state(norbs, pstat, state)

! accumulate the data
     prob(pstat,atom) = prob(pstat,atom) + one

     return
  end subroutine ctqmc_record_prob

!>>> reduce the gtau from all children processes
  subroutine ctqmc_reduce_gtau(gtau_mpi)
     use constants
     use context
     use plo_quantities, only : natom
     use mmpi

     implicit none

! external arguments
! impurity green's function
     real(dp), intent(out) :: gtau_mpi(ntime,norbs,norbs,natom)

! initialize gtau_mpi
     gtau_mpi = zero

! build gtau_mpi, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(gtau, gtau_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     gtau_mpi = gtau

# endif /* MPI */

! calculate the average
     gtau_mpi = gtau_mpi / real(nprocs)

     return
  end subroutine ctqmc_reduce_gtau

!>>> reduce the grnf from all children processes
  subroutine ctqmc_reduce_grnf(grnf_mpi)
     use constants
     use context
     use plo_quantities, only : natom
     use mmpi

     implicit none

! external arguments
! impurity green's function
     complex(dp), intent(out) :: grnf_mpi(mfreq,norbs,norbs,natom)

! initialize grnf_mpi
     grnf_mpi = zero

! build grnf_mpi, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(grnf, grnf_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     grnf_mpi = grnf

# endif /* MPI */

! calculate the average
     grnf_mpi = grnf_mpi / real(nprocs)

     return
  end subroutine ctqmc_reduce_grnf

!>>> reduce the hist from all children processes
  subroutine ctqmc_reduce_hist(hist_mpi)
     use constants
     use context
     use plo_quantities, only : natom
     use mmpi

     implicit none

! external arguments
! histogram for perturbation expansion series
     integer(int64), intent(out) :: hist_mpi(mkink,natom)

!     integer, intent(in) :: atom
! initialize hist_mpi
     hist_mpi = 0

! build hist_mpi, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(hist, hist_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     hist_mpi = hist

# endif /* MPI */

! calculate the average
     hist_mpi = hist_mpi / nprocs

     return
  end subroutine ctqmc_reduce_hist

!>>> reduce the nmat and nnmat from all children processes
  subroutine ctqmc_reduce_nmat(nmat_mpi, nnmat_mpi)
     use constants
     use context
     use plo_quantities,only : natom
     use mmpi

     implicit none

! external arguments

!     integer, intent(in) :: atom

! occupation number matrix
     real(dp), intent(out) :: nmat_mpi(norbs,natom)

! double occupation number matrix
     real(dp), intent(out) :: nnmat_mpi(norbs,norbs,natom)

! initialize nmat_mpi and nnmat_mpi
     nmat_mpi = zero
     nnmat_mpi = zero

! build nmat_mpi and nnmat_mpi, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(nmat, nmat_mpi)
     call mp_allreduce(nnmat, nnmat_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     nmat_mpi = nmat
     nnmat_mpi = nnmat

# endif /* MPI */

! calculate the average
     nmat_mpi = nmat_mpi / real(nprocs)
     nnmat_mpi = nnmat_mpi / real(nprocs)

     return
  end subroutine ctqmc_reduce_nmat

!>>> reduce the prob from all children processes
  subroutine ctqmc_reduce_prob(prob_mpi)
     use constants
     use context
     use plo_quantities, only : natom
     use mmpi

     implicit none

! external arguments

!     integer, intent(in) :: atom

! probability of atomic states
     real(dp), intent(out) :: prob_mpi(ncfgs,natom)

! initialize prob_mpi
     prob_mpi = zero

! build prob_mpi, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(prob, prob_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     prob_mpi = prob

# endif /* MPI */

! calculate the average
     prob_mpi = prob_mpi / real(nprocs)

     return
  end subroutine ctqmc_reduce_prob

!>>> symmetrize the nmat according to symm vector
  subroutine ctqmc_symm_nmat(symm, nmat, atom)
     use constants
     use control
     use plo_quantities, only : natom
     
     implicit none

! external arguments
! symmetry vector
     integer, intent(in) :: symm(norbs,natom)
     
     integer, intent(in) :: atom

! occupation number
     real(dp), intent(inout) :: nmat(norbs,natom)

! local variables
! loop index over bands
     integer  :: ibnd
     integer  :: jbnd

! dummy variables
     real(dp) :: raux

! histogram vector
! note: it is NOT the global one
     integer  :: hist(norbs)

! build histogram
     hist = 0
     do ibnd=1,norbs
         hist(symm(ibnd,atom)) = hist(symm(ibnd,atom)) + 1
     enddo ! over ibnd={1,norbs} loop

! perform symmetrization for those orbitals which symm index are identity
     if ( issun == 2 ) then
         do ibnd=1,norbs
             if ( hist(ibnd) > 0 ) then         ! need to enforce symmetry
                 raux = zero

                 do jbnd=1,norbs                ! gather the data
                     if ( symm(jbnd,atom) == ibnd ) then
                         raux = raux + nmat(jbnd,atom)
                     endif
                 enddo ! over jbnd={1,norbs} loop

                 raux = raux / real(hist(ibnd)) ! calculate average value

                 do jbnd=1,norbs                ! setup it
                     if ( symm(jbnd,atom) == ibnd ) then
                         nmat(jbnd,atom) = raux
                     endif
                 enddo ! over jbnd={1,norbs} loop
             endif
         enddo ! over ibnd={1,norbs} loop
     endif ! back if ( issun == 2 ) block

! symmetrize nmat over spin
     if ( isspn == 1 ) then
         do jbnd=1,nband
             raux = ( nmat(jbnd,atom) + nmat(jbnd+nband,atom) ) / two
             nmat(jbnd,atom) = raux
             nmat(jbnd+nband,atom) = raux
         enddo ! over jbnd={1,nband} loop
     endif ! back if ( isspn == 1 ) block

     return
  end subroutine ctqmc_symm_nmat

!>>> symmetrize the gtau according to symm vector
! only the diagonal elements are taken into considerations
  subroutine ctqmc_symm_gtau(symm, gtau, atom)
     use constants
     use control
     use plo_quantities,only : natom

     implicit none

! external arguments
! symmetry vector
     integer, intent(in) :: symm(norbs,natom)

     integer, intent(in) :: atom

! impurity green's function
     real(dp), intent(inout) :: gtau(ntime,norbs,norbs,natom)

! local variables
! loop index over bands
     integer  :: ibnd
     integer  :: jbnd

! loop index over imaginary-time points
     integer  :: ktau

! dummy variables
     real(dp) :: raux

! histogram vector
! note: it is NOT the global one
     integer  :: hist(norbs)

! build histogram
     hist = 0
     do ibnd=1,norbs
         hist(symm(ibnd,atom)) = hist(symm(ibnd,atom)) + 1
     enddo ! over ibnd={1,norbs} loop

! perform symmetrization for those orbitals which symm index are identity
     if ( issun == 2 ) then
         do ktau=1,ntime
             do ibnd=1,norbs
                 if ( hist(ibnd) > 0 ) then         ! need to enforce symmetry
                     raux = zero

                     do jbnd=1,norbs                ! gather the data
                         if ( symm(jbnd,atom) == ibnd ) then
                             raux = raux + gtau(ktau,jbnd,jbnd,atom)
                         endif
                     enddo ! over jbnd={1,norbs} loop

                     raux = raux / real(hist(ibnd)) ! calculate average value

                     do jbnd=1,norbs                ! setup it
                         if ( symm(jbnd,atom) == ibnd ) then
                             gtau(ktau,jbnd,jbnd,atom) = raux
                         endif
                     enddo ! over jbnd={1,norbs} loop
                 endif
             enddo ! over ibnd={1,norbs} loop
         enddo ! over ktau={1,ntime} loop
     endif ! back if ( issun == 2 ) block

! symmetrize gtau over spin
     if ( isspn == 1 ) then
         do ktau=1,ntime
             do jbnd=1,nband
                 raux = ( gtau(ktau,jbnd,jbnd,atom) + gtau(ktau,jbnd+nband,jbnd+nband,atom) ) / two
                 gtau(ktau,jbnd,jbnd,atom) = raux
                 gtau(ktau,jbnd+nband,jbnd+nband,atom) = raux
             enddo ! over jbnd={1,nband} loop
         enddo ! over ktau={1,ntime} loop
     endif ! back if ( isspn == 1 ) block

     return
  end subroutine ctqmc_symm_gtau

!>>> symmetrize the grnf according to symm vector
! only the diagonal elements are taken into considerations
  subroutine ctqmc_symm_grnf(symm, grnf, atom)
     use constants
     use control
     use plo_quantities, only : natom

     implicit none

! external arguments
! symmetry vector
     integer, intent(in) :: symm(norbs,natom)

     integer, intent(in) :: atom

! impurity green's function
     complex(dp), intent(inout) :: grnf(mfreq,norbs,norbs,natom)

! local variables
! loop index over bands
     integer :: ibnd
     integer :: jbnd

! loop index over matsubara frequencies
     integer :: kfrq

! dummy variables
     complex(dp) :: caux

! histogram vector
! note: it is NOT the global one
     integer :: hist(norbs)

! build histogram
     hist = 0
     do ibnd=1,norbs
         hist(symm(ibnd,atom)) = hist(symm(ibnd,atom)) + 1
     enddo ! over ibnd={1,norbs} loop

! perform symmetrization for those orbitals which symm index are identity
     if ( issun == 2 ) then
         do kfrq=1,mfreq
             do ibnd=1,norbs
                 if ( hist(ibnd) > 0 ) then         ! need to enforce symmetry
                     caux = czero

                     do jbnd=1,norbs                ! gather the data
                         if ( symm(jbnd,atom) == ibnd ) then
                             caux = caux + grnf(kfrq,jbnd,jbnd,atom)
                         endif
                     enddo ! over jbnd={1,norbs} loop

                     caux = caux / real(hist(ibnd)) ! calculate average value

                     do jbnd=1,norbs                ! setup it
                         if ( symm(jbnd,atom) == ibnd ) then
                             grnf(kfrq,jbnd,jbnd,atom) = caux
                         endif
                     enddo ! over jbnd={1,norbs} loop
                 endif
             enddo ! over ibnd={1,norbs} loop
         enddo ! over kfrq={1,mfreq} loop
     endif ! back if ( issun == 2 ) block

! symmetrize grnf over spin
     if ( isspn == 1 ) then
         do kfrq=1,mfreq
             do jbnd=1,nband
                 caux = ( grnf(kfrq,jbnd,jbnd,atom) + grnf(kfrq,jbnd+nband,jbnd+nband,atom) ) / two
                 grnf(kfrq,jbnd,jbnd,atom) = caux
                 grnf(kfrq,jbnd+nband,jbnd+nband,atom) = caux
             enddo ! over jbnd={1,nband} loop
         enddo ! over kfrq={1,mfreq} loop
     endif ! back if ( isspn == 1 ) block

     return
  end subroutine ctqmc_symm_grnf

!>>> smooth impurity self-energy function in low frequency region
  subroutine ctqmc_smth_sigf(sigf)
     use constants
     use control

     implicit none

! external arguments
! impurity self-energy function to be smoothen
     complex(dp), intent(inout) :: sigf(nfreq)

! local variables
! loop index
     integer  :: i
     integer  :: j
     integer  :: k

! smooth radius
     integer  :: lrad
     integer  :: srad

! imaginary part of self-energy function
     real(dp) :: ti

! real part of self-energy function
     real(dp) :: tr

! dummy variables for addition
     complex(dp) :: saux

! dummy self-energy function
     complex(dp) :: stmp(nfreq)

! determine smooth radius
     lrad = nfreq / 4  ! large radius
     srad = nfreq / 16 ! small radius

! |---------|---------|----------------------|
! 1         lrad      2*lrad                 nfreq
! deal with [1,lrad], head part
     do k=1,lrad
         stmp(k) = sigf(k)
     enddo ! over k={1,lrad} loop

! deal with [lrad+1,2*lrad], intermediate part
     do i=1,lrad
         k = lrad + i
         saux = czero
         do j=-srad,srad
             saux = saux + sigf(k+j)
         enddo ! over j={-srad,srad} loop
         stmp(k) = saux / real(2 * srad + 1)
         stmp(k) = ( (lrad - i) * sigf(k) + i * stmp(k) ) / real(lrad)
     enddo ! over i={1,lrad} loop

! deal with [nfreq-2*lrad+1,nfreq], tail part
     do k=nfreq-2*lrad+1,nfreq
         tr =  real( stmp(nfreq-2*lrad) )
         ti = aimag( stmp(nfreq-2*lrad) ) * real(nfreq - 2 * lrad) / real(k)
         stmp(k) = dcmplx(tr,ti)
     enddo ! over k={nfreq-2*lrad+1,nfreq} loop

! copy stmp to sigf
     do k=1,nfreq
         sigf(k) = stmp(k)
     enddo ! over k={1,nfreq} loop

     return
  end subroutine ctqmc_smth_sigf

!>>> build self-energy function by high frequency expansions, and then
! supplement the high frequency part of impurity green's function
  subroutine ctqmc_build_selfenergy1(atom)
     use constants
     use control
     use context
     use plo_quantities

     implicit none

     integer, intent(in   ) :: atom

! local variables
! loop index
     integer  :: i
     integer  :: j
     integer  :: k

! 0-order and 1-order expansion coefficients for self-energy function
     real(dp) :: u0(norbs)
     real(dp) :: u1(norbs)

! jump for self-energy function in the boundary between quantum Monte Carlo
! sampling and high frequency approximation
     real(dp) :: sdel(norbs)

! dummy imurity green's function: G^{-1}
     complex(dp) :: gaux(norbs,norbs)

! initialize variables
     u0 = zero
     u1 = zero

! calculate 0-order coefficient for self-energy function
     do i=1,norbs
         do j=1,norbs
             if ( i /= j ) then
                 u0(i) = u0(i) + uumat(j,i,atom) * nmat(j,atom)
             endif
         enddo ! over j={1,norbs} loop
     enddo ! over i={1,norbs} loop

! calculate 1-order coefficient for self-energy function
     do i=1,norbs
         do j=1,norbs
             do k=1,norbs
                 if ( k /= i .and. j /= i ) then
                     u1(i) = u1(i) + uumat(j,i,atom) * uumat(k,i,atom) * ( nnmat(j,k,atom) - nmat(j,atom) * nmat(k,atom) )
                 endif
             enddo ! over k={1,norbs} loop
         enddo ! over j={1,norbs} loop
     enddo ! over i={1,norbs} loop

! build self-energy function at low frequency region
!-------------------------------------------------------------------------
! filter grnf to suppress the fluctuation of its real part
!-------------------------------------------------------------------------
!<     do k=1,nfreq
!<         do i=1,norbs
!<             grnf(k,i,i) = dcmplx( zero, aimag( grnf(k,i,i) ) )
!<         enddo ! over i={1,norbs} loop
!<     enddo ! over k={1,nfreq} loop
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
     do k=1,nfreq
         gaux = grnf(k,:,:,atom)
         call ctqmc_zmat_inv(norbs, gaux)
         do i=1,norbs
             sig2(k,i,i,atom) = cmesh(k) - eimp(i,atom) - gaux(i,i) - hybf(k,i,i,atom) + dcmplx(dcc(1,atom),0.00_dp)
         enddo ! over i={1,norbs} loop
     enddo ! over k={1,nfreq} loop

!-------------------------------------------------------------------------
! filter sig2 to suppress the fluctuation of its imaginary part
!-------------------------------------------------------------------------
     do k=1,16
         do i=1,norbs
             call ctqmc_smth_sigf( sig2(1:nfreq,i,i,atom) ) ! smooth it 16 times
         enddo ! over i={1,norbs} loop
     enddo ! over k={1,16} loop
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! evaluate the high frequency part of self-energy function
     do i=1,norbs
         sdel(i) = -u1(i) / rmesh(nfreq) - aimag( sig2(nfreq,i,i,atom) )
     enddo ! over i={1,norbs} loop

     do k=nfreq+1,mfreq
         do i=1,norbs
             sig2(k,i,i,atom) = u0(i) - czi * u1(i) / rmesh(k) - dcmplx(zero, sdel(i)) + dcmplx(dcc(1,atom),0.00_dp)
         enddo ! over i={1,norbs} loop
     enddo ! over k={nfreq+1,mfreq} loop

! supply the high frequency part of impurity green's function
     do k=nfreq+1,mfreq
         gaux = czero
         do i=1,norbs
             gaux(i,i) = cmesh(k) - eimp(i,atom) - sig2(k,i,i,atom) - hybf(k,i,i,atom)
         enddo ! over i={1,norbs} loop
         call ctqmc_zmat_inv(norbs, gaux)
         grnf(k,:,:,atom) = gaux
     enddo ! over k={nfreq+1,mfreq} loop

     return
  end subroutine ctqmc_build_selfenergy1

!>>> build atomic green's function and self-energy function using improved
! Hubbard-I approximation, and then make interpolation for self-energy
! function between low frequency QMC data and high frequency Hubbard-I
! approximation data, the full impurity green's function can be obtained by
! using dyson's equation finally
  subroutine ctqmc_build_selfenergy2(atom)
     use constants
     use control
     use context
     use plo_quantities

     implicit none

!external
    integer, intent(in) :: atom

! local parameters
! maximum allowable number of non-zero elements in F matrix
     integer, parameter :: nzero = 1024

! local variables
! loop index
     integer  :: i
     integer  :: j
     integer  :: k
     integer  :: m
     integer  :: n

! dummy integer variables, used to build F matrix
     integer  :: start
     integer  :: value
     integer  :: permute

! dummy real variables, used to interpolate self-energy function
     real(dp) :: ob, oe
     real(dp) :: d0, d1

! dummy complex variables, used to interpolate self-energy function
     complex(dp) :: cb, ce
     complex(dp) :: sinf
     complex(dp) :: caux

! dummy atomic states: alpha, beta, gamma
     integer  :: sa(norbs)
     integer  :: sb(norbs)
     integer  :: sc(norbs)

! atomic basis sets
     integer  :: basis(ncfgs,norbs)

! F matrix, < alpha | f_{n} | beta >
     integer  :: fcounter(norbs)
     integer  :: fa(nzero,norbs)
     integer  :: fb(nzero,norbs)
     integer  :: fv(nzero,norbs)

! eigenvalues for local hmailtonian
     real(dp) :: eaux(ncfgs)

! dummy impurity green's function: G^{-1}
     complex(dp) :: gaux(norbs,norbs)

! atomic green's function and self-energy function in Hubbard-I approximation
     complex(dp) :: ghub(mfreq,norbs)
     complex(dp) :: shub(mfreq,norbs)

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

! evaluate atomic eigenvalues directly
     eaux = zero
     do i=1,ncfgs
         do j=1,norbs
             eaux(i) = eaux(i) + ( eimp(j,atom) ) * basis(i,j)
         enddo ! over j={1,norbs} loop
         do j=1,norbs-1
             do k=j+1,norbs
                 if ( basis(i,j) == 1 .and. basis(i,k) == 1 ) then
                     eaux(i) = eaux(i) + uumat(j,k,atom)
                 endif
             enddo ! over k={j+1,norbs} loop
         enddo ! over j={1,norbs-1} loop
     enddo ! over i={1,ncfgs} loop

! build F matrix < alpha | f_{n} | beta >
! note 1: to save the memory and accelerate the computation, we only store
! the non-zero element of F matrix
! note 2: it is crucial to check whether the number of non-zero elements
! exceed limit (nzero)
     fcounter = 0
     alpha_loop: do i=1,ncfgs
         sa = basis(i,:)
         beta_loop: do j=1,ncfgs
             sb = basis(j,:)

             orbital_loop: do m=1,norbs
                 sc = sb

                 if ( sc(m) == 1 ) then
                     permute = 1
                     do n=1,m-1
                         if ( sc(n) == 1 ) permute = -permute
                     enddo ! over n={1,m-1} loop
                     sc(m) = 0

                     value = 1
                     do n=1,norbs
                         if ( sa(n) /= sc(n) ) value = 0
                     enddo ! over n={1,norbs} loop
                     value = value * permute
                 else
                     value = 0
                 endif

                 if ( value /= 0 ) then
                     fcounter(m) = fcounter(m) + 1
                     if ( fcounter(m) > nzero ) then
                         call ctqmc_print_error('ctqmc_build_selfenergy2','non-zero elements exceed limit')
                     endif
                     fa(fcounter(m),m) = i
                     fb(fcounter(m),m) = j
                     fv(fcounter(m),m) = value
                 endif ! back if ( value /= 0 ) block
             enddo orbital_loop ! over m={1,norbs} loop

         enddo beta_loop ! over j={1,ncfgs} loop
     enddo alpha_loop ! over i={1,ncfgs} loop

! calculate atomic green's function using Hubbard-I approximation
     do i=1,norbs
         do k=1,mfreq
             caux = czero
             do m=1,fcounter(i)
                 ob = fv(m,i) * fv(m,i) * ( prob(fa(m,i),atom) + prob(fb(m,i),atom) )
                 cb = cmesh(k) + eaux(fa(m,i)) - eaux(fb(m,i))
                 caux = caux +  ob / cb
             enddo ! over m={1,fcounter(i)} loop
             ghub(k,i) = caux
         enddo ! over k={1,mfreq} loop
     enddo ! over i={1,norbs} loop

! calculate atomic self-energy function using dyson's equation
     do i=1,norbs
         do k=1,mfreq
             shub(k,i) = cmesh(k) - eimp(i,atom) - one / ghub(k,i)
         enddo ! over k={1,mfreq} loop
     enddo ! over i={1,norbs} loop

! build self-energy function at low frequency region
!-------------------------------------------------------------------------
! filter grnf to suppress the fluctuation of its real part
!-------------------------------------------------------------------------
!<     do k=1,nfreq
!<         do i=1,norbs
!<             ob =  real( grnf(k,i,i) )
!<             oe = aimag( grnf(k,i,i) )
!<             grnf(k,i,i) = dcmplx( zero, oe )
!<         enddo ! over i={1,norbs} loop
!<     enddo ! over k={1,nfreq} loop
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
     do k=1,nfreq
         gaux = grnf(k,:,:,atom)
         call ctqmc_zmat_inv(norbs, gaux)
         do i=1,norbs
             sig2(k,i,i,atom) = cmesh(k) - eimp(i,atom) - gaux(i,i) - hybf(k,i,i,atom)
         enddo ! over i={1,norbs} loop
     enddo ! over k={1,nfreq} loop
!-------------------------------------------------------------------------
! filter sig2 to suppress the fluctuation of its imaginary part
!-------------------------------------------------------------------------
     do k=1,16
         do i=1,norbs
             call ctqmc_smth_sigf( sig2(1:nfreq,i,i,atom) ) ! smooth it 16 times
         enddo ! over i={1,norbs} loop
     enddo ! over k={1,16} loop
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! interpolates self-energy function between low energy QMC data and high
! energy Hubbard-I approximation
     do i=1,norbs

! determine the base point, its value is calculated by using five points
         cb = czero
         do k=nfreq-4,nfreq
             cb = cb + sig2(k,i,i,atom)
         enddo ! over k={nfreq-4,nfreq} loop

         cb = cb / real(5)
         ob = rmesh(nfreq-2)

! for the imaginary part
! determine the intermediate region [nfreq+1,start] at first
         do k=nfreq+1,mfreq
             start = k
             d0 = aimag( shub(k,i) - cb ) / ( rmesh(k) - ob )
             d1 = aimag( shub(k,i) - shub(k-1,i) ) / ( rmesh(k) - rmesh(k-1) )
             if ( abs( d0 - d1 ) < 0.02_dp ) exit
         enddo ! over k={nfreq+1,mfreq} loop
         if ( start - nfreq < 32 ) start = nfreq + 32

         ce = shub(start,i)
         oe = rmesh(start)

! deal with the intermediate region, using linear interpolation
         do k=nfreq+1,start
             sig2(k,i,i,atom) = dcmplx( zero, aimag(cb) + aimag( ce - cb ) * ( rmesh(k) - ob ) / ( oe - ob ) )
         enddo ! over k={nfreq+1,start} loop

! deal with the tail region, using atomic self-energy function directly
         do k=start+1,mfreq
             sig2(k,i,i,atom) = dcmplx( zero, aimag( shub(k,i) ) )
         enddo ! over k={start+1,mfreq} loop

! for the real part
         sinf = shub(mfreq,i)
         do k=nfreq+1,mfreq
             sig2(k,i,i,atom) = sig2(k,i,i,atom) + real(sinf) + ( ob / rmesh(k) )**2 * real( cb - sinf )
         enddo ! over k={nfreq+1,mfreq} loop
!Doublecount the whole thing
         do k=1,mfreq
             sig2(k,i,i,atom) = sig2(k,i,i,atom) + dcmplx(dcc(1,atom),zero)
         enddo ! over k={nfreq+1,mfreq} loop

     enddo ! over i={1,norbs} loop

! calculate final impurity green's function using dyson's equation
     do k=1,mfreq
         gaux = czero
         do i=1,norbs
             gaux(i,i) = cmesh(k) - eimp(i,atom) - sig2(k,i,i,atom) - hybf(k,i,i,atom)
         enddo ! over i={1,norbs} loop
         call ctqmc_zmat_inv(norbs, gaux)
         grnf(k,:,:,atom) = gaux
     enddo ! over k={1,mfreq} loop

     return
  end subroutine ctqmc_build_selfenergy2
