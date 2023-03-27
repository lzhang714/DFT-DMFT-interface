!
! MK: modified, at least to account for multiple atoms, potentially also in other ways.
! Original header:
!
!-------------------------------------------------------------------------
! project : azalea
! program : ctqmc_insert_kink
!           ctqmc_remove_kink
!           ctqmc_lshift_kink
!           ctqmc_rshift_kink
!           ctqmc_reswap_kink
!           ctqmc_reflip_kink
!           ctqmc_reload_kink <<<---
!           cat_insert_matrix
!           cat_remove_matrix
!           cat_lshift_matrix
!           cat_rshift_matrix
!           cat_reswap_matrix
!           cat_reflip_matrix
!           cat_reload_matrix <<<---
!           cat_insert_detrat
!           cat_remove_detrat
!           cat_lshift_detrat
!           cat_rshift_detrat
!           cat_reswap_detrat
!           cat_reflip_detrat <<<---
! source  : ctqmc_update.f90
! type    : subroutines
! author  : li huang (email:huangli712@yahoo.com.cn)
! history : 09/16/2009 by li huang
!           09/18/2009 by li huang
!           09/20/2009 by li huang
!           09/24/2009 by li huang
!           09/26/2009 by li huang
!           09/30/2009 by li huang
!           10/02/2009 by li huang
!           10/25/2009 by li huang
!           10/29/2009 by li huang
!           11/02/2009 by li huang
!           11/08/2009 by li huang
!           11/17/2009 by li huang
!           11/20/2009 by li huang
!           11/24/2009 by li huang
!           11/27/2009 by li huang
!           11/30/2009 by li huang
!           12/09/2009 by li huang
!           12/18/2009 by li huang
!           12/26/2009 by li huang
!           01/05/2010 by li huang
!           02/27/2010 by li huang
!           03/22/2010 by li huang
!           06/09/2010 by li huang
! purpose : provide basic infrastructure (elementary updating subroutines)
!           for hybridization expansion version continuous time quantum
!           Monte Carlo (CTQMC) quantum impurity solver.
!           the following subroutines mainly deal with the \mathscr{M}
!           matrix: mmat, and \mathscr{G} matrix: gmat.
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
!>>> driver layer: updating perturbation expansion series              <<<
!-------------------------------------------------------------------------

!>>> insert new segment or anti-segment in the perturbation expansion series
  subroutine ctqmc_insert_kink(atom)
     use constants
     use control
     use context

     use spring

     implicit none

     integer, intent(in) :: atom

! local variables
! whether the new segment or anti-segment can be inserted diagrammatically
     logical  :: ladd

! whether it is an anti-segment
! anti = .true., anti-segment
! anti = .false., segment
     logical  :: anti

! whether the update operation is accepted
     logical  :: pass

! current flavor channel for both band and spin
     integer  :: flvr

! index address to insert new segment or anti-segment
! is and ie are for start and end points, respectively
     integer  :: is, ie

! transition probability for insert new segment or anti-segment
     real(dp) :: p

! \tau_s, start point of the new segment
     real(dp) :: tau_start

! \tau_e, end   point of the new segment
     real(dp) :: tau_end

! the possible maximum length of the new segment
     real(dp) :: tau_max

! ratio between old and new configurations, the local trace part
     real(dp) :: trace_ratio

! ratio between old and new configurations, the determinant part 
     real(dp) :: deter_ratio

! initialize logical variables
     ladd = .false.
     anti = .false.
     pass = .false.

! select the flavor channel randomly among 1 ~ norbs

     flvr = ceiling( spring_sfmt_stream() * norbs )

! get the perturbation expansion order ( number of existing segments or
! anti-segments ) for current flavor channel
     ckink = rank(flvr,atom)
     if ( ckink == mkink ) then
!<         call ctqmc_print_exception('ctqmc_insert_kink','can not insert any segments')
         insert_tcount = insert_tcount + one
         insert_reject = insert_reject + one
         RETURN
     endif

! randomly choose anti and tau_start, and then check whether tau_start is
! valid, if tau_start is valid, then determine tau_end, tau_max, is, and
! ie consistently, and set ladd to .true., if tau_start is not valid, then
! set ladd to .false.
     call ctqmc_make_flavor1(flvr, is, ie, anti, ladd, tau_start, tau_end, tau_max,atom)

! calculate the transition ratio between old and new configurations,
! for the local trace part
     if ( ladd .eqv. .true. ) then
         call cat_insert_ztrace(flvr, anti, tau_start, tau_end, trace_ratio,atom)
     else
         trace_ratio = zero
     endif

! calculate the transition ratio between old and new configurations,
! for the determinant part
     if ( ladd .eqv. .true. ) then
         call cat_insert_detrat(flvr, tau_start, tau_end, deter_ratio,atom)
     else
         deter_ratio = zero
     endif

! calculate the transition probability for insert new segment or anti-segment
     p = deter_ratio * trace_ratio * ( beta * tau_max / real( ckink + 1 ) )

! determine pass, using important sampling algorithm (metropolis algorithm)
     pass = ( min( one, abs(p) ) > spring_sfmt_stream() )

! if the update action is accepted
     if ( pass .eqv. .true. ) then

! update the mmat matrix and gmat matrix, respectively,
! cat_insert_segment() subroutine is invoked internally to update the
! perturbation expansion series
         call cat_insert_matrix(flvr, is, ie, tau_start, tau_end, deter_ratio,atom)

! update ckink for current flavor channel
         ckink = ckink + 1

! update stts for current flavor channel
         stts(flvr,atom) = cstat

! update rank for current flavor channel
         rank(flvr,atom) = rank(flvr,atom) + 1

     endif ! back if ( pass .eqv. .true. ) block

! update the insert statistics
     insert_tcount = insert_tcount + one
     if ( pass .eqv. .true. ) then
         insert_accept = insert_accept + one
     else
         insert_reject = insert_reject + one
     endif ! back if ( pass .eqv. .true. ) block

     return
  end subroutine ctqmc_insert_kink

!>>> remove old segment or anti-segment in the perturbation expansion series
  subroutine ctqmc_remove_kink(atom)
     use constants
     use control
     use context

     use spring

     implicit none

! local variables
! whether it is an anti-segment
! anti = .true., anti-segment
! anti = .false., segment

     integer, intent(in) :: atom

     logical  :: anti

! whether the update operation is accepted
     logical  :: pass

! current flavor channel for both band and spin
     integer  :: flvr

! index address to remove old segment or anti-segment
! is and ie are for start and end points, respectively
     integer  :: is, ie

! transition probability for remove old segment or anti-segment
     real(dp) :: p

! \tau_s, start point of the old segment
     real(dp) :: tau_start

! \tau_e, end   point of the old segment
     real(dp) :: tau_end

! the possible maximum length of the old segment
     real(dp) :: tau_max

! ratio between old and new configurations, the local trace part
     real(dp) :: trace_ratio

! ratio between old and new configurations, the determinant part 
     real(dp) :: deter_ratio

! initialize logical variables
     anti = .false.
     pass = .false.

! select the flavor channel randomly among 1 ~ norbs
     flvr = ceiling( spring_sfmt_stream() * norbs )

! get the perturbation expansion order ( number of existing segments or
! anti-segments ) for current flavor channel
     ckink = rank(flvr,atom)
     if ( ckink == 0 ) then
!<         call ctqmc_print_exception('ctqmc_remove_kink','can not remove any segments')
         remove_tcount = remove_tcount + one
         remove_reject = remove_reject + one
         RETURN
     endif

! at first determine anti and is randomly, then tau_start is obtained by
! is. and then ie, tau_end, and tau_max are evaluated carefully according
! to is and ie
     call ctqmc_make_flavor2(flvr, is, ie, anti, tau_start, tau_end, tau_max,atom)

! calculate the transition ratio between old and new configurations,
! for the local trace part
     call cat_remove_ztrace(flvr, anti, tau_start, tau_end, trace_ratio,atom)

! calculate the transition ratio between old and new configurations,
! for the determinant part
     call cat_remove_detrat(flvr, is, ie, deter_ratio,atom)

! calculate the transition probability for remove old segment or anti-segment
     p = deter_ratio * trace_ratio * real( ckink ) / ( beta * tau_max )

! determine pass, using important sampling algorithm (metropolis algorithm)
     pass = ( min( one, abs(p) ) > spring_sfmt_stream() )

! if update action is accepted
     if ( pass .eqv. .true. ) then

! update the mmat matrix and gmat matrix, respectively,
! cat_remove_segment() subroutine is invoked internally to update the
! perturbation expansion series
         call cat_remove_matrix(flvr, is, ie,atom)

! decrease ckink for current flavor channel
         ckink = ckink - 1

! update stts for current flavor channel
         stts(flvr,atom) = cstat

! update rank for current flavor channel
         rank(flvr,atom) = rank(flvr,atom) - 1

     endif ! back if ( pass .eqv. .true. ) block

! update the remove statistics
     remove_tcount = remove_tcount + one
     if ( pass .eqv. .true. ) then
         remove_accept = remove_accept + one
     else
         remove_reject = remove_reject + one
     endif ! back if ( pass .eqv. .true. ) block

     return
  end subroutine ctqmc_remove_kink

!>>> left shift old segment or anti-segment in the perturbation expansion series
  subroutine ctqmc_lshift_kink(atom)
     use constants
     use control
     use context

     use spring

     implicit none

     integer, intent(in) :: atom

! local variables
! whether the update operation winds around the circle
     logical  :: ring

! whether the update operation is accepted
     logical  :: pass

! current flavor channel for both band and spin
     integer  :: flvr

! index address to left shift old segment or anti-segment
! iso and isn are for old and new indices, respectively
     integer  :: iso, isn

! transition probability for left shift old segment or anti-segment
     real(dp) :: p

! \tau_s, start point of the old segment (old point)
     real(dp) :: tau_start1

! \tau_s, start point of the old segment (new point)
     real(dp) :: tau_start2

! ratio between old and new configurations, the local trace part
     real(dp) :: trace_ratio

! ratio between old and new configurations, the determinant part 
     real(dp) :: deter_ratio

! initialize logical variables
     ring = .false.
     pass = .false.

! select the flavor channel randomly among 1 ~ norbs
     flvr = ceiling( spring_sfmt_stream() * norbs )

! get the perturbation expansion order ( number of existing segments or
! anti-segments ) for current flavor channel
     ckink = rank(flvr,atom)
     if ( ckink == 0 ) then
!<         call ctqmc_print_exception('ctqmc_lshift_kink','can not lshift any segments')
         lshift_tcount = lshift_tcount + one
         lshift_reject = lshift_reject + one
         RETURN
     endif

! at first, we select iso randomly, and then obtain tau_start1. according
! to the existing segments, we determine tau_start2 and related index isn,
! finally ring is evaluated.
     call ctqmc_make_flavor3(flvr, iso, isn, ring, tau_start1, tau_start2,atom)

! calculate the transition ratio between old and new configurations,
! for the local trace part
     call cat_lshift_ztrace(flvr, ring, tau_start1, tau_start2, trace_ratio,atom)
 
! calculate the transition ratio between old and new configurations,
! for the determinant part
     call cat_lshift_detrat(flvr, iso, tau_start1, tau_start2, deter_ratio,atom)

! calculate the transition probability for left shift old segment or anti-segment
     p = deter_ratio * trace_ratio

! determine pass, using important sampling algorithm (metropolis algorithm)
     pass = ( min( one, abs(p) ) > spring_sfmt_stream() )

! if update action is accepted
     if ( pass .eqv. .true. ) then

! update the mmat matrix and gmat matrix, respectively, 
! cat_lshift_segment() subroutine is invoked internally to update the
! perturbation expansion series
         call cat_lshift_matrix(flvr, iso, isn, tau_start1, tau_start2, deter_ratio,atom)

! update stts for current flavor channel
         stts(flvr,atom) = cstat

     endif ! back if ( pass .eqv. .true. ) block

! update the lshift statistics
     lshift_tcount = lshift_tcount + one
     if ( pass .eqv. .true. ) then
         lshift_accept = lshift_accept + one
     else
         lshift_reject = lshift_reject + one
     endif ! back if ( pass .eqv. .true. ) block

     return
  end subroutine ctqmc_lshift_kink

!>>> right shift old segment or anti-segment in the perturbation expansion series
  subroutine ctqmc_rshift_kink(atom)
     use constants
     use control
     use context

     use spring

     implicit none

     integer, intent(in) :: atom

! local variables
! whether the update operation winds around the circle
     logical  :: ring

! whether the update operation is accepted
     logical  :: pass

! current flavor channel for both band and spin
     integer  :: flvr

! index address to right shift old segment or anti-segment
! ieo and ien are for old and new indices, respectively
     integer  :: ieo, ien

! transition probability for right shift old segment or anti-segment
     real(dp) :: p

! \tau_e, end   point of the old segment (old point)
     real(dp) :: tau_end1

! \tau_e, end   point of the old segment (new point)
     real(dp) :: tau_end2

! ratio between old and new configurations, the local trace part
     real(dp) :: trace_ratio

! ratio between old and new configurations, the determinant part 
     real(dp) :: deter_ratio

! initialize logical variables
     ring = .false.
     pass = .false.

! select the flavor channel randomly among 1 ~ norbs
     flvr = ceiling( spring_sfmt_stream() * norbs )

! get the perturbation expansion order ( number of existing segments or
! anti-segments ) for current flavor channel
     ckink = rank(flvr,atom)
     if ( ckink == 0 ) then
!<         call ctqmc_print_exception('ctqmc_rshift_kink','can not rshift any segments')
         rshift_tcount = rshift_tcount + one
         rshift_reject = rshift_reject + one
         RETURN
     endif

! at first, we select ieo randomly, and then obtain tau_end1. according
! to the existing segments, we determine tau_end2 and related index ien,
! finally ring is evaluated.
     call ctqmc_make_flavor4(flvr, ieo, ien, ring, tau_end1, tau_end2,atom)

! calculate the transition ratio between old and new configurations,
! for the local trace part
     call cat_rshift_ztrace(flvr, ring, tau_end1, tau_end2, trace_ratio,atom)
 
! calculate the transition ratio between old and new configurations,
! for the determinant part
     call cat_rshift_detrat(flvr, ieo, tau_end1, tau_end2, deter_ratio,atom)

! calculate the transition probability for right shift old segment or anti-segment
     p = deter_ratio * trace_ratio

! determine pass, using important sampling algorithm (metropolis algorithm)
     pass = ( min( one, abs(p) ) > spring_sfmt_stream() )

! if update action is accepted
     if ( pass .eqv. .true. ) then

! update the mmat matrix and gmat matrix, respectively,
! cat_rshift_segment() subroutine is invoked internally to update the
! perturbation expansion series
         call cat_rshift_matrix(flvr, ieo, ien, tau_end1, tau_end2, deter_ratio,atom)

! update stts for current flavor channel
         stts(flvr,atom) = cstat

     endif ! back if ( pass .eqv. .true. ) block

! update the rshift statistics
     rshift_tcount = rshift_tcount + one
     if ( pass .eqv. .true. ) then
         rshift_accept = rshift_accept + one
     else
         rshift_reject = rshift_reject + one
     endif ! back if ( pass .eqv. .true. ) block

     return
  end subroutine ctqmc_rshift_kink

!>>> perform a global update, exchange the segments and anti-segments
  subroutine ctqmc_reswap_kink(atom)
     use constants
     use control
     use context

     use spring

     implicit none
     
     integer, intent(in) :: atom

! local variables
! whether the update operation is accepted
     logical  :: pass

! current flavor channel for both band and spin
     integer  :: flvr

! transition probability for global swap
     real(dp) :: p

! ratio between old and new configurations, the local trace part
     real(dp) :: trace_ratio

! ratio between old and new configurations, the determinant part 
     real(dp) :: deter_ratio

! initialize logical variables
     pass = .false.

! select the flavor channel randomly among 1 ~ norbs
     flvr = ceiling( spring_sfmt_stream() * norbs )

! calculate the transition ratio between old and new configurations,
! for the local trace part
     call cat_reswap_ztrace(flvr, trace_ratio,atom)

! calculate the transition ratio between old and new configurations,
! for the determinant part
! note: here deter_ratio is -1 or 1 in fact
     call cat_reswap_detrat(flvr, deter_ratio,atom)

! calculate the transition probability for global swap
     p = deter_ratio * trace_ratio

! determine pass, using important sampling algorithm (metropolis algorithm)
     pass = ( min( one, abs(p) ) > spring_sfmt_stream() )

! if update action is accepted
     if ( pass .eqv. .true. ) then

! update the mmat matrix and gmat matrix, respectively,
! cat_reswap_segment() subroutine is invoked internally to update the
! perturbation expansion series
         call cat_reswap_matrix(flvr,atom)

! update stts for current flavor channel
         stts(flvr,atom) = cstat

     endif ! back if ( pass .eqv. .true. ) block

! update the reswap statistics
     reswap_tcount = reswap_tcount + one
     if ( pass .eqv. .true. ) then
         reswap_accept = reswap_accept + one
     else
         reswap_reject = reswap_reject + one
     endif ! back if ( pass .eqv. .true. ) block

     return
  end subroutine ctqmc_reswap_kink

!>>> perform a global update, exchange the states between spin up and spin
! down, it maybe useful for magnetic systems
  subroutine ctqmc_reflip_kink(cflip,atom)
     use constants
     use control
     use context

     use spring

     implicit none

! external arguments
! control flag
! if cflip = 1, flip inter-orbital spins randomly;
! if cflip = 2, flip intra-orbital spins one by one;
! if cflip = 3, flip intra-orbital spins globally
     integer, intent(in) :: cflip
     
     integer, intent(in) :: atom

! local variables
! whether the update operation is accepted
     logical  :: pass

! selected flavor pairs
     integer  :: fup, fdn

! loop index for flavor channel
     integer  :: flvr

! maximum rank order
     integer  :: kmax

! transition probability for global spin flip
     real(dp) :: p

! global flip determinant ratio
     real(dp) :: ratup
     real(dp) :: ratdn

! initialize logical variables
     pass = .false.

! initialize transition probability
     p = one

     if ( cflip == 1 ) then
! determine fup and fdn, and fup /= fdn
         fup = ceiling( spring_sfmt_stream() * norbs )
         do while ( .true. )
             fdn = ceiling( spring_sfmt_stream() * norbs )
             if ( fup /= fdn ) exit
         enddo ! over do while loop

! calculate the transition ratio between old and new configurations,
! for the determinant part, spin up case
         call cat_reflip_detrat(fup, fdn, ratup,atom)

! calculate the transition ratio between old and new configurations,
! for the determinant part, spin dn case
         call cat_reflip_detrat(fdn, fup, ratdn,atom)

! calculate the transition probability for global spin flip
         p = ratup * ratdn

! determine pass, using important sampling algorithm (metropolis algorithm)
         pass = ( min( one, abs(p) ) > spring_sfmt_stream() )

! if update action is accepted
         if ( pass .eqv. .true. ) then

! get maximum rank order in spin up and spin down states
             kmax = max( rank(fup,atom), rank(fdn,atom) )

! swap global variables between spin up and spin down states
             call cat_reflip_matrix(fup, fdn, kmax,atom)

         endif ! back if ( pass .eqv. .true. ) block

! update the reflip statistics
         reflip_tcount = reflip_tcount + one
         if ( pass .eqv. .true. ) then
             reflip_accept = reflip_accept + one
         else
             reflip_reject = reflip_reject + one
         endif ! back if ( pass .eqv. .true. ) block

     else if ( cflip == 2 ) then ! cflip = 2, local flip
         do flvr=1,nband

! get fup and fdn
             fup = flvr; fdn = flvr + nband

! calculate the transition ratio between old and new configurations,
! for the determinant part, spin up case
             call cat_reflip_detrat(fup, fdn, ratup,atom)

! calculate the transition ratio between old and new configurations,
! for the determinant part, spin dn case
             call cat_reflip_detrat(fdn, fup, ratdn,atom)

! calculate the transition probability for global spin flip
             p = ratup * ratdn

! determine pass, using important sampling algorithm (metropolis algorithm)
             pass = ( min( one, abs(p) ) > spring_sfmt_stream() )

! if update action is accepted
             if ( pass .eqv. .true. ) then

! get maximum rank order in spin up and spin down states
                 kmax = max( rank(fup,atom), rank(fdn,atom) )

! swap global variables between spin up and spin down states
                 call cat_reflip_matrix(fup, fdn, kmax,atom)

             endif ! back if ( pass .eqv. .true. ) block

! update the reflip statistics
             reflip_tcount = reflip_tcount + one
             if ( pass .eqv. .true. ) then
                 reflip_accept = reflip_accept + one
             else
                 reflip_reject = reflip_reject + one
             endif ! back if ( pass .eqv. .true. ) block

         enddo ! over flvr={1,nband} loop

     else if ( cflip == 3 ) then ! cflip = 3, global flip
         do flvr=1,nband

! get fup and fdn
             fup = flvr; fdn = flvr + nband

! calculate the transition ratio between old and new configurations,
! for the determinant part, spin up case
             call cat_reflip_detrat(fup, fdn, ratup,atom)

! calculate the transition ratio between old and new configurations,
! for the determinant part, spin dn case
             call cat_reflip_detrat(fdn, fup, ratdn,atom)

! calculate the transition probability for global spin flip
             p = p * ( ratup * ratdn )

         enddo ! over flvr={1,nband} loop

! determine pass, using important sampling algorithm (metropolis algorithm)
         pass = ( min( one, abs(p) ) > spring_sfmt_stream() )

! if update action is accepted
         if ( pass .eqv. .true. ) then

             do flvr=1,nband

! get fup and fdn
                 fup = flvr; fdn = flvr + nband

! get maximum rank order in spin up and spin down states
                 kmax = max( rank(fup,atom), rank(fdn,atom) )

! swap global variables between spin up and spin down states
                 call cat_reflip_matrix(fup, fdn, kmax,atom)

             enddo ! over flvr={1,nband} loop

         endif ! back if ( pass .eqv. .true. ) block

! update the reflip statistics
         reflip_tcount = reflip_tcount + one
         if ( pass .eqv. .true. ) then
             reflip_accept = reflip_accept + one
         else
             reflip_reject = reflip_reject + one
         endif ! back if ( pass .eqv. .true. ) block

     endif ! back if ( cflip == 1 ) block

     return
  end subroutine ctqmc_reflip_kink

!>>> global update all segments or anti-segments in the perturbation expansion series
  subroutine ctqmc_reload_kink(atom)
     use control
     use context

     implicit none
     
     integer, intent(in) :: atom
! local variables
! loop index for flavor channel
     integer :: flvr

     do flvr=1,norbs

! check the perturbation expansion order ( number of existing segments or
! anti-segments ) for current flavor channel
         if ( rank(flvr,atom) == 0 ) cycle

! regenerate the mmat matrix and gmat matrix from scratch
         call cat_reload_matrix(flvr,atom)

     enddo ! over flvr={1,norbs} loop

     return
  end subroutine ctqmc_reload_kink

!-------------------------------------------------------------------------
!>>> service layer: update M and G matrices                            <<<
!-------------------------------------------------------------------------

!>>> update the mmat matrix and gmat matrix for insert new segment
! or anti-segment
  subroutine cat_insert_matrix(flvr, is, ie, tau_start, tau_end, deter_ratio,atom)
     use constants
     use control
     use context

     implicit none

! external arguments
! current flavor channel
     integer, intent(in)  :: flvr
     
     integer, intent(in)  :: atom

! index address for insert new segment or anti-segment
     integer, intent(in)  :: is
     integer, intent(in)  :: ie

! imaginary time \tau_s for start point
     real(dp), intent(in) :: tau_start

! imaginary time \tau_e for end   point
     real(dp), intent(in) :: tau_end

! previous calculated determinant ratio
     real(dp), intent(in) :: deter_ratio

! local variables
! loop index over segments
     integer  :: i
     integer  :: j

! loop index over frequencies
     integer  :: k

! real(dp) dummy variables
     real(dp) :: p

! evaluate p at first
     p = one / deter_ratio

! shift lspace and rspace, and then supplement them with -1 at the end
     do i=ckink,ie,-1
         lspace(i+1, flvr,atom) = lspace(i, flvr,atom)
     enddo ! over i={ckink,ie,-1} loop
     lspace(ie, flvr,atom) = -one

     do j=ckink,is,-1
         rspace(j+1, flvr,atom) = rspace(j, flvr,atom)
     enddo ! over j={ckink,is,-1} loop
     rspace(is, flvr,atom) = -one

! scale lspace with p
     do i=1,ckink+1
         lspace(i, flvr,atom) = lspace(i, flvr,atom) * p
     enddo ! over i={1,ckink+1} loop

! shift mmat matrix
     do j=ckink,is,-1
         do i=ckink,ie,-1
             mmat(i+1, j+1, flvr,atom) = mmat(i, j, flvr,atom)
         enddo ! over i={ckink,ie,-1} loop
     enddo ! over j={ckink,is,-1} loop

     do j=ckink,is,-1
         do i=1,ie-1
             mmat(i, j+1, flvr,atom) = mmat(i, j, flvr,atom)
         enddo ! over i={1,ie-1} loop
     enddo ! over j={ckink,is,-1} loop
 
     do j=1,is-1
         do i=ckink,ie,-1
             mmat(i+1, j, flvr,atom) = mmat(i, j, flvr,atom)
         enddo ! over i={ckink,ie,-1} loop
     enddo ! over j={1,is-1} loop

! supplement mmat matrix with zero
     do i=1,ckink+1
         mmat(i, is, flvr,atom) = zero
     enddo ! over i={1,ckink+1} loop

     do j=1,ckink+1
         mmat(ie, j, flvr,atom) = zero
     enddo ! over j={1,ckink+1} loop

! finally evaluate mmat matrix
     do j=1,ckink+1
         do i=1,ckink+1
             mmat(i, j, flvr,atom) = mmat(i, j, flvr,atom) + lspace(i, flvr,atom) * rspace(j, flvr,atom)
         enddo ! over i={1,ckink+1} loop
     enddo ! over j={1,ckink+1} loop

! update the perturbation expansion series
     call cat_insert_segment(flvr, is, ie, tau_start, tau_end)

! update gmat matrix
     lsum(:, flvr,atom) = czero
     rsum(:, flvr,atom) = czero

     do i=1,ckink+1
         do k=1,nfreq
             lsum(k, flvr,atom) = lsum(k, flvr,atom) +         exp_e(k, index_e(i, flvr), flvr)   * lspace(i, flvr,atom)
             rsum(k, flvr,atom) = rsum(k, flvr,atom) + dconjg( exp_s(k, index_s(i, flvr), flvr) ) * rspace(i, flvr,atom)
         enddo ! over k={1,nfreq} loop
     enddo ! over i={1,ckink+1} loop

     p = one / beta ! we redefine p here
     do k=1,nfreq
         gmat(k, flvr, flvr,atom) = gmat(k, flvr, flvr,atom) - lsum(k, flvr,atom) * rsum(k, flvr,atom) * p
     enddo ! over k={1,nfreq} loop

! only for debug
!<     do i=1,ckink+1
!<         do j=1,ckink+1
!<             print *, 'M:', i, j, mmat(i, j, flvr)
!<         enddo ! over j={1,ckink+1} loop
!<     enddo ! over i={1,ckink+1} loop
!<
!<     print *, 'G1:', flvr, gmat(1, flvr, flvr)
!<     print *, 'G2:', flvr, gmat(2, flvr, flvr)
!<     print *, 'G3:', flvr, gmat(3, flvr, flvr)
!<     print *, 'Gn:', flvr, gmat(nfreq, flvr, flvr)

     return
  end subroutine cat_insert_matrix

!>>> update the mmat matrix and gmat matrix for remove old segment
! or anti-segment
  subroutine cat_remove_matrix(flvr, is, ie,atom)
     use constants
     use control
     use context

     implicit none

! external arguments
! current flavor channel
     integer, intent(in) :: flvr

     integer, intent(in) :: atom

! index address for remove old segment or anti-segment
     integer, intent(in) :: is
     integer, intent(in) :: ie

! local variables
! loop index over segments
     integer  :: i
     integer  :: j

! loop index over frequencies
     integer  :: k

! real(dp) dummy variables
     real(dp) :: p

! update gmat matrix
     lsum(:, flvr,atom) = czero
     rsum(:, flvr,atom) = czero

     do i=1,ckink
         do k=1,nfreq
             lsum(k, flvr,atom) = lsum(k, flvr,atom) +         exp_e(k, index_e(i, flvr), flvr)   * mmat(i, is, flvr,atom)
             rsum(k, flvr,atom) = rsum(k, flvr,atom) + dconjg( exp_s(k, index_s(i, flvr), flvr) ) * mmat(ie, i, flvr,atom)
         enddo ! over k={1,nfreq} loop
     enddo ! over i={1,ckink} loop

     p = one / ( mmat(ie, is, flvr,atom) * beta )
     do k=1,nfreq
         gmat(k, flvr, flvr,atom) = gmat(k, flvr, flvr,atom) + lsum(k, flvr,atom) * rsum(k, flvr,atom) * p
     enddo ! over k={1,nfreq} loop
     
! update mmat matrix
     p = one / mmat(ie, is, flvr,atom) ! we redefine p here
     do j=1,ckink
         do i=1,ckink
             if ( i /= ie .and. j /= is ) then
                 mmat(i, j, flvr,atom) = mmat(i, j, flvr,atom) - mmat(i, is, flvr,atom) * mmat(ie, j, flvr,atom) * p
             endif
         enddo ! over i={1,ckink} loop
     enddo ! over j={1,ckink} loop

     do j=is,ckink-1
         do i=ie,ckink-1
             mmat(i, j, flvr,atom) = mmat(i+1, j+1, flvr,atom)
         enddo ! over i={ie,ckink-1} loop
     enddo ! over j={is,ckink-1} loop

     do j=is,ckink-1
         do i=1,ie-1
             mmat(i, j, flvr,atom) = mmat(i, j+1, flvr,atom)
         enddo ! over i={1,ie-1} loop
     enddo ! over j={is,ckink-1} loop

     do j=1,is-1
         do i=ie,ckink-1
             mmat(i, j, flvr,atom) = mmat(i+1, j, flvr,atom)
         enddo ! over i={ie,ckink-1} loop
     enddo ! over j={1,is-1} loop

! update the perturbation expansion series
     call cat_remove_segment(flvr, is, ie)

! only for debug
!<     do i=1,ckink
!<         do j=1,ckink
!<             print *, 'M:', i, j, mmat(i, j, flvr)
!<         enddo ! over j={1,ckink} loop
!<     enddo ! over i={1,ckink} loop
!<
!<     print *, 'G1:', flvr, gmat(1, flvr, flvr)
!<     print *, 'G2:', flvr, gmat(2, flvr, flvr)
!<     print *, 'G3:', flvr, gmat(3, flvr, flvr)
!<     print *, 'Gn:', flvr, gmat(nfreq, flvr, flvr)

     return
  end subroutine cat_remove_matrix

!>>> update the mmat matrix and gmat matrix for left shift old segment
! or anti-segment
  subroutine cat_lshift_matrix(flvr, iso, isn, tau_start1, tau_start2, deter_ratio,atom)
     use constants
     use control
     use context

     implicit none

! external arguments
! current flavor channel
     integer, intent(in)  :: flvr

     integer, intent(in)  :: atom

! index address to left shift old segment or anti-segment
! iso and isn are for old and new indices, respectively
     integer, intent(in)  :: iso
     integer, intent(in)  :: isn

! imaginary time \tau_s for start point (the old one)
     real(dp), intent(in) :: tau_start1

! imaginary time \tau_s for start point (the new one)
     real(dp), intent(in) :: tau_start2

! previous calculated determinant ratio
     real(dp), intent(in) :: deter_ratio

! external arguments
! used to interpolate the hybridization function
     real(dp), external :: ctqmc_make_htau

! local variables
! loop index over segments
     integer  :: i
     integer  :: j

! loop index over frequencies
     integer  :: k

! used to store matrix element of mmat
     real(dp) :: md

! real(dp) dummy variables
     real(dp) :: xs
     real(dp) :: rs

! real(dp) dummy arrays, used to interpolate the hybridization function
     real(dp) :: lvec(mkink)
     real(dp) :: rvec(mkink)

! complex(dp) dummy arrays, used to calculate gmat matrix
     complex(dp) :: lexp(nfreq)
     complex(dp) :: gsum(nfreq)
     complex(dp) :: gdel(nfreq)

! evaluate lexp
     lexp = czero
     do k=1,nfreq
         xs = tau_start2 * rmesh(k)
         lexp(k) = - ( dcmplx( cos(xs), -sin(xs) ) - dconjg( exp_s(k, index_s(iso, flvr), flvr) ) ) / beta
     enddo ! over k={1,nfreq} loop

! evaluate gsum
     gsum = czero
     do i=1,ckink
         md = mmat(i, iso, flvr,atom)
         do k=1,nfreq
             gsum(k) = gsum(k) + exp_e(k, index_e(i, flvr), flvr) * md
         enddo ! over k={1,nfreq} loop
     enddo ! over i={1,ckink} loop

! evaluate gdel, \delta G for gmat matrix
     gdel = czero
     do k=1,nfreq
         gdel(k) = gsum(k) * lexp(k)
     enddo ! over k={1,nfreq} loop

! calculate rvec by cubic spline interpolation
     do i=1,ckink
         if ( tau_start1 < time_e(index_e(i, flvr), flvr) ) then
             rvec(i) = -ctqmc_make_htau(flvr, tau_start1 - time_e(index_e(i, flvr), flvr) + beta,atom)
         else
             rvec(i) =  ctqmc_make_htau(flvr, tau_start1 - time_e(index_e(i, flvr), flvr),atom)
         endif
     enddo ! over i={1,ckink} loop

! calculate lvec by cubic spline interpolation
     do j=1,ckink
         if ( tau_start2 < time_e(index_e(j, flvr), flvr) ) then
             lvec(j) = -ctqmc_make_htau(flvr, tau_start2 - time_e(index_e(j, flvr), flvr) + beta,atom)
         else
             lvec(j) =  ctqmc_make_htau(flvr, tau_start2 - time_e(index_e(j, flvr), flvr),atom)
         endif
     enddo ! over j={1,ckink} loop

! adjust rvec
     do i=1,ckink
         rvec(i) = lvec(i) - rvec(i)
     enddo ! over i={1,ckink} loop

! prepare rspace
     do i=1,ckink
         rs = zero
         do j=1,ckink
             rs = rs + rvec(j) * mmat(j, i, flvr,atom)
         enddo ! over j={1,ckink} loop
         rspace(i, flvr,atom) = rs / deter_ratio
     enddo ! over i={1,ckink} loop

! prepare lspace
     do i=1,ckink
         lspace(i, flvr,atom) = -mmat(i, iso, flvr,atom)
     enddo ! over i={1,ckink} loop

! calculate mmat matrix
     do j=1,ckink
         do i=1,ckink
             mmat(i, j, flvr,atom) = mmat(i, j, flvr,atom) + lspace(i, flvr,atom) * rspace(j, flvr,atom)
         enddo ! over i={1,ckink} loop
     enddo ! over j={1,ckink} loop

! shufle rows if time order changed because of move
     if ( isn /= iso ) then
         rs = rspace(iso, flvr,atom)
         do i=1,ckink
             rvec(i) = mmat(i, iso, flvr,atom)
         enddo ! over i={1,ckink} loop

         if ( isn < iso ) then
             do j=iso,isn+1,-1
                 do i=1,ckink
                     mmat(i, j, flvr,atom) = mmat(i, j-1, flvr,atom)
                 enddo ! over i={1,ckink} loop
                 rspace(j, flvr,atom) = rspace(j-1, flvr,atom)
             enddo ! over j={iso,isn+1,-1} loop
             do i=1,ckink
                 mmat(i, isn, flvr,atom) = rvec(i)
             enddo ! over i={1,ckink} loop
             rspace(isn, flvr,atom) = rs
         else
             do j=iso,isn-1
                 do i=1,ckink
                     mmat(i, j, flvr,atom) = mmat(i, j+1, flvr,atom)
                 enddo ! over i={1,ckink} loop
                 rspace(j, flvr,atom) = rspace(j+1, flvr,atom)
             enddo ! over j={iso,isn-1} loop
             do i=1,ckink
                 mmat(i, isn, flvr,atom) = rvec(i)
             enddo ! over i={1,ckink} loop
             rspace(isn, flvr,atom) = rs
         endif ! back if ( isn < iso ) block
     endif ! back if ( isn /= iso ) block

! update the perturbation expansion series
     call cat_lshift_segment(flvr, iso, isn, tau_start2)

! update gmat matrix
     lsum(:, flvr,atom) = czero
     rsum(:, flvr,atom) = czero

     do i=1,ckink
         do k=1,nfreq
             lsum(k, flvr,atom) = lsum(k, flvr,atom) +         exp_e(k, index_e(i, flvr), flvr)   * lspace(i, flvr,atom)
             rsum(k, flvr,atom) = rsum(k, flvr,atom) + dconjg( exp_s(k, index_s(i, flvr), flvr) ) * rspace(i, flvr,atom)
         enddo ! over k={1,nfreq} loop 
     enddo ! over i={1,ckink+1} loop

     do k=1,nfreq
         gmat(k, flvr, flvr,atom) = gmat(k, flvr, flvr,atom) - lsum(k, flvr,atom) * rsum(k, flvr,atom) / beta + gdel(k)
     enddo ! over k={1,nfreq} loop

! only for debug
!<     do i=1,ckink
!<         do j=1,ckink
!<             print *,'M:',i, j, mmat(i, j, flvr)
!<         enddo ! over j={1,ckink} loop
!<     enddo ! over i={1,ckink} loop
!<
!<     print *, 'G1:', flvr, gmat(1, flvr, flvr)
!<     print *, 'G2:', flvr, gmat(2, flvr, flvr)
!<     print *, 'G3:', flvr, gmat(3, flvr, flvr)
!<     print *, 'Gn:', flvr, gmat(nfreq, flvr, flvr)

     return
  end subroutine cat_lshift_matrix

!>>> update the mmat matrix and gmat matrix for right shift old segment
! or anti-segment
  subroutine cat_rshift_matrix(flvr, ieo, ien, tau_end1, tau_end2, deter_ratio,atom)
     use constants
     use control
     use context

     implicit none

! external arguments
! current flavor channel
     integer, intent(in)  :: flvr

     integer, intent(in)  :: atom

! index address to right shift old segment or anti-segment
! ieo and ien are for old and new indices, respectively
     integer, intent(in)  :: ieo
     integer, intent(in)  :: ien

! imaginary time \tau_e for end point (the old one)
     real(dp), intent(in) :: tau_end1

! imaginary time \tau_e for end point (the new one)
     real(dp), intent(in) :: tau_end2

! previous calculated determinant ratio
     real(dp), intent(in) :: deter_ratio

! external arguments
! used to interpolate the hybridization function
     real(dp), external :: ctqmc_make_htau

! local variables
! loop index over segments
     integer  :: i
     integer  :: j

! loop index over frequencies
     integer  :: k

! used to store matrix element of mmat
     real(dp) :: md

! real(dp) dummy variables
     real(dp) :: xe
     real(dp) :: ls

! real(dp) dummy arrays, used to interpolate the hybridization function
     real(dp) :: lvec(mkink)
     real(dp) :: rvec(mkink)

! complex(dp) dummy arrays, used to calculate gmat matrix
     complex(dp) :: rexp(nfreq)
     complex(dp) :: gsum(nfreq)
     complex(dp) :: gdel(nfreq)

! evaluate rexp
     rexp = czero
     do k=1,nfreq
         xe = tau_end2 * rmesh(k)
         rexp(k) = - ( dcmplx( cos(xe), sin(xe) ) - exp_e(k, index_e(ieo, flvr), flvr) ) / beta
     enddo ! over k={1,nfreq} loop

! evaluate gsum
     gsum = czero
     do i=1,ckink
         md = mmat(ieo, i, flvr,atom)
         do k=1,nfreq
             gsum(k) = gsum(k) + dconjg( exp_s(k, index_s(i, flvr), flvr) ) * md
         enddo ! over k={1,nfreq} loop
     enddo ! over i={1,ckink} loop

! evaluate gdel, \delta G for gmat matrix
     gdel = czero
     do k=1,nfreq
         gdel(k) = gsum(k) * rexp(k)
     enddo ! over k={1,nfreq} loop

! calculate lvec by cubic spline interpolation
     do i=1,ckink
         if ( time_s(index_s(i, flvr), flvr) < tau_end1 ) then
             lvec(i) = -ctqmc_make_htau(flvr, time_s(index_s(i, flvr), flvr) - tau_end1 + beta,atom)
         else
             lvec(i) =  ctqmc_make_htau(flvr, time_s(index_s(i, flvr), flvr) - tau_end1,atom)
         endif
     enddo ! over i={1,ckink} loop

! calculate rvec by cubic spline interpolation
     do j=1,ckink
         if ( time_s(index_s(j, flvr), flvr) < tau_end2 ) then
             rvec(j) = -ctqmc_make_htau(flvr, time_s(index_s(j, flvr), flvr) - tau_end2 + beta,atom)
         else
             rvec(j) =  ctqmc_make_htau(flvr, time_s(index_s(j, flvr), flvr) - tau_end2,atom)
         endif
     enddo ! over j={1,ckink} loop

! adjust lvec
     do i=1,ckink
         lvec(i) = rvec(i) - lvec(i)
     enddo ! over i={1,ckink} loop

! prepare lspace
     do i=1,ckink
         ls = zero
         do j=1,ckink
             ls = ls + mmat(i, j, flvr,atom) * lvec(j)
         enddo ! over j={1,ckink} loop
         lspace(i, flvr,atom) = ls / deter_ratio
     enddo ! over i={1,ckink} loop

! prepare rspace
     do i=1,ckink
         rspace(i, flvr,atom) = -mmat(ieo, i, flvr,atom)
     enddo ! over i={1,ckink} loop

! calculate mmat matrix
     do j=1,ckink
         do i=1,ckink
             mmat(i, j, flvr,atom) = mmat(i, j, flvr,atom) + lspace(i, flvr,atom) * rspace(j, flvr,atom)
         enddo ! over i={1,ckink} loop
     enddo ! over j={1,ckink} loop

! shufle columns if time order changed because of move
     if ( ien /= ieo ) then
         ls = lspace(ieo, flvr,atom)
         do i=1,ckink
             lvec(i) = mmat(ieo, i, flvr,atom)
         enddo ! over i={1,ckink} loop

         if ( ien < ieo ) then
             do j=ieo,ien+1,-1
                 do i=1,ckink
                     mmat(j, i, flvr,atom) = mmat(j-1, i, flvr,atom)
                 enddo ! over i={1,ckink} loop
                 lspace(j, flvr,atom) = lspace(j-1, flvr,atom)
             enddo ! over j={ieo,ien+1,-1} loop
             do i=1,ckink
                 mmat(ien, i, flvr,atom) = lvec(i)
             enddo ! over i={1,ckink} loop
             lspace(ien, flvr,atom) = ls
         else
             do j=ieo,ien-1
                 do i=1,ckink
                     mmat(j, i, flvr,atom) = mmat(j+1, i, flvr,atom)
                 enddo ! over i={1,ckink} loop
                 lspace(j, flvr,atom) = lspace(j+1, flvr,atom)
             enddo ! over j={ieo,ien-1} loop
             do i=1,ckink
                 mmat(ien, i, flvr,atom) = lvec(i)
             enddo ! over i={1,ckink} loop
             lspace(ien, flvr,atom) = ls
         endif ! back if ( ien < ieo ) block
     endif ! back if ( ien /= ieo ) block

! update the perturbation expansion series
     call cat_rshift_segment(flvr, ieo, ien, tau_end2)

! update gmat matrix
     lsum(:, flvr,atom) = czero
     rsum(:, flvr,atom) = czero

     do i=1,ckink
         do k=1,nfreq
             lsum(k, flvr,atom) = lsum(k, flvr,atom) +         exp_e(k, index_e(i, flvr), flvr)   * lspace(i, flvr,atom)
             rsum(k, flvr,atom) = rsum(k, flvr,atom) + dconjg( exp_s(k, index_s(i, flvr), flvr) ) * rspace(i, flvr,atom)
         enddo ! over k={1,nfreq} loop 
     enddo ! over i={1,ckink+1} loop

     do k=1,nfreq
         gmat(k, flvr, flvr,atom) = gmat(k, flvr, flvr,atom) - lsum(k, flvr,atom) * rsum(k, flvr,atom) / beta + gdel(k)
     enddo ! over k={1,nfreq} loop

! only for debug
!<     do i=1,ckink
!<         do j=1,ckink
!<             print *,'M:',i, j, mmat(i, j, flvr)
!<         enddo ! over j={1,ckink} loop
!<     enddo ! over i={1,ckink} loop
!<
!<     print *, 'G1:', flvr, gmat(1, flvr, flvr)
!<     print *, 'G2:', flvr, gmat(2, flvr, flvr)
!<     print *, 'G3:', flvr, gmat(3, flvr, flvr)
!<     print *, 'Gn:', flvr, gmat(nfreq, flvr, flvr)

     return
  end subroutine cat_rshift_matrix

!>>> global swap the time_s, time_e, exp_s, and exp_e, and then update
! mmat and gmat matrix. it is used to overcome the low acceptance ratio
! at high temperature region
  subroutine cat_reswap_matrix(flvr,atom)
     use constants
     use control
     use context

     implicit none

! external arguments
! current flavor channel
     integer, intent(in) :: flvr
     integer, intent(in) :: atom

! determine the new status of current flavor channel after swap operation
     select case ( stts(flvr,atom) )

         case (0)
             cstat = 3

         case (1)
             cstat = 2

         case (2)
             cstat = 1

         case (3)
             cstat = 0

     end select

! swap time_s, time_e, exp_s, and exp_e
     call cat_reswap_segment(flvr,atom)

! regenerate mmat and gmat matrix if need
     if ( ckink > 0 ) then
         call cat_reload_matrix(flvr,atom)
     endif ! back if ( ckink > 0 ) block

     return
  end subroutine cat_reswap_matrix

!>>> global flip the time_s, time_e, mmat matrix, gmat matrix, and other
! related global variables between spin up and spin down states. it is 
! used to avoid trapped by unphysical phase
  subroutine cat_reflip_matrix(fup, fdn, kmax,atom)
     use constants
     use control
     use context

     use stack

     implicit none

! external arguments
! current flavor channel
     integer, intent(in) :: fup
     integer, intent(in) :: fdn
     integer, intent(in) :: atom

! maximum rank order for current flavor channel
     integer, intent(in) :: kmax

! local variables
! maximum memory index accessed by index_s and index_e
     integer :: ismax
     integer :: iemax

! dummy copy for rank and stts
     integer :: Trank
     integer :: Tstts

! dummy copy for empty_s and empty_e
     type (istack) :: Tempty_s
     type (istack) :: Tempty_e

! allocate memory for Tempty_s and Tempty_e
     Tempty_s = istack_create(mkink)
     Tempty_e = istack_create(mkink)

! swap empty_s and empty_e
     call istack_copyer(empty_s(fup), Tempty_s)
     call istack_copyer(empty_e(fup), Tempty_e)

     call istack_copyer(empty_s(fdn), empty_s(fup))
     call istack_copyer(empty_e(fdn), empty_e(fup))

     call istack_copyer(Tempty_s, empty_s(fdn))
     call istack_copyer(Tempty_e, empty_e(fdn))

! deallocate memory for Tempty_s and Tempty_e
     call istack_destroy(Tempty_s)
     call istack_destroy(Tempty_e)

! swap rank
     Trank = rank(fup,atom)
     rank(fup,atom) = rank(fdn,atom)
     rank(fdn,atom) = Trank

! swap stts
     Tstts = stts(fup,atom)
     stts(fup,atom) = stts(fdn,atom)
     stts(fdn,atom) = Tstts

! swap gmat matrix when needed
     call zswap(nfreq, gmat(1:nfreq, fup, fup,atom), 1, gmat(1:nfreq, fdn, fdn,atom), 1)

     if ( kmax > 0 ) then

! determine ismax and iemax
         ismax = max( maxval( index_s(1:kmax, fup) ), maxval( index_s(1:kmax, fdn) ) )
         iemax = max( maxval( index_e(1:kmax, fup) ), maxval( index_e(1:kmax, fdn) ) )

! swap index_s and index_e
         call iswap(kmax, index_s(1:kmax, fup), index_s(1:kmax, fdn))
         call iswap(kmax, index_e(1:kmax, fup), index_e(1:kmax, fdn))

! swap time_s and time_e
         call dswap(ismax, time_s(1:ismax, fup), 1, time_s(1:ismax, fdn), 1)
         call dswap(iemax, time_e(1:iemax, fup), 1, time_e(1:iemax, fdn), 1)

! swap exp_s and exp_e
         call zswap(nfreq*ismax, exp_s(1:nfreq, 1:ismax, fup), 1, exp_s(1:nfreq, 1:ismax, fdn), 1)
         call zswap(nfreq*iemax, exp_e(1:nfreq, 1:iemax, fup), 1, exp_e(1:nfreq, 1:iemax, fdn), 1)

! update mmat and gmat matrix when needed
         if ( rank(fup,atom) > 0 ) call cat_reload_matrix(fup,atom)
         if ( rank(fdn,atom) > 0 ) call cat_reload_matrix(fdn,atom)

     endif ! back if ( kmax > 0 ) block

     return

  contains

!>>> extended BLAS subroutines, exchange two integer vectors
  pure subroutine iswap(n, ix, iy)
     implicit none

! external arguments
! dimension of integer vector
     integer, intent(in) :: n

! integer vector X
     integer, intent(inout) :: ix(n)

! integer vector Y
     integer, intent(inout) :: iy(n)

! local variables
! dummy integer vector
     integer :: it(n)

     it = ix
     ix = iy
     iy = it

     return
  end subroutine iswap

  end subroutine cat_reflip_matrix

!>>> global update the mmat matrix and gmat matrix from scratch
  subroutine cat_reload_matrix(flvr,atom)
     use constants
     use control
     use context

     implicit none

! external arguments
! current flavor channel
     integer, intent(in) :: flvr

     integer, intent(in) :: atom

! external functions
! used to interpolate the hybridization function
     real(dp), external :: ctqmc_make_htau

! local variables
! loop index over segments
     integer  :: i
     integer  :: j

! loop index over frequencies
     integer  :: k

! dummy perturbation expansion order
     integer  :: kaux

! used to store matrix element of mmat
     real(dp) :: maux

! imaginary time for start and end points
     real(dp) :: tau_start
     real(dp) :: tau_end

! complex(dp) dummy variables
     complex(dp) :: x_start
     complex(dp) :: x_end

! evaluate kaux
     kaux = rank(flvr,atom)

! reset mmat matrix
     mmat(1:kaux, 1:kaux, flvr,atom) = zero

! recalculate mmat from scratch
     do j=1,kaux
         tau_end = time_e(index_e(j, flvr), flvr)
         do i=1,kaux
             tau_start = time_s(index_s(i, flvr), flvr)
             if ( tau_start < tau_end ) then
                 mmat(i, j, flvr,atom) = -ctqmc_make_htau(flvr, tau_start - tau_end + beta,atom)
             else
                 mmat(i, j, flvr,atom) =  ctqmc_make_htau(flvr, tau_start - tau_end,atom)
             endif
         enddo ! over i={1,kaux} loop
     enddo ! over j={1,kaux} loop

! now we obtain dmat matrix, while what we need is its inversion
     call ctqmc_dmat_inv(kaux, mmat(1:kaux, 1:kaux, flvr,atom))

! reset gmat matrix
     gmat(:, flvr, flvr,atom) = czero

! recalculate gmat from scratch
     do j=1,kaux
         do i=1,kaux
             maux = -mmat(i, j, flvr,atom) / beta
             do k=1,nfreq
                 x_start = dconjg( exp_s(k, index_s(j, flvr), flvr) )
                 x_end   =         exp_e(k, index_e(i, flvr), flvr)
                 gmat(k, flvr, flvr,atom) = gmat(k, flvr, flvr,atom) + x_end * maux * x_start
             enddo ! over k={1,nfreq} loop
         enddo ! over i={1,kaux} loop
     enddo ! over j={1,kaux} loop

     return
  end subroutine cat_reload_matrix

!-------------------------------------------------------------------------
!>>> service layer: evaluate the determinant ratio                     <<<
!-------------------------------------------------------------------------

!>>> calculate the determinant ratio for insert new segment or anti-segment
  subroutine cat_insert_detrat(flvr, tau_start, tau_end, deter_ratio,atom)
     use constants
     use control
     use context

     implicit none

! external arguments
! current flavor channel
     integer, intent(in)   :: flvr

     integer, intent(in)   :: atom

! imaginary time \tau_s for start point
     real(dp), intent(in)  :: tau_start

! imaginary time \tau_e for end   point
     real(dp), intent(in)  :: tau_end

! the desired determinant ratio
     real(dp), intent(out) :: deter_ratio

! external arguments
! used to interpolate the hybridization function
     real(dp), external :: ctqmc_make_htau

! local variables
! loop index over segments
     integer  :: i
     integer  :: j

! real(dp) dummy variables
     real(dp) :: sl
     real(dp) :: sr

! real(dp) dummy arrays, used to interpolate the hybridization function
     real(dp) :: lvec(mkink)
     real(dp) :: rvec(mkink)

! calculate lvec by cubic spline interpolation
     do i=1,ckink
         if ( time_s(index_s(i, flvr), flvr) < tau_end   ) then
             lvec(i) = -ctqmc_make_htau(flvr, time_s(index_s(i, flvr), flvr) - tau_end + beta,atom)
         else
             lvec(i) =  ctqmc_make_htau(flvr, time_s(index_s(i, flvr), flvr) - tau_end,atom)
         endif
     enddo ! over i={1,ckink} loop

! calculate rvec by cubic spline interpolation
     do j=1,ckink
         if ( tau_start < time_e(index_e(j, flvr), flvr) ) then
             rvec(j) = -ctqmc_make_htau(flvr, tau_start - time_e(index_e(j, flvr), flvr) + beta,atom)
         else
             rvec(j) =  ctqmc_make_htau(flvr, tau_start - time_e(index_e(j, flvr), flvr),atom)
         endif
     enddo ! over j={1,ckink} loop

! calculate deter_ratio by cubic spline interpolation
     if ( tau_start > tau_end ) then
         deter_ratio =  ctqmc_make_htau(flvr, tau_start - tau_end,atom)
     else
         deter_ratio = -ctqmc_make_htau(flvr, tau_start - tau_end + beta,atom)
     endif

! calculate lspace and rspace
     do i=1,ckink
         sl = zero
         sr = zero

         do j=1,ckink
             sl = sl + mmat(i, j, flvr,atom) * lvec(j)
             sr = sr + rvec(j) * mmat(j, i, flvr,atom)
         enddo ! over j={1,ckink} loop

         lspace(i, flvr,atom) = sl
         rspace(i, flvr,atom) = sr
     enddo ! over i={1,ckink} loop

! calculate final determinant ratio
     do i=1,ckink
         deter_ratio = deter_ratio - rvec(i) * lspace(i, flvr,atom)
     enddo ! over i={1,ckink} loop

     return
  end subroutine cat_insert_detrat

!>>> calculate the determinant ratio for remove old segment or anti-segment
  subroutine cat_remove_detrat(flvr, is, ie, deter_ratio,atom)
     use context

     implicit none

! external arguments
! current flavor channel
     integer, intent(in)   :: flvr

     integer, intent(in)   :: atom

! index address to remove old segment or anti-segment
! is and ie are for start and end points, respectively
     integer, intent(in)   :: is
     integer, intent(in)   :: ie

! the desired determinant ratio
     real(dp), intent(out) :: deter_ratio

     deter_ratio = mmat(ie, is, flvr,atom)

     return
  end subroutine cat_remove_detrat

!>>> calculate the determinant ratio for left shift old segment or anti-segment
  subroutine cat_lshift_detrat(flvr, addr, tau_start1, tau_start2, deter_ratio,atom)
     use constants
     use control
     use context

     implicit none

! external arguments
! current flavor channel
     integer, intent(in)   :: flvr
     integer, intent(in)   :: atom

! index address for left  shift old segment or anti-segment (old index, = iso)
     integer, intent(in)   :: addr

! imaginary time \tau_s for start point (the old one)
     real(dp), intent(in)  :: tau_start1

! imaginary time \tau_s for start point (the new one)
     real(dp), intent(in)  :: tau_start2

! the desired determinant ratio
     real(dp), intent(out) :: deter_ratio

! external functions
! used to interpolate the hybridization function
     real(dp), external    :: ctqmc_make_htau

! local variables
! loop index over segments
     integer  :: i
     integer  :: j

! real(dp) dummy arrays, used to interpolate the hybridization function
     real(dp) :: lvec(mkink)
     real(dp) :: rvec(mkink)

! calculate rvec by cubic spline interpolation
     do i=1,ckink
         if ( tau_start1 < time_e(index_e(i, flvr), flvr) ) then
             rvec(i) = -ctqmc_make_htau(flvr, tau_start1 - time_e(index_e(i, flvr), flvr) + beta,atom)
         else
             rvec(i) =  ctqmc_make_htau(flvr, tau_start1 - time_e(index_e(i, flvr), flvr),atom)
         endif
     enddo ! over i={1,ckink} loop

! calculate lvec by cubic spline interpolation
     do j=1,ckink
         if ( tau_start2 < time_e(index_e(j, flvr), flvr) ) then
             lvec(j) = -ctqmc_make_htau(flvr, tau_start2 - time_e(index_e(j, flvr), flvr) + beta,atom)
         else
             lvec(j) =  ctqmc_make_htau(flvr, tau_start2 - time_e(index_e(j, flvr), flvr),atom)
         endif
     enddo ! over j={1,ckink} loop

! adjust rvec
     do i=1,ckink
         rvec(i) = lvec(i) - rvec(i)
     enddo ! over i={1,ckink} loop

! calculate final determinant ratio
     deter_ratio = one
     do i=1,ckink
         deter_ratio = deter_ratio + rvec(i) * mmat(i, addr, flvr,atom)
     enddo ! over i={1,ckink} loop

     return
  end subroutine cat_lshift_detrat

!>>> calculate the determinant ratio for right shift old segment or anti-segment
  subroutine cat_rshift_detrat(flvr, addr, tau_end1, tau_end2, deter_ratio,atom)
     use constants
     use control
     use context

     implicit none

! external arguments
! current flavor channel
     integer, intent(in)   :: flvr

     integer, intent(in)   :: atom

! index address for right shift old segment or anti-segment (old index = ieo)
     integer, intent(in)   :: addr

! imaginary time \tau_e for end   point (the old one)
     real(dp), intent(in)  :: tau_end1

! imaginary time \tau_e for end   point (the new one)
     real(dp), intent(in)  :: tau_end2

! the desired determinant ratio
     real(dp), intent(out) :: deter_ratio

! external functions
! used to interpolate the hybridization function
     real(dp), external    :: ctqmc_make_htau

! local variables
! loop index over segments
     integer  :: i
     integer  :: j

! real(dp) dummy arrays, used to interpolate the hybridization function
     real(dp) :: lvec(mkink)
     real(dp) :: rvec(mkink)

! calculate lvec by cubic spline interpolation
     do i=1,ckink
         if ( time_s(index_s(i, flvr), flvr) < tau_end1 ) then
             lvec(i) = -ctqmc_make_htau(flvr, time_s(index_s(i, flvr), flvr) - tau_end1 + beta,atom)
         else
             lvec(i) =  ctqmc_make_htau(flvr, time_s(index_s(i, flvr), flvr) - tau_end1,atom)
         endif
     enddo ! over i={1,ckink} loop

! calculate rvec by cubic spline interpolation
     do j=1,ckink
         if ( time_s(index_s(j, flvr), flvr) < tau_end2 ) then
             rvec(j) = -ctqmc_make_htau(flvr, time_s(index_s(j, flvr), flvr) - tau_end2 + beta,atom)
         else
             rvec(j) =  ctqmc_make_htau(flvr, time_s(index_s(j, flvr), flvr) - tau_end2,atom)
         endif
     enddo ! over j={1,ckink} loop

! adjust lvec
     do i=1,ckink
         lvec(i) = rvec(i) - lvec(i)
     enddo ! over i={1,ckink} loop

! calculate final determinant ratio
     deter_ratio = one
     do i=1,ckink
         deter_ratio = deter_ratio + mmat(addr, i, flvr,atom) * lvec(i) 
     enddo ! over i={1,ckink} loop

     return
  end subroutine cat_rshift_detrat

!>>> calculate the determinant ratio for global segment swap
  subroutine cat_reswap_detrat(flvr, ratio,atom)
     use constants
     use control
     use context

     implicit none

! external arguments
! current flavor channel
     integer, intent(in)   :: flvr
     
     integer, intent(in)   :: atom

! the desired determinant ratio
     real(dp), intent(out) :: ratio

! external functions
! used to interpolate the hybridization function
     real(dp), external :: ctqmc_make_htau

! local variables
! loop index over segments
     integer  :: i
     integer  :: j

! dummy perturbation expansion order
     integer  :: kaux

! status flag
     integer  :: istat

! imaginary time for start and end points
     real(dp) :: tau_start
     real(dp) :: tau_end

! dummy mmat matrix
     real(dp), allocatable :: Dmm(:,:)
     real(dp), allocatable :: Tmm(:,:)

! evaluate kaux
     kaux = rank(flvr,atom)

! check the status of kaux, if there does not exist any operators in flvr
! state ( kaux == 0 ), we need to return immediately and the ratio is one
     if ( kaux == 0 ) then
         ratio = one
         RETURN
     endif ! back if ( kaux == 0 ) block

! allocate memory
     allocate(Dmm(kaux,kaux), stat=istat)
     allocate(Tmm(kaux,kaux), stat=istat)
     if ( istat /= 0 ) then
         call ctqmc_print_error('cat_reswap_detrat','can not allocate enough memory')
     endif

! init Dmm and Tmm matrix
     Dmm = zero
     Tmm = zero

! recalculate Dmm from scratch
     do j=1,kaux
         tau_end = time_s(index_s(j, flvr), flvr)
         do i=1,kaux
             tau_start = time_e(index_e(i, flvr), flvr)
             if ( tau_start < tau_end ) then
                 Dmm(i, j) = -ctqmc_make_htau(flvr, tau_start - tau_end + beta,atom)
             else
                 Dmm(i, j) =  ctqmc_make_htau(flvr, tau_start - tau_end,atom)
             endif
         enddo ! over i={1,kaux} loop
     enddo ! over j={1,kaux} loop

! calculate Tmm: Tmm = Dmm . [M wDMmat
     call dgemm('N', 'N', kaux, kaux, kaux, one, Dmm, kaux, mmat(1:kaux, 1:kaux, flvr,atom), kaux, zero, Tmm, kaux)

! calculate the determinant of Tmm, it is the desired ratio
! note: in fact, here the ratio is -1 exactly!
     call ctqmc_dmat_det(kaux, Tmm, ratio)

! deallocate memory
     deallocate(Dmm)
     deallocate(Tmm)

     return
  end subroutine cat_reswap_detrat

!>>> calculate the determinant ratio for global spin flip
  subroutine cat_reflip_detrat(up, dn, ratio,atom)
     use constants
     use control
     use context

     implicit none

! external arguments
! band index for spin up state
     integer, intent(in)   :: up

     integer, intent(in)   :: atom

! band index for spin dn state
     integer, intent(in)   :: dn

! the desired determinant ratio
     real(dp), intent(out) :: ratio

! external functions
! used to interpolate the hybridization function
     real(dp), external :: ctqmc_make_htau

! local variables
! loop index over segments
     integer  :: i
     integer  :: j

! dummy perturbation expansion order
     integer  :: kaux

! status flag
     integer  :: istat

! imaginary time for start and end points
     real(dp) :: tau_start
     real(dp) :: tau_end

! dummy mmat matrix
     real(dp), allocatable :: Dmm(:,:)
     real(dp), allocatable :: Tmm(:,:)

! evaluate kaux
     kaux = rank(up,atom)

! check the status of kaux, if there does not exist any operators in up
! state ( kaux == 0 ), we need to return immediately and the ratio is one
     if ( kaux == 0 ) then
         ratio = one
         RETURN
     endif ! back if ( kaux == 0 ) block

! allocate memory
     allocate(Dmm(kaux,kaux), stat=istat)
     allocate(Tmm(kaux,kaux), stat=istat)
     if ( istat /= 0 ) then
         call ctqmc_print_error('cat_reflip_detrat','can not allocate enough memory')
     endif

! init Dmm and Tmm matrix
     Dmm = zero
     Tmm = zero

! recalculate Dmm from scratch
     do j=1,kaux
         tau_end = time_e(index_e(j, up), up)
         do i=1,kaux
             tau_start = time_s(index_s(i, up), up)
             if ( tau_start < tau_end ) then
                 Dmm(i, j) = -ctqmc_make_htau(dn, tau_start - tau_end + beta,atom)
             else
                 Dmm(i, j) =  ctqmc_make_htau(dn, tau_start - tau_end,atom)
             endif
         enddo ! over i={1,kaux} loop
     enddo ! over j={1,kaux} loop

! calculate Tmm: Tmm = Dmm . Mmat
     call dgemm('N', 'N', kaux, kaux, kaux, one, Dmm, kaux, mmat(1:kaux, 1:kaux, up,atom), kaux, zero, Tmm, kaux)

! calculate the determinant of Tmm, it is the desired ratio
     call ctqmc_dmat_det(kaux, Tmm, ratio)

! deallocate memory
     deallocate(Dmm)
     deallocate(Tmm)

     return
  end subroutine cat_reflip_detrat
