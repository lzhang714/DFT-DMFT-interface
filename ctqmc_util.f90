!
! MK: modified, at least to account for multiple atoms, potentially also in other ways.
! Original header:
!
!-------------------------------------------------------------------------
! project : azalea
! program : ctqmc_dmat_inv
!           ctqmc_zmat_inv
!           ctqmc_dmat_det
!           ctqmc_zmat_det
!           ctqmc_make_uumat
!           ctqmc_make_state
!           ctqmc_time_sorter
!           ctqmc_time_qsorter
!           ctqmc_time_qscorer
!           ctqmc_time_analyzer
! source  : ctqmc_util.f90
! type    : functions & subroutines
! author  : li huang (email:huangli712@yahoo.com.cn)
!         : michael karolak...
! history : 10/01/2008 by li huang
!           02/08/2009 by li huang
!           09/23/2009 by li huang
!           09/26/2009 by li huang
!           11/17/2009 by li huang
!           11/21/2009 by li huang
!           12/18/2009 by li huang
!           12/22/2009 by li huang
!           12/29/2009 by li huang
!           01/12/2010 by li huang
!           02/27/2010 by li huang
!           06/08/2010 by li huang
!           06/22/2010 by li huang
! purpose : to provide utility functions and subroutines for hybridization
!           expansion version continuous time quantum Monte Carlo (CTQMC)
!           quantum impurity solver
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

# define prefix '>>> used time: '

# define iolst1 prefix , mday, ' d ', mhou, ' h ', mmin, ' m in this iteration.'
# define iolst2 prefix , mhou, ' h ', mmin, ' m in this iteration.'
# define iolst3 prefix , mmin, ' m ', msec, ' s in this iteration.'
# define iolst4 prefix , msec, ' s in this iteration.'

# define iolst5 prefix , nday, ' d ', nhou, ' h ', nmin, ' m in total iteration.'
# define iolst6 prefix , nhou, ' h ', nmin, ' m in total iteration.'
# define iolst7 prefix , nmin, ' m ', nsec, ' s in total iteration.'
# define iolst8 prefix , nsec, ' s in total iteration.'

!>>> invert real(dp) matrix using lapack subroutines
  subroutine ctqmc_dmat_inv(ndim, dmat)
     use constants, only : dp

     implicit none

! external arguments
! dimension of dmat matrix
     integer, intent(in) :: ndim

! object matrix, on entry, it contains the original matrix, on exit,
! it is destroyed and replaced with the inversed matrix
     real(dp), intent(inout) :: dmat(ndim,ndim)

! local variables
! error flag
     integer  :: ierror

! working arrays for lapack subroutines
     integer  :: ipiv(ndim)
     real(dp) :: work(ndim)

! computes the LU factorization of a general m-by-n matrix, need lapack
! package, dgetrf subroutine
     call dgetrf(ndim, ndim, dmat, ndim, ipiv, ierror)
     if ( ierror /= 0 ) then
         call ctqmc_print_error('ctqmc_dmat_inv','error in lapack subroutine dgetrf')
     endif

! computes the inverse of an LU-factored general matrix, need lapack
! package, dgetri subroutine
     call dgetri(ndim, dmat, ndim, ipiv, work, ndim, ierror)
     if ( ierror /= 0 ) then
         call ctqmc_print_error('ctqmc_dmat_inv','error in lapack subroutine dgetri')
     endif

     return
  end subroutine ctqmc_dmat_inv

!>>> invert complex(dp) matrix using lapack subroutines
  subroutine ctqmc_zmat_inv(ndim, zmat)
     use constants, only : dp

     implicit none

! external arguments
! dimension of zmat matrix
     integer, intent(in) :: ndim

! object matrix, on entry, it contains the original matrix, on exit,
! it is destroyed and replaced with the inversed matrix
     complex(dp), intent(inout) :: zmat(ndim,ndim)

! local variables
! error flag
     integer     :: ierror

! working arrays for lapack subroutines
     integer     :: ipiv(ndim)
     complex(dp) :: work(ndim)

! computes the LU factorization of a general m-by-n matrix, need lapack
! package, zgetrf subroutine
     call zgetrf(ndim, ndim, zmat, ndim, ipiv, ierror)
     if ( ierror /= 0 ) then
         call ctqmc_print_error('ctqmc_zmat_inv','error in lapack subroutine zgetrf')
     endif

! computes the inverse of an LU-factored general matrix, need lapack
! package, zgetri subroutine
     call zgetri(ndim, zmat, ndim, ipiv, work, ndim, ierror)
     if ( ierror /= 0 ) then
         call ctqmc_print_error('ctqmc_zmat_inv','error in lapack subroutine zgetri')
     endif

     return
  end subroutine ctqmc_zmat_inv

!>>> calculate the determinant of a real(dp) matrix
  subroutine ctqmc_dmat_det(ndim, dmat, ddet)
     use constants, only : dp, one, cone

     implicit none

! external arguments
! dimension of dmat matrix
     integer, intent(in) :: ndim

! determinant of dmat matrix
     real(dp), intent(out) :: ddet

! object matrix, on entry, it contains the original matrix, on exit,
! it is destroyed and replaced with the L and U matrix
     real(dp), intent(inout) :: dmat(ndim,ndim)

! local variables
! loop index
     integer  :: i

! error flag
     integer  :: ierror

! size of working array work
     integer  :: lwork

! working arrays for lapack subroutines: dgetrf
     integer  :: ipiv(ndim)

! working arrays for lapack subroutines: dgeev
     real(dp) :: work(4*ndim)

! real and imaginary parts of the computed eigenvalues
     real(dp) :: wi(ndim)
     real(dp) :: wr(ndim)

! left and right eigenvectors
     real(dp) :: vl(ndim,ndim)
     real(dp) :: vr(ndim,ndim)

! dummy arrays, used to save dmat
     real(dp) :: amat(ndim,ndim)

! used to calculate determinant
     complex(dp) :: cres

! setup lwork
     lwork = 4 * ndim

! copy dmat to amat at first
     amat = dmat

!-------------------------------------------------------------------------
! method A:
!-------------------------------------------------------------------------
! computes the LU factorization of a general m-by-n matrix, need lapack
! package, dgetrf subroutine
     call dgetrf(ndim, ndim, dmat, ndim, ipiv, ierror)
     if ( ierror /= 0 ) then
         call ctqmc_print_exception('ctqmc_dmat_det','error in lapack subroutine dgetrf')
     endif

! calculate determinant
     ddet = one
     do i=1,ndim
         if ( ipiv(i) == i ) then
             ddet = ddet * ( +dmat(i,i) )
         else
             ddet = ddet * ( -dmat(i,i) )
         endif
     enddo ! over i={1,ndim} loop

! everything is ok!
     if ( ierror == 0 ) return

!-------------------------------------------------------------------------
! method B: as a backup
!-------------------------------------------------------------------------
! diagonalize amat to obtain its eigenvalues: wr and wi
     call dgeev('N', 'N', ndim, amat, ndim, wr, wi, vl, ndim, vr, ndim, work, lwork, ierror)
     if ( ierror /= 0 ) then
         call ctqmc_print_error('ctqmc_dmat_det','error in lapack subroutine dgeev')
     endif

! evaluate the final determinant
     cres = cone
     do i=1,ndim
         cres = cres * dcmplx( wr(i), wi(i) )
     enddo ! over i={1,ndim} loop
     ddet = cres

     return
  end subroutine ctqmc_dmat_det

!>>> calculate the determinant of a complex(dp) matrix
  subroutine ctqmc_zmat_det(ndim, zmat, zdet)
     use constants, only : dp, cone

     implicit none

! external arguments
! dimension of zmat matrix
     integer, intent(in) :: ndim

! determinant of zmat matrix
     real(dp), intent(out) :: zdet

! object matrix, on entry, it contains the original matrix, on exit,
! it is destroyed and replaced with the L and U matrix
     real(dp), intent(inout) :: zmat(ndim,ndim)

! local variables
! loop index
     integer :: i

! error flag
     integer :: ierror

! working arrays for lapack subroutines
     integer :: ipiv(ndim)

! computes the LU factorization of a general m-by-n matrix, need lapack
! package, zgetrf subroutine
     call zgetrf(ndim, ndim, zmat, ndim, ipiv, ierror)
     if ( ierror /= 0 ) then
         call ctqmc_print_error('ctqmc_zmat_det','error in lapack subroutine zgetrf')
     endif

! calculate determinant
     zdet = cone
     do i=1,ndim
         if ( ipiv(i) == i ) then
             zdet = zdet * ( +zmat(i,i) )
         else
             zdet = zdet * ( -zmat(i,i) )
         endif
     enddo ! over i={1,ndim} loop

     return
  end subroutine ctqmc_zmat_det

!
! MK: Additional Stuff for U matrix production
!     Adapted mostly from A. Poteryaevs code v 1.07
!

subroutine umatrix( umtrx, l, F )
!
!   Make <ij|V|kl> matrix of Coulomb interaction
!   (eq. (10) from JPCM 9, 767 (1997))
!   ATTENTION: The output matrix is in basis of complex spherical harmonics
!
  use constants
  
  implicit none 
  
  integer   , intent(in   ) :: l
  real(dp), intent(in   ) :: F(0:6)
  real(dp), intent(  out) :: umtrx(2*l+1,2*l+1,2*l+1,2*l+1)
  
  integer :: m1, m2, m3, m4, k

  umtrx = 0
  
  do m1 = -l, l
    do m2 = -l, l
      do m3 = -l, l
        do m4 = -l, l
          do k = 0, 2*l+1, 2
            umtrx(m1+l+1,m2+l+1,m3+l+1,m4+l+1) = umtrx(m1+l+1,m2+l+1,m3+l+1,m4+l+1) +   &
                                         F(k) * a( k, l, m1, m2, m3, m4 )
          end do
        end do
      end do
    end do
  end do

  contains

  real(dp) function Factorial( n )
    implicit none
    integer, intent(in   ) :: n
    
    integer :: l
    
    l = n

    Factorial = 1
    do while ( l > 0 )
      Factorial = Factorial * l
      l = l - 1
    end do
    
  end function Factorial

  real(dp) function ClebschGordan( j1,m1, j2,m2, j,m )
!
!  Return Clebsch-Gordan coefficient <j1j2m1m2|j1j2jm>
!  Last definition from Wolfram site
!  http://functions.wolfram.com/HypergeometricFunctions/ClebschGordan/06/01/
!
    integer, intent(in   ) :: j1, m1           !   Angular momentum and its projection
    integer, intent(in   ) :: j2, m2
    integer, intent(in   ) :: j,  m
  
    integer :: k, si
    real(dp) :: coef, a 
  
    ClebschGordan = 0
    
    if( abs(j1+j2) < j .or. j < abs(j1-j2) )  return          !  Triangle rule
    if( abs(m1+m2-m) > 0 )                    return          !  Selection rule
    if( mod(j1+j2+j,1) > 1.e-14 )             return          !  Implemented for integers only
    if( j1 < 0 .or. j2 < 0 )                  return 
    if( abs(m1) > j1 .or. abs(m2) > j2 )      return 
    
    coef = (2*j+1)*Factorial(1+j+j1+j2)/Factorial(j1+j2-j)/Factorial(j1-j2+j)/Factorial(-j1+j2+j)
    coef = coef * Factorial(j+m)*Factorial(j-m)*Factorial(j1+m1)*Factorial(j1-m1)
    coef = coef / Factorial(j2+m2) / Factorial(j2-m2)
    coef = sqrt(coef)

    if( mod(j1-m1,2) == 1 ) coef = -coef

    a  = 0
    si = 1
    do k = 0, min(j-m,j1-m1)
      a = a + si*Factorial(j1+j2-m-k)*Factorial(j+j2-m1-k) /        &
                 Factorial(k) / Factorial(j1-m1-k) /                &
                 Factorial(j-m-k)/Factorial(1+j+j1+j2-k)
      si = -si
    end do
    
    ClebschGordan = coef * a
    
  end function ClebschGordan

  real(dp) function a( k, l, m1, m2, m3, m4 )
    integer, intent(in   ) :: k, l
    integer, intent(in   ) :: m1, m2, m3, m4

    integer :: q

    a = 0
    do q = -k, k
      a = a + ClebschGordan( l,m3, k,q, l,m1 )*ClebschGordan( l,m2, k,q, l,m4 )
    end do
    a = a * ClebschGordan( l,0, k,0, l,0 )*ClebschGordan( l,0, k,0, l,0 )

  end function a

end subroutine umatrix

!
! \MK
!

!
!>>> to build general U interaction matrix: uumat
!
  subroutine ctqmc_make_uumat(full_u_mat,uumat,controltag,cfb_controltag)
     use constants
     use control
     use plo_quantities, only : natom, cfb_coeff
     use mmpi

     implicit none

! external arguments
! Coulomb interaction matrix
     real(dp), intent(out) :: uumat(norbs, norbs,natom)
     real(dp), intent(out) :: full_u_mat(nband,nband,nband,nband,natom)
     complex(dp) :: uu(nband,nband,nband,nband)
     complex(dp) :: ut(7,7)
     real(dp) :: u(nband,nband)
     integer, intent(in) :: controltag
     integer, intent(in) :: cfb_controltag

! local variables
! loop index
     integer  :: i
     integer  :: j
     integer  :: k
     integer  :: l
     integer  :: m1,m2,m3,m4
     
     integer :: ios1
     integer :: iat
     integer :: idum

     real(dp) :: sqrt2 = 1 / sqrt(2.0)

     character(20) :: filename

     logical  :: exists_umtrx


    select case(controltag)
!
!   simple Kanamori type matrix
!
    case(0)

       do iat=1,natom

! dummy array
    uu(:,:,:,:) = czero

          do i = 1, nband
           do j = 1, nband
             uu(i,j,i,j) = Uc(iat) - 2*Jz(iat)
             uu(i,j,j,i) = Jz(iat)
           end do
             uu(i,i,i,i) = Uc(iat)
          end do

   full_u_mat(:,:,:,:,iat) = uu(:,:,:,:)

!
!    Density-density part of matrix
!
        do i = 1, nband
         do j = 1, nband
           uumat(i,  j+nband,iat) = real(uu(i,j,i,j))
           uumat(i+nband,j,iat)   = real(uu(i,j,i,j))
           uumat(i,j,iat)     = real(uu(i,j,i,j)) -   &
                                            real(uu(i,j,j,i))
           uumat(i+nband,j+nband,iat) = real(uu(i,j,i,j)) -   &
                                            real(uu(i,j,j,i))
         end do
        end do

       enddo ! atom
!
! read dens-dens matrix from file as default.
!
    case(1)

   do iat=1,natom

    write(filename,*) iat
    filename = 'solver.uumat_'//trim(adjustl(filename))

    inquire (file = filename, exist = exists_umtrx)

    if(exists_umtrx .eqv. .false.) then
	stop 'solver.uumat_* file not found!'
    else

    open( 22,file=filename,form='formatted',status='old',iostat=ios1,&
              action='read',position='rewind' )

    end if

    do i=1,norbs
       read(22,*)(uumat(i,j,iat),j=1,norbs)
    enddo

    close(22)

    enddo
!
!   Compute the matrix from F0,F2,F4,F6
!   Dump it into solver.full_u_* then
!
    case(2)

! Check if the user is smart ;-)

    if(nband < 3) then
	stop 'For 2 bands Kanamori is exact, use case 0'
    end if

   do iat=1,natom

! dummy array
    uu(:,:,:,:) = czero

!     Define Slater integrals via U and J.
          Fslater(0,iat) = Uc(iat)

!     depending on shell, the slater integrals are different
!     d or f shells only
          select case (nband)
           case (3)
            Fslater(2,iat) = 5 * Jz(iat)
           case (5)
            Fslater(2,iat) = 14 * Jz(iat) / 1.63
            Fslater(4,iat) = 0.63 * Fslater(2,iat)
           case (7)
            Fslater(2,iat) = 6435 * Jz(iat)  /   &
                            ( 286 + 195 * 451 / 675 + 250 * 1001 / 2025 )
            Fslater(4,iat) = 451  * Fslater(2,iat) / 675
            Fslater(6,iat) = 1001 * Fslater(2,iat) / 2025
          end select

	l = ( nband - 1 ) / 2
!
!   compute matrix in basis of spherical harmonics
    call umatrix( full_u_mat(:,:,:,:,iat), l, Fslater(:,iat) )

!
!   transform to cubic harmonics (VASP order assumed)

  ut = czero

!
!   Define basis transformation
    select case (nband)
      case(3)      !  y  z  x
        ut(1,1) = sqrt2*czi
        ut(1,3) = sqrt2*czi
        ut(2,2) = 1
        ut(3,1) = sqrt2
        ut(3,3) =-sqrt2
      case(5)      !  xy  yz  z2  xz  x2-y2
        ut(1,1) = sqrt2*czi
        ut(1,5) =-sqrt2*czi
        ut(2,2) = sqrt2*czi
        ut(2,4) = sqrt2*czi
        ut(3,3) = 1
        ut(4,2) = sqrt2
        ut(4,4) =-sqrt2
        ut(5,1) = sqrt2
        ut(5,5) = sqrt2
      case(7)
        ut(1,1) = sqrt2*czi
        ut(1,7) = sqrt2*czi
        ut(2,2) = sqrt2*czi
        ut(2,6) =-sqrt2*czi
        ut(3,3) = sqrt2*czi
        ut(3,5) = sqrt2*czi
        ut(4,4) = 1
        ut(5,3) = sqrt2
        ut(5,5) =-sqrt2
        ut(6,2) = sqrt2
        ut(6,6) = sqrt2
        ut(7,1) = sqrt2
        ut(7,7) =-sqrt2
    end select

!
! Apply the transformation

          do i = 1, nband
           do j = 1, nband
            do k = 1, nband
             do l = 1, nband
               do m1 = 1, nband
                do m2 = 1, nband
                 do m3 = 1, nband
                  do m4 = 1, nband
                    uu(i,j,k,l) = uu(i,j,k,l) + conjg(ut(i,m1))*conjg(ut(j,m2)) *  &
                                           full_u_mat(m1,m2,m3,m4,iat) *  &
                                            ut(k,m3) * ut(l,m4)
                  end do
                 end do
                end do
               end do
             end do
            end do
           end do
          end do
!
!    Store it

    full_u_mat(:,:,:,:,iat) = uu(:,:,:,:)

!
!    Density-density part of matrix
        do i = 1, nband
         do j = 1, nband
           uumat(i,  j+nband,iat) = real(uu(i,j,i,j))
           uumat(i+nband,j,iat)   = real(uu(i,j,i,j))
           uumat(i,j,iat)     = real(uu(i,j,i,j)) -   &
                                            real(uu(i,j,j,i))
           uumat(i+nband,j+nband,iat) = real(uu(i,j,i,j)) -   &
                                            real(uu(i,j,j,i))
         end do
        end do

!
! Write the full matrix into a file

    write(filename,*) iat
    filename = 'solver.full_u_'//trim(adjustl(filename))
    open(mytmp, file=filename, form='formatted', status='unknown')

    do i = 1, nband
     do j = 1, nband
      do k = 1, nband
       do l = 1, nband
          write(mytmp,*) i, j, k, l, real(uu(i,j,k,l))
       enddo
      enddo
     enddo
    enddo

    close(mytmp)

    end do ! atom

    case(3)
!
! Read the full U matrix from a file. no rotation.

    full_u_mat = zero

   do iat=1,natom

    write(filename,*) iat
    filename = 'solver.full_u_'//trim(adjustl(filename))

    inquire (file = filename, exist = exists_umtrx)

    if(exists_umtrx .eqv. .false.) then

	stop 'solver.full_u_* file not found!'

    else

    open( 88,file=filename,form='formatted',status='old',iostat=ios1,&
              action='read',position='rewind' )

          do i = 1, nband
           do j = 1, nband
            do k = 1, nband
             do l = 1, nband
              read(88,*) idum,idum,idum,idum,full_u_mat(i,j,k,l,iat)
             enddo
           enddo
          enddo
         enddo

    close(88)

    end if
    
!
!   dummy array

    uu(:,:,:,:) = full_u_mat(:,:,:,:,iat)

!
!    Density-density part of matrix
        do i = 1, nband
         do j = 1, nband
           uumat(i,  j+nband,iat) = real(uu(i,j,i,j))
           uumat(i+nband,j,iat)   = real(uu(i,j,i,j))
           uumat(i,j,iat)     = real(uu(i,j,i,j)) -   &
                                            real(uu(i,j,j,i))
           uumat(i+nband,j+nband,iat) = real(uu(i,j,i,j)) -   &
                                            real(uu(i,j,j,i))
         end do
        end do

    end do! atom

    end select
!
!   CFB rotation or not. Full U_ijkl must be present!
!
    select case(cfb_controltag)

    case(0)
!
! Write the density-density matrix for checking only
!
    if(myid==master)then
     do iat = 1,natom
     write(*,*)'Density-density part of U tensor for atom',iat,':'
     write(*,*)''
      do i=1,norbs
    	write(*,'(10f10.5)')(uumat(i,j,iat),j=1,norbs)
      enddo
    write(*,*)''
     end do! atom
     write(*,*)''
    end if! master

      return

    case(1,2)

    if(myid==master)then
    write(*,*)'Rotated U matrix:'
    end if

    full_u_mat = zero
!
! Read the Full U matrix from file
!
   do iat=1,natom

    write(filename,*) iat
    filename = 'solver.full_u_'//trim(adjustl(filename))
    
    inquire (file = filename, exist = exists_umtrx)
    
    if(exists_umtrx .eqv. .false.) then
    
	stop 'solver.full_u_* file not found!'
    
    else

    open( 88,file=filename,form='formatted',status='old',iostat=ios1,&
              action='read',position='rewind' )

          do i = 1, nband
           do j = 1, nband
            do k = 1, nband
             do l = 1, nband
              read(88,*) idum,idum,idum,idum,full_u_mat(i,j,k,l,iat)
             enddo
           enddo
          enddo
         enddo

    close(88)
    
    end if
!
! Rotate the u matrix according to cfb_coeff
    uu(:,:,:,:) = czero

    u(:,:) = cfb_coeff(iat,:,:)

    u = transpose( u )

          do i = 1, nband
           do j = 1, nband
            do k = 1, nband
             do l = 1, nband
               do m1 = 1, nband
                do m2 = 1, nband
                 do m3 = 1, nband
                  do m4 = 1, nband
                    uu(i,j,k,l) = uu(i,j,k,l) + u(i,m1)*u(j,m2) *  &
                                           full_u_mat(m1,m2,m3,m4,iat) *  &
                                            u(k,m3) * u(l,m4)
                  end do
                 end do
                end do
               end do
             end do
            end do
           end do
          end do

!
!    Density-density part of matrix
        do i = 1, nband
         do j = 1, nband
           uumat(i,  j+nband,iat) = real(uu(i,j,i,j))
           uumat(i+nband,j,iat)   = real(uu(i,j,i,j))
           uumat(i,j,iat)     = real(uu(i,j,i,j)) -   &
                                            real(uu(i,j,j,i))
           uumat(i+nband,j+nband,iat) = real(uu(i,j,i,j)) -   &
                                            real(uu(i,j,j,i))
         end do
        end do

!
!   store the rotated matrix
    full_u_mat(:,:,:,:,iat) = uu(:,:,:,:)

    end do! atom

    end select

!
! Write the density-density matrix for checking
    if(myid==master)then
     do iat = 1,natom
     write(*,*)'Density-density part of U tensor for atom',iat,':'
     write(*,*)''
      do i=1,norbs
    	write(*,'(10f10.5)')(uumat(i,j,iat),j=1,norbs)
      enddo
      write(*,*)''
     end do! atom
     write(*,*)''
    end if! master

     return
  end subroutine ctqmc_make_uumat

!>>> convert current atomic state array into a decimal number (state index)
  subroutine ctqmc_make_state(norbs, pstat, state)
     implicit none

! external arguments
! index of atomic state
     integer, intent(out) :: pstat

! number of orbitals
     integer, intent(in)  :: norbs

! atomic state array
     integer, intent(in)  :: state(norbs)

! local variables
! loop index
     integer :: i

! init pstat
     pstat = 1

! evaluate pstat, for example, 0101 = 0*2^0 + 1*2^1 + 0*2^2 + 1*2^3 = 10
     do i=1,norbs
         if ( state(i) > 0 ) pstat = pstat + ishft(1, i-1)
     enddo ! over i={1,norbs} loop

     return
  end subroutine ctqmc_make_state

!>>> using bubble sort algorithm to sort a real dataset, the slowest algorithm
  subroutine ctqmc_time_sorter(nsize, list)
     use constants, only : dp

     implicit none

! external arguments
! grab the number of values from the calling code
     integer, intent(in) :: nsize

! dataset to be sorted
     real(dp), intent(inout) :: list(nsize)

! local variables
! dataset index
     integer  :: i = 0
     integer  :: j = 0

! dummy variables
     real(dp) :: swap

! basically we just loop through every element to compare it against
! every other element
! this loop increments i which is our starting point for the comparison
     sort_loop1: do i=nsize,1,-1
! this loop increments j which is the ending point for the comparison
         sort_loop2: do j=1,i-1
! swap the two elements here
             exchange: if ( list(j) > list(j+1) ) then
                 swap = list(j)
                 list(j) = list(j+1)
                 list(j+1) = swap
             endif exchange
         enddo sort_loop2 ! over j={1,i-1} loop
     enddo sort_loop1 ! over i={nsize,1,-1} loop

     return
  end subroutine ctqmc_time_sorter

!>>> sets up for the quick sort recursive method
  subroutine ctqmc_time_qsorter(nsize, list)
     use constants, only : dp

     implicit none

! external arguments
! grab the number of values from the calling code
     integer, intent(in) :: nsize

! dataset to be sorted
     real(dp), intent(inout) :: list(nsize)

! kicks off the recursive process
     call ctqmc_time_qscorer(1, nsize, nsize, list)

     return
  end subroutine ctqmc_time_qsorter

!>>> this is the actually recursive portion of the quicksort algorithm
  recursive &
  subroutine ctqmc_time_qscorer(pstart, pend, nsize, list)
     use constants, only : dp

     implicit none

! external arguments
! start point
     integer, intent(in) :: pstart

! end point
     integer, intent(in) :: pend

! size of array
     integer, intent(in) :: nsize

! dataset to be sorted
     real(dp), intent(inout) :: list(nsize)

! local variables
! used to find out list(left) > kaux and list(right) < kaux
     integer  :: left, right

! used to record list(pstart)
     real(dp) :: kaux

! used to swap data
     real(dp) :: taux

! setup left and right
     left = pstart
     right = pend + 1

! only in right > left, the data is to be sorted
     if ( right > left ) then

! record list(pstart) at first
         kaux = list(pstart)

         do while ( .true. )

! find out where list(left) < kaux
             do while ( .true. )
                 left = left + 1
                 if ( list(left)  > kaux .or. left  >= pend   ) exit
             enddo ! over do while loop

! find out where list(right) > kaux
             do while ( .true. )
                 right = right - 1
                 if ( list(right) < kaux .or. right <= pstart ) exit
             enddo ! over do while loop

! we should ensure right is larger than left
             if ( right <= left ) exit

! exchange data between list(left) and list(right)
             taux = list(left)
             list(left) = list(right)
             list(right) = taux

         enddo ! over do while loop

! exchange data between list(pstart) and list(right)
        list(pstart) = list(right)
        list(right) = kaux

! sort data from pstart to right-1
        call ctqmc_time_qscorer(pstart, right-1, nsize, list)

! sort data from right+1 to pend
        call ctqmc_time_qscorer(right+1, pend, nsize, list)

     endif ! back if ( right > left ) block

     return
  end subroutine ctqmc_time_qscorer

!>>> used to print the iteration timing information about continuous time
! quantum Monte Carlo quantum impurity solver.
  subroutine ctqmc_time_analyzer(time_iter, time_niter)
     use constants

     implicit none

! external arguments
! time used in this iteration
     real(dp), intent(in) :: time_iter

! time used in total iteration
     real(dp), intent(in) :: time_niter

! local variables
! number of days
     integer  :: mday, nday

! number of hours
     integer  :: mhou, nhou

! number of minutes
     integer  :: mmin, nmin

! number of seconds
     real(dp) :: msec, nsec

     mday = time_iter / 86400
     msec = time_iter - 86400 * mday
     mhou = msec / 3600
     msec = msec - 3600 * mhou
     mmin = msec / 60
     msec = msec - 60 * mmin

     nday = time_niter / 86400
     nsec = time_niter - 86400 * nday
     nhou = nsec / 3600
     nsec = nsec - 3600 * nhou
     nmin = nsec / 60
     nsec = nsec - 60 * nmin

     if      ( mday > 0 ) then
         write(mystd,'(4X,3(a,i2),a)')     iolst1

     else if ( mhou > 0 ) then
         write(mystd,'(4X,2(a,i2),a)')     iolst2

     else if ( mmin > 0 ) then
         write(mystd,'(4X,a,i2,a,f5.2,a)') iolst3

     else
         write(mystd,'(4X,a,f5.2,a)')      iolst4

     endif ! back if ( mday > 0 ) block

     if      ( nday > 0 ) then
         write(mystd,'(4X,3(a,i2),a)')     iolst5

     else if ( nhou > 0 ) then
         write(mystd,'(4X,2(a,i2),a)')     iolst6

     else if ( nmin > 0 ) then
         write(mystd,'(4X,a,i2,a,f5.2,a)') iolst7

     else
         write(mystd,'(4X,a,f5.2,a)')      iolst8

     endif ! back if ( nday > 0 ) block

     return
  end subroutine ctqmc_time_analyzer
