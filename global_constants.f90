!
! L.Zhang: modified, introduced int64...

  module constants

     implicit none

!=========================================================================
!>>> integer constants: numerical precision                            <<<
!=========================================================================

! single precision
     integer, public, parameter :: sp    = kind(1.0)

! double precision
     integer, public, parameter :: dp    = kind(1.0d0)

! precise precision
      integer, public, parameter :: int64  = selected_int_kind(18)

!=========================================================================
!>>> integer constants: file unit handler                              <<<
!=========================================================================

! standard console output
     integer, public, parameter :: mystd = 6

! log file out
     integer, public, parameter :: myout = 99

! common file output
     integer, public, parameter :: mytmp = 100

!=========================================================================
!>>> real constants: numerical constants                               <<<
!=========================================================================

! well-known $\pi$
     real(dp), public, parameter :: pi   = 3.141592653589793238462643383279_dp

! 0.0 in double precision form
     real(dp), public, parameter :: zero = 0.0_dp

! 1.0 in double precision form
     real(dp), public, parameter :: one  = 1.0_dp

! 2.0 in double precision form
     real(dp), public, parameter :: two  = 2.0_dp

! 0.5 in double precision form
     real(dp), public, parameter :: half = 0.5_dp

! $\epsilon$ in double precision form
     real(dp), public, parameter :: eps6 = 1.0E-6
     real(dp), public, parameter :: eps8 = 1.0E-8
     real(dp), public, parameter :: epst = 1.0E-10
     real(dp), public, parameter :: epss = 1.0E-12

!=========================================================================
!>>> real constants: physical constants                                <<<
!=========================================================================

! conversion factor from eV to kelvin
     real(dp), public, parameter :: ev2k = 11604.505008098_dp

!=========================================================================
!>>> complex constants: numerical constants                            <<<
!=========================================================================

! complex unit, i
     complex(dp), public, parameter :: czi   = dcmplx(0.0_dp, 1.0_dp)

! complex unit, 1
     complex(dp), public, parameter :: cone  = dcmplx(1.0_dp, 0.0_dp)

! complex unit, 0
     complex(dp), public, parameter :: czero = dcmplx(0.0_dp, 0.0_dp)

  end module constants
