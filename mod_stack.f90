
! purpose : the purpose of this module is to define a stack-type (LIFO)
!           data structure in fortran version

  module stack
     implicit none

!-------------------------------------------------------------------------
!::: declare global parameters                                         :::
!-------------------------------------------------------------------------

! stack size limit, default value
     integer, private, parameter :: limit = 1024

! mystd: device descriptor, console output
     integer, private, parameter :: mystd = 6

!-------------------------------------------------------------------------
!::: declare global data structure                                     :::
!-------------------------------------------------------------------------

! define integer type stack
     type istack

! top position of stack
         integer :: top

! size of allocatable array
         integer :: nsize

! allocatable array, which is used to store elements in stack
         integer, allocatable :: item(:)

     end type istack

!-------------------------------------------------------------------------
!::: declare accessibility for module routines                         :::
!-------------------------------------------------------------------------

     public :: istack_create
     public :: istack_clean
     public :: istack_destroy

     public :: istack_copyer
     public :: istack_setter
     public :: istack_getter

     public :: istack_push
     public :: istack_pop

     public :: istack_display

     public :: istack_gettop
     public :: istack_getrest
     public :: istack_getsize

     public :: istack_isfull
     public :: istack_isempty

  contains ! encapsulated functionality

!>>> create and initialize a integer type stack
  type (istack) &
  function istack_create(n) result (s)
     implicit none

! external arguments
! size of stack
     integer, optional, intent(in) :: n

! local variables
! status flag
     integer :: istat

! determine the capacity of stack
     if ( present (n) ) then
         s%nsize = n
     else
         s%nsize = limit
     endif

! setup the top position
     s%top = 0

! allocate memory for item array
     allocate(s%item(s%nsize), stat=istat)

     return
  end function istack_create

!>>> reset the integer type stack, clean all its elements
  subroutine istack_clean(s)
     implicit none

! external arguments
! integer type stack
     type (istack), intent(inout) :: s

! reset top position
     s%top = 0

     return
  end subroutine istack_clean

!>>> destroy and finalize a integer type stack
  subroutine istack_destroy(s)
     implicit none

! external arguments
! integer type stack
     type (istack), intent(inout) :: s

! deallocate memory
     if (allocated(s%item)) deallocate(s%item)

! reset top position
     s%top = 0

     return
  end subroutine istack_destroy

!>>> copy an istack object to another
  subroutine istack_copyer(sa, sb)
     implicit none

! external arguments
! integer type stack, input
     type (istack), intent(in) :: sa

! integer type stack, output
     type (istack), intent(inout) :: sb

! check nsize at first
     if ( sa%nsize /= sb%nsize ) then
         write(mystd,'(a)') 'istack: the sizes of two stacks are not equal'
         STOP
     endif

! sync the data
     sb%top = sa%top
     sb%item = sa%item

     return
  end subroutine istack_copyer

!>>> update the item's value of istack at special position
  subroutine istack_setter(s, item, pos)
     implicit none

! external arguments
! integer type stack
     type (istack), intent(inout) :: s

! elements to be setted
     integer, intent(in) :: item

! position of the element to be updated
     integer, intent(in) :: pos

     if ( pos < 1 .or. pos > s%nsize ) then
         write(mystd,'(a)') 'istack: the position is not correct'
         STOP
     else
         s%item(pos) = item
     endif

     return
  end subroutine istack_setter

!>>> return the item's value of istack at special position
  integer &
  function istack_getter(s, pos) result (item)
     implicit none

! external arguments
! integer type stack
     type (istack), intent(in) :: s

! position of the element
     integer, intent(in) :: pos

     if ( pos < 1 .or. pos > s%nsize ) then
         write(mystd,'(a)') 'istack: the position is not correct'
         STOP
     else
         item = s%item(pos)
     endif

     return
  end function istack_getter

!>>> push item on top of stack
  subroutine istack_push(s, item)
     implicit none

! external arguments
! integer type stack
     type (istack), intent(inout) :: s

! elements to be pushed in the stack
     integer, intent(in) :: item

     if ( s%top == s%nsize ) then
         write(mystd,'(a)') 'istack: the stack is full, can not push any item on it'
         STOP
     else
         s%top = s%top + 1
         s%item(s%top) = item
     endif

     return
  end subroutine istack_push

!>>> pop off item from the top of stack
  integer &
  function istack_pop(s) result (item)
     implicit none

! external arguments
! integer type stack
     type (istack), intent(inout) :: s

     if ( s%top == 0 ) then
         write(mystd,'(a)') 'istack: the stack is empty, can not pop off any item from it'
         STOP
     else
         item = s%item(s%top)
         s%top = s%top - 1
     endif

     return
  end function istack_pop

!>>> display the top item in the stack without pop it off
  integer &
  function istack_display(s) result (item)
     implicit none

! external arguments
! integer type stack
     type (istack), intent(in) :: s

     if ( s%top == 0 ) then
         write(mystd,'(a)') 'istack: the stack is empty, can not return the top item of it'
         STOP
     else
         item = s%item(s%top)
     endif

     return
  end function istack_display

!>>> return the top position of the stack, i.e, the number of items stored
! in the stack currently
  integer &
  function istack_gettop(s) result (t)
     implicit none

! external arguments
! integer type stack
     type (istack), intent(in) :: s

     t = s%top

     return
  end function istack_gettop

!>>> return the number of empty sites of the stack
  integer &
  function istack_getrest(s) result (r)
     implicit none

! external arguments
! integer type stack
     type (istack), intent(in) :: s

     r = s%nsize - s%top

     return
  end function istack_getrest

!>>> return the actual capacity of the stack
  integer &
  function istack_getsize(s) result (n)
     implicit none

! external arguments
! integer type stack
     type (istack), intent(in) :: s

     n = s%nsize

     return
  end function istack_getsize

!>>> check whether the stack is empty
  logical &
  function istack_isempty(s) result (b)
     implicit none

! external arguments
! integer type stack
     type (istack), intent(in) :: s

     b = ( s%top == 0 )

     return
  end function istack_isempty

!>>> check whether the stack is full of items
  logical &
  function istack_isfull(s) result (b)
     implicit none

! external arguments
! integer type stack
     type (istack), intent(in) :: s

     b = ( s%top == s%nsize )

     return
  end function istack_isfull

  end module stack
