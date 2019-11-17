module m_stretchylist

  implicit none
  integer,parameter :: dp=8
  integer,parameter :: LIST_SIZE=10
  integer,parameter :: LIST_STEP=10

  type :: stretchylist_t
    integer :: n
    integer,allocatable :: vals(:)
  end type stretchylist_t

  contains
  type(stretchylist_t) function stretchylist_init() result(new)
    new%n=0
    allocate(new%vals(LIST_SIZE))
  end function

  subroutine stretchylist_free(slist)
    type(stretchylist_t) :: slist
    deallocate(slist%vals)
  end subroutine

  subroutine stretchylist_append(slist,val)
    type(stretchylist_t) :: slist
    integer,intent(in) :: val
    integer,allocatable :: new_vals(:)
    slist%n = slist%n+1
    ! check if need to resize
    if (slist%n>size(slist%vals)) then
      allocate(new_vals(slist%n+LIST_STEP))
      new_vals(:size(slist%vals))=slist%vals
      call move_alloc(new_vals,slist%vals)
    end if
    ! add new element
    slist%vals(slist%n) = val
  end subroutine

end module m_stretchylist

program stretchylist_main

  use m_stretchylist
  implicit none

  type(stretchylist_t) :: slist
  integer :: ii

  slist = stretchylist_init()
  do ii=1,1000
    call stretchylist_append(slist,ii)
  end do
  write(*,'(10i6)') slist%vals(:slist%n)
  call stretchylist_free(slist)

end program stretchylist_main


