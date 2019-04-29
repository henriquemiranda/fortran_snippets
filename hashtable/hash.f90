module m_hashtable
  implicit none

  type :: hashbucket
    integer :: nitems
    integer,allocatable :: items(:,:)
  end type hashbucket

  type :: hashtable_t
    integer :: bucketsize
    integer :: bucketstep
    integer :: nbuckets
    type(hashbucket),allocatable :: buckets(:)
  end type hashtable_t

  contains
  type(hashtable_t) function hashtable_init(nbuckets,bucketsize,bucketstep) result(new)
    ! Create the hashtable structure
    integer,intent(in) :: bucketsize, bucketstep, nbuckets
    integer :: ibucket

    new%nbuckets = nbuckets
    new%bucketsize = bucketsize
    new%bucketstep = bucketstep
    allocate(new%buckets(nbuckets))
    do ibucket=1,nbuckets
      allocate(new%buckets(ibucket)%items(2,new%bucketsize))
      new%buckets(ibucket)%nitems = 0
    end do
  end function hashtable_init

  subroutine hashtable_add(self, key, val)
    ! Add a key to the hashtable structure
    type(hashtable_t) :: self
    integer,intent(in) :: key, val
    integer :: ihash, item, nitems
    integer,allocatable :: new_items(:,:)
    ! Compute bucket for this element
    ihash = mod(key,self%nbuckets)+1
    nitems = self%buckets(ihash)%nitems
    ! Check if the element already exists in the bucket
    do item=1,nitems
      if (self%buckets(ihash)%items(1,item) /= key) cycle
      ! Replace value
      self%buckets(ihash)%items(2,item) = val
      return
    end do
    ! Add the element to the bucket
    ! Check if the buckets are full
    if (size(self%buckets(ihash)%items,2)==nitems) then
      allocate(new_items(2,nitems+self%bucketstep))
      new_items(:,:nitems) = self%buckets(ihash)%items(:,:nitems)
      new_items(:,nitems+1:) = 0
      call move_alloc(new_items,self%buckets(ihash)%items)
    end if
    nitems = nitems + 1
    self%buckets(ihash)%items(:,nitems) = [key,val]
    self%buckets(ihash)%nitems = nitems
  end subroutine hashtable_add

  subroutine hashtable_print(self)
    ! Print the hashtable data
    type(hashtable_t) :: self
    integer :: ibucket,item
    do ibucket=1,self%nbuckets
      do item=1,self%buckets(ibucket)%nitems
        write(*,*) ibucket, self%buckets(ibucket)%items(:,item)
      end do
    end do
  end subroutine hashtable_print

  subroutine hashtable_get(self,key,val,ierr)
    ! Get the value of a key in the hashtable
    type(hashtable_t) :: self
    integer,intent(in)  :: key
    integer,intent(out) :: val,ierr
    integer :: item, ihash
    ierr = 0
    ihash = mod(key,self%nbuckets)+1
    do item=1,self%buckets(ihash)%nitems
      if (self%buckets(ihash)%items(1,item) /= key) cycle
      val = self%buckets(ihash)%items(2,item)
      return
    end do
    ierr = 1
  end subroutine hashtable_get

  subroutine hashtable_cleanup(self)
    ! Free up memory in all the buckets
    ! this should be done only after all the add operations are finished
    type(hashtable_t) ::  self
    integer :: ibucket
    do ibucket=1,self%nbuckets
    end do
  end subroutine hashtable_cleanup

  subroutine hashtable_free(self)
    type(hashtable_t) ::  self
    integer :: ibucket
    do ibucket=1,self%nbuckets
      deallocate(self%buckets(ibucket)%items)
    end do
    deallocate(self%buckets)
  end subroutine hashtable_free
end module m_hashtable

program hash_main
  use m_hashtable
  implicit none
  type(hashtable_t) :: hash
  integer :: val, ierr

  hash = hashtable_init(10,1,1)
  call hashtable_add(hash,1,3)
  call hashtable_add(hash,2,2)
  call hashtable_add(hash,3,2)
  call hashtable_add(hash,4,2)
  call hashtable_add(hash,5,2)
  call hashtable_add(hash,6,2)
  call hashtable_add(hash,7,2)
  call hashtable_add(hash,30,1235)
  call hashtable_add(hash,30,1236)
  call hashtable_add(hash,1000,123)
  call hashtable_add(hash,34,143)
  call hashtable_print(hash)
  call hashtable_get(hash,30,val,ierr)
  write(*,*) val, ierr
  call hashtable_free(hash)

end program hash_main
