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

  pure function compute_hash(self,key) result(ihash)
    type(hashtable_t),intent(in) :: self
    integer,intent(in) :: key
    integer :: ihash
    ihash = mod(key-1,self%nbuckets)+1
  end function compute_hash

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
    ! Check if the buckets are full
    if (size(self%buckets(ihash)%items,2)==nitems) then
      allocate(new_items(2,nitems+self%bucketstep))
      new_items(:,:nitems) = self%buckets(ihash)%items(:,:nitems)
      new_items(:,nitems+1:) = 0
      call move_alloc(new_items,self%buckets(ihash)%items)
    end if
    ! Add the element to the bucket
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

  subroutine hashtable_keys(self,keys)
    ! Get an array with all the keys in hashtable
    type(hashtable_t) ::  self
    integer :: item,ibucket,idx
    integer,allocatable :: keys(:)
    allocate(keys(sum(self%buckets(:)%nitems)))
    idx = 0
    do ibucket=1,self%nbuckets
      do item=1,self%buckets(ibucket)%nitems
        idx = idx + 1
        keys(idx) = self%buckets(ibucket)%items(1,item)
      end do
    end do
  end subroutine hashtable_keys

  subroutine hashtable_items(self,items)
    ! Get an array with all the keys in hashtable
    type(hashtable_t) ::  self
    integer :: item,ibucket,idx
    integer,allocatable :: items(:,:)
    allocate(items(2,sum(self%buckets(:)%nitems)))
    idx = 0
    do ibucket=1,self%nbuckets
      do item=1,self%buckets(ibucket)%nitems
        idx = idx + 1
        items(:,idx) = self%buckets(ibucket)%items(:,item)
      end do
    end do
  end subroutine hashtable_items

  subroutine hashtable_cleanup(self)
    ! Free up memory in all the buckets
    ! this should be done only after all the add operations are finished
    type(hashtable_t) ::  self
    integer :: ibucket,nitems
    integer,allocatable :: new_items(:,:)
    do ibucket=1,self%nbuckets
      ! get size of buckets and free unecessary memory
      nitems = self%buckets(ibucket)%nitems
      allocate(new_items(2,nitems))
      new_items = self%buckets(ibucket)%items(:,:nitems)
      call move_alloc(new_items,self%buckets(ibucket)%items)
    end do
  end subroutine hashtable_cleanup

  integer function hashtable_size(self) result(bsize)
    type(hashtable_t) ::  self
    integer :: ibucket
    ! Return the size of the hashtable in bytes
    bsize = storage_size(self%buckets)*self%nbuckets
    do ibucket=1,self%nbuckets
      bsize = bsize + storage_size(self%buckets(ibucket)%items)*size(self%buckets(ibucket)%items,2)*2
    end do
    bsize = bsize/8
  end function hashtable_size

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
  integer :: val, ierr, ii
  integer,allocatable :: keys(:), items(:,:)

  hash = hashtable_init(10,1,1)
  call hashtable_add(hash,1,3)
  call hashtable_add(hash,2,2)
  call hashtable_add(hash,7,2)
  call hashtable_add(hash,30,1235)
  call hashtable_add(hash,30,1236)
  call hashtable_add(hash,1000,123)
  call hashtable_add(hash,34,143)
  call hashtable_cleanup(hash)
  call hashtable_print(hash)
  call hashtable_get(hash,30,val,ierr)
  write(*,*) val, ierr

  call hashtable_keys(hash,keys)
  write(*,*) keys
  deallocate(keys)

  call hashtable_items(hash,items)
  do ii=1,size(items,2)
    write(*,'(i5,a,i5)') items(1,ii),'->',items(2,ii)
  end do
  deallocate(items)

  call hashtable_free(hash)

end program hash_main
