module m_octree

  implicit none

  integer,parameter :: POINTS_SIZE=10
  integer,parameter :: POINTS_STEP=10

  integer,parameter :: dp=8
  real(dp),parameter :: zero=0.0_dp,half=0.5_dp,one=1.0_dp

  real(dp) :: shifts(3,2,4)

  type :: octree_node_t
    integer :: id
    real(dp) :: hi(3)
    real(dp) :: lo(3)
    type(octree_node_t),pointer :: mother
    type(octree_node_t),pointer :: childs(:) ! if node this is allocated
    integer,allocatable :: ids(:)            ! if leaf this is allocated
  end type octree_node_t

  type :: octree_t
    integer :: max_npoints
    real(dp),pointer :: points(:,:)
    type(octree_node_t) :: first
  end type octree_t

!--------------------------------------------------------------------

  contains
  type(octree_t) function octree_init(points,max_npoints) result (new)
    integer,intent(in) :: max_npoints
    real(dp), target, intent(in) :: points(:,:)
    type(octree_node_t) :: mother
    real,parameter :: ieps = 0.1
    integer,allocatable :: ids(:)
    real(dp) :: hi(3), lo(3)
    integer :: ii, jj, kk, npoints, ioctant

    ! determine shifts
    do ii=0,1
      do jj=0,1
        do kk=0,1
          ioctant = ii*4+jj*2+kk+1
          shifts(:,1,ioctant) = [half*ii,half*jj,half*kk]
          shifts(:,2,ioctant) = shifts(:,1,ioctant) + [half,half,half]
        end do
      end do
    end do

    ! determine hi and lo
    do ii=1,3
      hi(ii) = maxval(points(ii,:))+ieps
      lo(ii) = minval(points(ii,:))-ieps
    end do
    hi = [1.0,1.0,1.0]
    lo = [0.0,0.0,0.0]

    !first octree contains all the points
    npoints = size(points,2)
    new%points => points
    new%max_npoints = max_npoints
    ids = [(ii,ii=1,npoints)]
    new%first = octree_node_build(new,mother,lo,hi,npoints,ids)

  end function octree_init

  type(octree_node_t) recursive function octree_node_build(octree,mother,lo,hi,nids,ids) result (new)
    type(octree_t),intent(in) :: octree
    type(octree_node_t),target,intent(in) :: mother
    integer,intent(in) :: nids
    integer,intent(in) :: ids(nids)
    integer :: octants(nids)
    integer :: id, counter, ioctant, ipoint
    integer :: new_ids(nids)
    real(dp) :: lo(3), hi(3), new_lo(3), new_hi(3)

    if (allocated(mother%ids)) then
      write(*,*) 'my mother is a leaf.'
    end if
    new%mother => mother
    new%lo = lo
    new%hi = hi
    ! check if this is a leaf node
    if (nids<octree%max_npoints) then
      allocate(new%ids(nids))
      new%ids = ids(:nids)
      return
    end if

    !determine the octants of each point
    call get_octants(lo,hi,nids,ids,octree%points,octants)

    allocate(new%childs(8))
    do ioctant=1,8
      ! How many points are in this octant?
      counter = 0
      do id=1,nids
        ipoint = ids(id)
        if (octants(id) /= ioctant) cycle
        counter = counter + 1
        new_ids(counter) = ipoint
      end do
      ! Build this octant
      call get_lo_hi(lo,hi,new_lo,new_hi,ioctant)
      new%childs(ioctant) = octree_node_build(octree,new,new_lo,new_hi,counter,new_ids)
    end do

  end function octree_node_build

  integer function octree_find(octree,point,dist) result(closest_id)
    ! Find the closest point in the box that contains it
    type(octree_t),target,intent(in) :: octree
    real(dp),intent(in) :: point(3)
    real(dp),intent(out) :: dist
    type(octree_node_t),pointer :: octn
    integer :: id, ipoint, ioctant, npoints
    real(dp) :: trial_dist

    dist = huge(dist)
    octn => octree%first
    do
      ! get octant of this point
      ioctant = get_octant_lohi(octn%lo,octn%hi,point)
      ! point to this node
      octn => octn%childs(ioctant)
      npoints = size(octn%ids)
      if (.not.allocated(octn%ids)) cycle
      ! if leaf node
      do id=1,npoints
        ipoint = octn%ids(id)
        trial_dist = dist_points(octree%points(:,ipoint),point)
        if (trial_dist > dist) cycle
        dist = trial_dist
        closest_id = ipoint
      end do
      return
    end do
  end function

  integer function octree_find_nearest(octree,point,dist) result(id)
    ! Find the ids of the points whose distance to point is smaller than max_dist
    ! counter is the number of points found so far
    type(octree_t),target,intent(in) :: octree
    real(dp),intent(in) :: point(3)
    real(dp),intent(inout) :: dist
    id = octn_find_nearest(octree,octree%first,point,dist)
  end function

  integer recursive function octn_find_nearest(octree,octn,point,min_dist) result(closest_id)
    ! find the nearest point by recursion
    type(octree_t),intent(in) :: octree
    type(octree_node_t),intent(in) :: octn
    real(dp),intent(in) :: point(3)
    real(dp),intent(inout) :: min_dist
    real(dp) :: dist
    integer :: id, ioctant, ipoint, trial_id
    closest_id = 0
    ! compute distance of point to this box (octant)
    dist = box_dist(octn%lo,octn%hi,point)
    ! if the distance is bigger than the closest point so far return
    if (dist>min_dist) return
    ! if this node is a leaf compare point by point
    if (allocated(octn%ids)) then
      do id=1,size(octn%ids)
        ipoint = octn%ids(id)
        dist = dist_points(octree%points(:,ipoint),point)
        if (dist>min_dist) cycle
        min_dist = dist
        closest_id = ipoint
      end do
      return
    end if
    ! if this is a node then find the nearest for all the childs
    do ioctant=1,8
      trial_id = octn_find_nearest(octree,octn%childs(ioctant),point,min_dist)
      if (trial_id==0) cycle
      closest_id = trial_id
    end do
  end function

  integer function octree_find_nearest_pbc(octree,point,dist,shift) result(id)
    ! Same as octree find but using periodic boundary conditions
    type(octree_t),target,intent(in) :: octree
    real(dp),intent(in) :: point(3)
    real(dp),intent(inout) :: dist
    real(dp),intent(out) :: shift(3)
    integer :: ii, jj, kk
    integer :: trial_id
    real(dp) :: hi(3), lo(3), de(3), trial_shift(3)
    real(dp) :: trial_dist
    ! find nearest in unshifted case
    id = octn_find_nearest(octree,octree%first,point,dist)
    ! find the shift that brings the point closer to the box
    lo = octree%first%lo
    hi = octree%first%hi
    de = hi-lo
    do ii=-1,1
      do jj=-1,1
        do kk=-1,1
          if (ii==0.and.jj==0.and.kk==0) cycle
          trial_shift = lo+[ii,jj,kk]*de
          trial_dist = box_dist(lo,hi,point+trial_shift)
          if (trial_dist>dist) cycle
          dist = trial_dist
          shift = trial_shift
        end do
      end do
    end do
    ! find the nearest point now
    trial_dist = dist
    trial_id = octn_find_nearest(octree,octree%first,point+shift,trial_dist)
    if (trial_dist>dist) return
    dist = trial_dist
    id = trial_id
  end function

  integer recursive function octn_free(octn) result(ierr)
    type(octree_node_t) :: octn
    integer :: ioctant
    ! if leaf deallocate ids
    if (allocated(octn%ids)) then
      deallocate(octn%ids)
    else
      do ioctant=1,8
        ierr = octn_free(octn%childs(ioctant))
      end do
    end if
  end function

  integer function octree_free(octree) result(ierr)
    ! Free octree datastructure
    type(octree_t),target,intent(in) :: octree
    ierr = octn_free(octree%first)
  end function

  pure real(dp) function dist_points(p1,p2) result(dist)
    real(dp),intent(in) :: p1(3),p2(3)
    dist = pow2(p1(1)-p2(1))+&
           pow2(p1(2)-p2(2))+&
           pow2(p1(3)-p2(3))
  end function

  pure logical function box_contains(lo,hi,po) result(inside)
    ! Find box that contains point
    real(dp),intent(in) :: lo(3), hi(3), po(3)
    inside = (po(1)>lo(1).and.po(1)<hi(1).and.&
              po(2)>lo(2).and.po(2)<hi(2).and.&
              po(3)>lo(3).and.po(3)<hi(3))
  end function

  pure real(dp) function box_dist(lo,hi,po) result(dist)
    ! Find the distance between point and the box
    real(dp),intent(in) :: lo(3), hi(3), po(3)
    dist = zero
    if (po(1)<lo(1)) dist = dist + pow2(po(1)-lo(1))
    if (po(1)>hi(1)) dist = dist + pow2(po(1)-hi(1))
    if (po(2)<lo(2)) dist = dist + pow2(po(2)-lo(2))
    if (po(2)>hi(2)) dist = dist + pow2(po(2)-hi(2))
    if (po(3)<lo(3)) dist = dist + pow2(po(3)-lo(3))
    if (po(3)>hi(3)) dist = dist + pow2(po(3)-hi(3))
  end function

  pure real(dp) function pow2(x) result(x2)
    real(dp),intent(in) :: x
    x2 = x*x
  end function

  pure integer function get_octant(mi,po) result(ioctant)
    real(dp),intent(in) :: po(3), mi(3)
    integer :: ii,jj,kk
    ii = 0; if (po(1)>=mi(1)) ii = 1
    jj = 0; if (po(2)>=mi(2)) jj = 1
    kk = 0; if (po(3)>=mi(3)) kk = 1
    ioctant = ii*4+jj*2+kk+1
  end function

  pure integer function get_octant_lohi(lo,hi,po) result(ioctant)
    real(dp),intent(in) :: lo(3),hi(3),po(3)
    real(dp) :: mi(3)
    mi = lo+half*(hi-lo)
    ioctant = get_octant(mi,po)
  end function

  pure subroutine get_octants(lo,hi,nids,ids,points,octants)
    ! From a list of points return the corresponding octant
    real(dp),intent(in) :: lo(3), hi(3)
    real(dp),intent(in) :: points(:,:)
    integer,intent(in) :: nids
    integer,intent(in) :: ids(nids)
    integer,intent(out) :: octants(nids)
    real(dp) :: mi(3)
    integer :: id, ipoint
    ! calculate midpoint
    mi = lo+half*(hi-lo)
    do id=1,nids
      ipoint = ids(id)
      octants(id) = get_octant(mi,points(:,ipoint))
    end do
  end subroutine

  pure subroutine get_lo_hi(lo_in,hi_in,lo_out,hi_out,ioctant)
    ! Subdivide a box in an octant
    integer,intent(in) :: ioctant
    real(dp),intent(in) :: lo_in(3), hi_in(3)
    real(dp),intent(out) :: lo_out(3), hi_out(3)
    real(dp) :: de(3)
    de = hi_in-lo_in
    lo_out = lo_in + shifts(:,1,ioctant)*de
    hi_out = lo_in + shifts(:,2,ioctant)*de
  end subroutine

end module m_octree

program octree_main

  use m_octree
  implicit none

  type(octree_t) :: oct
  integer,parameter :: n = 100**3
  integer :: id,ipoint
  real(dp) :: points(3,n)
  real(dp) :: start_time,stop_time
  real(dp) :: dist

  test_box_dist: block
    dist = box_dist([zero,zero,zero],[one,one,one],[-0.1_dp,0.1_dp,0.0_dp])
  end block test_box_dist

  test_octree_init: block
    call random_number(points)
    !do ipoint=1,n
    !  write(*,'(i3,3f12.6)') ipoint, points(:,ipoint)
    !end do

    call cpu_time(start_time)
    oct = octree_init(points,10)
    call cpu_time(stop_time)
    write(*,'(a30,f12.6,a)') 'octree init took:',stop_time-start_time,' [s]'
  end block test_octree_init

  test_octree_find: block
    call cpu_time(start_time)
    do ipoint=1,n
      id = octree_find(oct,points(:,ipoint),dist)
      if (id/=ipoint) call error('wrong point')
    end do
    call cpu_time(stop_time)
    write(*,'(a30,f12.6,a)') 'octree find took:',(stop_time-start_time), ' [s]'
  end block test_octree_find

  test_octree_find_nearest: block
    call cpu_time(start_time)
    do ipoint=1,n
      dist = zero
      ! find nearest point
      id = octree_find_nearest(oct,points(:,ipoint),dist)
      if (id/=ipoint) call error('wrong point')
    end do
    call cpu_time(stop_time)
    write(*,'(a30,f12.6,a)') 'octree find_nearest took:',(stop_time-start_time), ' [s]'
  end block test_octree_find_nearest

  test_octree_find_nearest_pbc: block
    real(dp) :: shift(3)
    call cpu_time(start_time)
    do ipoint=1,n
      dist = zero
      ! find nearest point
      id = octree_find_nearest_pbc(oct,points(:,ipoint),dist,shift)
      if (id/=ipoint) call error('wrong point')
    end do
    call cpu_time(stop_time)
    write(*,'(a30,f12.6,a)') 'octree find_nearest_pbc took:',(stop_time-start_time), ' [s]'
  end block test_octree_find_nearest_pbc

  test_octree_free: block
    integer :: ierr
    call cpu_time(start_time)
    ierr = octree_free(oct)
    call cpu_time(stop_time)
    write(*,'(a30,f12.6,a)') 'octree free took:',(stop_time-start_time), ' [s]'
  end block test_octree_free

  contains
  subroutine error(msg)
   character(len=*),intent(in) :: msg
   write(*,*) msg
   call exit(1)
  end subroutine
end program octree_main
