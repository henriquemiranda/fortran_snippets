!
! Author: Henrique Miranda
! Build an octree to find nearest negihbour points
!
module m_octree

  implicit none

  integer,parameter :: dp=8
  real(dp),parameter :: zero=0.0_dp,half=0.5_dp,one=1.0_dp,two=2.0_dp
  real(dp),parameter :: tol12=0.000000000001_dp
  real(dp),protected :: shifts(3,2,8)

  type :: octree_node_t
    type(octree_node_t),pointer :: childs(:) ! if node this is allocated
    integer,allocatable :: ids(:)            ! if leaf this is allocated
  end type octree_node_t

  type :: octree_t
    real(dp) :: hi(3)
    real(dp) :: lo(3)
    integer :: max_npoints
    real(dp),pointer :: points(:,:)
    type(octree_node_t) :: first
  end type octree_t

!--------------------------------------------------------------------

  contains
  type(octree_t) function octree_init(points,max_npoints,lo,hi) result (new)
    integer,intent(in) :: max_npoints
    real(dp), target, intent(in) :: points(:,:)
    real(dp), intent(in) :: lo(3), hi(3)
    real,parameter :: ieps = 0.1
    integer,allocatable :: ids(:)
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

    ! box dimensions
    new%lo = lo
    new%hi = hi

    !first octree contains all the points
    npoints = size(points,2)
    new%points => points
    new%max_npoints = max_npoints
    ids = [(ii,ii=1,npoints)]
    new%first = octree_node_build(new,new%lo,new%hi,npoints,ids)

  end function octree_init

  type(octree_node_t) recursive function octree_node_build(octree,lo,hi,nids,ids) result (new)
    type(octree_t),intent(in) :: octree
    integer,intent(in) :: nids
    integer,intent(in) :: ids(nids)
    integer :: id, counter, ioctant, ipoint
    integer :: octants(nids)
    integer :: new_ids(nids)
    real(dp) :: lo(3), hi(3), new_lo(3), new_hi(3)

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
      new%childs(ioctant) = octree_node_build(octree,new_lo,new_hi,counter,new_ids)
    end do

  end function octree_node_build

  integer function octree_find(octree,point,dist) result(closest_id)
    ! Find the closest point in the box that contains it
    type(octree_t),target,intent(in) :: octree
    real(dp),intent(in) :: point(3)
    real(dp),intent(out) :: dist
    type(octree_node_t),pointer :: octn
    integer :: id, ipoint, ioctant
    real(dp) :: hi(3),lo(3),hi_out(3),lo_out(3)
    real(dp) :: trial_dist

    closest_id = 0
    dist = huge(dist)
    octn => octree%first
    lo = octree%lo
    hi = octree%hi

    !check if the point is inside the initial box
    trial_dist = box_dist(lo,hi,point)
    if (trial_dist>0) then
      closest_id = -1
      return
    end if

    do
      ! if leaf node
      if (allocated(octn%ids)) then
        do id=1,size(octn%ids)
          ipoint = octn%ids(id)
          trial_dist = dist_points(octree%points(:,ipoint),point)
          if (trial_dist > dist) cycle
          dist = trial_dist
          closest_id = ipoint
        end do
        return
      end if
      ! get octant of this point
      ioctant = get_octant_lohi(lo,hi,point)
      ! point to this node
      octn => octn%childs(ioctant)
      ! get lo and hi
      call get_lo_hi(lo,hi,lo_out,hi_out,ioctant)
      lo = lo_out; hi = hi_out
    end do
  end function octree_find

  integer function octree_find_nearest(octree,point,dist) result(nearest_id)
    ! Find the ids of the points whose distance to point is smaller than max_dist
    ! counter is the number of points found so far
    type(octree_t),target,intent(in) :: octree
    real(dp),intent(in) :: point(3)
    real(dp),intent(inout) :: dist
    !check if the point is inside the initial box
    dist = box_dist(octree%lo,octree%hi,point)
    if (dist>0) then
      nearest_id = -1
      return
    end if
    nearest_id = octn_find_nearest(octree,octree%first,octree%lo,octree%hi,point,dist)
  end function octree_find_nearest

  integer recursive function octn_find_nearest(octree,octn,lo,hi,point,min_dist) result(closest_id)
    ! find the nearest point by recursion
    type(octree_t),intent(in) :: octree
    type(octree_node_t),intent(in) :: octn
    real(dp),intent(in) :: point(3)
    real(dp),intent(in) :: lo(3), hi(3)
    real(dp),intent(inout) :: min_dist
    real(dp) :: new_lo(3), new_hi(3)
    real(dp) :: dist
    integer :: id, ioctant, ipoint, trial_id
    closest_id = 0
    ! compute distance of point to this box (octant)
    dist = box_dist(lo,hi,point)
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
      call get_lo_hi(lo,hi,new_lo,new_hi,ioctant)
      trial_id = octn_find_nearest(octree,octn%childs(ioctant),new_lo,new_hi,point,min_dist)
      if (trial_id==0) cycle
      closest_id = trial_id
    end do
  end function octn_find_nearest

  integer function octree_find_nearest_pbc(octree,point,dist,shift) result(id)
    ! Same as octree find but using periodic boundary conditions
    type(octree_t),target,intent(in) :: octree
    real(dp),intent(in) :: point(3)
    real(dp),intent(inout) :: dist
    real(dp),intent(out) :: shift(3)
    logical :: found
    integer :: trial_id, first_id
    integer :: ii,jj,kk
    real(dp) :: trial_dist
    real(dp) :: first_shift(3), trial_shift(3)
    real(dp) :: po(3)
    ! bring the point inside the box
    po = modulo(point,one)
    ! compute shift
    first_shift = po-point
    first_id = octn_find_nearest(octree,octree%first,octree%lo,octree%hi,po,dist)
    id = first_id
    ! try unitary shifts
    found=.false.
    do ii=-1,1
      do jj=-1,1
        do kk=-1,1
          if (ii==0.and.jj==0.and.kk==0) cycle
          trial_shift = first_shift+[ii,jj,kk]
          ! compute shortest distance
          trial_dist = dist+tol12
          trial_id = octn_find_nearest(octree,octree%first,octree%lo,octree%hi,&
                                       point+trial_shift,trial_dist)
          ! if smaller than previous distance, store this shift and distance
          if (trial_dist>dist) cycle
          found = .true.
          dist  = trial_dist
          id    = trial_id
          shift = trial_shift
        end do
      end do
    end do
    if (.not.found) shift = first_shift
  end function octree_find_nearest_pbc

  integer recursive function octn_free(octn) result(ierr)
    type(octree_node_t) :: octn
    integer :: ioctant
    ! if leaf deallocate ids
    if (allocated(octn%ids)) then
      deallocate(octn%ids,stat=ierr)
    else
      do ioctant=1,8
        ierr = octn_free(octn%childs(ioctant))
      end do
    end if
  end function octn_free

  integer function octree_free(octree) result(ierr)
    ! Free octree datastructure
    type(octree_t),target,intent(in) :: octree
    ierr = octn_free(octree%first)
  end function octree_free

  pure real(dp) function dist_points(p1,p2) result(dist)
    real(dp),intent(in) :: p1(3),p2(3)
    dist = pow2(p1(1)-p2(1))+&
           pow2(p1(2)-p2(2))+&
           pow2(p1(3)-p2(3))
  end function dist_points

  pure logical function box_contains(lo,hi,po) result(inside)
    ! Find box that contains point
    real(dp),intent(in) :: lo(3), hi(3), po(3)
    inside = (po(1)>lo(1).and.po(1)<hi(1).and.&
              po(2)>lo(2).and.po(2)<hi(2).and.&
              po(3)>lo(3).and.po(3)<hi(3))
  end function box_contains

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
  end function box_dist

  pure real(dp) function pow2(x) result(x2)
    real(dp),intent(in) :: x
    x2 = x*x
  end function pow2

  pure integer function get_octant(mi,po) result(ioctant)
    real(dp),intent(in) :: po(3), mi(3)
    integer :: ii,jj,kk
    ii = 0; if (po(1)>=mi(1)) ii = 1
    jj = 0; if (po(2)>=mi(2)) jj = 1
    kk = 0; if (po(3)>=mi(3)) kk = 1
    ioctant = ii*4+jj*2+kk+1
  end function get_octant

  pure integer function get_octant_lohi(lo,hi,po) result(ioctant)
    real(dp),intent(in) :: lo(3),hi(3),po(3)
    real(dp) :: mi(3)
    mi = half*(hi+lo)
    ioctant = get_octant(mi,po)
  end function get_octant_lohi

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
    mi = half*(hi+lo)
    do id=1,nids
      ipoint = ids(id)
      octants(id) = get_octant(mi,points(:,ipoint))
    end do
  end subroutine get_octants

  pure subroutine get_lo_hi(lo_in,hi_in,lo_out,hi_out,ioctant)
    ! Subdivide a box in an octant
    integer,intent(in) :: ioctant
    real(dp),intent(in) :: lo_in(3), hi_in(3)
    real(dp),intent(out) :: lo_out(3), hi_out(3)
    real(dp) :: de(3)
    de = hi_in-lo_in
    lo_out = lo_in + shifts(:,1,ioctant)*de
    hi_out = lo_in + shifts(:,2,ioctant)*de
  end subroutine get_lo_hi

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

  test_octree_init: block
    call random_number(points)
    !do ipoint=1,n
    !  write(*,'(i3,3f12.6)') ipoint, points(:,ipoint)
    !end do

    call cpu_time(start_time)
    oct = octree_init(points,2**4,[-one,-one,-one],[two,two,two])
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
    integer :: ii,jj,kk
    integer :: naive_ipoint, iprobe
    real(dp) :: naive_dist, trial_dist
    real(dp) :: naive_shift(3), trial_shift(3), shift(3), point(3)
    call cpu_time(start_time)
    do ipoint=1,n
      dist = tol12
      ! find nearest point
      id = octree_find_nearest_pbc(oct,points(:,ipoint),dist,shift)
      if (id/=ipoint) call error('wrong point')
    end do
    call cpu_time(stop_time)
    write(*,'(a30,f12.6,a)') 'octree find_nearest_pbc took:',(stop_time-start_time), ' [s]'

    ! now compare with a naive approach
    call cpu_time(start_time)
    do iprobe=1,10
      naive_dist = huge(naive_dist)
      naive_ipoint = -1
      ! generate random point between -1 and 1
      call random_number(point)
      point = point*two-one
      ! find nearest point considering periodic boundary conditions
      id = octree_find(oct,point,dist)
      id = octree_find_nearest_pbc(oct,point,dist,shift)
      ! now find it naively
      trial_dist  = huge(trial_dist)
      naive_shift = huge(naive_shift)
      do ii=-2,2
        do jj=-2,2
          do kk=-2,2
            do ipoint=1,size(points,2)
              trial_shift = [ii,jj,kk]
              trial_dist = dist_points(points(:,ipoint),point+trial_shift)
              if (trial_dist>naive_dist) cycle
              naive_dist   = trial_dist
              naive_ipoint = ipoint
              naive_shift  = trial_shift
            end do
          end do
        end do
      end do
      if (id<1)                    call error('did not find point')
      if (id/=naive_ipoint)        call error('wrong point')
      if (any(shift/=naive_shift)) call error('wrong shift')
    end do
    call cpu_time(stop_time)
    write(*,'(a30,f12.6,a)') 'naive took:',(stop_time-start_time), ' [s]'
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
