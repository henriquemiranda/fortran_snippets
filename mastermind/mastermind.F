module m_itertools
    implicit none

    type product_type
        integer :: n ! total number of values
        integer :: reps ! number of repetitions
        integer,allocatable :: idx(:) ! indexes of values
        integer,allocatable :: values(:) ! values to permute
    end type

    contains
    subroutine product_init(pt,values,reps)
        type(product_type) :: pt
        integer :: values(:)
        integer,optional :: reps

        !values
        pt%n = size(values)
        allocate(pt%values(pt%n))
        pt%values = values

        !repetitions
        pt%reps = 1; if (present(reps)) pt%reps = reps
        allocate(pt%idx(pt%reps))
        pt%idx = 1
    end subroutine product_init

    subroutine product_iter(pt,values,done)
        type(product_type) :: pt
        integer :: values(pt%reps)
        integer :: i
        logical :: done

        done = .false.
        ! carry out if needed
        do i=1,pt%reps
            if (pt%idx(i)>pt%n) then
                pt%idx(i) = 1
                pt%idx(i+1) = pt%idx(i+1) + 1
            endif
        enddo
        ! return values
        do i=1,pt%reps
            values(i) = pt%values(pt%idx(i))
        enddo
        if (all(pt%idx == pt%n)) then
            done = .true.
            return
        endif
        pt%idx(1) = pt%idx(1) + 1
    end subroutine product_iter

    subroutine product_reset(pt)
        type(product_type) :: pt
        pt%idx = 1
    end subroutine product_reset

    subroutine product_free(pt)
        type(product_type) :: pt
        deallocate(pt%values)
        deallocate(pt%idx)
    end subroutine product_free
end module m_itertools

module m_mastermind
    use m_itertools
    implicit none

    type gamestate
        integer :: row_size
        integer :: ncolors
        integer :: max_plays ! maximum number of plays
        integer :: nplayed ! number of played moves
        integer,allocatable :: all_combinations(:,:)
        integer,allocatable :: played_combinations(:,:)
        integer,allocatable :: all_combinations_state(:) ! whether the combination is in the set (1) or not (0)
        integer,allocatable :: nblacks(:)
        integer,allocatable :: nwhites(:)
        integer,allocatable :: combination(:)
    end type

    contains
    subroutine init_mm(gs,ncolors,row_size,max_plays)
        type(gamestate) :: gs
        integer :: ncolors, row_size, max_plays

        type(product_type) :: pt
        integer :: n
        integer :: colors(ncolors)
        logical :: done
        real ::  tmp_combination(row_size)

        gs%ncolors = ncolors
        gs%row_size = row_size
        gs%max_plays = max_plays
        allocate(gs%all_combinations(gs%row_size,gs%ncolors**gs%row_size))
        allocate(gs%played_combinations(gs%row_size,gs%max_plays))
        allocate(gs%all_combinations_state(gs%ncolors**gs%row_size))
        allocate(gs%nblacks(gs%ncolors**gs%row_size))
        allocate(gs%nwhites(gs%ncolors**gs%row_size))
        gs%nblacks = 0
        gs%nwhites = 0
        gs%nplayed = 0
        gs%all_combinations_state = 1

        ! fill in all possible combinations
#if 0
        do n=1,gs%ncolors**gs%row_size
            gs%all_combinations(1,n) = mod(n-1,gs%ncolors)+1
            do i=2,gs%row_size
                gs%all_combinations(i,n) = mod((n-1)/gs%ncolors**(i-1),gs%ncolors)+1
            end do
        end do
#else
        colors = [(n, n = 1, gs%ncolors, 1)]
        call product_init(pt,colors,gs%row_size)
        done = .false.
        do n=1,gs%ncolors**gs%row_size
            call product_iter(pt,gs%all_combinations(:,n),done)
            if (done) exit
        end do
        call product_free(pt)
#endif
        ! combination
        allocate(gs%combination(gs%row_size))
        call random_number(tmp_combination)
        gs%combination = floor((gs%ncolors-1)*tmp_combination)+1
    end subroutine init_mm

    subroutine play_combination(gs,combination,nblacks,nwhites,removed_counter,won)
        type(gamestate) :: gs
        integer,intent(in) :: combination(gs%row_size)
        integer,intent(out) :: nblacks,nwhites
        integer,intent(out) :: removed_counter
        logical,intent(out) :: won

        won = .false.
        ! set combination in played combinations
        gs%nplayed = gs%nplayed + 1
        gs%played_combinations(:,gs%nplayed) = combination

        ! get result of move
        call evaluate_combination(gs,gs%combination,combination,nblacks,nwhites)
        if (nblacks==4) won = .true.

        ! remove  combinations
        call remove_combinations(gs,combination,nblacks,nwhites,.true.,removed_counter)
    end subroutine play_combination

    subroutine remove_combinations(gs,combination,nblacks,nwhites,remove,removed_counter)
        type(gamestate) :: gs
        integer :: combination(gs%row_size)
        integer,intent(in) :: nblacks, nwhites
        logical,intent(in) :: remove
        integer,intent(out) :: removed_counter

        integer :: n
        integer :: tmp_nblacks,tmp_nwhites

        ! remove all combinations that do not yield the same outcome
        removed_counter = 0
        do n=1,gs%ncolors**gs%row_size
            ! We skip combinations that we have already checked
            if (gs%all_combinations_state(n) == 0) cycle
            ! Test this combination against against one in S
            call evaluate_combination(gs,combination,gs%all_combinations(:,n),tmp_nblacks,tmp_nwhites)
            ! if the outcome is different then we remove
            !write(*,*) gs%all_combinations(:,n), '|', tmp_nblacks, tmp_nwhites, nblacks, nwhites
            if (tmp_nblacks/=nblacks.or.tmp_nwhites/=nwhites) then
                if (remove) gs%all_combinations_state(n) = 0
                removed_counter = removed_counter + 1
            end if
        end do
    end subroutine remove_combinations

    subroutine suggest_combination(gs,values)
        type(gamestate) :: gs
        integer :: values(gs%row_size)

        type(product_type) :: pt
        integer :: i, n, win_idx
        logical :: done
        integer :: combination(gs%row_size)
        integer :: outcome(2), pinvals(gs%row_size+1)
        integer :: active_counter, removed_counter, min_score, max_min_score

        ! if there is only one combination in S then play the winning hand
        win_idx = 0
        active_counter = 0
        do n=1,gs%ncolors**gs%row_size
            if (gs%all_combinations_state(n)==0) cycle
            active_counter = active_counter + 1
            win_idx = n
        enddo
        if (active_counter == 1) then
            values = gs%all_combinations(:,win_idx)
            return
        endif

        ! init generator of all possible outcomes of a move
        pinvals = [(n, n = 0, gs%row_size, 1)]
        call product_init(pt,pinvals,reps=2)

        ! look in all possible combinations
        max_min_score = 0
        all_combinations: do n=1,gs%ncolors**gs%row_size
            combination = gs%all_combinations(:,n)
            ! if this combination was played already then move on
            do i=1,gs%nplayed
                if (all(combination==gs%played_combinations(:,i))) cycle
            enddo
            ! if (gs%all_combinations_state(n)==0) cycle
            ! test all possible outcomes
            call product_reset(pt)
            done = .false.
            min_score = huge(min_score)
            do while (.not.done)
                call product_iter(pt,outcome,done)
                if (sum(outcome)>gs%row_size) cycle
                ! Acutally here we don't remove the combinations, but we get a count of how many would be removed
                call remove_combinations(gs,combination,outcome(1),outcome(2),.false.,removed_counter)
                if (min_score>removed_counter) min_score = removed_counter
                !write(*,*) combination,':', outcome, ':',removed_counter
            end do
            if (min_score>max_min_score) then
                max_min_score = min_score
                values = combination
            endif
        end do all_combinations

        ! free memory
        call product_free(pt)

    end subroutine suggest_combination

    subroutine evaluate_combination(gs,ref_combination,combination,nblacks,nwhites)
        type(gamestate) :: gs
        integer :: ref_combination(gs%row_size)
        integer :: combination(gs%row_size)

        logical :: skipb(gs%row_size),skipwi(gs%row_size),skipwj(gs%row_size)
        integer :: nblacks, nwhites
        integer :: i,j

        skipb = .false.
        skipwi = .false.
        skipwj = .false.
        nblacks=0
        nwhites=0
        !write(*,*) ref_combination
        !write(*,*) combination
        ! count blacks
        do i=1,gs%row_size
            if (skipb(i)) cycle
            if (combination(i)/=ref_combination(i)) cycle
            nblacks = nblacks + 1
            skipb(i) = .true.
        enddo
        ! count whites
        do i=1,gs%row_size
            if (skipb(i).or.skipwi(i)) cycle
            ! find color in ref_combination
            do j=1,gs%row_size
                if (skipb(j).or.skipwj(j).or.skipwi(i)) cycle
                if (combination(i)/=ref_combination(j)) cycle
                !write(*,*) i,j,combination(i),ref_combination(j),'-'
                nwhites = nwhites + 1
                skipwi(i)=.true.
                skipwj(j)=.true.
            end do
        end do
   end subroutine evaluate_combination

    subroutine gs_write(gs)
        type(gamestate) :: gs
        integer :: n, counter = 0
        counter = 0
        do n=1,gs%ncolors**gs%row_size
            if (gs%all_combinations_state(n)==0) cycle
            counter = counter + 1
            write(*,*) counter,n,gs%all_combinations(:,n)
        end do
    end subroutine gs_write

    subroutine free_mm(gs)
        type(gamestate) :: gs
        deallocate(gs%all_combinations)
    end subroutine free_mm

end module m_mastermind

!===============================================================
program master
    use m_itertools
    use m_mastermind
    implicit none

    type(gamestate) :: gs
    integer, parameter :: MAX_PLAYS=10
    integer, parameter :: ROW_SIZE=4
    integer, parameter :: NCOLORS=6
    integer, parameter :: reps = 3
    integer, parameter :: first_move(4) = [1,1,2,2]
    integer :: nblacks, nwhites
    integer :: values(ROW_SIZE)
    logical :: won
    integer :: removed_counter

    ! set intial state of the game
    call init_mm(gs,NCOLORS,ROW_SIZE,MAX_PLAYS)

    ! play once
    call play_combination(gs,first_move,nblacks,nwhites,removed_counter,won)
    write(*,*) gs%nplayed,'|',first_move,'|',nblacks,nwhites
    !call gs_write(gs)

    won = .false.
    do while (.not.won)
        call suggest_combination(gs,values)
        call play_combination(gs,values,nblacks,nwhites,removed_counter,won)
        write(*,*) gs%nplayed,'|',values,'|',nblacks,nwhites
        if (gs%nplayed == gs%max_plays) exit
        !call gs_write(gs)
    end do
    write(*,*)
    write(*,*) 0,'|',gs%combination

end program master
