! hash_wrappers.f90
! Wrappers for adding data to hash table "visits"
!
module hash_wrappers
    use globals, only: qp, Nmax, perc_p, surfaceQ, lattice_type
    use fhash_module__ints_ints
    use ints_module

contains
    ! ##############################################
    ! Hash table wrapper - set
    ! Adds step visit to input point;
    subroutine add_visit(this, pt, step)
        use lattice, only: coord
        implicit none
        type(fhash_type__ints_ints),intent(inout) :: this
        integer,intent(in)  :: pt(:), step   ! lattice point and step number
        type(ints_type)     :: key, value
        logical             :: success
        integer             :: c, v(coord/2) ! temp vector
        ! integer             :: n, q
        key%ints = pt
        call this%get(key, value, success)
        if (.not.success) then  ! unvisited;
            v = 0
            v(1) = step
        else    ! already exists; update
            v = value%ints
            c = count(v > 0)
            if ( c == size(v) ) then
                ! open(10, file='chain.out')
                ! do n = 1, step
                !     write(10, *) (chain(q,n), q = 1, size(chain, 1))
                ! enddo
                ! close(10)
                print*, 'point: ', pt
                print*, 'step: ', step
                print*, 'visits: ', v
                print*, 'count: ', c
                stop 'Error: already at max visits.'
            endif
            v(1 + c) = step
        endif

        value%ints = v
        call this%set(key, value)

    end subroutine add_visit

    ! ##############################################
    ! Hash table wrapper - get
    ! Removes last visit to input point; returns an error flag if something goes wrong
    ! if point has no visits remaining, removes point
    subroutine delete_last_visit(this, pt, err)
        implicit none
        type(fhash_type__ints_ints),intent(inout) :: this
        integer,intent(in)              :: pt(:)   ! lattice point and step number
        logical,optional,intent(out)    :: err
        type(ints_type)                 :: key, value
        integer                         :: c, rmv_success
        logical                         :: success

        key%ints = pt
        call this%get(key, value, success)
        if (.not.success) then
            ! key does not exist
            if (present(err)) err = .true.
            return
        else
            c = count(value%ints > 0)
            if (c == 0) then      ! point is unvisited, remove whole entry
                call this%remove(key, rmv_success)
                if (present(err) .and. rmv_success /= 1) err = .true.
                return
            elseif (c == 1) then  ! only single visit, can remove whole entry
                call this%remove(key, rmv_success)
            else                ! more than one visit, remove only last visit
                value%ints(c) = 0
                call this%set(key, value)
            endif
        endif
        if (present(err)) err = .false.

    end subroutine delete_last_visit

    ! ##############################################
    ! Hash table wrapper - clean
    ! Resets hash table
    subroutine reset_visits(this)
        implicit none
        type(fhash_type__ints_ints),intent(inout) :: this

        ! Hash table has already been declared in globals so clear should have no effect on an empty hash table
        call this%clear()
        call this%reserve(2*Nmax)

    end subroutine reset_visits

    ! ##############################################
    ! Hash table wrapper - clean
    ! Deallocates hash table
    subroutine clear_hash_table(this)
        implicit none
        type(fhash_type__ints_ints),intent(inout) :: this

        ! Hash table has already been declared in globals so clear should have no effect on an empty hash table
        call this%clear()

    end subroutine clear_hash_table

    ! ##############################################
    ! Hash table wrapper - get
    function get_visits(this, pt) result(steps)
        use lattice, only: coord
        implicit none
        type(fhash_type__ints_ints),intent(inout) :: this
        integer,intent(in)  :: pt(:)
        type(ints_type)     :: key, value
        logical             :: success
        integer             :: steps(coord/2)

        ! step = 0     ! error case
        key%ints = pt
        call this%get(key, value, success)
        if (.not.success) then
            steps = 0
        else
            steps = value%ints
        endif

    end function get_visits

    ! ##############################################
    ! Hash table wrapper - get
    function count_visits(this, pt) result(c)
        implicit none
        type(fhash_type__ints_ints),intent(inout) :: this
        integer,intent(in)  :: pt(:)
        type(ints_type)     :: key, value
        logical             :: success
        integer             :: c

        c = 0
        key%ints = pt
        call this%get(key, value, success)
        if (success) then
            c = count(value%ints > 0)
        endif

    end function count_visits

    ! ##############################################
    ! Hash table wrapper - get
    function get_last_visit(this, pt) result(step)

        type(fhash_type__ints_ints),intent(inout) :: this
        integer,intent(in)  :: pt(:)
        type(ints_type)     :: key, value
        logical             :: success
        integer             :: c, step

        step = 0     ! error case
        key%ints = pt
        call this%get(key, value, success)
        if (success) then
            c = count(value%ints > 0)
            if (c > 0) then
                step = value%ints(c)
            endif
        endif

    end function get_last_visit

    ! ##############################################
    ! routines for list of holes/imperfections
    ! ##############################################

    ! ##############################################
    ! Hash table wrapper - set
    ! Adds hole to input point;
    subroutine add_hole(this, pt, step)
        use lattice, only: coord
        implicit none
        type(fhash_type__ints_ints),intent(inout) :: this
        integer,intent(in)  :: pt(:), step   ! lattice point and step number
        type(ints_type)     :: key, value
        integer             :: i

        key%ints = pt
        value%ints = [(0, i = 1, coord/2)]   ! unnecessary but might be useful later
        call this%set(key, value)
        ! OLD - ANNEALED DISORDER METHOD DOES NOT WORK
        ! history_holes(:,this%key_count()) = [ step, pt ]

    end subroutine add_hole

    ! ##############################################
    ! Hash table wrapper - clean
    ! Resets hash table
    subroutine reset_holes(this, d, probq)
        implicit none
        type(fhash_type__ints_ints),intent(inout) :: this
        integer,intent(in)  :: d
        real(qp),intent(out)  :: probq
        integer             :: q, n
        real,allocatable    :: r(:)
        integer,allocatable :: b(:)
        real,parameter      :: pi = 3.141592653589793
        ! real                :: v(4)

        if (d < 2) stop 'Error: incorrect dimension for generating holes'
        n = 0
        ! volume of unit sphere
        ! v = [2., pi, 4*pi/3, pi**2/2]
        if (lattice_type == 'squ' .or. lattice_type == 'cub') then
            ! volume of system - octohedron inside cube
            ! + 1 to compensate for integer division
            if (d == 3) then
                if (surfaceQ) then
                    n = (2*Nmax + 1)**(d - 1) / 2 + 1
                else
                    n = (2*Nmax + 1)**d / 6 + 1
                end if
            end if
            if (d == 2) then
                if (surfaceQ) then
                    n = (2*Nmax + 1)**(d - 1)
                else
                    n = (2*Nmax + 1)**d / 2 + 1
                end if
            end if
        else if (lattice_type == 'tri') then
            ! volume of system - hexagon with side length Nmax
            n = 6*Nmax*(Nmax - 1)/2
            ! transforms to the following in rectilinear format
            !     /--|
            !    /| /|
            !   / |/ |
            !  /--/--/
            !  | /| /
            !  |/ |/
            !  |--/
            ! so upper left and lower right quadrants are same as squ lattice, other quadrants include the whole square
        else
            ! hex lattice?
            n = 0
        end if

        allocate(r(d))
        allocate(b(d))

        ! calculate number of holes
        ! normal approximation of binomial distribution for large n
        ! Box-Muller transformation for uniform to normal
        call RANDOM_NUMBER(r)
        q = floor(n*perc_p + sqrt(n*perc_p*(1 - perc_p)) * sqrt(-2*log(r(1))) * cos(2*pi*r(2))) + 1
        probq = exp(-(q - n*perc_p)**2 / (2*n*perc_p*(1 - perc_p))) / sqrt(2*pi*n*perc_p*(1-perc_p))

        call this%clear()
        call this%reserve(q)

        ! generate points and add to hash table
        do while (count_holes(this) < q)
            ! ! first pass to avoid extra checks
            ! call RANDOM_NUMBER(r)
            ! if (lattice_type == 'squ' .or. lattice_type == 'cub') then
            !     if (surfaceQ) then
            !         b(1:d-1) = floor( (2*Nmax + 1)*r(1:d-1) ) - Nmax
            !         b(d) = floor( (Nmax + 1)*r(d) )
            !     else
            !         b = floor( (2*Nmax + 1)*r ) - Nmax
            !     endif
            ! else if (lattice_type == 'tri') then
            !     b = floor( (2*Nmax + 1)*r ) - Nmax
            ! else
            !     ! hex lattice?
            ! end if

            ! simple start
            b = 0

            ! Check with metric and exclude origin and duplicates
            if (lattice_type == 'squ' .or. lattice_type == 'cub') then
                do while (get_holes(this, b) .or. sum(abs(b)) > Nmax .or. all(b == 0) .or. maxval(abs(b)) > Nmax)
                    call RANDOM_NUMBER(r)
                    if (surfaceQ) then
                        b(1:d-1) = floor( (2*Nmax + 1)*r(1:d-1) ) - Nmax
                        b(d) = floor( (Nmax + 1)*r(d) )
                    else
                        b = floor( (2*Nmax + 1)*r ) - Nmax
                    endif
                enddo
            else if (lattice_type == 'tri') then
                ! second condition is to account for rectilinear format
                do while (get_holes(this, b) .or. (product(b) < 0 .and. sum(abs(b)) > Nmax) .or. all(b == 0) &
                        .or. maxval(abs(b)) > Nmax)
                    call RANDOM_NUMBER(r)
                    b = floor( (2*Nmax + 1)*r ) - Nmax
                enddo
            else
                ! hex lattice?
            end if

            call add_hole(this, b, 0)  ! value does not matter
        enddo

        deallocate(r, b)

    end subroutine reset_holes

    ! ##############################################
    ! Hash table wrapper - get
    function get_holes(this, pt) result(success)
        implicit none
        type(fhash_type__ints_ints),intent(inout) :: this
        integer,intent(in)  :: pt(:)
        type(ints_type)     :: key, value
        logical             :: success

        key%ints = pt
        ! print*,pt
        call this%get(key, value, success)

    end function get_holes

    ! ##############################################
    ! Hash table wrapper - get
    function count_holes(this) result(c)
        implicit none
        type(fhash_type__ints_ints),intent(inout) :: this
        integer             :: c

        c = this%key_count()

    end function count_holes

    ! ##############################################
    ! Hash table wrapper - convert hash table to list
    subroutine list_holes(this, l)
        implicit none
        type(fhash_type__ints_ints),intent(inout) :: this
        integer,intent(in)  :: l
        integer,allocatable :: holes_list(:,:)
        integer             :: i, q, status
        type(ints_type)     :: key, value
        type(fhash_type_iterator__ints_ints)    :: it

        q = count_holes(this)
        ! if (l == 3) q = 2*floor( perc_p*(2*Nmax + 1)**l / 6 ) + 1
        ! if (l == 2) q = 2*floor( perc_p*(2*Nmax + 1)**l / 2 ) + 1
        if (.not. allocated(holes_list)) allocate(holes_list(l,q))
        holes_list = 0

        ! print*,'Defects: ',count_holes()
        ! print*,shape(holes_list)

        call it%begin(this)
        status = 0
        i = 0
        ! do while (status == 0)
        do i = 1, q
            call it%next(key, value, status)
            ! i = i + 1
            if (status /= 0) then
                print*,'Error: fhash iterator failed ',i
                exit
            endif

            holes_list(:,i) = key%ints
        enddo

        print*,'Holes converted to list; list dicarded.'

    end subroutine list_holes

end module hash_wrappers
