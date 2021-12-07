    ! lattice.f90
    ! Set up variables to describe the lattice structure
    !
module lattice
    use globals, only: qp
    implicit none
    !private
    !public init_lattice, coord, directions, origin, lattdim

    ! declare lattice quantities here to be global variables
    integer             :: coord, lattdim
    integer,allocatable :: directions(:,:), origin(:)
    real(qp)            :: mu_cc


contains
    ! ##############################################
    ! Initialise lattice structure
    subroutine init_lattice(latt, trailsQ)
        implicit none
        character(len = 3), intent(in)  :: latt
        logical, intent(in)             :: trailsQ
        integer             :: i
        real(qp)            :: mu_vals(6)

        select case (latt)
        case('squ')
            coord = 4
            lattdim = 2
            allocate(directions(lattdim,coord))
            directions(:,:) = reshape([1,0,0,1,-1,0,0,-1], [lattdim,coord])
            if (trailsQ) then
                mu_cc = 2.72062_qp
            else
                mu_cc = 2.63815853_qp
            end if
        case('cub')
            coord = 6
            lattdim = 3
            allocate(directions(lattdim,coord))
            directions(:,:) = reshape([1,0,0,0,1,0,0,0,1,-1,0,0,0,-1,0,0,0,-1], [lattdim,coord])
            mu_cc = 4.65987_qp
        ! case('hc4')
        !     coord = 8
        !     lattdim = 4
        !     allocate(directions(lattdim,coord))
        !     directions(:,:) = reshape([   1,0,0,0,&
        !                                 0,1,0,0,&
        !                                 0,0,1,0,&
        !                                 0,0,0,1,&
        !                                 -1,0,0,0,&
        !                                 0,-1,0,0,&
        !                                 0,0,-1,0,&
        !                                 0,0,0,-1],&
        !         [lattdim,coord])
        !     mu_cc = 6.77404_qp
        case('tri')
            coord = 6
            lattdim = 2
            allocate(directions(lattdim,coord))
            directions(:,:) = reshape([1,0,1,1,0,1,-1,0,-1,-1,0,-1], [lattdim,coord])
            if (trailsQ) then
                mu_cc = 4.524_qp
            else
                mu_cc = 4.15079_qp
            end if
        case('tet')
            coord = 4
            lattdim = 3
            allocate(directions(lattdim,coord))
            directions(:,:) = reshape([1,1,1,1,-1,-1,-1,1,-1,-1,-1,1], [lattdim,coord])
            mu_cc = 2d0
        case('man')
            coord = 4
            lattdim = 2
            allocate(directions(lattdim,coord))
            directions(:,:) = reshape([1,0,0,1,-1,0,0,-1], [lattdim,coord])
            mu_cc = 1.733535_qp
        case('hex')
            ! NB: uses 'brickworks' lattice. If surface exists, it is horizontal such that
            ! origin is not in surface. I.e. x is point on surface, origin lies between them.
            !      |
            !   x__|__x
            !   |     |
            coord = 3
            lattdim = 2
            allocate(directions(lattdim,coord))
            directions(:,:) = reshape([1,0,-1,0,0,1], [lattdim,coord])
            mu_cc = sqrt(2.0_qp + sqrt(2.0_qp))
        case default
            ! general hyoercubic lattice
            if (latt(1:2) == 'hc') then
                ! cD=
                read(latt(3:3),*) lattdim
                coord = 2*lattdim
                allocate(directions(lattdim,coord))
                directions = 0
                do i = 1, coord/2
                    directions(i,i) = 1
                enddo
                directions(:,coord/2+1:coord) = -directions(:,1:coord/2)

                if (lattdim < 7) then
                    mu_vals = [1.0_qp,2.63815853_qp, 4.65987_qp, 6.77404_qp, 8.83843_qp, 10.87809_qp]
                    mu_cc = mu_vals(lattdim)
                else
                    mu_cc = 2.06584_qp*lattdim - 1.50565_qp
                endif
            else
                print *, "Lattice *",latt,"* undefined."
            endif
        end select

        allocate(origin(lattdim))
        origin = 0

    end subroutine init_lattice

    ! ##############################################
    ! Clean up lattice structures
    subroutine clean_lattice()
        implicit none

        if (allocated(directions)) deallocate(directions)
        if (allocated(origin)) deallocate(origin)

    end subroutine clean_lattice

end module
