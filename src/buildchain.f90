    ! buildchain.f90
    ! Routines for building chains and finding atmosphere
    !
module buildchain
    use hash_wrappers
    use globals, only: Nmax, collisionsQ, groovesQ, surfaceQ, superQ, lattice_type, percQ, perc_p, boxQ, Lmax
    use lattice, only: lattdim, coord, directions, origin
    use fhash_module__ints_ints, only: fhash_type__ints_ints
    implicit none
    private
    public neighbours, super_neighbours, percolation


contains
    ! ##############################################
    ! Neighbour check (atmosphere)
    function neighbours(chain,visits) result(atmo)
        implicit none
        integer, intent(in) :: chain(:,:)       ! input
        type(fhash_type__ints_ints), intent(inout) :: visits
        integer             :: atmo(coord)          ! output
        integer             :: n, k
        integer             :: next(lattdim), this(lattdim), prev(lattdim), p(coord/2), q(coord/2)
        integer             :: jprev(lattdim), jnext(lattdim)
        integer             :: i, j             ! loop variables
        integer             :: e(2), t(2)      ! grooves: directions of prev and next for existing and test branches

        n = size(chain, dim = 2)
        this = chain(:,n)
        if (collisionsQ) then
            q = get_visits(visits, this)
        else
            q = 0
        end if
        atmo = (/ (i, i = 1, coord) /)
        if (n > 1) then
            prev = chain(:, n-1)
        else
            if (surfaceQ .and. lattice_type /= 'hex') then
                if (lattice_type == 'tri') then
                    atmo(4:5) = 0
                else    ! square and cubic lattices
                    atmo(coord) = 0
                end if
            end if
            return
        end if

dirs:   do i = 1, coord
            if (lattice_type == 'hex') then
                if (i == coord) then
                    if (modulo(sum(this), 2) == 1) then
                        next = this - directions(:,i)
                    else
                        next = this + directions(:,i)
                    end if
                else
                    next = this + directions(:,i)
                end if
            else
                next = this + directions(:,i)
            end if   ! lattice_type=='hex'

            ! Prevent retracing last step
            if (all( next == prev )) then
                atmo(i) = 0
                cycle dirs
            end if

            ! prevent crossing adsorption surface
            if (surfaceQ) then
                if (next(lattdim) < 0) then
                    atmo(i) = 0
                    cycle dirs
                end if
            end if   ! surfaceQ

            ! confine to box
            ! NB: box size is number of lattice SITES, whereas chain length is number of STEPS, and chains start from origin.
            ! Hence the equality in cases below.
            if (boxQ) then
                if (lattice_type == 'hex') then
                    !
                else
                    k = mod(i - 1, coord/2) + 1
                    if (Lmax <= abs(next(k) - maxval(chain(k,:)))) then     ! exiting box to the left/bottom/down
                        atmo(i) = 0
                        cycle dirs
                    end if
                    if (Lmax <= abs(next(k) - minval(chain(k,:)))) then     ! exiting box to the right/top/up
                        atmo(i) = 0
                        cycle dirs
                    end if
                end if
            end if

            ! count visits
            p = get_visits(visits,next)
            ! check if already in chain
            if (count(p > 0) > 0) then
                if (collisionsQ) then
                    ! prevent loops
                    if ( n == Nmax .and. all(next == origin)) then
                        atmo(i) = 0
                        cycle dirs
                    end if

                    ! skip if maximally visited
                    ! NB for hex lattice this means max visits is 1
                    if (count(p > 0) >= coord/2) then
                        atmo(i) = 0
                        cycle dirs
                    end if
                    ! prevent retracing earlier step
                    if (count(q > 0) > 1) then
                        ! loop over past visits
past:                   do j = 1, count(q > 0) - 1     ! Obviously, don't check current point
                            if (q(j) == 1) then     ! special check for origin
                                if (all(next == chain(:,2)) .and. all(this == origin)) then
                                    atmo(i) = 0
                                    cycle dirs
                                end if
                            else
                                if ((all(next == chain(:,q(j)+1)) .or. all(next == chain(:,q(j)-1))) &
                                    .and. all(this == chain(:,q(j))) ) then
                                    atmo(i) = 0
                                    cycle dirs
                                end if
                            end if
                        end do past
                    end if
                else
                    atmo(i) = 0   ! No collisions allowed!
                    cycle dirs
                end if   ! collisionsQ
            end if   ! count(p > 0) > 0
        end do dirs

        ! Grooves: prevent crossings if grooves are enabled
        if (groovesQ) then
            if (count(atmo>0)>0) then
ilp:            do i = 1, coord
                    if (atmo(i) == 0) cycle ilp
                    ! for remaining possible steps
                    next = this + directions(:,atmo(i))
                    ! test branch
                    t(1) = get_dir(this, prev)
                    t(2) = atmo(i)
                    if (t(1) > t(2)) t = t(2:1:-1)  ! sort
jlp:                do j = 1, count(q > 0) - 1
                        if (q(j) == 1) cycle jlp
                        jnext = chain(:,q(j)+1)
                        jprev = chain(:,q(j)-1)
                        ! existing branch
                        e(1) = get_dir(this,jprev)
                        e(2) = get_dir(this,jnext)
                        if (e(1) > e(2)) e = e(2:1:-1)  ! sort
                        if (e(1) < t(1) .and. t(1) < e(2) .and. e(2) < t(2)) then
                            ! Branches cross: reject
                            atmo(i) = 0
                            cycle ilp
                        elseif (t(1) < e(1) .and. e(1) < t(2) .and. t(2) < e(2)) then
                            ! Branches cross: reject
                            atmo(i) = 0
                            cycle ilp
                        end if
                    end do jlp
                end do ilp
            end if
        end if

        ! failsafe
        if (count(atmo > 0) == coord .and. n > 1) then
            print *, this
            print *, next
            print *, atmo
            print *, n
            open(10, file = 'chain.out')
            do i = 1, n
                write(10, *) (chain(j,i), j = 1, lattdim)
            end do
            close(10)
            stop 'Invalid atmo. (buildchain)'
        end if

    end function neighbours

    ! ##############################################
    ! Neighbour check for superSAWs (atmosphere)
    function super_neighbours(chain,visits) result(atmo)
        implicit none
        integer,intent(in)  :: chain(:,:)       ! input
        type(fhash_type__ints_ints),intent(inout) :: visits
        integer             :: atmo(coord)!,atmo2(coord)          ! output
        integer             :: n
        integer             :: next(lattdim),next2next(lattdim),this(lattdim),prev(lattdim),p(coord/2),q(coord/2)
        ! integer             :: jprev(lattdim),jnext(lattdim)
        integer             :: i,j             ! loop variables
        ! integer             :: e(2),t(2)      ! grooves: directions of prev and next for existing and test branches

        n=size(chain,dim=2)
        this=chain(:,n)
        if (collisionsQ .or. .not. superQ) then
            atmo=0
            return
        end if
        atmo = (/ (i,i=1,coord) /)
        if (n>1) then
            prev=chain(:,n-1)
        else
            return
        end if

dir1:   do i=1,coord
            if (lattice_type=='hex') then
                if (i==coord) then
                    if (modulo(sum(this),2)==1) then
                        next=this-directions(:,i)
                    else
                        next=this+directions(:,i)
                    end if
                else
                    next=this+directions(:,i)
                end if
            else
                next=this+directions(:,i)
            end if   ! lattice_type=='hex'

            ! Prevent retracing last step
            if (all( next == prev )) then
                atmo(i)=0
                cycle dir1
            end if

            ! prevent crossing adsorption surface
            if (surfaceQ) then
                if (next(lattdim)<0) then
                    atmo(i)=0
                    cycle dir1
                end if
            end if   ! surfaceQ

            ! count visits
            p=get_visits(visits,next)
            ! check if already in chain
            if ( count(p>0) > 0 ) then
                ! shouldn't ever get to this point if next part works properly
                atmo(i)=0
                cycle dir1
            else
                ! atmo2 = (/ (j,j=1,coord) /)

dir2:           do j=1,coord
                    if (lattice_type=='hex') then
                        if (j==coord) then
                            if (modulo(sum(next),2)==1) then
                                next2next=next-directions(:,j)
                            else
                                next2next=next+directions(:,j)
                            end if
                        else
                            next2next=next+directions(:,j)
                        end if
                    else
                        next2next=next+directions(:,j)
                    end if   ! lattice_type=='hex'

                    ! ignore previous
                    if (all(next2next==this)) then
                        cycle dir2
                    end if

                    ! surface check???
                    if (surfaceQ) then
                        if (next2next(lattdim)<0) cycle dir2
                    end if   ! surfaceQ

                    q=get_visits(visits,next2next)
                    if (count(q>0) > 0) then
                        atmo(i)=0
                        cycle dir1
                    end if

                end do dir2
            end if   ! count(p>0) > 0
        end do dir1

        ! failsafe
        if (count(atmo>0)==coord .and. n>1) then
            print*,this
            print*,next
            print*,atmo
            print*,n
            open(10,file='chain.out')
            do i=1,n
                write(10,*) (chain(j,i),j=1,lattdim)
            end do
            close(10)
            stop 'Invalid atmo. (buildchain)'
        end if

    end function super_neighbours

    ! ##############################################
    ! check atmo for imperfections
    subroutine percolation(last, holes, atmo)
        integer,intent(in)      :: last(:)!,l
        type(fhash_type__ints_ints),intent(inout) :: holes
        integer,intent(inout)   :: atmo(:)
        integer                 :: i, pt(lattdim)

        if (count(atmo > 0) == 0) return

        do i = 1, coord
            if (atmo(i) == 0) cycle

            ! check previous holes
            pt = last + directions(:,atmo(i))
            if (get_holes(holes,pt)) then
                atmo(i) = 0
            end if
        end do

    end subroutine percolation

    ! ##############################################
    ! get direction from point x to point y
    ! could possibly use where() but that can be fiddly with arrays
    function get_dir(x,y) result(d)
        integer,intent(in)  :: x(:),y(:)
        integer             :: d,i

        d=0
        do i=1,coord
            if (all(y-x==directions(:,i))) then
                d=i
                exit
            end if
        end do

        ! failsafe
        if (d==0) then
            print*,'from: ',x
            print*,'to: ',y
            print*,'vector: ',y-x
            print*,'dir. number: ',i
            print*,'dir. vector: ',directions(:,i)
            write(*,'(I4,I4)') (directions(:,d),d=1,coord)
            stop 'Direction cannot be matched'
        end if

    end function get_dir

end module
