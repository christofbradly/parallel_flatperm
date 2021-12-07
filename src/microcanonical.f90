    ! microcanonical.f90
    ! Functions for calculating microcanonical parameters (i.e. number of interactions) of chain
    ! Outputs should be integers
    !
module microcanonical
    ! Don't include 'chain'; avoids conflict with argument names below
    use hash_wrappers
    use globals, only: lattice_type, params, fixedQ, fixed_param, fixed2Q, fixed2_param, numpar, collisionsQ, qp, ext
    use lattice, only: lattdim, coord, directions, origin
    implicit none
    private
    public update_micro, get_local_weight

contains
    ! ##############################################
    ! Get microcanonical weights of chain and return them as a list of indices
    subroutine update_micro(chain, visits, micro)
        implicit none
        integer,intent(in)      :: chain(:,:)       ! input
        type(fhash_type__ints_ints),intent(inout)   :: visits
        integer,intent(inout)   :: micro(3)
        character               :: para
        integer                 :: i, c, ci, di, p,q
        logical                 :: already_checkedQ

        already_checkedQ = .false.
        ! micro(1) = chain_length(chain) + 1
        micro(1) = size(chain, dim = 2)
        ! cx=0

        ! if (scan(params, 'x') /= 0) coll_i = scan(params, 'c')
        ! if (scan(params,'c')/=0) cx=collisions_all(chain)

        do i = 2, numpar
            para = params(i:i)
            select case (para)
            ! case('x')   ! must be done before 'c' in order to use previous step
            !     if (prev_coll > 0) then
            !         micro(i) = micro(i) + crossings_end(chain) ! i.e. 'x'
            !         micro(coll_i) = micro(coll_i) - 1   ! only testing current step so only one new crossing possible
            !         prev_coll = 0     ! probably redundant
            !     endif
            case('c', 'd')
                if (.not.already_checkedQ) then
                    c = collisions_end(chain, visits)   ! returns 0, 1, 2 for number of visits OTHER THAN current endpoint and origin
                    if (c > 0) then
                        ci = scan(params, 'c')
                        di = scan(params, 'd')
                        if (di == 0) then       ! 'd' not present or fixed
                            if (fixedQ .and. fixed_param == 'd') then
                                ! if triply-visited sites are counted separately only count doubly-visited site:
                                ! if (c == 1) m = 1       ! new double encountered
                                ! if (c == 2) m = -1      ! double converted to triple
                                ! else m = 0
                                micro(ci) = micro(ci) + (-1)**(c + 1)        ! one line
                            else
                                ! micro(i) = micro(i) + c         ! count number of visits
                                if (c == 1) micro(i) = micro(i) + 1         ! count ONLY double visits; equivalent to allowing but not weighting triple visits
                                ! micro(i) = micro(i) + c*(c + 1)/2   ! count no. of 2-body interactions
                            end if
                        else if (ci == 0) then      ! 'c' not present or fixed
                            if (c == 2) micro(i) = micro(i) + 1
                        else                ! both 'c' and 'd' present
                            if (c == 1) then        ! double visit
                                micro(ci) = micro(ci) + 1
                            elseif (c == 2) then    ! triple visit
                                micro(di) = micro(di) + 1
                                micro(ci) = micro(ci) - 1   ! replace 2-body interaction w/ 3-body interaction
                            else
                                ! c == 0; no collision
                            end if
                        endif
                        already_checkedQ = .true.     ! flag to avoid calling collisions_end twice
                    end if
                end if
            case('s')
                micro(i) = micro(i) + stiffs_end(chain)
                ! micro(i) = stiffs_all(chain) + 1
            case('m')
                micro(i) = micro(i) + near_neighbours_end(chain, visits)
                ! micro(i) = near_neighbours_all(chain) + 1
            case('a')
                micro(i) = micro(i) + adsorbs_end(chain)
            case('h')
                micro(i) = chain(lattdim,micro(1)) + 1
                ! micro(i) = chain(lattdim,floor(alpha*(micro(1) - 1)) + 1) + 1
            case('g')   ! generalisation of 'h' to any direction w/o a surface
                micro(i) = maxval(abs(chain(:,micro(1)))) + 1
            case('l')
                micro(i) = box_size(chain)
            case('e')
                if (end_on_bounding_box(chain)) then
                    micro(i) = 2        ! walk ends on bounding box
                else
                    micro(i) = 1        ! walk ends in bulk
                end if
            end select
        enddo

        if (any(micro(1:numpar) .le. 0) .or. any(micro(2:numpar) > ext(2:numpar))) then
            open(10, file = 'chain.out')
            do p = 1, micro(1)
                write(10, *) (chain(q,p), q = 1, lattdim)
            enddo
            close(10)
            print*, 'micro:', micro
            print*, 'ext:  ', ext
            print*, 'Invalid microcanonical weights (step).'
            stop -1
        endif

    end subroutine update_micro

    ! ##############################################
    ! Get weight that is only used locally, not as part of flatPERM algorithm
    subroutine get_local_weight(chain, visits, micro, m, m2)
        implicit none
        integer,intent(in)      :: chain(:,:)       ! input
        type(fhash_type__ints_ints),intent(inout)   :: visits
        integer,intent(in)      :: micro(3)
        integer,intent(out)     :: m, m2
        integer                 :: c

        m = 0
        m2 = 0
        if (.not.fixedQ) return     ! safety check

        select case (fixed_param)
        case('c')
            c = collisions_end(chain, visits)   ! returns 0, 1, 2 for number of visits OTHER THAN current endpoint and origin

            if (scan(params, 'd') > 0) then
                ! if triply-visited sites are counted separately only count doubly-visited site:
                ! if (c == 1) m = 1       ! new double encountered
                ! if (c == 2) m = -1      ! double converted to triple
                ! else m = 0
                if (c > 0) m = (-1)**(c + 1)        ! one line
            else if (fixed2_param == 'd') then      ! do this all in one go, probably fails if c and d are reversed
                if (c > 0) m = (-1)**(c + 1)
                if (c == 2) m2 = 1
                return      ! don't need to do d below
            else
                m = c
                ! m = c*(c + 1)/2     ! count no. of 2-body interactions ([0,1,2] --> [0,1,3])
            end if
        case('d')
            if (collisions_end(chain, visits) == 2) m = 1   ! specifically count triply-visited sites here
        ! case('x')
        !     stop 'I probably do not know how to fix crossings weight'
        case('s')
            m = stiffs_end(chain)
        case('m')
            m = near_neighbours_end(chain, visits)
        case('a')
            m = adsorbs_end(chain)
        case('h')
            m = chain(lattdim,micro(1)) + 1
            ! m = chain(lattdim,floor(alpha*(micro(1)-1)) + 1) + 1
        case default
            stop 'Undefined fixed parameter: '//fixed_param
        end select

        if (.not.fixed2Q) return     ! safety check

        select case (fixed2_param)
        case('c')
            c = collisions_end(chain, visits)   ! returns 0, 1, 2 for number of visits OTHER THAN current endpoint and origin

            if (scan(params, 'd') > 0) then
                ! if triply-visited sites are counted separately only count doubly-visited site:
                ! if (c == 1) m = 1       ! new double encountered
                ! if (c == 2) m = -1      ! double converted to triple
                ! else m = 0
                if (c > 0) m = (-1)**(c + 1)        ! one line
            else
                m = c
                ! m = c*(c + 1)/2     ! count no. of 2-body interactions ([0,1,2] --> [0,1,3])
            end if
        case('d')
            if (collisions_end(chain, visits) == 2) m2 = 1   ! specifically count triply-visited sites here
        ! case('x')
        !     stop 'I probably do not know how to fix crossings weight'
        case('s')
            m2 = stiffs_end(chain)
        case('m')
            m2 = near_neighbours_end(chain, visits)
        case('a')
            m2 = adsorbs_end(chain)
        case('h')
            m2 = chain(lattdim,micro(1)) + 1
            ! m = chain(lattdim,floor(alpha*(micro(1)-1)) + 1) + 1
        case default
            stop 'Undefined fixed parameter: '//fixed_param
        end select

    end subroutine get_local_weight

    ! ##############################################
    ! Length of chain
    function chain_length(chain) result(m)
        implicit none
        integer,intent(in)  :: chain(:,:)       ! input
        integer             :: m                ! output

        m = size(chain, dim = 2) - 1

    end function chain_length

    ! ##############################################
    ! Number of corner collisions (for VISAWs)
    ! aka multiply-visited sites (for ISATs)
    function collisions_end(chain, visits) result(c)
        implicit none
        integer,intent(in)  :: chain(:,:)       ! input
        type(fhash_type__ints_ints)   :: visits
        integer             :: last(lattdim), c

        last = chain(:,size(chain, dim = 2))

        ! visits method
        c = count_visits(visits, last) - 1  ! do not count current point
        if (all(last == origin)) c = c - 1 ! do not count origin;

    end function collisions_end

    ! ##############################################
    ! Number of crossings
    ! NB this routine is testing the previous step but input is full chain!
    ! Currently only valid for 'squ' and 'tri' lattices
    function crossings_end(chain) result(m)
        implicit none
        integer,intent(in)  :: chain(:,:)       ! input
        type(fhash_type__ints_ints)   :: visits
        integer             :: m         ! output
        integer             :: n, j, t(2), e(2)
        integer             :: this(lattdim), prev(lattdim), next(lattdim)
        integer             :: q(coord/2)
        logical             :: good_cross = .false.

        n = size(chain, dim = 2)
        m = 0

        this = chain(:,n-1)
        next = chain(:,n)
        prev = chain(:,n-2)
        t(1) = get_dir(this,prev)
        t(2) = get_dir(this,next)
        if (t(1) > t(2)) t = t(2:1:-1)
        ! Find points in chain that are visited more than once
        q = get_visits(visits, this)
lpi:    do j = 1, count(q > 0) - 1!2,n-2
            if (q(j) == 1) cycle lpi
                ! jnext=chain(:,j+1)
                ! jprev=chain(:,j-1)
                e(1) = get_dir(this,chain(:,q(j)-1))
                e(2) = get_dir(this,chain(:,q(j)+1))
                if (e(1) > e(2)) e = e(2:1:-1)
                if (e(1) < t(1) .and. t(1) < e(2) .and. e(2) < t(2)) then
                    ! Branches cross: increment
                    good_cross = .true.
                    cycle lpi
                elseif (t(1) < e(1) .and. e(1) < t(2) .and. t(2) < e(2)) then
                    ! Branches cross: increment
                    good_cross = .true.
                    cycle lpi
                else
                    ! Branches do not cross
                    good_cross = .false.
                endif
        enddo lpi

        if (good_cross) m = 1


    end function crossings_end

    ! ##############################################
    ! get direction from point x to point y
    ! could possibly use where() but that can be fiddly with arrays
    function get_dir(x, y) result(d)
        integer,intent(in)  :: x(:), y(:)
        integer             :: d, i

        d = 0
        do i = 1,coord
            if (all(y - x == directions(:,i))) then
                d = i
                exit
            endif
        enddo

        ! failsafe
        if (d == 0) then
            print*, 'from: ', x
            print*, 'to: ', y
            print*, 'vector: ', y - x
            print*, 'dir. number: ', i
            print*, 'dir. vector: ', directions(:,i)
            write(*, '(I4,I4)') (directions(:,d), d = 1 , coord)
            stop 'Direction cannot be matched'
        endif

    end function get_dir


    ! ##############################################
    ! Number of monomers where preceeding step is parallel to next step
    function stiffs_end(chain) result(m)
        implicit none
        integer,intent(in)  :: chain(:,:)       ! input
        integer             :: m                ! output
        integer             :: n
        integer             :: stiff(lattdim)

        m = 0
        n = size(chain, dim = 2)
        if (n < 3) return

        stiff = 2*chain(:,n - 1)-chain(:,n - 2)-chain(:,n)

        if (all( origin == stiff )) m = 1

    end function stiffs_end

    ! ##############################################
    ! Number of adjacent non-consecutive monomers -- check only end; O(1) time
    ! Update 8/9/2017: New routine uses hash table 'visits'. Decrements are handled
    ! correctly now for all lattice types and topologies. Adjacent points can
    ! produce multiple interaction counts if either or both points are multiply visited.
    function near_neighbours_end(chain, visits) result(m)
        implicit none
        integer,intent(in)  :: chain(:,:)       ! input
        type(fhash_type__ints_ints)   :: visits
        integer             :: m                ! output
        integer             :: i, k, n
        integer             :: this(lattdim), prev(lattdim), next(lattdim)
        integer             :: p(coord/2), q(coord/2), this_vis, next_vis
        logical             :: bond

        m = 0
        n = size(chain, dim = 2)
        this = chain(:,n)
        if (n < 2) then
            return
        else
            prev = chain(:,n-1)
        endif
        q = get_visits(visits, this)
        this_vis = count(q > 0)
        ! Decrement for newly formed collisions
        if (collisionsQ) then
            if (all(this == origin)) then
                m = m - (this_vis - 2)
            else
                if ( this_vis > 1 .and. all(prev /= origin) ) then
                    m = m - (this_vis-2 + count(get_visits(visits,prev) > 0))
                else
                    m = m - (this_vis - 1)
                endif
            endif
        endif

dirs:   do i = 1,coord
            bond = .false.
            ! pick adjacent point
            if (lattice_type == 'hex') then
                if (i == coord) then
                    if (modulo(sum(this),2) == 1) then
                        next = this - directions(:,i)
                    else
                        next = this + directions(:,i)
                    endif
                else
                    next = this + directions(:,i)
                endif
            else
                next = this + directions(:,i)
            endif   ! lattice_type=='hex'

            ! easy discard
            if (all(next == prev)) cycle dirs

            p = get_visits(visits,next)
            next_vis = count(p > 0)
            ! check visits of adjacent point
            if (next_vis > 0) then
                if (collisionsQ) then
klp:                do k = 1, next_vis     ! loop over adjacent point visits
                        ! check if a bond exists
                        if (all(this == chain(:,p(k)+1))) then
                            bond = .true.
                        elseif (p(k) > 1) then
                            ! special case for origin
                            if (all(this == chain(:,p(k)-1))) then
                                bond = .true.
                            endif
                        endif
                        if (bond) exit klp
                    enddo klp
                endif   ! implicit else: bond cannot exist for SAWs

                if (.not.bond) then
                    ! special case: no interaction with origin step
                    if (all(next == origin)) then
                        m = m + next_vis - 1
                    else
                        m = m + next_vis
                    endif
                endif
            endif   ! implicit else: nothing at adjacent point

        enddo dirs

        return

    end function near_neighbours_end

    ! ##############################################
    ! Number of monomers incident on surface
    ! Ignores origin
    function adsorbs_end(chain) result(m)
        implicit none
        integer,intent(in)  :: chain(:,:)       ! input
        integer             :: m                ! output
        integer             :: n, last(lattdim)

        m = 0
        n = size(chain, dim = 2)
        if (n < 2) return

        last = chain(:,n)
        if (lattice_type.eq.'hex') then
            if ( last(lattdim) == 0 ) then
                ! Horizontal surface as defined by Batchelor Bennett-Wood Owczarek EPJB 1998 (Fig. 3)
                if (modulo(last(1), 2) == 1) m = 1
            endif
        else
            if ( last(lattdim) == 0 ) m = 1
        endif
    end function adsorbs_end

    ! ##############################################
    ! Size of smallest containing square box
    function box_size(chain) result(m)
        implicit none
        integer,intent(in)  :: chain(:,:)       ! input
        integer             :: m                ! output
        integer             :: i, n, w(lattdim)

        m = 1
        n = size(chain, dim = 2)
        if (n < 2) return

        do i = 1, lattdim
            w(i) = maxval(chain(i,:)) - minval(chain(i,:)) + 1
        end do

        m = maxval(w)
    end function box_size

    ! ##############################################
    ! Test if the last point of chain lies on the bounding box
    function end_on_bounding_box(chain) result(m)
        implicit none
        integer,intent(in)  :: chain(:,:)       ! input
        logical             :: m                ! output
        integer             :: i, n

        m = .false.
        n = size(chain, dim = 2)
        if (n < 2) return

        ! corners of bounding box
        if (any(chain(:,n) == minval(chain, dim = 2))) then
            m = .true.
            return
        end if
        if (any(chain(:,n) == maxval(chain, dim = 2))) then
            m = .true.
            return
        end if
        if (any(chain(:,1) == minval(chain, dim = 2))) then
            m = .true.
            return
        end if
        if (any(chain(:,1) == maxval(chain, dim = 2))) then
            m = .true.
            return
        end if
    end function end_on_bounding_box

! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
! # # # # unused  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    ! ! ##############################################
    ! ! Number of adjacent non-consecutive monomers -- check only end; O(N) time
    ! ! Update 20/7/2017: New routine only counts 'gaps' as one interaction. Cannot have
    ! ! interactions w/ origin UNLESS there is another step crossing the origin; the latter
    ! ! step can produce interactions. Routine also correctly DECREMENTS count if
    ! ! new collision/crossing removes a gap.
    ! function near_neighbours_end(chain) result(m)
    !     implicit none
    !     integer,intent(in)  :: chain(:,:)       ! input
    !     integer             :: m                ! output
    !     integer             :: latticescope,i,j,n,these
    !     integer             :: this(lattdim),prev(lattdim),dir(lattdim),atm(coord)
    !     Logical             :: skip_dbl
    !
    !     n=size(chain,dim=2)
    !     m=0
    !     latticescope=2      ! Jump to prevent detections between adjacent steps
    !     if (any(['squ','cub','hex']==lattice_type)) latticescope=3
    !     if ( n<=latticescope ) then
    !         return
    !     endif
    !
    !     this=chain(:,n)
    !     prev=chain(:,n-1)
    !     ! Setup neighbours
    !     do i=1,coord
    !         dir=directions(:,i)
    !         if (all(prev-this==dir)) then   ! avoid previous step
    !             atm(i)=0
    !         else
    !             atm(i)=i
    !         endif
    !     enddo
    !
    !     these=1     ! track no. of points coincident to current point
    !     ! Decrement for new-formed collisions -- once only; do not count origin
    !     if (collisionsQ) then
    !         do i=2,n-1
    !             if (all( chain(:,i) == this )) then
    !                 if ( these == 1 ) then
    !                     m=m-1
    !                 endif
    !                 these=these+1
    !             endif
    !         enddo
    !     endif
    !
    !     ! Check neighbours -- increment once only
    ! nbs:    do j=1,coord
    !         if (atm(j)==0) cycle nbs
    !         skip_dbl=.true.   ! flag to prevent edge cases when collision is at/adjacent to origin
    !         ! loop over chain
    ! pts:        do i=2,n-2
    !             if (i==2) then  ! check to toggle skip_dbl
    !                 if (collisionsQ) then
    !                     if (all(this==origin)) skip_dbl=.false.
    !                 endif
    !             endif
    !
    !             if (all( chain(:,i) == this + directions(:,atm(j)) )) then  ! standard adjacency test
    !                 if (.not.(all( chain(:,i-1) == this ) .or. all( chain(:,i+1) == this ))) then  ! consecutivity test
    !                     if (these==1) then  ! do not increment if bond counted previously
    !                         if (.not.skip_dbl .and.all( chain(:,i)==chain(:,2) )) then
    !                             cycle pts
    !                         else
    !                             m=m+1
    !                             exit pts    ! increment once only
    !                         endif
    !                     endif
    !                 endif
    !             endif
    !         enddo pts
    !     enddo nbs
    !
    !     return
    !
    ! end function near_neighbours_end

!     ! ##############################################
!     ! Number of corner collisions (for ISAWs)
!     ! aka multiply-visited sites (for ISATs)
!     ! Currently only valid for 'squ' lattice
!     ! *** ALSO INCORRECT W/ REGARDS TO CROSSINGS ***
!     function collisions_all(chain) result(m)
!         implicit none
!         integer,intent(in)  :: chain(:,:)       ! input
!         integer             :: m(2)         ! output
!         integer             :: n
!         integer             :: stiff(lattdim),c,x,i,j,k,h
!
!         n=size(chain,dim=2)
!         c=0
!         x=0
!         h=0
!         if ( endcollQ ) h=1
!
!         ! Find points in chain that are visited more than once
! lpi:    do i=1,n
! lpj:        do j=i+1,n
!                 if (all( chain(:,i) == chain(:,j) )) then
!                     c=c+1
!                     ! Determine which multiply visited sites are crossings
!                     if ( crossingsQ ) then
!                         ! only needs to test first step at multiply visited site
!                         stiff=chain(:,i+1)-chain(:,i-1)
! dirs:                   do k=1,coord
!                             if (all( stiff == 2*directions(:,k) )) then
!                                 x=x+1
!                                 exit dirs
!                             endif
!                         enddo dirs
!                     endif
!                 endif
!             enddo lpj
!         enddo lpi
!
!         m=[c-x,x]
!
!     end function collisions_all
!
!     ! ##############################################
!     ! Number of adjacent non-consecutive monomers -- entire chain; O(N^2) time
!     ! # # # # # # # # # # # # #
!     ! # # # # unused  # # # # #
!     ! # # # # # # # # # # # # #
!     function near_neighbours_all(chain) result(m)
!         implicit none
!         integer,intent(in)  :: chain(:,:)       ! input
!         integer             :: m                ! output
!         integer             :: latticescope,i,j,k,h,a,b,n
!         integer             :: nnsteps(2,size(chain,2)),diff(lattdim),posa(coord),posb(coord)
!
!         n=size(chain,dim=2)
!         latticescope=2      ! Jump to prevent detections between adjacent steps
!         if (any(['squ','cub','hex']==lattice_type)) latticescope=3
!         if ( n<=latticescope ) then
!             m=0
!             return
!         endif
!
!         ! find pairs of adjacent (non-consecutive) sites
!         nnsteps=0
!         k=0
! ilp:    do i=1,n-1      ! this could be parallelised
! jlp:        do j=i+latticescope,n
!                 diff=chain(:,i)-chain(:,j)
! dirs:           do h=1,coord
!                     if (all( diff == directions(:,h) )) then
!                         k=k+1
!                         if (collisionsQ) nnsteps(:,k)=[i,j]
!                         exit dirs
!                     endif
!                 enddo dirs
!             enddo jlp
!         enddo ilp
!
!         ! remove those pairs that are actually two adjacent collisions/crossings
!         if ( k /= 0 .and. collisionsQ ) then
!             posa=0
!             posb=0
!             m=k
! pairs:      do i=1,k    ! this could be parallelised
!                 a=0
!                 b=0
!                 ! find positions in chain of detected pair sites
! pos:            do j=1,n
!                     if (all( chain(:,j) == chain(:,nnsteps(1,i)) )) then
!                         a=a+1
!                         posa(a)=j
!                     elseif (all( chain(:,j) == chain(:,nnsteps(2,i)) )) then
!                         b=b+1
!                         posb(b)=j
!                     endif
!                 enddo pos
!                 !
!                 if ( a==0 .or. b==0 ) cycle pairs
!
!                 ! Test neighbourhoods of sites for consecutivity
! consec:         do j=1,a
!                     if ( any(posa(j)-1 == posb(1:b)) .or. any(posa(j)+1 == posb(1:b)) ) then
!                         m=m-1
!                     endif
!                 enddo consec
!             enddo pairs
!         else
!             m=k
!         endif
!
!     end function near_neighbours_all
!
!     ! ##############################################
!     ! Number of monomers where preceeding step is parallel to next step
!     ! # # # # # # # # # # # # #
!     ! # # # # unused  # # # # #
!     ! # # # # # # # # # # # # #
!     function stiffs_all(chain) result(m)
!         implicit none
!         integer,intent(in)  :: chain(:,:)       ! input
!         integer             :: m                ! output
!         integer             :: n,i
!         integer             :: stiff(lattdim)
!
!         m=0
!         n=size(chain,dim=2)
!         if (n<3) return
!
!         do i=2,n-1
!             stiff=2*chain(:,i)-chain(:,i-1)-chain(:,i+1)
!             if (all( origin == stiff )) then
!                 m=m+1
!             endif
!         enddo
!
!         ! This method might make a copy?
!         !stiffs=0
!         !where (2*chain(:,2:n-1)-chain(:,1:n-2)-chain(:,3:n) == origin(:))
!         !    stiffs = 1
!         !end where
!         !m=sum(stiffs)
!
!     end function stiffs_all

end module
