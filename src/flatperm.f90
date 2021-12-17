    ! flatperm.f90
    ! Main subroutines for running flatPERM algorithm
    !
module flatperm
    use globals
    ! use iso_fortran_env
    use mt_stream
    use hash_wrappers
    use lattice, only: coord, directions, lattdim, origin
    use microcanonical, only: update_micro, get_local_weight
    implicit none
    private
    public run_flatperm

contains
    ! ##############################################
    ! Run the algorithm for a single tour
    subroutine run_flatperm(t, mts)
        use buildchain, only: neighbours, super_neighbours, percolation
        implicit none
        integer(8),intent(in)   :: t
        type(mt_state),intent(inout)    :: mts
        type(fhash_type__ints_ints)     :: visits, holes
        integer     :: chain(lattdim,ext(1)), history_enr(numpar,ext(1)), history_atmo(coord,ext(1)), &
            history_steps(ext(1)), history_fixed(ext(1)), history_fixed2(ext(1))
        real(qp)                :: history_W(ext(1))
        integer                 :: micro(3), i_hist, fixm_total, fix2m_total
        real(qp)                :: W
        integer                 :: atmo(coord), prevenr, copies

        ! reset history
        i_hist = 0
        history_enr = 0
        if (fixedQ) history_fixed = 0
        if (fixedQ) fixm_total = 0
        if (fixed2Q) history_fixed2 = 0
        if (fixed2Q) fix2m_total = 0
        history_atmo = 0
        history_steps = 0
        history_W = 0.0_qp
        ! initialise state
        chain = 0
        micro = 1
        call reset_visits(visits)
        if (percQ) then
            call reset_holes(holes, lattdim, W)
            !$omp atomic update
                normq = normq + W
        else
            W = 1.0_qp
        end if
        call add_visit(visits, origin, micro(1))
        ! add trivial state to data arrays
        ! W = 1.0_qp
        !$omp atomic update
            Z(micro(1),micro(2),micro(3)) = Z(micro(1),micro(2),micro(3)) + W
        !$omp atomic update
            samples(micro(1),micro(2),micro(3)) = samples(micro(1),micro(2),micro(3)) + 1_8
        !$omp atomic update
            effsamples(micro(1),micro(2),micro(3)) = effsamples(micro(1),micro(2),micro(3)) + 1.0_qp


        ! loop within a tour
tour:   do while (.not. interruptQ)

            ! get atmosphere
            if (superQ) then
                atmo = super_neighbours(chain(:,1:micro(1)),visits)
            else
                atmo = neighbours(chain(:,1:micro(1)),visits)
                if (percQ) call percolation(chain(:,micro(1)), holes, atmo)
            endif

            ! Prune/enrich: Get copies
            call get_copies(micro, atmo, t, W, mts, copies)
            ! calculate R&R weight
            W = W*count(atmo > 0)/mu

            ! Prune/enrich: Commit copies
            if ( copies == 0 ) then
                !$omp atomic update
                    prunes(micro(1),micro(2),micro(3)) = prunes(micro(1),micro(2),micro(3)) + 1_8

                ! Check if tour is complete
                if ( i_hist == 0 ) then
                    prevenr = 0
                    exit tour
                else
                    ! get previous enrich point
                    prevenr = history_enr(1,i_hist)
                endif
                ! shrink chain to previous enrichment
shrink:         do while (micro(1) > prevenr)
                    call unstep(chain, visits, micro, i_hist, history_enr)
                enddo shrink
                ! call update_micro(chain(:,1:micro(1)),0,micro)    ! updated later
                if (any( micro .le. 0 )) then
                    write(*,*) 'Invalid microcanonical weights (unstep).'
                    write(*,*) micro
                    stop
                endif
                ! restore some state variables to earlier state
                W = history_W(i_hist)
                if (fixedQ) fixm_total = history_fixed(i_hist)
                if (fixed2Q) fix2m_total = history_fixed2(i_hist)
            else
                ! record copies
                if (count(atmo > 0) == 0 ) stop 'Cannot Enrich if atmosphere is empty. (record)'
                if ( copies > 1 ) then
                    !$omp atomic update
                        enriches(micro(1),micro(2),micro(3)) = enriches(micro(1),micro(2),micro(3)) + int(copies - 1, 8)
                endif
                ! update history
                i_hist = i_hist + 1
                history_enr(:,i_hist) = micro(1:numpar)
                if (fixedQ) history_fixed(i_hist) = fixm_total
                if (fixed2Q) history_fixed2(i_hist) = fix2m_total
                history_atmo(:,i_hist) = shuffle_atmo(copies, atmo, mts)
                history_W(i_hist) = W
            endif
            ! commit next step
            call step(W, chain, visits, micro, i_hist, fixm_total, fix2m_total, history_enr, history_atmo, &
                history_steps, history_fixed, history_fixed2, history_W)

        end do tour

        ! deallocate hash tables
        call clear_hash_table(visits)
        if (percQ) call clear_hash_table(holes)

    end subroutine run_flatperm

    ! ##############################################
    ! Shuffle the atmosphere - Knuth shuffle
    function shuffle_atmo(c, vec_in, mts) result(vec_out)
        implicit none
        integer,intent(in)      :: c
        integer,intent(in)      :: vec_in(:)
        type(mt_state),intent(inout)    :: mts
        integer                 :: vec_out(coord)
        integer                 :: a, i, j, temp(count(vec_in > 0))
        real(dp)                :: b

        temp = pack(vec_in, vec_in > 0)
        do i=1,size(temp)-1
            ! Generate random integer in range [1,...,size of atmosphere]
            ! b=genrand64_real1()   ! MT RNG
            b = genrand_double1(mts)    ! MT RNG for multiple streams
            ! call random_integer(i,size(temp),j)   ! wrapper for gcc RNG
            j = i + floor((size(temp) + 1 - i)*b)
            a = temp(j)
            temp(j) = temp(i)
            temp(i) = a
        enddo

        vec_out = 0
        vec_out(1:min(c, size(temp))) = temp(1:min(c, size(temp)))

    end function shuffle_atmo


    ! ##############################################
    ! Commit next step and update weights/thermo
    subroutine step(weight, chain, visits, micro, i_hist, fixm_total, fix2m_total, history_enr, history_atmo, &
                    history_steps, history_fixed, history_fixed2, history_W)
        use radius
        use lattice
        implicit none
        real(qp),intent(inout)  :: weight, history_W(:)
        integer,intent(inout)   :: chain(:,:), history_enr(:,:), history_atmo(:,:), &
            history_steps(:), history_fixed(:), history_fixed2(:)
        integer,intent(inout)   :: micro(3), i_hist, fixm_total, fix2m_total
        type(fhash_type__ints_ints)   :: visits
        integer                 :: n_ind, p, q, m, m2
        real(qp)                :: add_term

        !if ( micro(1) > 1 ) then
            ! These checks are redundant in principle, but worth having just in case
            if ( i_hist == 0 ) stop 'Cannot Step if history is empty.'
            if (all( history_atmo(:,i_hist) == 0 )) stop 'Cannot Step if atmosphere is empty.'
        !endif

        ! Add step to chain
        call add_to_chain(micro(1))    ! Takes current  length
        ! Update micro
        call update_micro(chain(:,1:micro(1)+1), visits, micro)

        ! Update weight if fixed local weight is set
        if (fixedQ) then
            call get_local_weight(chain(:,1:micro(1)), visits, micro, m, m2)
            if ( m /= 0 ) weight = weight * fixed_weight**m
            if ( m2 /= 0 ) weight = weight * fixed2_weight**m2
            fixm_total = fixm_total + m
            add_term = weight * fixm_total
            !$omp atomic update
                FixMom1(micro(1),micro(2),micro(3)) = FixMom1(micro(1),micro(2),micro(3)) + add_term
            add_term = weight * fixm_total**2
            !$omp atomic update
                FixMom2(micro(1),micro(2),micro(3)) = FixMom2(micro(1),micro(2),micro(3)) + add_term
            if (fixed2Q) then
                fix2m_total = fix2m_total + m2
                add_term = weight * fix2m_total
                !$omp atomic update
                    Fix2Mom1(micro(1),micro(2),micro(3)) = Fix2Mom1(micro(1),micro(2),micro(3)) + add_term
                add_term = weight * fix2m_total**2
                !$omp atomic update
                    Fix2Mom2(micro(1),micro(2),micro(3)) = Fix2Mom2(micro(1),micro(2),micro(3)) + add_term
            endif
        endif

        ! Update thermo/radius data arrays
        !$omp atomic update
            Z(micro(1),micro(2),micro(3)) = Z(micro(1),micro(2),micro(3)) + weight
        !$omp atomic update
            samples(micro(1),micro(2),micro(3)) = samples(micro(1),micro(2),micro(3)) + 1_8
        if ( i_hist == 0 ) then
            n_ind=0
        else
            ! Number of 'independent' steps
            n_ind = history_enr(1,i_hist) - 1
        endif
        add_term = 1.0_qp - real(n_ind, kind = qp)/(micro(1) - 1)
        !$omp atomic update
            effsamples(micro(1),micro(2),micro(3)) = effsamples(micro(1),micro(2),micro(3)) + add_term
        if (surfaceQ) then
            p = lattdim - 1
            q = lattdim
            !transverse
            add_term = weight*end_to_end(chain(1:p,1:micro(1)))
            !$omp atomic update
                Re2W(micro(1),micro(2),micro(3)) = Re2W(micro(1),micro(2),micro(3)) + add_term
            ! perpendicular
            add_term = weight*end_to_end(chain(q:q,1:micro(1)))
            !$omp atomic update
                Re2Wp(micro(1),micro(2),micro(3)) = Re2Wp(micro(1),micro(2),micro(3)) + add_term
        else
            p = lattdim
            ! Normal
            add_term = weight*end_to_end(chain(1:p,1:micro(1)))
            !$omp atomic update
                Re2W(micro(1),micro(2),micro(3)) = Re2W(micro(1),micro(2),micro(3)) + add_term
        endif

        ! anisotropy
        if (anisoQ) then
            add_term = weight*anisotropy_step(history_steps(2:micro(1)), coord)
            !$omp atomic update
                aniso(micro(1),micro(2),micro(3)) = aniso(micro(1),micro(2),micro(3)) + add_term
            if (surfaceQ .and. lattice_type=='cub') then
                add_term = weight*anisotropy_step_2(history_steps(2:micro(1)), coord)
                !$omp atomic update
                    aniso2(micro(1),micro(2),micro(3)) = aniso2(micro(1),micro(2),micro(3)) + add_term
            endif
        endif

        ! Save walks
        if (bestQ .and. micro(1) > Nmax) then
            if (weight > bestWeights(micro(2),micro(3))) then
                !$omp atomic write
                    bestWeights(micro(2),micro(3)) = weight
                !$omp critical(write_bestWalks)
                    bestWalks(:,:,micro(2),micro(3)) = chain
                !$omp end critical(write_bestWalks)
            end if
        end if

    contains
        ! Append a point to the chain
        subroutine add_to_chain(L)
            use lattice
            implicit none
            integer,intent(in)      :: L    ! L is length of chain NOT INCLUDING ADDED STEP
            integer                 :: i, j, d(lattdim)

            ! Choose step from atmosphere (already shuffled)
            ! minloc should find the first 0 in (shuffled) atmo, so this chooses the last direction in the list
            ! Caution: minloc is a tricky function
            if (count(history_atmo(:,i_hist) > 0) == size(history_atmo, 1)) then
                i = size(history_atmo,1)
            else
                i = minloc(history_atmo(:,i_hist),1) - 1
            endif
            ! i==0 means atmo is all zeroes, i.e. no valid directions to choose
            if ( i == 0 .and. micro(1) > 1) stop 'Cannot step if atmo is empty.'

            j = history_atmo(i,i_hist)
            d = directions(:,j)

            ! choose random point from atmosphere for next step
            if (lattice_type == 'hex') then
                ! if (i==coord) then
                if (history_atmo(i,i_hist) == coord) then
                    if (modulo(sum(chain(:,L)),2) == 1) then
                        chain(:,L+1) = chain(:,L) - d
                    else
                        chain(:,L+1) = chain(:,L) + d
                    endif
                else
                    chain(:,L+1) = chain(:,L) + d
                endif
            else
                chain(:,L+1) = chain(:,L) + d
            endif
            ! update visits hash table
            call add_visit(visits, chain(:,L+1), L+1)
            ! record step taken
            history_steps(L+1) = j

            ! remove chosen point from atmosphere
            history_atmo(i,i_hist) = 0
            ! If atmosphere is empty, clean history
            if (all( history_atmo(:,i_hist) == 0 )) then
                history_W(i_hist) = 0.0_qp
                history_enr(:,i_hist) = 0
                if (fixedQ) history_fixed(i_hist) = 0
                if (fixed2Q) history_fixed2(i_hist) = 0
                i_hist = i_hist - 1
            endif

        end subroutine add_to_chain

    end subroutine step

    ! ##############################################
    ! Remove step from chain
    ! NB: history is already cleaned by subroutine step()
    subroutine unstep(chain, visits, micro, i_hist, history_enr)
        implicit none
        integer,intent(inout)       :: chain(:,:), history_enr(:,:)
        integer,intent(inout)       :: micro(3), i_hist
        type(fhash_type__ints_ints),intent(inout)   :: visits
        integer     :: i, old_n
        logical     :: err

        ! Don't need to recalculate all microcanonical weights until shrinking is finished
        old_n = micro(1)
        micro(1:numpar) = history_enr(:,i_hist)
        
        ! remove visits; need loop in case jump back is more than one step
        do i = old_n, micro(1) + 1, -1
            ! remove visit from hash table
            call delete_last_visit(visits,chain(:,i), err)
        enddo
        ! remove points from chain
        chain(:,micro(1)+1:old_n) = 0

    end subroutine unstep

    ! ##############################################
    ! Get number of prune/enrich copies with probabilistic ratio
    subroutine get_copies(micro, atm, S, weight, mts, c)
        implicit none
        integer,intent(in)      :: micro(3), atm(:)
        integer(8),intent(in)   :: S
        integer,intent(out)     :: c
        real(qp),intent(inout)  :: weight
        type(mt_state),intent(inout)    :: mts
        real(qp)                :: ratio, Zval, effval, eff_incr
        real(dp)                :: b
        integer                 :: a, i
        integer(8)              :: thread_delay
        !$ integer              :: omp_get_thread_num

        a = count(atm > 0)

        if (fixedQ .and. numpar == 2 .and. scan(params, fixed_param) /= 0) then
            !$omp critical(read_array_Zval)
                Zval = sum(Z(micro(1),1:ext(2),1))
            !$omp end critical(read_array_Zval)
            effval = 0._qp
            do i = 1, ext(2)
                !$omp atomic read
                    eff_incr = effsamples(micro(1),i,1)
                effval = effval + eff_incr
            enddo
        else
            !$omp atomic read
                Zval = Z(micro(1),micro(2),micro(3))
            !$omp atomic read
                effval = effsamples(micro(1),micro(2),micro(3))
        endif

        ! delay factor should also account for OMP parallel loop scheduling
        thread_delay = 0_8
        !$ if (.not.resumeQ .or. OldTours < delay*(Nmax+1)) thread_delay = omp_get_thread_num()*NewTours/ntasks

        if ( micro(1) - 1 < Nmax .and. a > 0 .and. delay*(micro(1) - 1) < S - thread_delay ) then
            ratio = (weight*(S - thread_delay - floor(delay*(micro(1) - 1), 8))**2) / Zval / effval
            ! ratio = weight*(S-floor(delay*(micro(1)-1),8)) / Zval
            if ( ratio < 1.0_qp ) then
                ! prune first
                ! b=genrand64_real1()
                b = genrand_double1(mts)
                if ( b < ratio ) then
                    c = 1
                    weight = weight/ratio
                else
                    c = 0
                    weight = 0.0_qp
                endif
            else
                ! enriches
                ! if ( probenrQ ) then
                !     b=genrand64_real1()
                !     if ( b>ratio-floor(ratio,8) ) then
                !         c=floor(ratio,8)+1
                !     else
                !         c=floor(ratio,8)
                !     endif
                !     copies=min(a, c)
                ! else
                    ! standard method
                    c = int(min(int(a, 8), floor(ratio, 8)), 4)
                ! endif
                weight = weight/c
            endif
        else
            ! Ensure safe fail state
            c = 0
            weight = 0.0_qp
        endif

    end subroutine get_copies

end module
