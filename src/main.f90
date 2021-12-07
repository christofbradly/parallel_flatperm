!  main.f90
!****************************************************************************
!  PROGRAM: main
!  PURPOSE:  Entry point for the flatPERM console application.
!****************************************************************************

program main
    use globals
    use mt_stream
    use frontend
    use save_data
    use flatperm, only: run_flatperm
    implicit none

    external    :: sig_handler
    intrinsic   :: signal
    integer, parameter   :: SIGINT = 2, SIGTERM = 15, SIGCONT = 18
    integer     :: save_hours = 24, v(8), tID
    integer(8)  :: s, count_rate, now
    integer(8),allocatable  :: s_complete(:)
    integer     :: seed_size, i, j
    integer,allocatable :: seed_put(:)
    logical     :: saveQ, emailQ
    real(dp)    :: wstart, wsavetime
    logical,allocatable :: thread_wait(:)
    type(mt_state)  :: mts_main
    type(mt_state),save  :: mts
    !$omp threadprivate(mts)

    ! Signal handling
    interruptQ = .false.
    call signal(SIGINT, sig_handler)
    call signal(SIGTERM, sig_handler)
    call signal(SIGCONT, sig_handler)

    ! Get input parameters
    call start_program()

    ! Create or load save file
    ! Setup structures and other variables
    if ( resumeQ ) then
        stop 'Resume no longer works with parallel implementation.'
        call data_interface_hdf5('load')
        OldSamples = sum(samples)
        OldSamplesN = sum(samples(Nmax+1,:,:))
    else
        call data_interface_hdf5('new')
        OldSamples = 0_8
        OldSamplesN = 0_8
    end if


    ! Parallelisation
    call omp_set_num_threads(ntasks)
    allocate(s_complete(ntasks))
    s_complete = 0_8
    allocate(thread_wait(ntasks))
    thread_wait = .false.
    saveQ = .false.

    ! OMP locks are not necessary
    ! if (bestQ) then
    !     allocate(best_lock(ext(2),ext(3)))
    !     do i = 1, ext(2)
    !         do j = 1, ext(3)
    !             call omp_init_lock(best_lock(i,j))
    !         end do
    !     end do
    ! end if
    ! if (fixedQ .and. numpar == 2 .and. scan(params,fixed_param) /= 0) then
    !     allocate(Zval_lock(ext(1)))
    !     allocate(effval_lock(ext(1)))
    !     do i = 1, ext(1)
    !         call omp_init_lock(Zval_lock(i))
    !         call omp_init_lock(effval_lock(i))
    !     end do
    ! end if

    ! start the clock
    call SYSTEM_CLOCK(count = sim_start_time)
    last_save_time=sim_start_time
    allocate(wtime(ntasks)); wtime = 0._dp

    ! RNG
    ! call init_genrand64(input_seed)   ! serial MT implementation
    call set_mt19937()
    call new(mts_main)
    call init(mts_main,input_seed)
    if (percQ) then
        call RANDOM_SEED(size = seed_size)   ! dont need a sophisticated RNG for this part
        allocate(seed_put(seed_size))
        seed_put = int(modulo(input_seed + 373_8 * instance, int(huge(1), kind = 8)))
        call RANDOM_SEED(put = seed_put)
    end if

    ! Run flatPERM simulation
    print*; print*, 'Initialising ', ntasks, ' threads...'; print*
    !$omp parallel default(shared) &
    !$omp private(tID, wstart, wsavetime, s)

    tID = omp_get_thread_num() + 1

    !$omp critical
    call create_stream(mts_main, mts, tID)
    !$omp end critical
    !$omp barrier
    !$omp single
    print*, 'All threads initialised.'
    print*, 'Starting flatPERM algorithm.'
    !$omp end single

    ! start thread clock
    wstart = omp_get_wtime()
    wsavetime = 0._dp

    !$omp do schedule(static)
    do s = OldTours + 1, OldTours + NewTours
        call run_flatperm(s, mts)

        ! Regular saving only needed for resuming, which does not work with parallel implementation
        ! check for regular save point
        ! if (tID == 1) then
        !     call SYSTEM_CLOCK(count=now,count_rate=count_rate)
        !     if (now-last_save_time > count_rate*save_hours*3600) then
        !         last_save_time=now
        !         call date_and_time(values=v)
        !         print*; write(*,'(1X,6(A,I2.2),A)') 'Save point triggered! (', &
        !             v(5), ':', v(6), ':', v(7), ', ', v(3), '-', v(2),'-', v(1)-2000, ')'
        !         ! !$omp atomic write
        !         !     saveQ = .true.
        !         ! !$omp atomic write
        !         !     thread_wait(tID) = .true.
        !         ! pause clock
        !         wsavetime = omp_get_wtime()
        !     end if
        ! end if

        ! if save is triggered pause threads until save is complete
        ! if (saveQ) then
        !     if (tID == 1) then
        !         write(*,'(1X,A,I2,A,I2,A)') &
        !             'Thread', tID, ' waiting for others... (', count(.not.thread_wait), ' remaining)'
        !         emailQ = .false.
        !         do while (.not.all(thread_wait))
        !             call sleep(1)
        !             !$omp cancel do if(interruptQ)
        !             if (.not.emailQ .and. omp_get_wtime() - wsavetime > 3600._dp) then
        !                 call email_notice()
        !                 emailQ = .true.
        !             end if
        !         end do
        !         if (emailQ) emailQ = .false.
        !         write(*,'(1X,A,I2,A)') 'All threads paused, thread', tID, ' proceeding with save...'
        !         call data_interface_hdf5('save', OldTours + sum(s_complete))
        !         write(*,'(1X,A)') 'Save complete! releasing other threads...'
        !         !$omp atomic write
        !             saveQ = .false.
        !         do i = 1, ntasks
        !             !$omp atomic write
        !                 thread_wait(i) = .false.
        !         end do
        !         ! restart clock
        !         wstart = wstart + omp_get_wtime() - wsavetime
        !     else
        !         ! pause clock
        !         wsavetime = omp_get_wtime()
        !         write(*,'(5X,A,I2,A)') '... thread', tID, ' waiting...'
        !         !$omp atomic write
        !             thread_wait(tID) = .true.
        !         call sleep(1)
        !         do while (thread_wait(tID))
        !             call sleep(1)
        !             !$omp cancel do if(interruptQ)
        !         end do
        !         write(*,'(5X,A,I2,A)') '... thread', tID, ' released, resuming loop...'
        !         ! restart clock
        !         wstart = wstart + omp_get_wtime() - wsavetime
        !     end if
        ! end if

        ! record completed tours
        s_complete(tID) = s_complete(tID) + 1_8

        !$omp cancel do if(interruptQ)
    end do
    !$omp end do

    wtime(tID) = omp_get_wtime() - wstart

    !$omp end parallel

    ! Export data to file
    call data_interface_hdf5('save', OldTours + sum(s_complete))
    call data_interface_hdf5('close')
    ! if (resumeQ .and. .not. interruptQ) call data_repack()

    ! Clean up and output
    call output_sim(sum(s_complete))
    deallocate(s_complete, thread_wait, wtime, datadir)
    ! !$ if (allocated(Zval_lock)) deallocate(Zval_lock)
    ! !$ if (allocated(effval_lock)) deallocate(effval_lock)
    ! !$ if (bestQ) deallocate(best_lock)

    if (interruptQ) call exit(2)
    call exit(0)

contains

    subroutine email_notice()
        implicit none
        character(len=256)  :: cmd, cf
        character(len=8)    :: ci
        character(len=12)   :: fmt
        character           :: a
        integer             :: e, c

        print*, "Master still waiting (1hr) for other threads; sending email notice..."

        write(fmt,'(A,I0, A)') '(', ntasks, 'L1)'
        write(ci, *) instance
        ci = adjustl(ci)
        write(cf, fmt) thread_wait
        cf = adjustl(cf)

        cmd="fpscripts/mail_signal.sh "//trim(datadir)//" "//trim(ci)//" "//trim(cf)
        call execute_command_line(trim(cmd), wait = .false., exitstat = e, cmdstat = c, cmdmsg = a)
        print*, e, c, a
    end subroutine email_notice

end program main

! ##############################################
! Signal handlers
subroutine sig_handler
    use globals, only: interruptQ
    implicit none

    interruptQ = .true.
    print*, 'Signal detected...'

end subroutine sig_handler
