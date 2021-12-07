    ! frontend.f90
    ! Subroutines for setting up simulation/major console interactions
    !
module frontend
    use globals
    implicit none
    private
    public start_program, output_sim!,data_repack, output_tour

    character(len=3), dimension(9)  :: valid_latts = ['squ','cub','hc4','hc5','hc6','tri','hex','man','tet']

contains
    ! ##############################################
    ! Parse command line options
    subroutine start_program
        implicit none
        integer             :: v(8)
        character(len=64)   :: opt

        call date_and_time(values = v)
        write(*,*)
        write(*,*) '# # # # # # # # # # # # # # # # # # # #'
        write(*,*) '# #    flatPERM, now in parallel    # #'
        write(*,*) '# # # # # # # # # # # # # # # # # # # #'
        write(*,*)
        write(*, '(1X,6(A,I2.2))') 'Program started at ',v(5),':',v(6),':',v(7),', ',v(3),'-',v(2),'-',v(1)-2000

        ! failsafe: set all flags to false
        collisionsQ = .false.
        endcollQ    = .false.
        bestQ       = .false.
        probenrQ    = .false.
        surfaceQ    = .false.
        anisoQ      = .false.
        fixedQ      = .false.
        fixed2Q     = .false.
        groovesQ    = .false.
        superQ      = .false.
        percQ       = .false.
        boxQ        = .false.

        ! No options
        if (command_argument_count()>0) then
            ! Check for particular options
            call get_command_argument(1,opt)

            select case (opt)
            case ('--new')
                if ( 1 == command_argument_count() ) stop 'No command line options given.'
                call input_from_cmd_line()
            case ('--resume')
                if ( 1 == command_argument_count() ) stop 'No resume file specified.'
                call resume_setup()
                return
            case ('-h','--help')
                call print_help()
                stop
            case default
                print*
                print *, '  Unrecognized command-line option: ', opt
                print*
                call print_help()
                stop
            end select
        else
            call print_help()
            stop
        endif

    end subroutine start_program

    ! ##############################################
    ! Read simulation parameters from command line
    subroutine input_from_cmd_line()
        ! use iso_varying_string
        implicit none
        integer             :: i
        ! character(len=32)   :: cN,cR,cp,cf
        ! character(len=1)    :: cW,cS
        character(len=64)   :: opt, str
        logical             :: skipnext
        integer(8)          :: tours_per_thread

        ! call get_command(str)
        ! print*,str

        skipnext = .false.
        do i = 2, command_argument_count(), 2
            call get_command_argument(i, opt)
            select case (opt)
            case ('-dir')
                call get_command_argument(i + 1, str)
                datadir = str
            case ('-i')
                call get_command_argument(i + 1, str)
                read(str, *) instance
            case ('-tasks')
                call get_command_argument(i + 1, str)
                read(str, *) ntasks
            case ('-l')
                call get_command_argument(i + 1, lattice_type)
                ! read(str, *) lattice_type
            case ('-p')
                call get_command_argument(i + 1, params)
                ! read(str, *) params
            case ('-n')
                call get_command_argument(i + 1, str)
                read(str, *) Nmax
            case ('-s')
                call get_command_argument(i + 1, str)
                read(str, *) input_seed
            case ('-t')
                call get_command_argument(i + 1, str)
                read(str, *) tours_per_thread
            case ('-best')
                call get_command_argument(i + 1, str)
                read(str, *) bestQ
            case ('-perc')
                call get_command_argument(i + 1, str)
                read(str, *) percQ
            case ('-perc_p')
                call get_command_argument(i + 1, str)
                read(str, *) perc_p
                ! read(str, *) cp
            case ('-trails')
                call get_command_argument(i + 1, str)
                read(str, *) collisionsQ
            case ('-grooves')
                call get_command_argument(i + 1, str)
                read(str, *) groovesQ
            case ('-super')
                call get_command_argument(i + 1, str)
                read(str, *) superQ
            case ('-surface')
                call get_command_argument(i + 1, str)
                read(str, *) surfaceQ
                ! if (str=='t') surfaceQ=.true.
            case ('-box')
                call get_command_argument(i + 1, str)
                read(str, *) boxQ
            case ('-max_box')
                call get_command_argument(i + 1, str)
                read(str, *) Lmax
            case ('-fixed')
                call get_command_argument(i + 1, str)
                read(str, *) fixedQ
                ! if (str=='t') fixedQ=.true.
            case ('-fix_p')
                call get_command_argument(i + 1, str)
                read(str, *) fixed_param
            case ('-fix_w')
                call get_command_argument(i + 1, str)
                read(str, *) fixed_weight
            case ('-fixed2')
                call get_command_argument(i + 1, str)
                read(str, *) fixed2Q
                ! if (str=='t') fixedQ=.true.
            case ('-fix2_p')
                call get_command_argument(i + 1, str)
                read(str, *) fixed2_param
            case ('-fix2_w')
                call get_command_argument(i + 1, str)
                read(str, *) fixed2_weight
            case default
                stop 'Invalid command line inputs.'
            end select
        enddo

        OldTours = 0
        NewTours = ntasks * tours_per_thread
        call validate_inputs()  ! If return is successful then inputs are valid.

        if (input_seed == 0) call time_seed()
        compress_strength = 1
        numpar = len(trim(params))
        if (numpar == 1) bestQ = .false.

        print *, 'This is instance ', instance
        write(*, *) 'Data will be saved to folder: ', datadir

    end subroutine input_from_cmd_line

    ! ##############################################
    ! Check command line inputs
    subroutine validate_inputs()
        implicit none
     !$ integer     :: omp_get_max_threads

        ! Check string inputs are valid
        if (all(valid_latts /= lattice_type)) stop &
            'Invalid input parameter; lattice_type must be string.'

        if (scan(params, 'n') /= 1) stop &
            'Invalid input parameter; first index must be n.'

        if (scan(params, 'ncxsmah') == 0) stop &
            'Invalid input parameter; params must be string.'

        ! Check integer inputs
        if (Nmax .le. 0) stop &
            'Invalid input parameter; Nmax must be positive integer.'

        if (instance .le. 0) stop &
            'Invalid input parameter; instance must be positive integer.'

        !$ if (ntasks .le. 0) stop &
        !$     'Invalid input parameter; ntasks must be positive integer.'
        !$ if (ntasks > omp_get_max_threads()) stop &
        !$     'Invalid input parameter; ntasks exceeds available threads.'

        if (input_seed < 0) stop &
            'Invalid input parameter; seed must be non-negative integer (0 for time seed).'

        ! Set new tours
        if (NewTours .le. 0) stop &
            'Invalid input parameter; NewTours must be non-negative integer.'

        ! fixed local weight -- special cases
        if (fixedQ) then
            if (lattice_type == 'hex' .and. fixed_param == 'a') then
                fixed_weight = 1.0_qp + sqrt(2.0_qp)
                print*, 'Fixed weight set to kappa_c = ', fixed_weight
            endif
            if (scan(params, fixed_param) /= 0) then
                print *, 'Note: counting fixed parameter: ', fixed_param
                print *, 'PERM ratio not flattening over ', fixed_param
            endif
        else
            fixed_param = '-'
            fixed2_param = '-'
        endif

        if (fixed2Q) then
            if (.not.fixedQ) stop 'Error: multiple fixed parameters ill-defined'
            print *, 'Note: counting fixed parameter: ', fixed2_param
            print *, 'PERM ratio not flattening over ', fixed2_param
        else
            fixed2_param = '-'
        end if


        ! ISATs
        if (scan(params//fixed_param//fixed2_param, 'cd') > 0) collisionsQ = .true.
        if (.not.collisionsQ) groovesQ = .false.

        ! Adsorption
        if (scan(params//fixed_param//fixed2_param, 'a') > 0 ) surfaceQ = .true.

        ! Anisotropy
        if (scan(params//fixed_param//fixed2_param, 'sb') > 0 ) anisoQ = .true.

        ! percolation/impurities
        if (percQ) then
            if (superQ) stop 'Impurities incompatible with neighbour-avoiding walks'
            if (perc_p < 0 .or. perc_p > 1) stop 'Bad impurity fraction.'
        endif

        ! pulling
        if (scan(params, 'h') > 0) then
            if (.not.surfaceQ) stop 'Cannot measure height without surface.'
        endif

        ! confined to box
        ! if (scan(params, 'l') > 0 .or. fixed_param == 'l') then
        !     boxQ = .true.
        ! endif
        if (boxQ) then
            if (scan(lattice_type, 'squ') == 0 .or. scan(lattice_type, 'cub') == 0) &
                stop 'Confinement unsupported on lattice: '//lattice_type
        end if

        ! Various exceptions
        if ( lattice_type == 'hex' ) then
            if (scan(params, 'scx') /= 0 .or. fixed_param == 'c' .or. fixed_param == 's' &
                .or. collisionsQ .or. groovesQ) then
                stop 'Cannot have collisions or stiffness on hex lattice.'
            endif
        endif

        if (scan(params//fixed_param//fixed2_param, 'd') > 0) then
            if (.not.collisionsQ) stop 'Need trails for triply-visited sites'
            if ( lattice_type /= 'cub' .and. lattice_type /= 'tri') &
                stop 'Can only have triply-visited sites on cub or tri lattices.'
        endif

        if (groovesQ) then
            if (scan(params, 'x') /= 0) then
                stop 'Grooves cannot have crossings.'
            elseif (lattice_type == 'cub') then
                stop 'Grooves not implemented for cubic lattice.'
            endif
        elseif (scan(params, 'x') /= 0 .and. lattice_type == 'cub') then
            stop 'Counting crossings not implemented for cubic lattice.'
        endif

        if (scan(params, 'x') /= 0 .and. scan(params, 'x') < scan(params, 'c')) &
            stop 'Collisions must be tested before crossings.'

        if (superQ .and. scan(params, 'cxm') /= 0) &
            stop 'Cannot have monomer-monomer interactions on superSAWs'

        if (superQ .and. (collisionsQ .or. groovesQ)) &
            stop 'Cannot have superSATs'

    end subroutine validate_inputs

    ! ##############################################
    ! Setup for resuming from ifle
    subroutine resume_setup()
        implicit none
        character(len=64)   :: str

        call get_command_argument(2,str)
        read(str, *) instance

        call get_command_argument(3,str)
        read(str, *) input_seed
        if (input_seed==0) call time_seed()
        resumeQ=.true.

        call get_command_argument(4,str)
        read(str, *) NewTours

        call get_command_argument(5,str)
        read(str, *) ntasks

        call get_command_argument(6,str)
        datadir=trim(str)
        print *, 'Resuming from folder: ',trim(datadir),', instance ',instance



        !call load_parameters_hdf5()

    end subroutine resume_setup

    ! ##############################################
    ! Seed RNG with cpu clock
    ! input_seed is int(4) for parallel MT and int (8) for serial MT
    subroutine time_seed()
        implicit none
        integer      :: count, count_rate, count_max

        call SYSTEM_CLOCK(count, count_rate, count_max)
        input_seed = count
        print*, 'RNG seeded with clock time: ', input_seed
    end subroutine time_seed

    ! ##############################################
    ! Help console
    subroutine print_help()
        implicit none
        print *, 'flatPERM options (choose one only):'
        print *, '  [none], -h, --help      Print this message.'
        print *,
        print *, '  --new <options>         Run flatPERM with options taken from command line'
        print *, '          Required:'
        print *, '              -dir        name of directory for this run (within /savedata/)'
        print *, '              -i          number of this instance, for array jobs'
        print *, '              -tasks      number of OMP threads to use'
        print *, '              -l          lattice_type; one of squ, cub, hex, man, tet, tri'
        print *, '              -p          params; must start with n; see below for options'
        print *, '              -n          Nmax'
        print *, '              -s          RNG seed; enter 0 for time seed'
        print *, '              -t          NewTours'
        print *, '          Optional:'
        print *, '              -best       record samples with highest weights'
        print *, '              -trails     force SATs; ignored if already required'
        print *, '              -grooves    force grooves; ignored if not trails'
        print *, '              -super      force neighbour-avoiding walks; ignored if not trails'
        print *, '              -surface    force surface; ignored if already required'
        print *, '              -fixed      force fixed weight'
        print *, '              -fix_p      set fixed parameter; ignored if not fixed'
        print *, '              -fix_w      set fixed weight; ignored if not fixed'
        print *, '              -fixed2     EXPERIMENTAL! force second fixed weight'
        print *, '              -fix2_p     EXPERIMENTAL! set second fixed parameter; ignored if not fixed'
        print *, '              -fix2_w     EXPERIMENTAL! set second fixed weight; ignored if not fixed'
        print *, '              -perc       force lattice defects'
        print *, '              -perc_p     set lattice defect fraction; ignored if not perc'
        print *, '              -box        force confining box; squ lattice only'
        print *, '              -max_box    set maximum box size; ignored if not box'
        print *, '          Valid parameters for "param" argument:'
        print *, '              (see subroutine update_micro() for details and possible other options)'
        print *, '                  n       length of chain'
        print *, '                  m       (non-consecutive) near-neighbours (ISAW)'
        print *, '                  c       number of multiply-visited or doubly visited sites (trails only)'
        print *, '                  d       number of triply-visited sites (ISAT cub or tri lattice)'
        print *, '                  a       number of contacts with surface'
        print *, '                  h       height of endpoint above surface'
        print *, '                  s       number of straight segments (stiffness)'
        print *, '                  l       size of bounding box'
        print*
        print *, '  --resume <options>      DEFUNCT! Resume flatPERM from previous run. Options must be entered as:'
        print *, '                          <seed> <NewTours> <filename w/o .h5 extension>'
        print *,
        print *, 'Program created by:'
        print *
        print *, '  Chris Bradly (2021)'
        print *, '  School of Mathematics & Statistics'
        print *, '  University of Melbourne'
        print *
    end subroutine print_help

    ! ##############################################
    ! Output - no longer used; functionality moved to main.f90
    ! subroutine output_tour(t, saveflag, rpkflag)
    !     implicit none
    !     integer(8), intent(in)  :: t
    !     logical, intent(out)    :: saveflag, rpkflag
    !     integer                 :: v(8), hours = 24, savecount = 0
    !     real(dp)                :: now, last_save_time = 0.0 ! initialised here to save value across subroutine calls
    !
    !     saveflag = .false.
    !     rpkflag = .false.
    !     if (t == 1) last_save_time = sim_start_time
    !
    !     call cpu_time(now)
    !     if (now - last_save_time > hours*3600.0_dp) then
    !         saveflag = .true.
    !         savecount = savecount + 1
    !         call date_and_time(values = v)
    !         write(*, '(3X,A,I8)')   'Finished tour: ', t
    !         write(*, '(1X,6(A,I2.2))') 'Updating save file; ', &
    !             v(5), ':', v(6), ':', v(7), ', ', v(3), '-', v(2), '-', v(1) - 2000
    !         if (modulo(savecount, 2) == 0) rpkflag = .true.
    !         last_save_time = now
    !     endif
    !
    ! end subroutine output_tour

    ! ##############################################
    ! Output some info about the run
    subroutine output_sim(t)
        implicit none
        integer(8),intent(in)  :: t
        integer             :: v(8)!,i,j
        ! character           :: a
        character(len = 64) :: fmtI, fmtR, fmtC1, fmtC2
        integer(8)          :: s, sN, sample_rate, now, count_rate
        real(dp)            :: run_time, tour_rate
        real(dp)            :: e, eN

        s = sum(samples)
        sN = sum(samples(size(samples, 1),:,:))
        e = real(sum(effsamples), dp)  ! this down conversion of real kind shouldn't ever be a problem, but in principle it is a bit unsafe.
        eN = real(sum(effsamples(size(effsamples, 1),:,:)), dp)  ! this down conversion of real kind shouldn't ever be a problem, but in principle it is a bit unsafe.
        ! call cpu_time(now)
        call SYSTEM_CLOCK(count = now, count_rate = count_rate)
        run_time = real(now - sim_start_time, dp)/count_rate
        if (run_time < 0.1E-10) then
            tour_rate = 0
            sample_rate = 0
        else
            tour_rate = (t - OldTours)/run_time
            sample_rate = int((s - OldSamples)/run_time, 8)
        endif

        write(fmtI, '( "I",I2,"," )') 15
        write(fmtR, '( "",I1,"X, ES8.2," )') 7
        ! write(fmtC1,'( "A",I2,"," )') 24
        ! write(fmtC2,'( "A",I2,"," )') 16
        write(fmtC1, '( "A",I2,"," )')
        write(fmtC2, '( "A",I2,"," )')

        call date_and_time(values = v)
        print*
        write(*, '(1X,6(A,I2.2))') 'Program finished at ', &
            v(5), ':', v(6), ':', v(7), ', ', v(3), '-', v(2),'-', v(1) - 2000

        ! write(*, '(5X,A24,I12,A16,I12,A)')       'Total samples:',s,' (max length:',sN,')'
        write(*, '(5X,A,10X,'//trim(fmtI)//')') 'Total samples:', s
        write(*, '(9X,A,8X,'//trim(fmtI)//'A)')  '(max length:', sN, ')'
        write(*, '(5X,A,8X,'//trim(fmtR)//')')  'Total effective:', e
        write(*, '(9X,A,8X,'//trim(fmtR)//'A)')  '(max length:', eN, ')'

        write(*, '(5X,A,14X,'//trim(fmtR)//'A)') 'Tour rate:',  tour_rate, ' tours/sec'
        write(*, '(5X,A,12X,'//trim(fmtR)//'A)') '(per thread:',  tour_rate/ntasks, ' tours/sec)'
        write(*, '(5X,A,12X,'//trim(fmtI)//'A)') 'Sample rate:',  sample_rate,' samples/sec'

        write(*, '(1X,A,6X,4X,'//trim(fmtR)//'A,F6.1,A)') 'Algorithm runtime:',  run_time,' sec (', run_time/3600, ' hours)'
        write(*, '(5X,A,6X,'//trim(fmtR)//'A,'//'I2,'//'A)') 'Total thread time:', sum(wtime),' sec (' ,ntasks, ' threads)'
        write(*, '(5X,A,4X,'//trim(fmtR)//'A)') 'Average thread time:', sum(wtime)/ntasks, ' sec'

        call dealloc_arrays()
        print*

    contains
        ! ##############################################
        ! Deallocate data arrays -- free up space even if program is not exited.
        subroutine dealloc_arrays()
            implicit none

            ! deallocate history and chain arrays
            ! deallocate(history_enr,history_atmo,history_W,history_steps,chain)
            ! deallocate data arrays
            deallocate(samples,prunes,enriches,Z)
            ! conditionally utilised arrays
            if (allocated(effsamples)) deallocate(effsamples)
            if (allocated(Re2W)) deallocate(Re2W)
            if (allocated(Re2Wp)) deallocate(Re2Wp)
            if (allocated(Rm2W)) deallocate(Rm2W)
            if (allocated(Rg2W)) deallocate(Rg2W)
            if (allocated(aniso)) deallocate(aniso)
            if (allocated(aniso2)) deallocate(aniso2)
            if (allocated(FixMom1)) deallocate(FixMom1)
            if (allocated(FixMom2)) deallocate(FixMom2)
            if (allocated(Fix2Mom1)) deallocate(Fix2Mom1)
            if (allocated(Fix2Mom2)) deallocate(Fix2Mom2)
            if (allocated(bestWalks)) deallocate(bestWalks)
            if (allocated(bestWeights)) deallocate(bestWeights)
            print*
            print*, 'Arrays deallocated.'

        end subroutine dealloc_arrays

    end subroutine output_sim

end module
