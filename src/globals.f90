    ! globals.f90
    ! Declare global variabes and data arrays
    !

module globals
    use omp_lib
    use, intrinsic :: iso_fortran_env
    ! use fhash_module__ints_ints,  only: fhash_type__ints_ints
    implicit none

    ! signal handling
    logical     :: interruptQ

    ! Standardised definition for reals precision
    integer, parameter  :: dp=selected_real_kind(p=10, r=100)
    !! Extended precision:   p>15,  r>308 (,  kind=10)
    integer, parameter  :: qp=selected_real_kind(p=16, r=1000)
    !! Double precision
    ! integer, parameter  :: qp=selected_real_kind(p=10, r=100)
    !! Quad precision:       p>18,  r>308 (,  kind=16)
    ! integer, parameter  :: qp=selected_real_kind(p=20, r=1000)
    ! integer, parameter  :: qp=16


    ! time
    integer(8)              :: last_save_time, sim_start_time
    real(dp), allocatable   :: wtime(:)

    ! Declare flags
    logical     :: collisionsQ, endcollQ, bestQ, probenrQ, surfaceQ, anisoQ
    logical     :: fixedQ, fixed2Q, groovesQ, superQ, percQ, boxQ
    integer     :: compress_strength

    ! Declare input parameters
    integer             :: Nmax, Lmax, numpar, instance, max_instances!, max_holes
    integer             :: input_seed
    integer(8)          :: OldTours, NewTours, OldSamples, OldSamplesN
    real(qp)            :: mu, fixed_weight, fixed2_weight, normq
    character(len=:), allocatable  :: datadir
    logical             :: resumeQ=.false.
    character(len=7)    :: params
    character(len=1)    :: fixed_param, fixed2_param
    character(len=3)    :: lattice_type
    real(qp)            :: delay=0.0_qp!0.1_qp
    real(dp)            :: perc_p!, exp_holes

    ! parallel
    integer             :: ntasks, break_rank
    ! integer(kind=omp_lock_kind), allocatable  :: Zval_lock(:)
    ! integer(kind=omp_lock_kind), allocatable  :: effval_lock(:)
    ! integer(kind=omp_lock_kind), allocatable  :: best_lock(:, :)

    ! Declare data arrays
        ! This method leaves room for up to 7 microcanonical variables without needing to rewrite the program
    integer, allocatable     :: bestWalks(:, :, :, :)
    real(qp), allocatable    :: bestWeights(:, :)
    integer(8), allocatable  :: samples(:, :, :), prunes(:, :, :), enriches(:, :, :)
    real(qp), allocatable    :: Re2W(:, :, :), Rm2W(:, :, :), Rg2W(:, :, :)
    real(qp), allocatable    :: Re2Wp(:, :, :), aniso(:, :, :), aniso2(:, :, :)
    real(qp), allocatable    :: FixMom1(:, :, :), FixMom2(:, :, :), Fix2Mom1(:, :, :), Fix2Mom2(:, :, :)
    real(qp), allocatable    :: Z(:, :, :)
    real(qp), allocatable    :: effsamples(:, :, :)
    integer     :: ext(3)

contains
    ! ##############################################
    ! Random integer in range [m, ..., n]
    subroutine random_integer(m, n, i)
        implicit none
        integer, intent(in)  :: m, n
        integer, intent(out) :: i
        real(dp)            :: a

        call random_number(a)
        i = m + floor((n + 1 - m) * a)
    end subroutine random_integer

    ! Wrappers for MKL RNG -- no longer used
!    ! ##############################################
!    ! Wrapper for any type of RNG method
!    ! Random real in range [0.0, 1.0]
!    subroutine random_real(b)
!        !use mkl_vsl_type
!        !use mkl_vsl
!        implicit none
!
!        real(dp), intent(out) :: b
!        real                :: r(1)
!        integer             :: errcode
!
!        if ( brng == 0 ) then
!            call random_number(b)
!        else
!            errcode=vsrnguniform( VSL_RNG_METHOD_UNIFORM_STD,  stream,  1,  r,  0.0,  1.0 )
!            if (errcode < 0) then
!                print *,  'MKL VSL error (at call); code: ', errcode
!                print *,  'Refer to module MKL_VSL_TYPE in mkl_vsl.f90'
!                stop
!            endif
!            b=r(1)
!        endif
!
!    end subroutine random_real
!
!    ! ##############################################
!    ! Seed the RNG
!    subroutine setup_RNG()
!        implicit none
!        integer                 :: i_seed, i, errcode
!        integer, allocatable     :: a_seed(:)
!
!        ! brng=0 uses default
!        if ( brng == 0 ) then
!            call random_seed(size=i_seed)
!            allocate(a_seed(i_seed))
!            a_seed=input_seed + 37*(/ (i, i=0, i_seed-1) /)
!            call random_seed(put=a_seed)
!        else
!            call init_genrand64(input_seed)
!        endif
!
!    end subroutine setup_RNG

    ! ##############################################
    ! Pause
    subroutine dbg_pause(msg)
        implicit none
        character(len=*), intent(in) :: msg

        write(*, *) msg
        write(*, *) 'Press Enter to continue.'
        read(*, *)
    end subroutine dbg_pause

end module
