    ! save_data.f90
    ! Frontend specifically for save/load data arrays to file
    !

module save_data
    use hdf5
    use globals
    use lattice
    implicit none
    private
    public data_interface_hdf5!,save_data_basic

    integer(HID_T)  :: file_id=0

contains
    ! ##############################################
    ! Export data to hdf5 file
    subroutine data_interface_hdf5(action,s)
        implicit none
        character(len=*),intent(in)     :: action
        integer(8),intent(in),optional  :: s
        character(len=256)              :: filename!,homedir
        ! real(dp)                        :: save_access_time,now
        integer                         :: err,config_flags
        logical                         :: status
        character(len=2)                :: c_inst


        write(c_inst,'(I2.2)') instance
        filename='savedata/'//trim(datadir)//'/'//c_inst//'.h5'
        print*
        print *, '  HDF5 action: ', action

        select case (action)
        case ('new')
            ! create new file
            call h5open_f(err)
            call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, err)

            ! zipping
            call h5zfilter_avail_f(H5Z_FILTER_DEFLATE_F,status,err)
            if (.not. status) stop 'gzip not registered.'
            call h5zget_filter_info_f(H5Z_FILTER_DEFLATE_F, config_flags, err)
            if (0==iand(config_flags,H5Z_FILTER_ENCODE_ENABLED_F)) stop &
                'gzip encode disabled.'
            if (0==iand(config_flags,H5Z_FILTER_DECODE_ENABLED_F)) stop &
                'gzip decode disabled.'

            ! Initialise lattice and arrays
            call init_lattice(lattice_type, collisionsQ)
            ! if (boxQ) then
            !     Lmax = Nmax
            !     if (collisionsQ) then
            !         Nmax = (coord/2)*Lmax**lattdim - 1
            !     else
            !         Nmax = Lmax**lattdim - 1
            !     end if
            ! end if

            ! if (mu.le.0.0_qp) then
                mu=real(mu_cc,qp)
                print*, 'Renormalisation mu set to connective constant: ', mu
            ! endif
            ! if (percQ) exp_holes=perc_p*Nmax*(coord-2)
            ! Create new attributes
            call new_parameters_hdf5()
            call init_data_arrays()
            call create_datasets_hdf5()
            ! close file
            call h5fclose_f(file_id,err)

        case ('load')
            ! open existing file
            call h5open_f(err)
            call h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, err)
            ! zipping
            call h5zfilter_avail_f(H5Z_FILTER_DEFLATE_F,status,err)
            if (.not. status) stop 'gzip not registered.'
            call h5zget_filter_info_f(H5Z_FILTER_DEFLATE_F, config_flags, err)
            if (0==iand(config_flags,H5Z_FILTER_ENCODE_ENABLED_F)) stop &
                'gzip encode disabled.'
            if (0==iand(config_flags,H5Z_FILTER_DECODE_ENABLED_F)) stop &
                'gzip decode disabled.'
            ! Load attributes
            call load_parameters_hdf5()
            ! Load arrays and initialise lattice
            call init_lattice(lattice_type, collisionsQ)
            print*, 'Renormalisation mu set to connective constant: ', mu

            ! if (percQ) exp_holes=perc_p*Nmax*(coord-2)

            call init_data_arrays()
            call load_datasets_hdf5()
            ! close file
            call h5fclose_f(file_id,err)
            ! print*, 'Current samples: ', sum(samples)

        case ('save')
            call h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, err)
            ! call cpu_time(save_access_time)
            call save_datasets_hdf5()
            call update_parameters_hdf5(s)
            ! call cpu_time(now)
            ! save_time=save_time+(now-save_access_time)
            ! close file
            call h5fclose_f(file_id,err)
        case ('close')
            ! close file
            ! call h5fclose_f(file_id,err)
            ! close HDF5 interface
            call h5close_f(err)

            call clean_lattice()
        case default
            call h5close_f(err)
            print *, 'Unknown call to HDF5.'
            stop
        end select

    end subroutine data_interface_hdf5

    ! ##############################################
    ! Save parameters for new file
    subroutine new_parameters_hdf5()
        implicit none
        integer(HID_T)                  :: attrspace_id,attr_id,atype_id
        integer(HSIZE_T)                :: num_attr(1)=1,attr_dims(1)=1
        integer                         :: err,flags(11)
        integer(SIZE_T)                 :: str_len
        character(len=8)                :: date_str
        integer(HID_T)      :: EXT_REAL_TYPE
        EXT_REAL_TYPE = h5kind_to_type(qp,H5_REAL_KIND)
        ! Type testing
        ! print*,qp,EXT_REAL_TYPE,H5T_NATIVE_DOUBLE,H5T_IEEE_F64LE
        ! print*,kind(qp),kind(EXT_REAL_TYPE),kind(H5T_NATIVE_DOUBLE),kind(h5kind_to_type(qp,H5_REAL_KIND)),HID_T,SIZE_T
        ! call h5fclose_f(file_id,err)
        ! stop
        ! Create attribute dataspace
        call h5screate_simple_f(1, num_attr, attrspace_id, err)
        ! Numeric parameters
        call h5acreate_f(file_id, 'Nmax', H5T_NATIVE_INTEGER, attrspace_id, attr_id, err)
        call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, Nmax, attr_dims, err)
        call h5aclose_f(attr_id, err)

        if (boxQ) then
            call h5acreate_f(file_id, 'Lmax', H5T_NATIVE_INTEGER, attrspace_id, attr_id, err)
            call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, Lmax, attr_dims, err)
            call h5aclose_f(attr_id, err)
        end if

        call h5acreate_f(file_id, 'Tours', H5T_STD_I64LE, attrspace_id, attr_id, err)
        call h5awrite_f(attr_id, H5T_STD_I64LE, 0, attr_dims, err)
        call h5aclose_f(attr_id, err)

        call h5acreate_f(file_id, 'seed', H5T_STD_I64LE, attrspace_id, attr_id, err)
        call h5awrite_f(attr_id, H5T_STD_I64LE, input_seed, attr_dims, err)
        call h5aclose_f(attr_id, err)

        call h5acreate_f(file_id, 'mu', EXT_REAL_TYPE, attrspace_id, attr_id, err)
        call h5awrite_f(attr_id, EXT_REAL_TYPE, mu, attr_dims, err)
        call h5aclose_f(attr_id, err)

        call h5acreate_f(file_id, 'instance', H5T_NATIVE_INTEGER, attrspace_id, attr_id, err)
        call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, instance, attr_dims, err)
        call h5aclose_f(attr_id, err)

        if (percQ) then
            call h5acreate_f(file_id, 'perc_p', H5T_NATIVE_DOUBLE, attrspace_id, attr_id, err)
            call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, perc_p, attr_dims, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(file_id, 'normq', EXT_REAL_TYPE, attrspace_id, attr_id, err)
            call h5awrite_f(attr_id, EXT_REAL_TYPE, normq, attr_dims, err)
            call h5aclose_f(attr_id, err)
        endif

        ! String parameters
        CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, err)

        str_len=int(len(lattice_type),kind=SIZE_T)
        CALL h5tset_size_f(atype_id, str_len, err)
        call h5acreate_f(file_id, 'lattice_type', atype_id, attrspace_id, attr_id, err)
        call h5awrite_f(attr_id, atype_id, lattice_type, attr_dims, err)
        call h5aclose_f(attr_id, err)

        str_len=int(numpar,kind=SIZE_T)
        CALL h5tset_size_f(atype_id, str_len, err)
        call h5acreate_f(file_id, 'params', atype_id, attrspace_id, attr_id, err)
        call h5awrite_f(attr_id, atype_id, trim(params), attr_dims, err)
        call h5aclose_f(attr_id, err)

        call date_and_time(date=date_str)
        str_len=int(len(date_str),kind=SIZE_T)
        CALL h5tset_size_f(atype_id, str_len, err)
        call h5acreate_f(file_id, 'date', atype_id, attrspace_id, attr_id, err)
        call h5awrite_f(attr_id, atype_id, date_str, attr_dims, err)
        call h5aclose_f(attr_id, err)

        call h5sclose_f(attrspace_id, err)

        ! Logical parameters
        flags=merge([1,1,1,1,1,1,1,1,1,1,1],[0,0,0,0,0,0,0,0,0,0,0],&
            [collisionsQ,surfaceQ,bestQ,probenrQ,fixedQ,fixed2Q,groovesQ,superQ,percQ,anisoQ,boxQ])
        num_attr=size(flags)
        attr_dims=size(flags)
        call h5screate_simple_f(1, num_attr, attrspace_id, err)
                attr_dims=size(flags)
        call h5acreate_f(file_id, 'flags', H5T_NATIVE_INTEGER, attrspace_id, attr_id, err)
        call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, flags, attr_dims, err)
        call h5aclose_f(attr_id, err)
        call h5sclose_f(attrspace_id, err)

        ! Fixed weight
        if (fixedQ) then
            num_attr=1
            attr_dims=1
            call h5screate_simple_f(1, num_attr, attrspace_id, err)
            call h5acreate_f(file_id, 'fixedW', EXT_REAL_TYPE, attrspace_id, attr_id, err)
            call h5awrite_f(attr_id, EXT_REAL_TYPE, fixed_weight, attr_dims, err)
            call h5aclose_f(attr_id, err)

            str_len=int(1,kind=SIZE_T)
            CALL h5tset_size_f(atype_id, str_len, err)
            call h5acreate_f(file_id, 'fixed_param', atype_id, attrspace_id, attr_id, err)
            call h5awrite_f(attr_id, atype_id, trim(fixed_param), attr_dims, err)
            call h5aclose_f(attr_id, err)
            call h5sclose_f(attrspace_id, err)

            print*,'Running with (uncounted) fixed parameter ', fixed_param,'; weight = ', fixed_weight

            if (fixed2Q) then
                num_attr=1
                attr_dims=1
                call h5screate_simple_f(1, num_attr, attrspace_id, err)
                call h5acreate_f(file_id, 'fixed2W', EXT_REAL_TYPE, attrspace_id, attr_id, err)
                call h5awrite_f(attr_id, EXT_REAL_TYPE, fixed2_weight, attr_dims, err)
                call h5aclose_f(attr_id, err)

                str_len=int(1,kind=SIZE_T)
                CALL h5tset_size_f(atype_id, str_len, err)
                call h5acreate_f(file_id, 'fixed2_param', atype_id, attrspace_id, attr_id, err)
                call h5awrite_f(attr_id, atype_id, trim(fixed2_param), attr_dims, err)
                call h5aclose_f(attr_id, err)
                call h5sclose_f(attrspace_id, err)

                print*,'Running with (uncounted) fixed parameter ', fixed2_param,'; weight = ', fixed2_weight
            endif
        endif

    end subroutine new_parameters_hdf5

    ! ##############################################
    ! Extract simulation parameters from resumed file
    subroutine load_parameters_hdf5()
        implicit none
        integer(HID_T)                  :: attr_id,atype_id
        integer(HSIZE_T)                :: num_attr(1)=1,attr_dims(1)=1
        integer                         :: err,flags(11),file_inst
        integer(SIZE_T)                 :: str_len
        integer(HID_T)      :: EXT_REAL_TYPE
        EXT_REAL_TYPE = h5kind_to_type(qp,H5_REAL_KIND)

        ! Create attribute dataspace
        ! call h5sopen_f(1, num_attr, attrspace_id, err)
        ! Numeric parameters
        call h5aopen_f(file_id, 'instance', attr_id, err)
        call h5aread_f(attr_id, H5T_NATIVE_INTEGER, file_inst, attr_dims, err)
        call h5aclose_f(attr_id, err)
        if (instance/=file_inst) then
            print*, 'ERROR: Instance mismatch. file_inst: ', file_inst,'; slurm instance: ', instance
            stop
        endif

        call h5aopen_f(file_id, 'Nmax', attr_id, err)
        call h5aread_f(attr_id, H5T_NATIVE_INTEGER, Nmax, attr_dims, err)
        call h5aclose_f(attr_id, err)

        call h5aopen_f(file_id, 'Tours', attr_id, err)
        call h5aread_f(attr_id, H5T_STD_I64LE, OldTours, attr_dims, err)
        write(*,'(1X,A,I8,A,I8,A)', advance='no') 'Already have ', OldTours,' tours, running ', NewTours,' more.'
        ! read(*,*) NewTours
        ! call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, OldTours+NewTours, attr_dims, err)
        call h5aclose_f(attr_id, err)

        ! RNG reset with clock time in frontend.f90/resume_setup()
        ! call h5aopen_f(file_id, 'seed', attr_id, err)
        ! call h5aread_f(attr_id, H5T_STD_I64LE, input_seed, attr_dims, err)
        ! call h5aclose_f(attr_id, err)

        call h5aopen_f(file_id, 'mu', attr_id, err)
        call h5aread_f(attr_id, EXT_REAL_TYPE, mu, attr_dims, err)
        call h5aclose_f(attr_id, err)

        call h5aexists_f(file_id,'perc_p', percQ,err)
        if (percQ) then
            call h5aopen_f(file_id, 'perc_p', attr_id, err)
            call h5aread_f(attr_id, H5T_NATIVE_DOUBLE, perc_p, attr_dims, err)
            call h5aclose_f(attr_id, err)
            call h5aopen_f(file_id, 'normq', attr_id, err)
            call h5aread_f(attr_id, EXT_REAL_TYPE, normq, attr_dims, err)
            call h5aclose_f(attr_id, err)
        else
            perc_p=0.0_dp
            normq = 0.0_qp
        endif

        ! String parameters
        CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, err)

        str_len=int(3,kind=SIZE_T)
        CALL h5tset_size_f(atype_id, str_len, err)
        call h5aopen_f(file_id, 'lattice_type', attr_id, err)
        call h5aread_f(attr_id, atype_id, lattice_type, attr_dims, err)
        call h5aclose_f(attr_id, err)

        str_len=int(7,kind=SIZE_T)
        CALL h5tset_size_f(atype_id, str_len, err)
        call h5aopen_f(file_id, 'params', attr_id, err)
        call h5aread_f(attr_id, atype_id, params, attr_dims, err)
        call h5aclose_f(attr_id, err)
        numpar=len(trim(params))

        ! Logical parameters
        num_attr=size(flags)
        attr_dims=size(flags)
        flags=0
        ! call h5screate_simple_f(1, num_attr, attrspace_id, err)
        call h5aopen_f(file_id, 'flags', attr_id, err)
        call h5aread_f(attr_id, H5T_NATIVE_INTEGER, flags, attr_dims, err)
        call h5aclose_f(attr_id, err)
        ! print*,flags

        collisionsQ = i2l(flags(1))
        surfaceQ    = i2l(flags(2))
        bestQ       = i2l(flags(3))
        probenrQ    = i2l(flags(4))
        fixedQ      = i2l(flags(5))
        fixed2Q      = i2l(flags(6))
        groovesQ    = i2l(flags(7))
        superQ      = i2l(flags(8))
        anisoQ      = i2l(flags(9))
        percQ       = i2l(flags(10))
        boxQ        = i2l(flags(11))

        if (boxQ) then
            num_attr=1
            attr_dims=1
            call h5aopen_f(file_id, 'Lmax', attr_id, err)
            call h5aread_f(attr_id, H5T_NATIVE_INTEGER, Lmax, attr_dims, err)
            call h5aclose_f(attr_id, err)
        end if

        ! Fixed weight
        if (fixedQ) then
            num_attr=1
            attr_dims=1
            call h5aopen_f(file_id, 'fixedW', attr_id, err)
            call h5aread_f(attr_id, EXT_REAL_TYPE, fixed_weight, attr_dims, err)
            call h5aclose_f(attr_id, err)

            str_len=int(1,kind=SIZE_T)
            CALL h5tset_size_f(atype_id, str_len, err)
            call h5aopen_f(file_id, 'fixed_param', attr_id, err)
            call h5aread_f(attr_id, atype_id, fixed_param, attr_dims, err)
            call h5aclose_f(attr_id, err)
            print*,
            print*,'Resuming with fixed parameter ', fixed_param,'; local weight = ', fixed_weight

            if (fixed2Q) then
                num_attr=1
                attr_dims=1
                call h5aopen_f(file_id, 'fixed2W', attr_id, err)
                call h5aread_f(attr_id, EXT_REAL_TYPE, fixed2_weight, attr_dims, err)
                call h5aclose_f(attr_id, err)

                str_len=int(1,kind=SIZE_T)
                CALL h5tset_size_f(atype_id, str_len, err)
                call h5aopen_f(file_id, 'fixed2_param', attr_id, err)
                call h5aread_f(attr_id, atype_id, fixed2_param, attr_dims, err)
                call h5aclose_f(attr_id, err)
                print*,
                print*,'Resuming with fixed parameter ', fixed2_param,'; local weight = ', fixed2_weight
            endif
        endif

    contains
        ! integer to logical conversion
        pure function i2l(i) result(l)
            implicit none
            integer, intent(in) :: i
            logical             :: l
            l = .true.
            if ( i == 0 ) l = .false.
        end function i2l

    end subroutine load_parameters_hdf5

    ! ##############################################
    ! Extract simulation parameters from resumed file
    subroutine update_parameters_hdf5(s)
        implicit none
        integer(8),intent(in)           :: s
        integer(HID_T)                  :: attr_id
        integer(HSIZE_T)                :: attr_dims(1)=1
        integer                         :: err
        integer(HID_T)      :: EXT_REAL_TYPE
        EXT_REAL_TYPE = h5kind_to_type(qp,H5_REAL_KIND)

        call h5aopen_f(file_id, 'Tours', attr_id, err)
        call h5awrite_f(attr_id, H5T_STD_I64LE, s, attr_dims, err)
        call h5aclose_f(attr_id, err)

        if (percQ) then
            call h5aopen_f(file_id, 'normq', attr_id, err)
            call h5awrite_f(attr_id, EXT_REAL_TYPE, normq, attr_dims, err)
            call h5aclose_f(attr_id, err)
        endif

    end subroutine update_parameters_hdf5

    ! ##############################################
    ! Create new datasets in .h5 file
    subroutine create_datasets_hdf5()
        implicit none
        integer(HSIZE_T)    :: chunk_dims(numpar)
        ! chunk_dims = int(Nmax + 1, kind = HSIZE_T)
        chunk_dims = int(ext(1:numpar), kind = HSIZE_T)
        do while ((16_8*chunk_dims(1)**int(numpar,8)) .ge. (38_8*10_8**8_8))
            ! Large chunks are faster and smaller but 4GB is hard limit.
            chunk_dims = chunk_dims/2
        enddo

        call create_array_int_hdf5('Sn')
        call create_array_int_hdf5('Enr')
        call create_array_int_hdf5('Pru')
        call create_array_real_hdf5('Z')
        call create_array_real_hdf5('Se')
        if (bestQ) then
            call create_array_walks_hdf5('bestWalks')
            call create_array_weights_hdf5('bestWeights')
        endif
        ! if (percQ) call create_array_holes_hdf5('holes')

        ! if (fixedQ .and. numpar == 2) then
        !     call create_array_real_hdf5('FixMom1')
        !     call create_array_real_hdf5('FixMom2')
        ! endif
        if (fixedQ) then
            call create_array_real_hdf5('FixMom1')
            call create_array_real_hdf5('FixMom2')
            if (fixed2Q) then
                call create_array_real_hdf5('Fix2Mom1')
                call create_array_real_hdf5('Fix2Mom2')
            endif
        endif

        if (surfaceQ) then
            call create_array_real_hdf5('Re2W')
            call create_array_real_hdf5('Re2Wp')
        else
            call create_array_real_hdf5('Re2W')
            ! call create_array_real_hdf5('Rm2W')
            ! call create_array_real_hdf5('Rg2W')
        endif

        if (anisoQ) then
            call create_array_real_hdf5('aniso')
            if (surfaceQ .and. lattice_type == 'cub') call create_array_real_hdf5('aniso2')
        endif

    contains
        subroutine create_array_real_hdf5(name)
            implicit none
            character(len=*),intent(in) :: name
            integer(HID_T)              :: dataset_id,dataspace_id,prop_id
            integer(HSIZE_T)            :: file_dims(numpar)
            integer                     :: err,rank
            integer(HID_T)      :: EXT_REAL_TYPE
            EXT_REAL_TYPE = h5kind_to_type(qp,H5_REAL_KIND)

            file_dims=int(ext(1:numpar),kind=HSIZE_T)

            rank=numpar
            ! create dataspace
            call h5screate_simple_f(rank, file_dims , dataspace_id, err)
            ! Chunk and enable compression
            call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,err)
            call h5pset_chunk_f(prop_id,rank,chunk_dims,err)
            call h5pset_deflate_f(prop_id,compress_strength,err)
            ! create dataset
            call h5dcreate_f(file_id, name, EXT_REAL_TYPE, dataspace_id,&
                                dataset_id, err, prop_id)
            ! close
            call h5dclose_f(dataset_id, err)
            call h5sclose_f(dataspace_id, err)
            call h5pclose_f(prop_id,err)
        end subroutine create_array_real_hdf5

        subroutine create_array_int_hdf5(name)
            implicit none
            character(len=*),intent(in) :: name
            integer(HID_T)              :: dataset_id,dataspace_id,prop_id
            integer(HSIZE_T)            :: file_dims(numpar)
            integer                     :: err,rank
            ! file_dims=int(Nmax+1,kind=HSIZE_T)
            file_dims=int(ext(1:numpar),kind=HSIZE_T)
            rank=numpar
            ! create dataspace
            call h5screate_simple_f(rank, file_dims , dataspace_id, err)
            ! Chunk and enable compression
            call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,err)
            call h5pset_chunk_f(prop_id,rank,chunk_dims,err)
            call h5pset_deflate_f(prop_id,compress_strength,err)
            ! create dataset
            call h5dcreate_f(file_id, name, H5T_STD_I64LE, dataspace_id,&
                                dataset_id, err, prop_id)
            ! close
            call h5dclose_f(dataset_id, err)
            call h5sclose_f(dataspace_id, err)
            call h5pclose_f(prop_id,err)
        end subroutine create_array_int_hdf5

        subroutine create_array_walks_hdf5(name)
            implicit none
            character(len=*), intent(in)    :: name
            integer(HID_T)                  :: dataset_id, dataspace_id, prop_id
            integer(HSIZE_T), allocatable   :: file_dims(:)
            integer                         :: err, rank

            rank = numpar + 1
            allocate(file_dims(rank))
            file_dims(1) = lattdim
            file_dims(2:numpar+1) = int(ext(1:numpar), kind = HSIZE_T)
            ! create dataspace
            call h5screate_simple_f(rank, file_dims, dataspace_id, err)
            ! Chunk and enable compression
            call h5pcreate_f(H5P_DATASET_CREATE_F, prop_id, err)
            call h5pset_chunk_f(prop_id, rank, file_dims, err)
            call h5pset_deflate_f(prop_id, compress_strength, err)
            ! create dataset
            call h5dcreate_f(file_id, name, H5T_STD_I32LE, dataspace_id, dataset_id, err, prop_id)
            ! close
            call h5dclose_f(dataset_id, err)
            call h5sclose_f(dataspace_id, err)
            call h5pclose_f(prop_id, err)
            print*,'walks array created'
        end subroutine create_array_walks_hdf5

        subroutine create_array_weights_hdf5(name)
            implicit none
            character(len=*), intent(in)    :: name
            integer(HID_T)                  :: dataset_id, dataspace_id, prop_id
            integer(HSIZE_T), allocatable   :: file_dims(:)
            integer                         :: err, rank
            integer(HID_T)                  :: EXT_REAL_TYPE
            EXT_REAL_TYPE = h5kind_to_type(qp,H5_REAL_KIND)

            rank = numpar - 1
            allocate(file_dims(rank))
            file_dims = int(ext(2:numpar), kind = HSIZE_T)
            ! create dataspace
            call h5screate_simple_f(rank, file_dims, dataspace_id, err)
            ! Chunk and enable compression
            call h5pcreate_f(H5P_DATASET_CREATE_F, prop_id, err)
            call h5pset_chunk_f(prop_id, rank, file_dims, err)
            call h5pset_deflate_f(prop_id, compress_strength, err)
            ! create dataset
            call h5dcreate_f(file_id, name, EXT_REAL_TYPE, dataspace_id, dataset_id, err, prop_id)
            ! close
            call h5dclose_f(dataset_id, err)
            call h5sclose_f(dataspace_id, err)
            call h5pclose_f(prop_id, err)
            print*,'weights array created'
        end subroutine create_array_weights_hdf5

    end subroutine create_datasets_hdf5

    ! ##############################################
    ! Load exiting datasets from .h5 file
    subroutine load_datasets_hdf5()
        implicit none

        call load_array_int_hdf5('Sn', samples)
        call load_array_int_hdf5('Enr', enriches)
        call load_array_int_hdf5('Pru', prunes)
        call load_array_real_hdf5('Z', Z)
        call load_array_real_hdf5('Se', effsamples)
        if (bestQ) then
            call load_array_walks_hdf5('bestWalks', bestWalks)
            call load_array_weights_hdf5('bestWeights', bestWeights)
        endif

        ! if (fixedQ .and. numpar==2) then
        !     call load_array_real_hdf5('FixMom1', FixMom1)
        !     call load_array_real_hdf5('FixMom2', FixMom2)
        ! endif
        if (fixedQ) then
            call load_array_real_hdf5('FixMom1', FixMom1)
            call load_array_real_hdf5('FixMom2', FixMom2)
            if (fixed2Q) then
                call load_array_real_hdf5('Fix2Mom1', Fix2Mom1)
                call load_array_real_hdf5('Fix2Mom2', Fix2Mom2)
            endif
        endif

        if (surfaceQ) then
            call load_array_real_hdf5('Re2W', Re2W)
            call load_array_real_hdf5('Re2Wp', Re2Wp)
        else
            call load_array_real_hdf5('Re2W', Re2W)
            ! call load_array_real_hdf5('Rm2W', Rm2W)
            ! call load_array_real_hdf5('Rg2W', Rg2W)
        endif

        if (anisoQ) then
            call load_array_real_hdf5('aniso', aniso)
            if (surfaceQ .and. lattice_type=='cub') call load_array_real_hdf5('aniso2', aniso2)
        endif

    contains
        subroutine load_array_real_hdf5(name,array)
            implicit none
            character(len=*),intent(in) :: name
            real(qp),intent(inout)      :: array(:,:,:)
            integer(HID_T)              :: dataset_id
            integer(HSIZE_T)            :: file_dims(numpar)
            integer                     :: err
            integer(HID_T)      :: EXT_REAL_TYPE
            EXT_REAL_TYPE = h5kind_to_type(qp,H5_REAL_KIND)
            file_dims=int(ext(1:numpar),kind=HSIZE_T)

            call h5dopen_f(file_id, name, dataset_id, err)
            call h5dread_f(dataset_id, EXT_REAL_TYPE, array, file_dims, err)
            call h5dclose_f(dataset_id, err)
        end subroutine load_array_real_hdf5

        subroutine load_array_int_hdf5(name,array)
            implicit none
            character(len=*),intent(in) :: name
            integer(8),intent(inout)    :: array(:,:,:)
            integer(HID_T)              :: dataset_id
            integer(HSIZE_T)            :: file_dims(numpar)
            integer                     :: err
            ! file_dims=int(Nmax+1,kind=HSIZE_T)
            file_dims=int(ext(1:numpar),kind=HSIZE_T)

            call h5dopen_f(file_id, name, dataset_id, err)
            call h5dread_f(dataset_id, H5T_STD_I64LE, array, file_dims, err)
            call h5dclose_f(dataset_id, err)
        end subroutine load_array_int_hdf5

        subroutine load_array_walks_hdf5(name,array)
            implicit none
            character(len=*),intent(in) :: name
            integer,intent(inout)       :: array(:,:,:,:)
            integer(HID_T)              :: dataset_id
            integer(HSIZE_T)            :: file_dims(numpar+1)
            integer                     :: err
            file_dims(1)=lattdim
            file_dims(2:numpar+1)=int(ext(1:numpar),kind=HSIZE_T)

            call h5dopen_f(file_id, name, dataset_id, err)
            call h5dread_f(dataset_id, H5T_STD_I32LE, array, file_dims, err)
            call h5dclose_f(dataset_id, err)
        end subroutine load_array_walks_hdf5

        subroutine load_array_weights_hdf5(name,array)
            implicit none
            character(len=*),intent(in) :: name
            real(qp),intent(inout)      :: array(:,:)
            integer(HID_T)              :: dataset_id
            integer(HSIZE_T)            :: file_dims(numpar-1)
            integer                     :: err
            integer(HID_T)      :: EXT_REAL_TYPE
            EXT_REAL_TYPE = h5kind_to_type(qp,H5_REAL_KIND)
            file_dims=int(ext(2:numpar),kind=HSIZE_T)

            call h5dopen_f(file_id, name, dataset_id, err)
            call h5dread_f(dataset_id, EXT_REAL_TYPE, array, file_dims, err)
            call h5dclose_f(dataset_id, err)
        end subroutine load_array_weights_hdf5

    end subroutine load_datasets_hdf5

    ! ##############################################
    ! Save datasets to .h5 file
    subroutine save_datasets_hdf5()
        implicit none

        write(*, '(A)', advance = 'no') '   - Arrays saved:'

        call save_array_int_hdf5('Sn', samples)
        call save_array_int_hdf5('Enr', enriches)
        call save_array_int_hdf5('Pru', prunes)
        call save_array_real_hdf5('Z', Z)
        call save_array_real_hdf5('Se', effsamples)
        if (bestQ) then
            call save_array_walks_hdf5('bestWalks', bestWalks)
            call save_array_weights_hdf5('bestWeights', bestWeights)
        endif
        ! if (percQ .and. .not. resumeQ) then
        !     call create_array_holes_hdf5('holes', int(shape(holes_list), kind=HSIZE_T))
        !     call save_array_holes_hdf5('holes', holes_list)
        ! endif

        ! if (fixedQ .and. numpar==2) then
        !     call save_array_real_hdf5('FixMom1', FixMom1)
        !     call save_array_real_hdf5('FixMom2', FixMom2)
        ! endif
        if (fixedQ) then
            call save_array_real_hdf5('FixMom1', FixMom1)
            call save_array_real_hdf5('FixMom2', FixMom2)
            if (fixed2Q) then
                call save_array_real_hdf5('Fix2Mom1', Fix2Mom1)
                call save_array_real_hdf5('Fix2Mom2', Fix2Mom2)
            endif
        endif

        if (surfaceQ) then
            call save_array_real_hdf5('Re2W', Re2W)
            call save_array_real_hdf5('Re2Wp', Re2Wp)
        else
            call save_array_real_hdf5('Re2W', Re2W)
            ! call save_array_real_hdf5('Rm2W', Rm2W)
            ! call save_array_real_hdf5('Rg2W', Rg2W)
        endif

        if (anisoQ) then
            call save_array_real_hdf5('aniso', aniso)
            if (surfaceQ .and. lattice_type=='cub') call save_array_real_hdf5('aniso2', aniso2)
        endif

        write(*,*)

    contains
        subroutine save_array_real_hdf5(name,array)
            implicit none
            character(len=*),intent(in) :: name
            real(qp),intent(in)         :: array(:,:,:)
            integer(HID_T)              :: dataset_id
            integer(HSIZE_T)            :: file_dims(numpar)
            integer                     :: err
            integer(HID_T)      :: EXT_REAL_TYPE
            EXT_REAL_TYPE = h5kind_to_type(qp,H5_REAL_KIND)
            file_dims=int(ext(1:numpar),kind=HSIZE_T)

            call h5dopen_f(file_id, name, dataset_id, err)
            call h5dwrite_f(dataset_id, EXT_REAL_TYPE, array, file_dims, err)
            call h5dclose_f(dataset_id, err)
            write(*, '(A,A)', advance = 'no') ' ', name
        end subroutine save_array_real_hdf5

        subroutine save_array_int_hdf5(name,array)
            implicit none
            character(len=*),intent(in) :: name
            integer(8),intent(in)       :: array(:,:,:)
            integer(HID_T)              :: dataset_id
            integer(HSIZE_T)            :: file_dims(numpar)
            integer                     :: err
            ! file_dims=int(Nmax+1,kind=HSIZE_T)
            file_dims=int(ext(1:numpar),kind=HSIZE_T)

            call h5dopen_f(file_id, name, dataset_id, err)
            call h5dwrite_f(dataset_id, H5T_STD_I64LE, array, file_dims, err)
            call h5dclose_f(dataset_id, err)
            write(*,'(A,A)', advance='no') ' ', name
        end subroutine save_array_int_hdf5

        subroutine create_array_holes_hdf5(name,file_dims)
            implicit none
            character(len=*),intent(in) :: name
            integer(HID_T)              :: dataset_id,dataspace_id,prop_id
            integer(HSIZE_T),intent(in) :: file_dims(2)
            integer                     :: err,rank

            rank=2
            print*,'file_dims:  ', file_dims
            ! create dataspace
            call h5screate_simple_f(rank, file_dims , dataspace_id, err)
            ! Chunk and enable compression
            call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,err)
            ! call h5pset_chunk_f(prop_id,rank,chunk_dims,err)
            ! call h5pset_deflate_f(prop_id,compress_strength,err)
            ! create dataset
            call h5dcreate_f(file_id, name, H5T_NATIVE_INTEGER, dataspace_id,&
                                dataset_id, err, prop_id)
            ! close
            call h5dclose_f(dataset_id, err)
            call h5sclose_f(dataspace_id, err)
            call h5pclose_f(prop_id,err)
        end subroutine create_array_holes_hdf5

        subroutine save_array_holes_hdf5(name,array)
            implicit none
            character(len=*),intent(in) :: name
            integer,intent(in)          :: array(:,:)
            integer(HID_T)              :: dataset_id
            integer(HSIZE_T)            :: file_dims(2)
            integer                     :: err

            file_dims=int(shape(array), kind=HSIZE_T)

            call h5dopen_f(file_id, name, dataset_id, err)
            call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, array, file_dims, err)
            call h5dclose_f(dataset_id, err)
            write(*,'(A,A)', advance='no') ' ', name
        end subroutine save_array_holes_hdf5

        subroutine save_array_walks_hdf5(name,array)
            implicit none
            character(len=*),intent(in) :: name
            integer,intent(in)          :: array(:,:,:,:)
            integer(HID_T)              :: dataset_id
            integer(HSIZE_T)            :: file_dims(numpar+1)
            integer                     :: err
            file_dims(1)=lattdim
            file_dims(2:numpar+1)=int(ext(1:numpar),kind=HSIZE_T)

            call h5dopen_f(file_id, name, dataset_id, err)
            call h5dwrite_f(dataset_id, H5T_STD_I32LE, array, file_dims, err)
            call h5dclose_f(dataset_id, err)
            write(*,'(A,A)', advance='no') ' ', name
        end subroutine save_array_walks_hdf5

        subroutine save_array_weights_hdf5(name,array)
            implicit none
            character(len=*),intent(in) :: name
            real(qp),intent(in)         :: array(:,:)
            integer(HID_T)              :: dataset_id
            integer(HSIZE_T)            :: file_dims(numpar-1)
            integer                     :: err
            integer(HID_T)      :: EXT_REAL_TYPE
            EXT_REAL_TYPE = h5kind_to_type(qp,H5_REAL_KIND)
            file_dims=int(ext(2:numpar),kind=HSIZE_T)

            call h5dopen_f(file_id, name, dataset_id, err)
            call h5dwrite_f(dataset_id, EXT_REAL_TYPE, array, file_dims, err)
            call h5dclose_f(dataset_id, err)
            write(*,'(A,A)', advance='no') ' ', name
        end subroutine save_array_weights_hdf5

    end subroutine save_datasets_hdf5

    ! ##############################################
    ! Allocate and initialise data arrays for saving sample and thermo data
    subroutine init_data_arrays()
        use lattice
        implicit none
        integer     :: err, i

        ext = 1
        ! ext(1:numpar)=Nmax+1
        do i = 1, numpar
            if ((lattice_type == 'cub' .or. lattice_type == 'tri') &     !.or. scan(params,'d')==i
                  .and. scan(params, 'm') == i ) then
                ext(i) = 2*Nmax + 1
            elseif (lattice_type == 'hex' .and. fixedQ .and. scan(params, 'a') == i .and. Nmax > 4000) then
                ext(i) = floor(Nmax/8.0,4)
            ! elseif (scan(params,'h') == i) then
            !     ext(i) = floor(alpha*Nmax) + 1
            elseif (scan(params, 'c') == i) then
                if (scan(params, 'd') > 0) then
                    ext(i) = Nmax/2 + 1
                else
                    ext(i) = 2*Nmax + 1
                endif
            elseif (scan(params, 'd') == i) then
                ext(i) = Nmax/3 + 1
            elseif (scan(params, 'l') == i) then
                if (boxQ) then
                    ext(i) = Lmax
                else
                    if (collisionsQ) then
                        ext(i) = (Nmax + 1)*(coord/2)
                    else
                        ext(i) = Nmax + 1
                    end if
                end if
            elseif (scan(params, 'e') == i) then
                ext(i) = 2
            else
                ext(i) = Nmax + 1
            endif
        enddo

        ! allocate best walks arrays
        if (bestQ) then
            allocate(bestWalks(lattdim,ext(1),ext(2),ext(3)), stat = err)
            if (err /= 0) stop 'Allocation failed.'
            allocate(bestWeights(ext(2),ext(3)), stat = err)
            if (err /= 0) stop 'Allocation failed.'
        endif

        ! if (fixedQ .and. numpar == 2) then
        !     allocate(FixMom1(ext(1),ext(2),ext(3)), stat = err)
        !     if (err /= 0) stop 'Allocation failed.'
        !     allocate(FixMom2(ext(1),ext(2),ext(3)), stat = err)
        !     if (err /= 0) stop 'Allocation failed.'
        ! endif
        if (fixedQ) then
            allocate(FixMom1(ext(1),ext(2),ext(3)), stat = err)
            if (err /= 0) stop 'Allocation failed.'
            allocate(FixMom2(ext(1),ext(2),ext(3)), stat = err)
            if (err /= 0) stop 'Allocation failed.'
            if (fixed2Q) then
                allocate(Fix2Mom1(ext(1),ext(2),ext(3)), stat = err)
                if (err /= 0) stop 'Allocation failed.'
                allocate(Fix2Mom2(ext(1),ext(2),ext(3)), stat = err)
                if (err /= 0) stop 'Allocation failed.'
            endif
        endif

        ! allocate data arrays
        allocate(samples(ext(1),ext(2),ext(3)), stat = err)
        if (err /= 0) stop 'Allocation failed.'
        allocate(prunes(ext(1),ext(2),ext(3)), stat = err)
        if (err /= 0) stop 'Allocation failed.'
        allocate(enriches(ext(1),ext(2),ext(3)), stat = err)
        if (err /= 0) stop 'Allocation failed.'
        allocate(Z(ext(1),ext(2),ext(3)), stat = err)
        if (err /= 0) stop 'Allocation failed.'
        allocate(effsamples(ext(1),ext(2),ext(3)), stat = err)
        if (err /= 0) stop 'Allocation failed.'

        if (surfaceQ) then
            allocate(Re2W(ext(1),ext(2),ext(3)), stat = err)
            if (err /= 0) stop 'Allocation failed.'
            allocate(Re2Wp(ext(1),ext(2),ext(3)), stat = err)
            if (err /= 0) stop 'Allocation failed.'
        else
            allocate(Re2W(ext(1),ext(2),ext(3)), stat = err)
            if (err /= 0) stop 'Allocation failed.'
            ! allocate(Rm2W(ext(1),ext(2),ext(3)), stat = err)
            ! if (err /= 0) stop 'Allocation failed.'
            ! allocate(Rg2W(ext(1),ext(2),ext(3)), stat = err)
            ! if (err /= 0) stop 'Allocation failed.'
        endif

        if (anisoQ) then
            allocate(aniso(ext(1),ext(2),ext(3)), stat = err)
            if (err /= 0) stop 'Allocation failed.'
            if (surfaceQ .and. lattice_type == 'cub') then
                allocate(aniso2(ext(1),ext(2),ext(3)), stat = err)
                if (err /= 0) stop 'Allocation failed.'
            endif
        endif

        ! initialise data arrays
        samples = 0
        prunes = 0
        enriches = 0
        Z = 0.0_qp
        effsamples = 0.0_qp
        if (bestQ) then
            bestWalks = 0
            bestWeights = 0.0_qp
        endif

        ! if (fixedQ .and. numpar == 2) then
        !     FixMom1 = 0.0_qp
        !     FixMom2 = 0.0_qp
        ! endif
        if (fixedQ) then
            FixMom1 = 0.0_qp
            FixMom2 = 0.0_qp
            if (fixed2Q) then
                Fix2Mom1 = 0.0_qp
                Fix2Mom2 = 0.0_qp
            endif
        endif

        if (surfaceQ) then
            Re2W = 0.0_qp
            Re2Wp = 0.0_qp
        else
            Re2W = 0.0_qp
            ! Rm2W = 0.0_qp
            ! Rg2W = 0.0_qp
        endif

        if (anisoQ) then
            aniso = 0.0_qp
            if (surfaceQ .and. lattice_type == 'cub') aniso2 = 0.0_qp
        endif

        print*, '  - Arrays initialised.'

    end subroutine init_data_arrays

end module save_data
