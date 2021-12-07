    ! radius.f90
    ! Functions for calculating various radius-squared quantities
    ! The outputs should be real numbers for later multiplication with weight
    !
module radius
    use globals, only: dp, qp
    ! use lattice, only: lattdim
    implicit none
    private
    public end_to_end, monomer_monomer, centre_of_mass, radius_of_gyration, anisotropy_step, anisotropy_step_2

contains

    ! ##############################################
    ! end_to_end
    pure function end_to_end(chain) result(r2)
        implicit none
        integer,intent(in)      :: chain(:,:)       ! input
        real(qp)                :: r2               ! output
        integer                 :: last(size(chain,dim=1))

        last=chain(:,size(chain,2))
        r2 = real(dot_product(last,last),qp)

    end function end_to_end

    ! ##############################################
    ! monomer_monomer
    pure function monomer_monomer(chain) result(r2)
        implicit none
        integer,intent(in)      :: chain(:,:)       ! input
        real(qp)                :: r2               ! output
        ! integer                 :: i,dots(size(chain,2)),point(size(chain,dim=1))
        integer                 :: dots!,i

        ! dots=0
        ! do i=2,size(chain,2)    ! Origin gives no contribution
        !     dots=dots+dot_product(chain(:,i),chain(:,i))
        ! enddo

        dots=sum(chain**2)
        r2=real(dots,qp)/size(chain,2)

    end function monomer_monomer

    ! ##############################################
    ! centre_of_mass
    pure function centre_of_mass(chain) result(r2)
        implicit none
        integer,intent(in)      :: chain(:,:)       ! input
        real(qp)                :: r2               ! output
        integer                 :: Rcom(size(chain,dim=1))

        Rcom = sum(chain, 2)
        r2 = real(dot_product(Rcom,Rcom),qp)/size(chain,2)**2

    end function centre_of_mass

    ! ##############################################
    ! radius_of_gyration - not usually called explicitly
    pure function radius_of_gyration(chain) result(r2)
        implicit none
        integer,intent(in)      :: chain(:,:)       ! input
        real(qp)                :: r2               ! output

        r2 = monomer_monomer(chain) - centre_of_mass(chain)

    end function radius_of_gyration

    ! ##############################################
    ! Convert square-representation of chain to proper hexagonal chain
    pure function hex_scale(chain) result(hexchain)
        implicit none
        integer,intent(in)      :: chain(:,:)       ! input
        real(qp)                :: hexchain(size(chain,1),size(chain,2))   ! output
        integer                 :: i
        real(qp)                :: Q(size(chain,1),size(chain,1)),P(size(chain,1),size(chain,2))

        do i=1,size(chain,2)
            P(:,i)=[ 0.0_qp, -0.5_qp*modulo(i-1,2) ]
        enddo
        Q=reshape( [ sqrt(0.75_qp), 0.0_qp, 0.0_qp, 1.5_qp ], [ size(chain,1),size(chain,1) ] )
        hexchain=matmul(chain,Q)+P

    end function hex_scale

    ! ##############################################
    ! anisotropy parameter that counts steps
    pure function anisotropy_step(dirs,c) result(r2)
        implicit none
        integer,intent(in)      :: dirs(:),c      ! input
        integer                 :: a(max(c/2,2)),i
        real(qp)                :: r2               ! output

        if (c/=3) then
            do i=1,c/2
                a(i)=count(dirs==i)+count(dirs==i+c/2)
            enddo
        else    ! hex lattice
            a(1)=count(dirs==1)+count(dirs==2)
            a(2)=count(dirs==3)
        endif
            ! hex lattice
        ! a = abs(dirs)
        r2 = 1.0_qp - minval(1.0_qp*a)/maxval(a, mask=(a>0))

    end function anisotropy_step

    ! ##############################################
    ! anisotropy parameter that counts steps ONLY IN X-Y PLANE
    pure function anisotropy_step_2(dirs,c) result(r2)
        implicit none
        integer,intent(in)      :: dirs(:),c      ! input
        integer                 :: a(2),i
        real(qp)                :: r2               ! output

        if (c/=3) then
            do i=1,2
                a(i)=count(dirs==i)+count(dirs==i+c/2)
            enddo
        else    ! hex lattice
            a(1)=count(dirs==1)+count(dirs==2)
            a(2)=count(dirs==3)
        endif
            ! hex lattice
        ! a = abs(dirs)
        r2 = 1.0_qp - minval(1.0_qp*a)/maxval(a, mask=(a>0))

    end function anisotropy_step_2

end module
