! Define the module for the key type.
! Override the hash_value and == operator interface.
module ints_module

    implicit none

    type ints_type
        integer, allocatable :: ints(:)
    end type

    interface hash_value
        module procedure hash_value_ints
    end interface

    interface operator (==)
        module procedure ints_equal
    end interface

contains

    ! ##############################################
    ! Hashing function
    function hash_value_ints(key) result(hash)
        type(ints_type), intent(in) :: key
        integer :: hash
        integer :: i

        hash = 0
        do i = 1, size(key%ints)
            hash = xor(hash, key%ints(i) + 1640531527 + ishft(hash, 6) + ishft(hash, -2))
        enddo
    end function hash_value_ints

    ! ##############################################
    ! overloaded equality operator - integers
    function ints_equal(lhs, rhs)
        type(ints_type), intent(in) :: lhs, rhs
        logical :: ints_equal
        ! integer :: i

        if (size(lhs%ints) /= size(rhs%ints)) then
            ints_equal = .false.
            return
        endif

        ! do i = 1, size(lhs%ints)
        !     if (lhs%ints(i) /= rhs%ints(i)) then
        !         ints_equal = .false.
        !         return
        !     endif
        ! enddo
        ! ints_equal = .true.

        ints_equal=all(lhs%ints==rhs%ints)

    end function ints_equal

end module ints_module

! Define the macros needed by fhash and include fhash.f90
#define KEY_TYPE type(ints_type)
! #define VALUE_TYPE real(real64)
#define VALUE_TYPE type(ints_type)
#define KEY_USE use ints_module
#define VALUE_USE use, intrinsic :: iso_fortran_env
! #define SHORTNAME ints_double
#define SHORTNAME ints_ints
#include "fhash.f90"
