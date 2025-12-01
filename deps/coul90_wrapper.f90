!----------------------------------------------------------------------
! Fortran wrapper for COUL90 with ISO_C_BINDING for Julia interoperability
!----------------------------------------------------------------------
module coul90_wrapper
    use iso_c_binding
    use coulfunc
    implicit none

contains

    !------------------------------------------------------------------
    ! C-callable wrapper for COUL90
    ! Computes Coulomb wave functions F, G and their derivatives
    !
    ! Arguments:
    !   x      - argument (rho = k*r)
    !   eta    - Sommerfeld parameter
    !   l      - angular momentum quantum number
    !   fc     - output: F_l(eta, x)
    !   gc     - output: G_l(eta, x)
    !   fcp    - output: F'_l(eta, x)
    !   gcp    - output: G'_l(eta, x)
    !   ifail  - output: error flag (0 = success)
    !------------------------------------------------------------------
    subroutine coul90_single(x, eta, l, fc, gc, fcp, gcp, ifail) bind(c, name='coul90_single')
        real(c_double), intent(in), value :: x
        real(c_double), intent(in), value :: eta
        integer(c_int), intent(in), value :: l
        real(c_double), intent(out) :: fc
        real(c_double), intent(out) :: gc
        real(c_double), intent(out) :: fcp
        real(c_double), intent(out) :: gcp
        integer(c_int), intent(out) :: ifail

        ! Local arrays for COUL90
        real*8 :: fc_arr(0:l), gc_arr(0:l), fcp_arr(0:l), gcp_arr(0:l)
        real*8 :: xlmin
        integer :: lrange, kfn

        xlmin = 0.0d0
        lrange = l
        kfn = 0  ! Coulomb functions

        call COUL90(x, eta, xlmin, lrange, fc_arr, gc_arr, fcp_arr, gcp_arr, kfn, ifail)

        fc = fc_arr(l)
        gc = gc_arr(l)
        fcp = fcp_arr(l)
        gcp = gcp_arr(l)

    end subroutine coul90_single

    !------------------------------------------------------------------
    ! C-callable wrapper for computing multiple l values
    !
    ! Arguments:
    !   x      - argument (rho = k*r)
    !   eta    - Sommerfeld parameter
    !   lmax   - maximum angular momentum
    !   fc     - output array: F_l(eta, x) for l = 0, ..., lmax
    !   gc     - output array: G_l(eta, x) for l = 0, ..., lmax
    !   fcp    - output array: F'_l(eta, x) for l = 0, ..., lmax
    !   gcp    - output array: G'_l(eta, x) for l = 0, ..., lmax
    !   ifail  - output: error flag (0 = success)
    !------------------------------------------------------------------
    subroutine coul90_array(x, eta, lmax, fc, gc, fcp, gcp, ifail) bind(c, name='coul90_array')
        real(c_double), intent(in), value :: x
        real(c_double), intent(in), value :: eta
        integer(c_int), intent(in), value :: lmax
        real(c_double), intent(out) :: fc(0:lmax)
        real(c_double), intent(out) :: gc(0:lmax)
        real(c_double), intent(out) :: fcp(0:lmax)
        real(c_double), intent(out) :: gcp(0:lmax)
        integer(c_int), intent(out) :: ifail

        real*8 :: xlmin
        integer :: lrange, kfn

        xlmin = 0.0d0
        lrange = lmax
        kfn = 0  ! Coulomb functions

        call COUL90(x, eta, xlmin, lrange, fc, gc, fcp, gcp, kfn, ifail)

    end subroutine coul90_array

    !------------------------------------------------------------------
    ! C-callable wrapper for Coulomb phase shifts
    !------------------------------------------------------------------
    subroutine coulph_wrapper(eta, cph, lmax) bind(c, name='coulph_wrapper')
        real(c_double), intent(in), value :: eta
        real(c_double), intent(out) :: cph(0:lmax)
        integer(c_int), intent(in), value :: lmax

        call coulph(eta, cph, lmax)

    end subroutine coulph_wrapper

end module coul90_wrapper
