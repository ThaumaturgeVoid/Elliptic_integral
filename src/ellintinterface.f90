module module_calc_ellint
    interface
        subroutine ellint(k, F, E) bind(C, name="ellint")
            use, intrinsic :: iso_c_binding
            real(c_double), intent(in) :: k
            real(c_double), intent(out) :: F,E
        end subroutine ellint
    end interface
end module module_calc_ellint

