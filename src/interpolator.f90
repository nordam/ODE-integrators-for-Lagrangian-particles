module interpolator_module

use bspline_module, only: bspline_3d
use parameters, only: SP, WP

implicit none

private
public :: interpolator


type interpolator
    type(bspline_3d) :: fvx, fvy
    contains
        private
        procedure, public  :: eval => evaluate
        procedure, public  :: init => initialise
end type interpolator

contains

    subroutine initialise(this, xc, yc, tc, gvx, gvy, order)
        implicit none
        class(interpolator), intent(inout)     :: this
        real(WP), intent(in), dimension(:)     :: xc, yc, tc
        real(WP), intent(in), dimension(:,:,:) :: gvx, gvy
        integer,  intent(in)                   :: order
        integer                                :: iflag
        call this%fvx%initialize(xc, yc, tc, gvx, order, order, order, iflag)
        if (iflag /= 0) then
            print*, this % fvx % status_message(iflag)
        endif
        call this%fvy%initialize(xc, yc, tc, gvy, order, order, order, iflag)
        if (iflag /= 0) then
            print*, this % fvy % status_message(iflag)
        endif
    end subroutine

    function evaluate(this, X, t) result(V)
        implicit none
        ! inputs
        class(interpolator),    intent(inout) :: this
        real(WP), dimension(2), intent(in)    :: X
        real(WP),               intent(in)    :: t
        ! output
        real(WP), dimension(2)                :: V
        ! local variables
        integer                               :: iflag
        integer, parameter                    :: d = 0

        call this%fvx%evaluate(X(1), X(2), t, d, d, d, V(1), iflag)
        call this%fvy%evaluate(X(1), X(2), t, d, d, d, V(2), iflag)
    end function

end module interpolator_module

! vim: ai ts=4 sts=4 et sw=4 tw=79 fenc=utf-8
