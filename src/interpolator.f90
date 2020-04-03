module interpolator_module

use bspline_module, only: bspline_3d
use parameters, only: SP, WP

implicit none

private
public :: interpolator

interface
    function rhs_interface(X, t) result( val )
        import :: WP
        real(WP), dimension(:)       :: X
        real(WP)                     :: t
        real(WP), dimension(size(X)) :: val
    end function
end interface

type interpolator
    private
    type(bspline_3d) :: fvx, fvy
    procedure(rhs_interface), nopass, pointer, public  :: f
    contains
        private
        procedure, public  :: init => initialize
        procedure, public  :: destroy => destroy
end type interpolator

class(interpolator), pointer :: this

contains

    subroutine destroy(this_)
        implicit none
        class(interpolator), intent(inout)     :: this_
        call this_ % fvx % destroy()
        call this_ % fvy % destroy()
    end subroutine

    subroutine initialize(this_, xc, yc, tc, gvx, gvy, order)
        implicit none
        class(interpolator), intent(inout)     :: this_
        real(WP), intent(in), dimension(:)     :: xc, yc, tc
        real(WP), intent(in), dimension(:,:,:) :: gvx, gvy
        integer,  intent(in)                   :: order
        integer                                :: iflag
        call this_ % fvx % initialize(xc, yc, tc, gvx, order, order, order, iflag, extrap = .true.)
        if (iflag /= 0) then
            print*, this_ % fvx % status_message(iflag)
        endif
        call this_ % fvy % initialize(xc, yc, tc, gvy, order, order, order, iflag, extrap = .true.)
        if (iflag /= 0) then
            print*, this_ % fvy % status_message(iflag)
        endif
        ! Assoctiate this_ % f with the evaluate function
        this_ % f => evaluate
        ! assoctiate globally accessible pointer this with this_ instance
        this => this_
    end subroutine

    function evaluate(X, t) result(V)
        implicit none
        !!!! inputs !!!!
        real(WP), dimension(2), intent(in)    :: X
        real(WP),               intent(in)    :: t

        !!!! output !!!!
        ! velocity vector to return
        real(WP), dimension(2)                :: V

        !!!! local variables !!!!
        ! status flag
        integer                               :: iflag
        ! order of derivative to return (0 returns the function value itself)
        integer, parameter                    :: d = 0

        ! Interpolate x-component of velocity
        call this%fvx%evaluate(X(1), X(2), t, d, d, d, V(1), iflag)
        if (iflag /= 0) then
            print*, this%fvx%status_message(iflag)
            print*, 'X = ', X, 't = ', t
        endif

        ! Interpolate y-component of velocity
        call this%fvy%evaluate(X(1), X(2), t, d, d, d, V(2), iflag)
        if (iflag /= 0) then
            print*, this%fvx%status_message(iflag)
            print*, 'X = ', X, 't = ', t
        endif
    end function

end module interpolator_module

! vim: ai ts=4 sts=4 et sw=4 tw=79 fenc=utf-8
