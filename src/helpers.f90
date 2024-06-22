module helpers 
    use odepack_mod 
    IMPLICIT NONE 
    type, public :: Transport_Coefficients
    real(dp), DIMENSION(:), ALLOCATABLE :: D
    real(dp), DIMENSION(:), ALLOCATABLE :: X
    real(dp) :: V
    REAL(DP) :: U
    end type Transport_Coefficients
    type, public :: BoundaryConditions
    real(dp) :: T_EDGE
    real(dp) :: N_EDGE
    real(dp) :: Q_IN
    real(dp) :: G_IN
    end type BoundaryConditions
    REAL(DP), PARAMETER :: pi=3.14159265359_dp
CONTAINS 
SUBROUTINE ALLOCATE_TRANSPORT_COEFFS(coeffs, d1, d2, d3, x1, x2, x3, isep, x)
    TYPE(Transport_Coefficients), INTENT(INOUT) :: coeffs
    INTEGER, INTENT(IN) :: isep
    REAL(DP), DIMENSION(:), INTENT(IN) :: x
    REAL(DP), INTENT(IN) :: d1, d2, d3, x1, x2, x3
    INTEGER :: io, n
    REAL(DP) :: A1, A2
    n = size(x)
    ALLOCATE(coeffs%D(n))
    ALLOCATE(coeffs%X(n))

    io = isep - (n - isep) 
    A1 = abs(d1 - d2)
    A2 = abs(d2 - d3)

    coeffs%D(1:io) = d1
    coeffs%D(io:isep) = (A1/2.0_dp)*(sin(((x(io:isep) -x(io))  /(x(isep)-x(io)) * pi) +pi/2.0_dp)  +1.0_dp) +d2
    coeffs%D(isep:n)  = (A2/2.0_dp)*(sin(((x(isep:n)  -x(isep))/(x(n)-x(isep))*pi)    -pi/2.0_dp)  +1.0_dp) +d2

    coeffs%X(1:io) = x1
    coeffs%X(io:isep) = (abs(x1 - x2)/2.0_dp) *(sin(((x(io:isep) -x(io))  /(x(isep)-x(io))*pi)+pi/2.0_dp)+1.0_dp)+x2
    coeffs%X(isep:n)  = (abs(x2 - x3)/2.0_dp) *(sin(((x(isep:n)  -x(isep))/(x(n)-x(isep)) *pi)-pi/2.0_dp)+1.0_dp)+x2
END SUBROUTINE ALLOCATE_TRANSPORT_COEFFS
subroutine linspace(from, to, array)
    real(dp), intent(in) :: from, to
    real(dp), intent(out) :: array(:)
    real(dp) :: range
    integer :: n, i
    n = size(array)
    range = to - from

    if (n == 0) return

    if (n == 1) then
      array(1) = from
      return
    end if

    do i=1, n
      array(i) = from + range * (i - 1) / (n - 1)
    end do
  end subroutine
SUBROUTINE CENTRALDIFF(y, x, dydx)
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: DYDX
    REAL(DP), DIMENSION(:), INTENT(IN) :: y, x
    INTEGER :: i
    dydx = 0.0_dp
    do i=2, size(y)-1
        dydx(i) = (y(i+1) - y(i-1)) / ((x(i+1) - x(i-1)) / 2.0)
    end do 
    ! dydx(1) = (y(2) - y(1)) / (x(2) - x(1))
    ! dydx(size(y)) = (y(size(y)) - y(size(y)-1)) / (x(size(y)) - x(size(y)-1))
    ! dydx(1) = y(2)
    ! dydx(size(y)) = y(size(y)-1)
END SUBROUTINE CENTRALDIFF
SUBROUTINE TIMESTEP_INTERELM(self, neq, t, u, du, ierr)
    class(lsoda_class), intent(inout) :: self
        integer, intent(in) :: neq
        real(dp), intent(in) :: t
        real(dp) :: x(neq)
        real(dp), intent(in) :: u(neq)
        real(dp), intent(out) :: du(neq)
        integer, intent(out) :: ierr
        TYPE(Transport_Coefficients) :: coeffs
        TYPE(BoundaryConditions) :: bcs
        integer :: i, nx, isep
        real(dp) :: dudx(neq), GU(neq)

        call linspace(0.9_dp, 1.05_dp, x)
        isep = MINLOC(abs(x - 1.0), 1)

        ! call ALLOCATE_TRANSPORT_COEFFS(coeffs, 1.0_dp, 0.4_dp, 1.0_dp, 1.0_dp, 0.4_dp, 1.0_dp, isep, x)
        call ALLOCATE_TRANSPORT_COEFFS(coeffs, 1.0_dp, 0.01_dp, 1.0_dp, 1.0_dp, 0.01_dp, 1.0_dp, isep, x)
        coeffs%V = 0.01_dp
        coeffs%U = 0.01_dp

        ! SET BOUNDARY CONDITIONS 
        bcs%T_EDGE = 1.0_dp
        bcs%N_EDGE = 1.0_dp
        bcs%Q_IN = -300.0_dp
        bcs%G_IN = -300.0_dp

        ! u = max(u, bcs%N_EDGE)
        nx = size(u)

        ! u = MAX(u, bcs%N_EDGE)
        call CENTRALDIFF(u, x, dudx)

        dudx(1) = (bcs%G_IN - coeffs%V*u(1)) / coeffs%D(1)

        dudx(nx) = -(u(nx-1) - bcs%n_edge) / (x(nx) - x(nx-1))
        GU = coeffs%D * dudx + coeffs%V * u
        GU(1) = bcs%G_IN
        call CENTRALDIFF(GU, x, du)

        du(1) = (GU(2) - GU(1)) / (x(2) - x(1))
        du(nx) = 0.0_dp 
END SUBROUTINE TIMESTEP_INTERELM
SUBROUTINE TIMESTEP_INTERELM_LMODE(self, neq, t, u, du, ierr)
    class(lsoda_class), intent(inout) :: self
        integer, intent(in) :: neq
        real(dp), intent(in) :: t
        real(dp) :: x(neq)
        real(dp), intent(in) :: u(neq)
        real(dp), intent(out) :: du(neq)
        integer, intent(out) :: ierr
        TYPE(Transport_Coefficients) :: coeffs
        TYPE(BoundaryConditions) :: bcs
        integer :: i, nx, isep
        real(dp) :: dudx(neq), GU(neq)

        call linspace(0.9_dp, 1.05_dp, x)
        isep = MINLOC(abs(x - 1.0), 1)

        call ALLOCATE_TRANSPORT_COEFFS(coeffs, 1.0_dp, 0.4_dp, 1.0_dp, 1.0_dp, 0.4_dp, 1.0_dp, isep, x)
        ! call ALLOCATE_TRANSPORT_COEFFS(coeffs, 1.0_dp, 0.05_dp, 1.0_dp, 1.0_dp, 0.05_dp, 1.0_dp, isep, x)
        coeffs%V = 0.01_dp
        coeffs%U = 0.01_dp

        ! SET BOUNDARY CONDITIONS 
        bcs%T_EDGE = 1.0_dp
        bcs%N_EDGE = 1.0_dp
        bcs%Q_IN = -300.0_dp
        bcs%G_IN = -300.0_dp

        nx = size(u)

        ! u = MAX(u, bcs%N_EDGE)
        call CENTRALDIFF(u, x, dudx)

        dudx(1) = (bcs%G_IN - coeffs%V*u(1)) / coeffs%D(1)
        dudx(nx) = -(u(nx-1) - bcs%n_edge) / (x(nx) - x(nx-1))
        GU = coeffs%D * dudx + coeffs%V * u
        GU(1) = bcs%G_IN
        call CENTRALDIFF(GU, x, du)

        du(1) = (GU(2) - GU(1)) / (x(2) - x(1))
        du(nx) = 0.0_dp 
END SUBROUTINE TIMESTEP_INTERELM_LMODE
SUBROUTINE TIMESTEP_INTRAELM(self, neq, t, u, du, ierr)
    class(lsoda_class), intent(inout) :: self
        integer, intent(in) :: neq
        real(dp), intent(in) :: t
        real(dp) :: x(neq)
        real(dp), intent(in) :: u(neq)
        real(dp), intent(out) :: du(neq)
        integer, intent(out) :: ierr
        TYPE(Transport_Coefficients) :: coeffs
        TYPE(BoundaryConditions) :: bcs
        integer :: i, nx, isep
        real(dp) :: dudx(neq), GU(neq)

        call linspace(0.9_dp, 1.05_dp, x)
        isep = MINLOC(abs(x - 1.0), 1)

        call ALLOCATE_TRANSPORT_COEFFS(coeffs, 1.0_dp, 0.4_dp, 0.4_dp, 1.0_dp, 0.4_dp, 0.4_dp, isep, x)
        coeffs%V = 0.01_dp
        coeffs%U = 0.01_dp

        ! SET BOUNDARY CONDITIONS 
        bcs%T_EDGE = 1.0_dp
        bcs%N_EDGE = 1.0_dp
        bcs%Q_IN = -300.0_dp
        bcs%G_IN = -300.0_dp

        nx = size(u)

        ! u = MAX(u, bcs%N_EDGE)
        call CENTRALDIFF(u, x, dudx)

        dudx(1) = (bcs%G_IN - coeffs%V*u(1)) / coeffs%D(1)
        dudx(nx) = -(u(nx-1) - bcs%n_edge) / (x(nx) - x(nx-1))

        GU = coeffs%D * dudx + coeffs%V * u
        GU(1) = bcs%G_IN

        call CENTRALDIFF(GU, x, du)

        du(1) = (GU(2) - GU(1)) / (x(2) - x(1))
        du(nx) = 0.0_dp 
END SUBROUTINE TIMESTEP_INTRAELM  
end module helpers 