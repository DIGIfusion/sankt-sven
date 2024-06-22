program main 
    use odepack_mod
    use helpers
    IMPLICIT NONE 
    type(lsoda_class) :: ls_interelm, ls_intraelm, ls_interelm_lmode
    integer, parameter :: neq=100
    integer :: itask, istate_interelm, istate_intraelm, istate_interelm_lmode
    real(dp) :: y(neq), t, tout, rtol, atol(1), x(neq)
    real(dp) :: t_eval(1001)
    real(dp) :: ysol(neq,1001)
    logical :: in_elm, in_l, in_h
    integer :: i, j, isep
    
    call linspace(0.9_dp, 1.05_dp, x)
    isep = MINLOC(abs(x - 1.0), 1)
    print *, "neq = ", neq
    call ls_interelm%initialize(TIMESTEP_INTERELM, neq, istate=istate_interelm)
    call ls_intraelm%initialize(TIMESTEP_INTRAELM, neq, istate=istate_intraelm)
    call ls_interelm_lmode%initialize(TIMESTEP_INTERELM_LMODE, neq, istate=istate_interelm_lmode)

    if (istate_interelm < 0 .or. istate_intraelm < 0 .or.  istate_interelm_lmode < 0) then
        print*,istate_interelm, istate_intraelm, istate_interelm_lmode
        error stop
    endif
    
    call linspace(0.0_dp, 1.0_dp, t_eval)

    rtol = 1.0e-5_dp
    atol = 1.0e-5_dp
    itask = 1

      
    istate_interelm = 1
    istate_intraelm = 1 
    istate_interelm_lmode = 1

    ! INITIALIZE PRESSURES?
    call linspace(100.0_dp, 1.0_dp, y)

    in_elm = .False.
    in_L   = .False.
    in_H   = .True.

    t = 0.0_dp
    do i = 1,size(t_eval)
        y(size(y)) = 1.0_dp
        y = max(y, 1.0_dp)
        tout = t_eval(i)
        ! if (i == 1) then 
        !    call ls_interelm_lmode%integrate(y, t, tout, rtol, atol, itask, istate_interelm_lmode)
        print *, "TOUT = ", tout, "nesep", y(isep)
        if (i == 1 .or. in_l) then 
            call ls_interelm_lmode%integrate(y, t, tout, rtol, atol, itask, istate_interelm_lmode)
        else if (in_h) then 
            call ls_interelm%integrate(y, t, tout, rtol, atol, itask, istate_interelm)
        else if (in_elm) then 
            call ls_intraelm%integrate(y, t, tout, rtol, atol, itask, istate_intraelm)
        end if

        if ((y(isep) > 30.0_dp) .and. (.not. in_elm) .and. (.not. in_l)) then 
        ! else if ((y(1) / y(isep) > 2.0_dp) .and. (.not. in_elm) .and. (.not. in_l)) then 
            ! ELM is triggered 
            in_elm = .True.
            in_L   = .False.
            in_H   = .False.
            ! print *, "INTRA --- ELM", x(isep), y(isep), isep, i
            print *, "IN ELM"
            ! call ls_intraelm%integrate(y, t, tout, rtol, atol, itask, istate_intraelm)
        else if (y(isep) < 25.0_dp .and. in_elm) then 
            ! BUILIDING UP 
            in_elm = .False. 
            in_L   = .False. 
            in_H   = .True.
            ! print *, "INTER +++ ELM", x(isep), y(isep), isep, i
            print *, "Finished Crashed, HMODE"
        else 
            ! print *, "ELSE", x(isep), y(isep), isep, i
        end if    


            
        if (istate_interelm < 0 .or. istate_intraelm < 0 .or.  istate_interelm_lmode < 0) then
            print*,istate_interelm, istate_intraelm, istate_interelm_lmode
            error stop
        endif
        ysol(:,i) = y(:)
    end do 
    ! Write ysol to file named fort.12
    open(unit=12, file='fort.12', status='unknown')
    ! write grid 
    write(12,*) size(t_eval), size(x), isep
    write(12,*) x
    write(12, *)
    do i = 1, size(t_eval)
        write(12,*) ysol(:,i)
    end do 

    ! do i = 1, size(t_eval)
    !     print *, ysol(:,i)
    ! 
    ! end do
CONTAINS 
    
end program main 