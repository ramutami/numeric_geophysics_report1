program burgers
    ! 非粘性burgersを一次元風上差分と二段階lax_wendroffで。
    ! -1<x<1で周期境界をやると、膨張波が衝撃波匂いついた後の挙動が見えずらいので、-2<x<2でやる。
    implicit none

    real(8), parameter :: m_pi = 4*atan(1.0d0)
    integer :: meshnum_x 
    real(8) :: nu !delta_t/delta_x
    real(8) :: epsilon !人工粘性
    integer :: numstep !dx*numstemp秒だけ時間を進める。
    integer :: outputstep

    real(8) :: xmax,xmin
    real(8) :: dt,dx !meshnumx->dx->dt with respect to nu -> tfinal with respect to numstep と決まっていく。
    real(8) , allocatable :: u_exact(:),u_upwind(:),u_lax(:),x(:)
    character(len=128) :: filename

    integer :: i,istep,ix,ifile
    real(8) :: time
     
    !======parameters==========
        xmax = 2.0d0
        xmin = -2.0d0
        meshnum_x = 200
        nu = 0.2d0
        epsilon = 0.2d0
        numstep = 100000
        outputstep = 2
    !==========================

    !=======初期条件(矩形波)=======
        allocate(x(meshnum_x),u_upwind(meshnum_x),u_lax(meshnum_x),u_exact(meshnum_x))
        call setup(x,u_upwind,dx)
        call setup(x,u_lax,dx)
        call setup(x,u_exact,dx)
        dt = nu*dx
    !==========================

    !====計算and出力===========
        ifile = 0
        do istep = 0 , numstep
            time = istep*dt

            !周期で回ってきて、衝撃波が膨張波の端っこに追いついたらおしまい！
            if (xs(time) >= xmax+(-1.0d0/3.0d0-xmin)) then
                exit 
            end if

            !解析解の計算
            call update_u_exact(u_exact,x,time)

            !出力
            if (mod(istep,outputstep)==0) then
                write(*,*) istep , "out of" , numstep , "time:",time,"filenum",ifile
                write(filename,fmt='("./dataout/burgersmesh200nu0d2",i6.6,".dat")') ifile
                open(10,file=filename,status='replace',action='write')
                write(10,*) "time:x:u_exact(x):u_upwind(x):u_lax(x)"
                do ix = 1,size(x)
                    write(10,*) time,x(ix),u_exact(ix),u_upwind(ix),u_lax(ix)
                end do
                close(10)
                ifile = ifile + 1
            end if

            !更新
            call push_upwind(u_upwind,dt,dx)
            call push_lax(u_lax,dt,dx,epsilon)
            
        end do
    !==========================


    stop
    contains

    subroutine setup(x, u, dx)
        implicit none
        real(8), intent(inout) :: x(:)
        real(8), intent(inout) :: u(:)
        real(8), intent(out)   :: dx

        integer :: ix, nx

        nx = size(x)

        ! 矩形波; -2 < x < 2
        dx = (xmax-xmin)/real(nx, kind=8)

        !有限体積方的な考え方
        !nx個の格子を作る
        ! |---------o---------|--------o---------|-- ...
        ! x=-2    (ix=1)     1*dx     (ix=2)    2*dx
        ! ...--|---------o---------|--------o---------|
        !    x=1-2dx  (ix=nx-1)   x=1-dx  (ix=nx)     x=2  

        do ix = 1, nx
            x(ix) = xmin + dx*ix - dx/2.0d0 
        end do

        do ix = 1, nx
            if( x(ix) < -1.0d0/3.0d0 .or. x(ix) > +1.0d0/3.0d0 ) then
                u(ix) = 0.0d0
            else
                u(ix) = 1.0d0
            end if
        end do

    end subroutine setup

    !====解析解関係======
        real(8) function u_exact_func(x,t)
            implicit none
            real(8) , intent(in) :: t,x
            real(8) :: t_crit

            t_crit = 4.0d0/3.0d0

            if (t <= 0.0d0) then
                if (x < -1.0d0/3.0d0 .or. x > 1.0d0/3.0d0) then
                    u_exact_func = 0.0d0
                else
                    u_exact_func = 1.0d0
                end if
                return
            end if

            if (t <= t_crit) then
                if (x <= -1.0d0/3.0d0) then
                    u_exact_func = 0.0d0
                else if ( -1.0d0/3.0d0 < x .and. x <= -1.0d0/3.0d0 + t) then
                    u_exact_func = (x + 1.0d0/3.0d0)/t
                else if ( -1.0d0/3.0d0 + t < x .and. x <= 1.0d0/3.0d0 + 0.5d0*t ) then
                    u_exact_func = 1.0d0
                else 
                    u_exact_func = 0.0d0
                end if
            else 
                if (x <= -1.0d0/3.0d0) then
                    u_exact_func = 0.0d0
                else if (-1.0d0/3.0d0 < x .and. x <= xs(t) ) then 
                    u_exact_func = (x + 1.0d0/3.0d0)/t
                else
                    u_exact_func = 0.0d0
                end if
            end if

        end function u_exact_func
        
        ! 衝撃波の位置
        real(8) function xs(t)
            implicit none
            real(8),intent(in) :: t
            real(8) :: t_crit
            t_crit = 4.0d0/3.0d0

            if (t<= t_crit) then
                xs = 1.0d0/3.0d0 + 0.5d0 * t
            else
                xs = -1.0d0/3.0d0 + (2.0d0/sqrt(3.0d0))*sqrt(t)
            end if
        end function xs

        !周期境界を考慮したu_exactの更新
        subroutine update_u_exact(u,x,t)
            implicit none
            real(8), intent(inout) :: u(:),x(:)
            real(8), intent(in)    :: t
            real(8) :: xcoord

            integer :: ix

            if (xs(t) < xmax) then !衝撃波が周期で回り込む前!
                do ix = 1, size(u)
                    u(ix) = u_exact_func(x(ix),t)
                end do
            else 
                do ix = 1,size(u)
                    xcoord = x(ix)
                    if (xcoord<= -1.0d0/3.0d0) then
                        u(ix) = u_exact_func(xcoord+4.0d0,t)
                    else
                        u(ix) = u_exact_func(xcoord,t)
                    end if
                end do
            end if

        end subroutine update_u_exact
    !============

    !風上
    subroutine push_upwind(u, dt, dx)
        !非粘性
        implicit none
        real(8), intent(inout) :: u(:)
        real(8), intent(in)    :: dt
        real(8), intent(in)    :: dx

        integer :: n, ix, lbx, ubx
        real(8) :: flux(size(u))

        lbx = 1
        ubx = size(u) 

        ! 数値流束 f(ix) = f_(i+1/2) = (u_i^n)^2 / 2
        do ix = lbx, ubx
            flux(ix) = (u(ix)**2.0_8 /2.0_8)
        end do

        ! 更新
        do ix = lbx+1, ubx
            u(ix) = u(ix) - dt/dx * (flux(ix) - flux(ix-1))
        end do

        ! 境界条件
        u(lbx) = u(lbx) - dt/dx * (flux(lbx) - flux(ubx))

    end subroutine push_upwind

    !lax
    subroutine push_lax(u,dt,dx,epsilon)
        implicit none
        real(8), intent(inout) :: u(:)
        real(8), intent(in) :: dt
        real(8), intent(in) :: dx
        real(8), intent(in) :: epsilon

        integer :: n, ix, lbx, ubx
        real(8) :: flux(size(u))
        real(8) :: gamma !人工粘性
        real(8) :: uh

        lbx = 1
        ubx = size(u) 

        ! 数値流束 f(ix) = f_(i+1/2) = f(u_{i+1/2}^{n+1/2})-epsilon*|u_i+1 - ui|
        ! uh = u_{i+1/2}^{n+1/2}と書く。
        do ix = lbx, ubx-1
            uh = 0.50d0*(u(ix+1) + u(ix)) - 0.25d0*dt/dx*(u(ix+1)**2 - u(ix)**2)
            gamma = epsilon * abs(u(ix+1) - u(ix))
            flux(ix) = 0.50d0*uh**2 - gamma*dx/dt*(u(ix+1) - u(ix))
        end do
        !境界
        uh = 0.50d0*(u(lbx) + u(ubx)) - 0.25d0*dt/dx*(u(lbx)**2 - u(ubx)**2)
        gamma = epsilon * abs(u(lbx) - u(ubx))
        flux(ubx) = 0.50d0*uh**2 - gamma*dx/dt*(u(lbx) - u(ubx))

        !　更新
        do ix = lbx+1,ubx
            u(ix) = u(ix) - dt/dx * (flux(ix) - flux(ix-1))
        end do
        !境界
        u(lbx) = u(lbx) - dt/dx * (flux(lbx)-flux(ubx))

    end subroutine push_lax


end program burgers