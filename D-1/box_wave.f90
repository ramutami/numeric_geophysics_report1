! 線形移流方程式を初期条件u(x)=sin(2\pi x)、領域0<x<1、周期境界条件で、一次元風上差分とLax_wendroffで解く
! t=1まで数値解を計算し、解析解と数値解の二乗誤差を計算し、Delta_x依存性を見る
!　というのを、\Delta_xに関するループで回す。



program main
    implicit none
    ! 解析解u_act(:),u_n_FTBS(:),u_n_Lax(:),u_ntemp(:)(計算用に一時保存しておくやつ),解析解u_analytic(:)を用意する。
    ! t_iに関してループを回し、各時刻におけるt;x,u_n_FTBS(:),u_n_Lax(:),u_analytic(:)を出力する。
    ! gnuplotで時刻毎の画像ファイルを出力して、ffmpegで動画ファイルにする。
    real(8) , allocatable :: u_act(:),u_n_FTBS(:),u_n_Lax(:),u_ntemp(:),u_n_FTCS(:),u_n_QUICKEST(:),u_n_kawakuwa(:)
    real(8) :: c_speed !移流速度(c_speed>0にする。)
    real(8) :: delta_x,delta_t,cfl_ratio !x,tのグリッド幅及びCFL条件に関する値(c_speed*delta_t/delta_x = cfl_ratio)
    real(8) :: t_min,t_max !一応、変えられるように
    real(8) :: x_min,x_max !一応、変えられるように
    integer :: meshnum_x,meshnum_t
    real(8) :: x,t
    integer :: i,j
    real(8) :: pi = 4.0_8*atan(1.0_8)
    character(len=128) :: filename

    !===========パラメタ調整エリア===================
        delta_x = 0.01_8
        cfl_ratio = 1.0_8/100.0_8
        c_speed = 0.5_8
        x_min = 0.0_8
        x_max = 1.0_8
        t_min = 0.0_8
        t_max = 10.0_8
    !=============================================

    !===========諸々のパラメタの計算=================
        !delta_tはCFL条件に従って、(cfl_ratio)*delta_x/cとすることにする。(ratio_grid<1)
        delta_t = cfl_ratio*delta_x/c_speed
        meshnum_x = (x_max-x_min)/real(delta_x,8)
        meshnum_t = (t_max-t_min)/real(delta_t,8)
        allocate(u_act(0:meshnum_x-1),u_n_FTBS(0:meshnum_x-1),u_n_Lax(0:meshnum_x-1),&
        & u_ntemp(0:meshnum_x-1),u_n_FTCS(0:meshnum_x-1),u_n_QUICKEST(0:meshnum_x-1),u_n_kawakuwa(0:meshnum_x-1))
    !=============================================

    !==========初期条件の設定======================
        do j = 0,meshnum_x-1
            x = x_min + real(j,8)*delta_x
            u_n_FTBS(j) = u_init(x)
            u_n_Lax(j) = u_init(x)
            u_n_FTCS(j) = u_init(x)
            u_n_QUICKEST(j) = u_init(x)
            u_n_kawakuwa(j) = u_init(x)
        end do
    !===========================================

    !==========主要計算部分(ただし、Delta_xは固定)====
        
        do i = 0,meshnum_t

            !書き出し
            if (mod(i,10)==0) then
                write(filename,"('./dataout/box',I6.6,'.dat')") i/10
                if (i==0) then
                    write(*,*) "first file name:",filename
                else if (i==meshnum_t) then
                    write(*,*) "last file name:",filename
                end if
                open(10,file=filename,status='replace',action='write')
                write(10,*) "#t,x,u_act(x,t),u_upwind(x,t),u_lax(x,t),u_FTCS(x,t),u_QUICKEST(x,t),u_KAWAMURA-KUWAHARA(x,t)"
                t = t_min + real(i,8) * delta_t
                write(*,*) t , "out of",t_max
                do j = 0,meshnum_x-1
                    x = x_min + real(j,8)*delta_x
                    !解析解
                    u_act(j) = u_actual(x,t,c_speed)
                    write(10,*) t,x,u_act(j),u_n_FTBS(j),u_n_Lax(j),u_n_FTCS(j),u_n_QUICKEST(j),u_n_kawakuwa(j)
                end do
            end if

            !一次元風上差分
            u_ntemp = u_n_FTBS
            call u_FTBS_update(u_ntemp,u_n_FTBS,c_speed,delta_t,delta_x)

            !lax_wendroff
            u_ntemp = u_n_Lax
            call u_Lax_update(u_ntemp,u_n_Lax,c_speed,delta_t,delta_x)

            !FTCS
            u_ntemp = u_n_FTCS
            call u_FTCS_update(u_ntemp,u_n_FTCS,c_speed,delta_t,delta_x)

            !QUICKEST
            u_ntemp = u_n_QUICKEST
            call u_QUICKEST_update(u_ntemp,u_n_QUICKEST,c_speed,delta_t,delta_x)

            !kawamura-kuwahara
            u_ntemp = u_n_kawakuwa
            call u_kawakuwa_update(u_ntemp,u_n_kawakuwa,c_speed,delta_t,delta_x)
 

            close(10)
        end do


    !=============================================





    stop
    contains
    ! 時刻t,位置xでの解析解u(x,t)を返す関数u_actual(x,t)を定義する
    ! u^n(i=1~meshnumxベクトル)をうけて、u^n+1(i=1~meshnumxベクトル)を一次元風上差分で周期境界条件のもと返す関数u_FTBS(u_n)を定義する。
    ! u^n(i=1~meshnumxベクトル)をうけて、u^n+1(i=1~meshnumxベクトル)をLax_wendroffで周期境界条件のもと返す関数u_LaxWendroff(u_n)を定義する。

    !==========解析解=============
        real(8) function u_init(x)
            !初期条件u_0(x)
            implicit none
            real(8), intent(in) :: x
            real(8) :: x_r

            x_r = modulo(x,1.0_8)
            if (x_r < 0.5_8)   then
                u_init = 1.0_8
            else
                u_init = -1.0_8
            end if
        end function u_init

        real(8) function u_actual(x,t,c)
            ! 解析解(u(x,t)=u_0(x=x-ct))
            implicit none
            real(8),intent(in)::x,t,c
            u_actual = u_init(x-t*c)
        end function u_actual
    !==========================

    !===一次元風上差分による更新========
        subroutine u_FTBS_update(u_old,u_new,c,delta_t,delta_x)
            implicit none
            real(8) , allocatable :: u_old(:)
            real(8) , allocatable :: u_new(:)
            real(8) , intent(in) :: c,delta_t,delta_x
            integer :: i
            if(.not. (size(u_old) == size(u_new))) then
                write(*,*) "array size miss match at subroutine u_FTBS_update"
            end if

            do i=lbound(u_old,1)+1,ubound(u_old,1)
                u_new(i) = u_old(i)-(c*delta_t/delta_x)*(u_old(i)-u_old(i-1))
            end do
            !周期境界
            u_new(lbound(u_new,1)) = u_old(lbound(u_old,1))-(c*delta_t/delta_x)*(u_old(lbound(u_old,1))-u_old(ubound(u_old,1)))

        end subroutine u_FTBS_update
    !===============================

    !====Lax-wendroffによる更新========
        subroutine u_Lax_update(u_old,u_new,c,delta_t,delta_x)
            implicit none
            real(8) , allocatable :: u_old(:)
            real(8) , allocatable :: u_new(:)
            real(8) , intent(in) :: c,delta_t,delta_x
            real(8) :: courant_num
            integer :: i

            courant_num = c*delta_t/delta_x

            if(.not. (size(u_old) == size(u_new))) then
                write(*,*) "array size miss match at subroutine u_LAX_update"
            end if

            do i = lbound(u_old,1)+1,ubound(u_old,1 )-1
                u_new(i) = u_old(i) - 0.5_8*courant_num*(u_old(i+1)-u_old(i-1)) + 0.5_8 * (courant_num**2.0_8) * (u_old(i+1)-2.0_8*u_old(i)+u_old(i-1))
            end do

            !周期境界（二つ！）
            u_new(lbound(u_new,1)) = u_old(lbound(u_old,1)) - 0.5_8*courant_num*(u_old(lbound(u_old,1)+1)-u_old(ubound(u_old,1))) + 0.5_8 * (courant_num**2.0_8) * (u_old(lbound(u_old,1)+1)-2.0_8*u_old(lbound(u_old,1))+u_old(ubound(u_old,1)))
            u_new(ubound(u_new,1)) = u_old(ubound(u_old,1)) - 0.5_8*courant_num*(u_old(lbound(u_old,1))-u_old(ubound(u_old,1)-1)) + 0.5_8 * (courant_num**2.0_8) * (u_old(lbound(u_old,1))-2.0_8*u_old(ubound(u_old,1))+u_old(ubound(u_old,1)-1))

        end subroutine u_Lax_update
    !================================

    !=======FTCSスキーム==============
        subroutine u_FTCS_update(u_old,u_new,c,delta_t,delta_x)
            implicit none
            real(8) , allocatable :: u_old(:)
            real(8) , allocatable :: u_new(:)
            real(8) , intent(in) :: c,delta_t,delta_x
            real(8) :: courant_num
            integer :: i

            courant_num = c*delta_t/delta_x

            if(.not. (size(u_old) == size(u_new))) then
                write(*,*) "array size miss match at subroutine u_FTCS_update"
            end if

            do i = lbound(u_old,1)+1,ubound(u_old,1 )-1
                u_new(i) = u_old(i) - 0.5_8*courant_num*(u_old(i+1)-u_old(i-1))
            end do

            !周期境界（二つ！）
            u_new(lbound(u_new,1)) = u_old(lbound(u_old,1)) - 0.5_8*courant_num*(u_old(lbound(u_old,1)+1)-u_old(ubound(u_old,1)))
            u_new(ubound(u_new,1)) = u_old(ubound(u_old,1)) - 0.5_8*courant_num*(u_old(lbound(u_old,1))-u_old(ubound(u_old,1)-1)) 

        end subroutine u_FTCS_update
    !================================

    !=======QUICKEST==============
        subroutine u_QUICKEST_update(u_old,u_new,c,delta_t,delta_x)
            implicit none
            real(8) , allocatable :: u_old(:)
            real(8) , allocatable :: u_new(:)
            real(8) , intent(in) :: c,delta_t,delta_x
            integer :: i,ip2,ip1,im2,im1
            real(8) :: c_partial_u,courant_num,courant_num2

            courant_num = c*delta_t/delta_x
            courant_num2 = courant_num*courant_num

            if(.not. (size(u_old) == size(u_new))) then
                write(*,*) "array size miss match at subroutine u_FTCS_update"
            end if

            do i = lbound(u_old,1),ubound(u_old,1 )

                !境界処理４つ！
                if (i==lbound(u_new,1)) then
                    ip1 = i+1
                    ip2 = i+2
                    im1 = ubound(u_old,1)
                    im2 = ubound(u_old,1)-1
                else if (i==lbound(u_new,1)+1) then
                    ip1 = i+1
                    ip2 = i+2
                    im1 = i-1
                    im2 = ubound(u_old,1)
                else if (i==ubound(u_new,1)) then
                    ip1 = lbound(u_old,1)
                    ip2 = lbound(u_old,1)+1
                    im1 = i-1
                    im2 = i-2
                else if (i==ubound(u_new,1)-1) then
                    ip1 = i+1
                    ip2 = lbound(u_old,1)
                    im1 = i-1
                    im2 = i-2
                else
                    ip1 = i+1
                    ip2 = i+2
                    im1 = i-1
                    im2 = i-2
                end if

                u_new(i) = courant_num*(courant_num2-1)*u_old(im2)/6.0_8 &
                & - courant_num*(courant_num2-courant_num-2.0_8)*u_old(im1)/2.0_8 &
                & + (1.0_8 + 0.5_8*courant_num*(courant_num2-2.0_8*courant_num-1.0_8))*u_old(i) &
                & - courant_num*(courant_num2-3.0_8*courant_num+2.0_8)*u_old(ip1)/6.0_8
            end do
        end subroutine u_QUICKEST_update

        !https://hydroeurope.upc.edu/wp-content/uploads/2025/04/QUICKEST-Numerical-scheme.pdf
    !================================

    !=======Kawamura-kuwahara==============
        subroutine u_kawakuwa_update(u_old,u_new,c,delta_t,delta_x)
            implicit none
            real(8) , allocatable :: u_old(:)
            real(8) , allocatable :: u_new(:)
            real(8) , intent(in) :: c,delta_t,delta_x
            integer :: i,ip2,ip1,im2,im1
            real(8) :: courant_num

            courant_num = c*delta_t/delta_x

            if(.not. (size(u_old) == size(u_new))) then
                write(*,*) "array size miss match at subroutine u_FTCS_update"
            end if

            do i = lbound(u_old,1),ubound(u_old,1 )

                !境界処理４つ！
                if (i==lbound(u_new,1)) then
                    ip1 = i+1
                    ip2 = i+2
                    im1 = ubound(u_old,1)
                    im2 = ubound(u_old,1)-1
                else if (i==lbound(u_new,1)+1) then
                    ip1 = i+1
                    ip2 = i+2
                    im1 = i-1
                    im2 = ubound(u_old,1)
                else if (i==ubound(u_new,1)) then
                    ip1 = lbound(u_old,1)
                    ip2 = lbound(u_old,1)+1
                    im1 = i-1
                    im2 = i-2
                else if (i==ubound(u_new,1)-1) then
                    ip1 = i+1
                    ip2 = lbound(u_old,1)
                    im1 = i-1
                    im2 = i-2
                else
                    ip1 = i+1
                    ip2 = i+2
                    im1 = i-1
                    im2 = i-2
                end if

                u_new(i) = u_old(i) - courant_num*(-1.0_8*u_old(ip2)+8.0_8*u_old(ip1)-8.0_8*u_old(im1)+u_old(im2))/12.0_8 &
                & - 0.25_8*courant_num*(1.0_8*u_old(ip2)-4.0_8*u_old(ip1)+6.0_8*u_old(i)-4.0_8*u_old(im1)+u_old(im2))
            end do
        end subroutine u_kawakuwa_update
        !https://pbcglab.jp/cgi-bin/wiki/?%E6%B2%B3%E6%9D%91%E3%83%BB%E6%A1%91%E5%8E%9F%E3%82%B9%E3%82%AD%E3%83%BC%E3%83%A0
    !================================


end program main
