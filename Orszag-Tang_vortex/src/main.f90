program main
    use mhd
    implicit none
    ! rho: 密度
    ! p:   圧力
    ! u:   x方向速度
    ! v:   y方向速度
    ! w:   z方向速度
    ! mx:  x方向運動量
    ! my:  y方向運動量
    ! mz:  z方向運動量
    ! bx:  x方向磁場
    ! by:  y方向磁場
    ! bz:  z方向磁場
    ! e:   エネルギー
    ! psi: 磁場発散抑制のための人工ポテンシャル (Dedner et al., 2002)
    ! q: 保存変数ベクトル
    ! -- 1: rho
    ! -- 2: mx
    ! -- 3: my
    ! -- 4: mz
    ! -- 5: bx
    ! -- 6: by
    ! -- 7: bz
    ! -- 8: e
    ! -- 9: psi
    real(8) :: rho(XN, YN), p(XN, YN), u(XN, YN), v(XN, YN), w(XN, YN)
    real(8) :: mx(XN, YN), my(XN, YN), mz(XN, YN)
    real(8) :: bx(XN, YN), by(XN, YN), bz(XN, YN), psi(XN, YN), e(XN, YN)
    real(8) :: q(XN, YN, VN), q_n(XN, YN, VN)
    real(8) :: x, y, t
    integer :: i, j, k, n, step
    real(8) :: a, ca, cax, cay, cfx, cfy, lxmax, lymax
    real(8) :: dt, ch, cd
    character(64) :: var(VN), dir, name(VN), footer
    ! 初期値
    !$omp parallel do private(i, x, y)
    do j = 1, YN
        y = (j-1)*dy
        do i = 1, XN
            x = (i-1)*dx
            rho(i, j) = gam*gam
            p(i, j) = gam
            u(i, j) = -sin(y)
            v(i, j) = sin(x)
            w(i, j) = 0.0d0
            bx(i, j) = -sin(y)
            by(i, j) = sin(2.0d0*x)
            bz(i, j) = 0.0d0
            e(i, j) = p(i, j)/(gam-1.0d0) &
                     +0.5d0*rho(i, j)*(u(i, j)*u(i, j)+v(i, j)*v(i, j)+w(i, j)*w(i, j)) &
                     +0.5d0*(bx(i, j)*bx(i, j)+by(i, j)*by(i, j)+bz(i, j)*bz(i, j))
            psi(i, j) = 0.0d0
            q(i, j, 1) = rho(i, j)
            q(i, j, 2) = rho(i, j)*u(i, j)
            q(i, j, 3) = rho(i, j)*v(i, j)
            q(i, j, 4) = rho(i, j)*w(i, j)
            q(i, j, 5) = bx(i, j)
            q(i, j, 6) = by(i, j)
            q(i, j, 7) = bz(i, j)
            q(i, j, 8) = e(i, j)
            q(i, j, 9) = psi(i, j)
        end do
    end do
    ! 計算開始
    step = 0; n = 0
    t = 0.0d0
    do while (t <= TL)
    !do while (n <= 0)
        lxmax = 0.0d0; lymax = 0.0d0
        !$omp parallel do private(i, a, ca, cax, cay, cfx, cfy)
        do j = 1, YN
            do i = 1, XN
                ! 密度
                rho(i, j) = q(i, j, 1)
                ! 運動量
                mx(i, j) = q(i, j, 2)
                my(i, j) = q(i, j, 3)
                mz(i, j) = q(i, j, 4)
                ! 磁場
                bx(i, j) = q(i, j, 5)
                by(i, j) = q(i, j, 6)
                bz(i, j) = q(i, j, 7)
                ! エネルギー
                e(i, j) = q(i, j, 8)
                ! 磁場発散抑制のための人工ポテンシャル (Dedner et al., 2002)
                psi(i, j) = q(i, j, 9)
                ! 速度
                u(i, j) = mx(i, j)/rho(i, j)
                v(i, j) = my(i, j)/rho(i, j)
                w(i, j) = mz(i, j)/rho(i, j)
                ! (熱的) 圧力
                p(i, j) = (gam-1.0d0)*(e(i, j) &
                         -0.5d0*rho(i, j)*(u(i, j)*u(i, j)+v(i, j)*v(i, j)+w(i, j)*w(i, j)) &
                         -0.5d0*(bx(i, j)*bx(i, j)+by(i, j)*by(i, j)+bz(i, j)*bz(i, j)))
                ! 音速
                a = sqrt(gam*p(i, j)/rho(i, j))
                ! Alfvén速度
                ca = sqrt((bx(i, j)*bx(i, j)+by(i, j)*by(i, j)+bz(i, j)*bz(i, j))/rho(i, j))
                cax = sqrt(bx(i, j)*bx(i, j)/rho(i, j))
                cay = sqrt(by(i, j)*by(i, j)/rho(i, j))
                ! 速進磁気音波速度
                cfx = sqrt(0.5d0*(a*a+ca*ca+sqrt((a*a+ca*ca)**2-4.0d0*(a*cax)**2)))
                cfy = sqrt(0.5d0*(a*a+ca*ca+sqrt((a*a+ca*ca)**2-4.0d0*(a*cay)**2)))
                !$omp critical
                lxmax = max(abs(u(i, j))+cfx, lxmax)
                lymax = max(abs(v(i, j))+cfy, lymax)
                !$omp end critical
            end do
        end do
        ! CFL数をもとにdtを設定
        dt = CFL*min(dx/lxmax, dy/lymax)
        ! 磁場発散抑制のためのパラメータ (Dedner et al., 2002)
        ch = CFL*min(dx, dy)/dt
        cd = exp(-dt*ch/cr)
        ! 時間発展
        !call euler(q_n, q, dt)
        !call ssprk2(q_n, q, dt)
        call ssprk3(q_n, q, dt)
        q_n(:, :, 9) = cd*q_n(:, :, 9)
        q = q_n
        ! ファイル出力
        if (int(t*PN/TL+1) .ne. int((t-dt)*PN/TL+1)) then
            write(*, '("n: ", i3, ", t: ", f4.2)') n, t
            var(1) = "r"; var(2) = "p"
            var(3) = "u"; var(4) = "v"; var(5) = "w"
            var(6) = "bx"; var(7) = "by"; var(8) = "bz"
            var(9) = "ps"
            dir = "../data/"
            write(footer, '("_", i3.3, ".csv")') n
            do k = 1, VN
                name(k) = trim(dir)//trim(var(k))//trim(footer)
                print *, name(k)
                open(10+k, file = name(k))
            end do
            do j = 1, YN
                do i = 1, XN-1
                    write(11, '(e12.4, ",")', advance = "no") rho(i, j)
                    write(12, '(e12.4, ",")', advance = "no") p(i, j)
                    write(13, '(e12.4, ",")', advance = "no") u(i, j)
                    write(14, '(e12.4, ",")', advance = "no") v(i, j)
                    write(15, '(e12.4, ",")', advance = "no") w(i, j)
                    write(16, '(e12.4, ",")', advance = "no") bx(i, j)
                    write(17, '(e12.4, ",")', advance = "no") by(i, j)
                    write(18, '(e12.4, ",")', advance = "no") bz(i, j)
                    write(19, '(e12.4, ",")', advance = "no") psi(i, j)
                end do
                write(11, '(e12.4)') rho(XN, j)
                write(12, '(e12.4)') p(XN, j)
                write(13, '(e12.4)') u(XN, j)
                write(14, '(e12.4)') v(XN, j)
                write(15, '(e12.4)') w(XN, j)
                write(16, '(e12.4)') bx(XN, j)
                write(17, '(e12.4)') by(XN, j)
                write(18, '(e12.4)') bz(XN, j)
                write(19, '(e12.4)') psi(XN, j)
            end do
            do k = 1, VN
                close(10+k)
            end do
            n = n+1
        end if
        t = t+dt
        step = step+1
    end do
    ! 計算終了
    ! 計算失敗の検知
    if (isnan(t-t)) then
        write(*, '("Computaion failed: step = ", i3)') step
    end if
    ! y = πでの値
    print *, "  x    rho      p      u      v      w     bx     by     bz    psi"
    do i = 1, XN, 5
        j = YN/2+1
        write(*, '(f4.2, 9("  ", f5.2))') &
                dx*(i-1), rho(i, j), p(i, j), u(i, j), v(i, j), w(i, j), bx(i, j), by(i, j), bz(i, j), psi(i, j)
    end do
end program main
