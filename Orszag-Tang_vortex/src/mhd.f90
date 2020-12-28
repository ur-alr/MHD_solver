module mhd
    implicit none
    real(8), parameter :: pi = 3.141592653589793d0
    ! グリッド数
    integer, parameter :: XN = 200, YN = 200
    ! ファイル出力回数
    integer, parameter :: PN = 100
    ! 保存変数の個数
    integer, parameter :: VN = 9
    ! 計算領域のサイズ
    real(8), parameter :: XL = 2.0d0*pi, YL = 2.0d0*pi
    real(8), parameter :: TL = pi
    ! グリッドサイズ
    real(8), parameter :: dx = XL/XN, dy = YL/YN
    ! CFL数
    real(8), parameter :: CFL = 0.4d0
    ! 比熱比
    real(8), parameter :: gam = 5.0d0/3.0d0
    ! 磁場発散抑制のためのパラメータ (Dedner et al., 2002)
    real(8), parameter :: cr = 0.18d0
    contains
    ! MINMOD関数
    real(8) function minmod(a, b)
        implicit none
        real(8), intent(in) :: a, b
        minmod = 0.5d0*(sign(1.0d0, a)+sign(1.0d0, b))*min(abs(a), abs(b))
        return
    end function minmod
    ! 中央値
    real(8) function median(a, b, c)
        implicit none
        real(8), intent(in) :: a, b, c
        median = a+minmod(b-a, c-a)
        return
    end function median
    ! x方向フラックス
    subroutine mhdfx(fx, rho, u, v, w, bx, by, bz, e, psi, pt, ch)
        implicit none
        real(8), intent(out) :: fx(VN)
        real(8), intent(in) :: rho, u, v, w, bx, by, bz, e, psi, pt, ch
        fx(1) = rho*u
        fx(2) = rho*u*u+pt-bx*bx
        fx(3) = rho*v*u-bx*by
        fx(4) = rho*w*u-bx*bz
        fx(5) = psi
        fx(6) = by*u-bx*v
        fx(7) = bz*u-bx*w
        fx(8) = (e+pt)*u-bx*(u*bx+v*by+w*bz)
        fx(9) = ch*ch*bx
    end subroutine mhdfx
    ! y方向フラックス
    subroutine mhdfy(fy, rho, u, v, w, bx, by, bz, e, psi, pt, ch)
        implicit none
        real(8), intent(out) :: fy(VN)
        real(8), intent(in) :: rho, u, v, w, bx, by, bz, e, psi, pt, ch
        fy(1) = rho*v
        fy(2) = rho*u*v-by*bx
        fy(3) = rho*v*v+pt-by*by
        fy(4) = rho*w*v-by*bz
        fy(5) = bx*v-by*u
        fy(6) = psi
        fy(7) = bz*v-by*w
        fy(8) = (e+pt)*v-by*(u*bx+v*by+w*bz)
        fy(9) = ch*ch*by
    end subroutine mhdfy
    ! HLLDリーマンソルバ (Miyoshi and Kusano, 2005)
    subroutine hlldx(hlldfx, ql, qr, ch)
        implicit none
        real(8), intent(out) :: hlldfx(VN)
        real(8), intent(in) :: ql(VN), qr(VN), ch
        real(8) :: rhol, rhor, mxl, myl, mzl, mxr, myr, mzr
        real(8) :: bxl, byl, bzl, bxr, byr, bzr, el, er, psil, psir
        real(8) :: ul, vl, wl, ur, vr, wr, pl, pr, ptl, ptr
        real(8) :: al, ar, cal, car, caxl, caxr, cml, cmr, sfl, sfr
        real(8) :: sm, um, bxm, psim, ptm
        real(8) :: rhoml, rhomr, vol, vor, wol, wor, byol, byor, bzol, bzor, eol, eor
        real(8) :: srrhoml, srrhomr, sal, sar, vi, wi, byi, bzi, eil, eir
        ! 密度
        rhol = ql(1)
        rhor = qr(1)
        ! 運動量
        mxl = ql(2)
        myl = ql(3)
        mzl = ql(4)
        mxr = qr(2)
        myr = qr(3)
        mzr = qr(4)
        ! 磁場
        bxl = ql(5)
        byl = ql(6)
        bzl = ql(7)
        bxr = qr(5)
        byr = qr(6)
        bzr = qr(7)
        ! エネルギー
        el = ql(8)
        er = qr(8)
        ! 磁場発散抑制のための人工ポテンシャル (Dedner et al., 2002)
        psil = ql(9)
        psir = qr(9)
        ! 速度
        ul = mxl/rhol
        vl = myl/rhol
        wl = mzl/rhol
        ur = mxr/rhor
        vr = myr/rhor
        wr = mzr/rhor
        ! (熱的) 圧力
        pl = (gam-1.0d0)*(el-0.5d0*rhol*(ul*ul+vl*vl+wl*wl)-0.5d0*(bxl*bxl+byl*byl+bzl*bzl))
        pr = (gam-1.0d0)*(er-0.5d0*rhor*(ur*ur+vr*vr+wr*wr)-0.5d0*(bxr*bxr+byr*byr+bzr*bzr))
        ptl = pl+0.5d0*(bxl*bxl+byl*byl+bzl*bzl)
        ! 総圧力 (熱的圧力 + 磁気圧)
        ptr = pr+0.5d0*(bxr*bxr+byr*byr+bzr*bzr)
        ! 音速
        al = sqrt(gam*pl/rhol)
        ar = sqrt(gam*pr/rhor)
        ! Alfvén速度
        cal = sqrt((bxl*bxl+byl*byl+bzl*bzl)/rhol)
        car = sqrt((bxr*bxr+byr*byr+bzr*bzr)/rhor)
        caxl = sqrt(bxl*bxl/rhol)
        caxr = sqrt(bxr*bxr/rhor)
        ! 速進磁気音波速度
        cml = sqrt(0.5d0*(al*al+cal*cal+sqrt((al*al+cal*cal)**2-4.0d0*(al*caxl)**2)))
        cmr = sqrt(0.5d0*(ar*ar+car*car+sqrt((ar*ar+car*car)**2-4.0d0*(ar*caxr)**2)))
        sfl = min(ul-cml, ur-cmr)
        sfr = max(ul+cml, ur+cmr)
        sm = ((sfr-ur)*rhor*ur-(sfl-ul)*rhol*ul-ptr+ptl)/((sfr-ur)*rhor-(sfl-ul)*rhol)
        um = sm
        bxm = bxl+0.5d0*(bxr-bxl)-0.5d0/ch*(psir-psil)
        psim = psil+0.5d0*(psir-psil)-0.5d0*ch*(bxr-bxl)
        ! リーマンファン外側
        ptm = ((sfr-ur)*rhor*ptl-(sfl-ul)*rhol*ptr+rhol*rhor*(sfr-ur)*(sfl-ul)*(ur-ul))/((sfr-ur)*rhor-(sfl-ul)*rhol)
        rhoml = rhol*(sfl-ul)/(sfl-sm)
        rhomr = rhor*(sfr-ur)/(sfr-sm)
        vol = vl-bxm*byl*(sm-ul)/(rhol*(sfl-ul)*(sfl-sm)-bxm*bxm)
        vor = vr-bxm*byr*(sm-ur)/(rhor*(sfr-ur)*(sfr-sm)-bxm*bxm)
        wol = wl-bxm*bzl*(sm-ul)/(rhol*(sfl-ul)*(sfl-sm)-bxm*bxm)
        wor = wr-bxm*bzr*(sm-ur)/(rhor*(sfr-ur)*(sfr-sm)-bxm*bxm)
        byol = byl*(rhol*(sfl-ul)*(sfl-ul)-bxm*bxm)/(rhol*(sfl-ul)*(sfl-sm)-bxm*bxm)
        byor = byr*(rhor*(sfr-ur)*(sfr-ur)-bxm*bxm)/(rhor*(sfr-ur)*(sfr-sm)-bxm*bxm)
        bzol = bzl*(rhol*(sfl-ul)*(sfl-ul)-bxm*bxm)/(rhol*(sfl-ul)*(sfl-sm)-bxm*bxm)
        bzor = bzr*(rhor*(sfr-ur)*(sfr-ur)-bxm*bxm)/(rhor*(sfr-ur)*(sfr-sm)-bxm*bxm)
        eol = ((sfl-ul)*el-ptl*ul+ptm*sm+bxm*(ul*bxl+vl*byl+wl*bzl-um*bxm-vol*byol-wol*bzol))/(sfl-sm)
        eor = ((sfr-ur)*er-ptr*ur+ptm*sm+bxm*(ur*bxr+vr*byr+wr*bzr-um*bxm-vor*byor-wor*bzor))/(sfr-sm)
        srrhoml = sqrt(rhoml)
        srrhomr = sqrt(rhomr)
        sal = sm-sqrt(bxm*bxm/rhoml)
        sar = sm+sqrt(bxm*bxm/rhomr)
        ! リーマンファン内側
        vi = (srrhoml*vol+srrhomr*vor+(byor-byol)*sign(1.0d0, bxm))/(srrhoml+srrhomr)
        wi = (srrhoml*wol+srrhomr*wor+(bzor-bzol)*sign(1.0d0, bxm))/(srrhoml+srrhomr)
        byi = (srrhoml*byor+srrhomr*byol+srrhoml*srrhomr*(vor-vol)*sign(1.0d0, bxm))/(srrhoml+srrhomr)
        bzi = (srrhoml*bzor+srrhomr*bzol+srrhoml*srrhomr*(wor-wol)*sign(1.0d0, bxm))/(srrhoml+srrhomr)
        eil = eol-srrhoml*(vol*byol+wol*bzol-vi*byi-wi*bzi)*sign(1.0d0, bxm)
        eir = eor+srrhomr*(vor*byor+wor*bzor-vi*byi-wi*bzi)*sign(1.0d0, bxm)
        if (sfl > 0) then
            call mhdfx(hlldfx, rhol, ul, vl, wl, bxm, byl, bzl, el, psim, ptl, ch)
        else if (sal > 0) then
            call mhdfx(hlldfx, rhoml, um, vol, wol, bxm, byol, bzol, eol, psim, ptm, ch)
        else if (sm > 0) then
            call mhdfx(hlldfx, rhoml, um, vi, wi, bxm, byi, bzi, eil, psim, ptm, ch)
        else if (sar > 0) then
            call mhdfx(hlldfx, rhomr, um, vi, wi, bxm, byi, bzi, eir, psim, ptm, ch)
        else if (sfr > 0) then
            call mhdfx(hlldfx, rhomr, um, vor, wor, bxm, byor, bzor, eor, psim, ptm, ch)
        else
            call mhdfx(hlldfx, rhor, ur, vr, wr, bxm, byr, bzr, er, psim, ptr, ch)
        end if
    end subroutine hlldx
    subroutine hlldy(hlldfy, ql, qr, ch)
        implicit none
        real(8), intent(out) :: hlldfy(VN)
        real(8), intent(in) :: ql(VN), qr(VN), ch
        real(8) :: rhol, rhor, mxl, myl, mzl, mxr, myr, mzr
        real(8) :: bxl, byl, bzl, bxr, byr, bzr, el, er, psil, psir
        real(8) :: ul, vl, wl, ur, vr, wr, pl, pr, ptl, ptr
        real(8) :: al, ar, cal, car, cayl, cayr, cml, cmr, sfl, sfr
        real(8) :: sm, vm, bym, psim, ptm
        real(8) :: rhoml, rhomr, wol, wor, uol, uor, bzol, bzor, bxol, bxor, eol, eor
        real(8) :: srrhoml, srrhomr, sal, sar, wi, ui, bzi, bxi, eil, eir
        ! 密度
        rhol = ql(1)
        rhor = qr(1)
        ! 運動量
        mxl = ql(2)
        myl = ql(3)
        mzl = ql(4)
        mxr = qr(2)
        myr = qr(3)
        mzr = qr(4)
        ! 磁場
        bxl = ql(5)
        byl = ql(6)
        bzl = ql(7)
        bxr = qr(5)
        byr = qr(6)
        bzr = qr(7)
        ! エネルギー
        el = ql(8)
        er = qr(8)
        ! 磁場発散抑制のための人工ポテンシャル (Dedner et al., 2002)
        psil = ql(9)
        psir = qr(9)
        ! 速度
        ul = mxl/rhol
        vl = myl/rhol
        wl = mzl/rhol
        ur = mxr/rhor
        vr = myr/rhor
        wr = mzr/rhor
        ! (熱的) 圧力
        pl = (gam-1.0d0)*(el-0.5d0*rhol*(ul*ul+vl*vl+wl*wl)-0.5d0*(bxl*bxl+byl*byl+bzl*bzl))
        pr = (gam-1.0d0)*(er-0.5d0*rhor*(ur*ur+vr*vr+wr*wr)-0.5d0*(bxr*bxr+byr*byr+bzr*bzr))
        ! 総圧力 (熱的圧力 + 磁気圧)
        ptl = pl+0.5d0*(bxl*bxl+byl*byl+bzl*bzl)
        ptr = pr+0.5d0*(bxr*bxr+byr*byr+bzr*bzr)
        ! 音速
        al = sqrt(gam*pl/rhol)
        ar = sqrt(gam*pr/rhor)
        ! Alfvén速度
        cal = sqrt((bxl*bxl+byl*byl+bzl*bzl)/rhol)
        car = sqrt((bxr*bxr+byr*byr+bzr*bzr)/rhor)
        cayl = sqrt(byl*byl/rhol)
        cayr = sqrt(byr*byr/rhor)
        ! 速進磁気音波速度
        cml = sqrt(0.5d0*(al*al+cal*cal+sqrt((al*al+cal*cal)**2-4.0d0*(al*cayl)**2)))
        cmr = sqrt(0.5d0*(ar*ar+car*car+sqrt((ar*ar+car*car)**2-4.0d0*(ar*cayr)**2)))
        sfl = min(vl-cml, vr-cmr)
        sfr = max(vl+cml, vr+cmr)
        sm = ((sfr-vr)*rhor*vr-(sfl-vl)*rhol*vl-ptr+ptl)/((sfr-vr)*rhor-(sfl-vl)*rhol)
        vm = sm
        bym = byl+0.5d0*(byr-byl)-0.5d0/ch*(psir-psil)
        psim = psil+0.5d0*(psir-psil)-0.5d0*ch*(byr-byl)
        ! リーマンファン外側
        ptm = ((sfr-vr)*rhor*ptl-(sfl-vl)*rhol*ptr+rhol*rhor*(sfr-vr)*(sfl-vl)*(vr-vl))/((sfr-vr)*rhor-(sfl-vl)*rhol)
        rhoml = rhol*(sfl-vl)/(sfl-sm)
        rhomr = rhor*(sfr-vr)/(sfr-sm)
        wol = wl-bym*bzl*(sm-vl)/(rhol*(sfl-vl)*(sfl-sm)-bym*bym)
        wor = wr-bym*bzr*(sm-vr)/(rhor*(sfr-vr)*(sfr-sm)-bym*bym)
        uol = ul-bym*bxl*(sm-vl)/(rhol*(sfl-vl)*(sfl-sm)-bym*bym)
        uor = ur-bym*bxr*(sm-vr)/(rhor*(sfr-vr)*(sfr-sm)-bym*bym)
        bzol = bzl*(rhol*(sfl-vl)*(sfl-vl)-bym*bym)/(rhol*(sfl-vl)*(sfl-sm)-bym*bym)
        bzor = bzr*(rhor*(sfr-vr)*(sfr-vr)-bym*bym)/(rhor*(sfr-vr)*(sfr-sm)-bym*bym)
        bxol = bxl*(rhol*(sfl-vl)*(sfl-vl)-bym*bym)/(rhol*(sfl-vl)*(sfl-sm)-bym*bym)
        bxor = bxr*(rhor*(sfr-vr)*(sfr-vr)-bym*bym)/(rhor*(sfr-vr)*(sfr-sm)-bym*bym)
        eol = ((sfl-vl)*el-ptl*vl+ptm*sm+bym*(ul*bxl+vl*byl+wl*bzl-uol*bxol-vm*bym-wol*bzol))/(sfl-sm)
        eor = ((sfr-vr)*er-ptr*vr+ptm*sm+bym*(ur*bxr+vr*byr+wr*bzr-uor*bxor-vm*bym-wor*bzor))/(sfr-sm)
        srrhoml = sqrt(rhoml)
        srrhomr = sqrt(rhomr)
        sal = sm-sqrt(bym*bym/rhoml)
        sar = sm+sqrt(bym*bym/rhomr)
        ! リーマンファン内側
        wi = (srrhoml*wol+srrhomr*wor+(bzor-bzol)*sign(1.0d0, bym))/(srrhoml+srrhomr)
        ui = (srrhoml*uol+srrhomr*uor+(bxor-bxol)*sign(1.0d0, bym))/(srrhoml+srrhomr)
        bzi = (srrhoml*bzor+srrhomr*bzol+srrhoml*srrhomr*(wor-wol)*sign(1.0d0, bym))/(srrhoml+srrhomr)
        bxi = (srrhoml*bxor+srrhomr*bxol+srrhoml*srrhomr*(uor-uol)*sign(1.0d0, bym))/(srrhoml+srrhomr)
        eil = eol-srrhoml*(uol*bxol+wol*bzol-ui*bxi-wi*bzi)*sign(1.0d0, bym)
        eir = eor+srrhomr*(uor*bxor+wor*bzor-ui*bxi-wi*bzi)*sign(1.0d0, bym)
        if (sfl > 0) then
            call mhdfy(hlldfy, rhol, ul, vl, wl, bxl, bym, bzl, el, psim, ptl, ch)
        else if (sal > 0) then
            call mhdfy(hlldfy, rhoml, uol, vm, wol, bxol, bym, bzol, eol, psim, ptm, ch)
        else if (sm > 0) then
            call mhdfy(hlldfy, rhoml, ui, vm, wi, bxi, bym, bzi, eil, psim, ptm, ch)
        else if (sar > 0) then
            call mhdfy(hlldfy, rhomr, ui, vm, wi, bxi, bym, bzi, eir, psim, ptm, ch)
        else if (sfr > 0) then
            call mhdfy(hlldfy, rhomr, uor, vm, wor, bxor, bym, bzor, eor, psim, ptm, ch)
        else
            call mhdfy(hlldfy, rhor, ur, vr, wr, bxr, bym, bzr, er, psim, ptr, ch)
        end if
    end subroutine hlldy
    ! 1次精度風上差分
    subroutine upwindx(ql, qr, q, i, j)
        implicit none
        real(8), intent(out) :: ql(VN), qr(VN)
        real(8), intent(in) :: q(XN, YN, VN)
        integer, intent(in) :: i, j
        integer :: il, k
        il = mod(i-2+XN, XN)+1
        do k = 1, VN
            ql(k) = q(il, j, k)
            qr(k) = q(i, j, k)
        end do
    end subroutine upwindx
    subroutine upwindy(ql, qr, q, i, j)
        implicit none
        real(8), intent(out) :: ql(VN), qr(VN)
        real(8), intent(in) :: q(XN, YN, VN)
        integer, intent(in) :: i, j
        integer :: jl, k
        jl = mod(j-2+YN, YN)+1
        do k = 1, VN
            ql(k) = q(i, jl, k)
            qr(k) = q(i, j, k)
        end do
    end subroutine upwindy
    ! 2次精度MUSCL (minmod) (van Leer, 1979)
    subroutine musclx(ql, qr, q, i, j)
        implicit none
        real(8), intent(out) :: ql(VN), qr(VN)
        real(8), intent(in) :: q(XN, YN, VN)
        integer, intent(in) :: i, j
        integer :: ill, il, ir, k
        ill = mod(i-3+XN, XN)+1
        il = mod(i-2+XN, XN)+1
        ir = mod(i, XN)+1
        do k = 1, VN
            ql(k) = q(il, j, k)+0.5d0*minmod(q(i, j, k)-q(il, j, k), q(il, j, k)-q(ill, j, k))
            qr(k) = q(i, j, k)-0.5d0*minmod(q(ir, j, k)-q(i, j, k), q(i, j, k)-q(il, j, k))
        end do
    end subroutine musclx
    subroutine muscly(ql, qr, q, i, j)
        implicit none
        real(8), intent(out) :: ql(VN), qr(VN)
        real(8), intent(in) :: q(XN, YN, VN)
        integer, intent(in) :: i, j
        integer :: jll, jl, jr, k
        jll = mod(j-3+YN, YN)+1
        jl = mod(j-2+YN, YN)+1
        jr = mod(j, YN)+1
        do k = 1, VN
            ql(k) = q(i, jl, k)+0.5d0*minmod(q(i, j, k)-q(i, jl, k), q(i, jl, k)-q(i, jll, k))
            qr(k) = q(i, j, k)-0.5d0*minmod(q(i, jr, k)-q(i, j, k), q(i, j, k)-q(i, jl, k))
        end do
    end subroutine muscly
    ! 5次精度MP5 (Suresh and Huynh, 1997)
    subroutine mp5x(ql, qr, q, i, j)
        implicit none
        real(8), intent(out) :: ql(VN), qr(VN)
        real(8), intent(in) :: q(XN, YN, VN)
        integer, intent(in) :: i, j
        integer :: illl, ill, il, ir, irr, k
        real(8) :: dll, dl, d, dr, dmml, dmm, dmmr
        real(8) :: qull, qulr, qav, qmd, qlcl, qlcr
        real(8) :: qminl, qmaxl, qminr, qmaxr, q5l, q5r
        illl = mod(i-4+XN, XN)+1
        ill = mod(i-3+XN, XN)+1
        il = mod(i-2+XN, XN)+1
        ir = mod(i, XN)+1
        irr = mod(i+1, XN)+1
        do k = 1, VN
            dll = q(illl, j, k)-2.0d0*q(ill, j, k)+q(il, j, k)
            dl = q(ill, j, k)-2.0d0*q(il, j, k)+q(i, j, k)
            d = q(il, j, k)-2.0d0*q(i, j, k)+q(ir, j, k)
            dr = q(i, j, k)-2.0d0*q(ir, j, k)+q(irr, j, k)
            dmml = minmod(dll, dl)
            dmm = minmod(dl, d)
            dmmr = minmod(d, dr)
            qull = q(il, j, k)+2.0d0*(q(il, j, k)-q(ill, j, k))
            qulr = q(i, j, k)+2.0d0*(q(i, j, k)-q(ir, j, k))
            !qull = q(il, j, k)+4.0d0*(q(il, j, k)-q(ill, j, k))
            !qulr = q(i, j, k)+4.0d0*(q(i, j, k)-q(ir, j, k))
            qav = 0.5d0*(q(il, j, k)+q(i, j, k))
            qmd = qav-0.5d0*dmm
            qlcl = q(il, j, k)+0.5d0*(q(il, j, k)-q(ill, j, k))+4.0d0/3.0d0*dmml
            qlcr = q(i, j, k)+0.5d0*(q(i, j, k)-q(ir, j, k))+4.0d0/3.0d0*dmmr
            qminl = max(min(q(il, j, k), q(i, j, k), qmd), min(q(il, j, k), qull, qlcl))
            qmaxl = min(max(q(il, j, k), q(i, j, k), qmd), max(q(il, j, k), qull, qlcl))
            qminr = max(min(q(i, j, k), q(il, j, k), qmd), min(q(i, j, k), qulr, qlcr))
            qmaxr = min(max(q(i, j, k), q(il, j, k), qmd), max(q(i, j, k), qulr, qlcr))
            q5l = (2.0d0*q(illl, j, k)-13.0d0*q(ill, j, k)+47.0d0*q(il, j, k)+27.0d0*q(i, j, k)-3.0d0*q(ir, j, k))/60.0d0
            q5r = (2.0d0*q(irr, j, k)-13.0d0*q(ir, j, k)+47.0d0*q(i, j, k)+27.0d0*q(il, j, k)-3.0d0*q(ill, j, k))/60.0d0
            ql(k) = median(q5l, qminl, qmaxl)
            qr(k) = median(q5r, qminr, qmaxr)
        end do
    end subroutine mp5x
    subroutine mp5y(ql, qr, q, i, j)
        implicit none
        real(8), intent(out) :: ql(VN), qr(VN)
        real(8), intent(in) :: q(XN, YN, VN)
        integer, intent(in) :: i, j
        integer :: jlll, jll, jl, jr, jrr, k
        real(8) :: dll, dl, d, dr, dmml, dmm, dmmr
        real(8) :: qull, qulr, qav, qmd, qlcl, qlcr
        real(8) :: qminl, qmaxl, qminr, qmaxr, q5l, q5r
        jlll = mod(j-4+YN, YN)+1
        jll = mod(j-3+YN, YN)+1
        jl = mod(j-2+YN, YN)+1
        jr = mod(j, YN)+1
        jrr = mod(j+1, YN)+1
        do k = 1, VN
            dll = q(i, jlll, k)-2.0d0*q(i, jll, k)+q(i, jl, k)
            dl = q(i, jll, k)-2.0d0*q(i, jl, k)+q(i, j, k)
            d = q(i, jl, k)-2.0d0*q(i, j, k)+q(i, jr, k)
            dr = q(i, j, k)-2.0d0*q(i, jr, k)+q(i, jrr, k)
            dmml = minmod(dll, dl)
            dmm = minmod(dl, d)
            dmmr = minmod(d, dr)
            qull = q(i, jl, k)+2.0d0*(q(i, jl, k)-q(i, jll, k))
            qulr = q(i, j, k)+2.0d0*(q(i, j, k)-q(i, jr, k))
            !qull = q(i, jl, k)+4.0d0*(q(i, jl, k)-q(i, jll, k))
            !qulr = q(i, j, k)+4.0d0*(q(i, j, k)-q(i, jr, k))
            qav = 0.5d0*(q(i, jl, k)+q(i, j, k))
            qmd = qav-0.5d0*dmm
            qlcl = q(i, jl, k)+0.5d0*(q(i, jl, k)-q(i, jll, k))+4.0d0/3.0d0*dmml
            qlcr = q(i, j, k)+0.5d0*(q(i, j, k)-q(i, jr, k))+4.0d0/3.0d0*dmmr
            qminl = max(min(q(i, jl, k), q(i, j, k), qmd), min(q(i, jl, k), qull, qlcl))
            qmaxl = min(max(q(i, jl, k), q(i, j, k), qmd), max(q(i, jl, k), qull, qlcl))
            qminr = max(min(q(i, j, k), q(i, jl, k), qmd), min(q(i, j, k), qulr, qlcr))
            qmaxr = min(max(q(i, j, k), q(i, jl, k), qmd), max(q(i, j, k), qulr, qlcr))
            q5l = (2.0d0*q(i, jlll, k)-13.0d0*q(i, jll, k)+47.0d0*q(i, jl, k)+27.0d0*q(i, j, k)-3.0d0*q(i, jr, k))/60.0d0
            q5r = (2.0d0*q(i, jrr, k)-13.0d0*q(i, jr, k)+47.0d0*q(i, j, k)+27.0d0*q(i, jl, k)-3.0d0*q(i, jll, k))/60.0d0
            ql(k) = median(q5l, qminl, qmaxl)
            qr(k) = median(q5r, qminr, qmaxr)    
        end do
    end subroutine mp5y
    ! dq/dtの計算
    subroutine dqdt(dqdtij, q, ch)
        implicit none
        real(8), intent(out) :: dqdtij(XN, YN, VN)
        real(8), intent(in) :: q(XN, YN, VN), ch
        real(8) :: ql(VN), qr(VN), fxi(VN), fyj(VN)
        real(8) :: fx(XN, YN, VN), fy(XN, YN, VN)
        integer :: i, j, ir, jr, k
        ! x方向
        !$omp parallel do private(i, ql, qr)
        do j = 1, YN
            do i = 1, XN
                !call upwindx(ql, qr, q, i, j)
                !call musclx(ql, qr, q, i, j)
                call mp5x(ql, qr, q, i, j)
                call hlldx(fx(i, j, :), ql, qr, ch)
            end do
        end do
        ! y方向
        !$omp parallel do private(i, ql, qr)
        do j = 1, YN
            do i = 1, XN
                !call upwindy(ql, qr, q, i, j)
                !call muscly(ql, qr, q, i, j)
                call mp5y(ql, qr, q, i, j)
                call hlldy(fy(i, j, :), ql, qr, ch)
            end do
        end do
        ! 保存則の計算
        !$omp parallel do private(i, ir, jr, k)
        do j = 1, YN
            jr = mod(j, YN)+1
            do i = 1, XN
                ir = mod(i, XN)+1
                do k = 1, VN
                    dqdtij(i, j, k) = -(fx(ir, j, k)-fx(i, j, k))/dx-(fy(i, jr, k)-fy(i, j, k))/dy
                end do
            end do
        end do
    end subroutine dqdt
    ! 1次精度Euler法
    subroutine euler(q1, q, dt)
        implicit none
        real(8), intent(out) :: q1(XN, YN, VN)
        real(8), intent(in) :: q(XN, YN, VN), dt
        real(8) :: ch, dqdtij(XN, YN, VN)
        ch = CFL*min(dx, dy)/dt
        call dqdt(dqdtij, q, ch)
        q1 = q+dqdtij*dt
    end subroutine euler
    ! 2次精度Runge-Kutta
    subroutine ssprk2(q2, q, dt)
        implicit none
        real(8), intent(out) :: q2(XN, YN, VN)
        real(8), intent(in) :: q(XN, YN, VN), dt
        real(8) :: ch, dqdtij(XN, YN, VN), q1(XN, YN, VN)
        ch = CFL*min(dx, dy)/dt
        call dqdt(dqdtij, q, ch)
        q1 = q+dqdtij*dt
        call dqdt(dqdtij, q1, ch)
        q2 = 0.5d0*q+0.5d0*(q1+dqdtij*dt)
    end subroutine ssprk2
    ! 3次精度Runge-Kutta
    subroutine ssprk3(q3, q, dt)
        implicit none
        real(8), intent(out) :: q3(XN, YN, VN)
        real(8), intent(in) :: q(XN, YN, VN), dt
        real(8) :: ch, dqdtij(XN, YN, VN), q1(XN, YN, VN), q2(XN, YN, VN)
        ch = CFL*min(dx, dy)/dt
        call dqdt(dqdtij, q, ch)
        q1 = q+dqdtij*dt
        call dqdt(dqdtij, q1, ch)
        q2 = 0.75d0*q+0.25d0*(q1+dqdtij*dt)
        call dqdt(dqdtij, q2, ch)
        q3 = 1.0d0/3.0d0*q+2.0d0/3.0d0*(q2+dqdtij*dt)
    end subroutine ssprk3
end module mhd
