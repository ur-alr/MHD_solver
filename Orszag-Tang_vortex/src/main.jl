using Printf
using PyCall
using CSV
using DataFrames
using Formatting
plt = pyimport("matplotlib.pyplot")
cm = pyimport("matplotlib.cm")
anim = pyimport("matplotlib.animation")

# グリッド数
const XN = 200
const YN = 200
# ファイル出力回数
const PN = 100
# 保存変数の個数
const VN = 9
# 計算領域のサイズ
const XL = 2.0*pi
const YL = 2.0*pi
const TL = pi
# グリッドサイズ
const dx = XL/XN
const dy = YL/YN
# CFL数
const CFL = 0.4
# 比熱比
const gam = 5.0/3.0
# 磁場発散抑制のためのパラメータ (Dedner et al., 2002)
const cr = 0.18

# MINMOD関数
function minmod(a, b)
    return 0.5*(sign(a)+sign(b))*min(abs(a), abs(b))
end
# 中央値
function median(a, b, c)
    return a+minmod(b-a, c-a)
end
# x方向フラックス
function mhdfx(rho, u, v, w, bx, by, bz, e, psi, pt, ch)
    return [
        rho*u,
        rho*u*u+pt-bx*bx,
        rho*v*u-bx*by,
        rho*w*u-bx*bz,
        psi,
        by*u-bx*v,
        bz*u-bx*w,
        (e+pt)*u-bx*(u*bx+v*by+w*bz),
        ch*ch*bx
    ]
end
# y方向フラックス
function mhdfy(rho, u, v, w, bx, by, bz, e, psi, pt, ch)
    return [
        rho*v
        rho*u*v-by*bx
        rho*v*v+pt-by*by
        rho*w*v-by*bz
        bx*v-by*u
        psi
        bz*v-by*w
        (e+pt)*v-by*(u*bx+v*by+w*bz)
        ch*ch*by
    ]
end
# HLLDリーマンソルバ (Miyoshi and Kusano, 2005)
function hlldx(ql, qr, ch)
    # 密度
    rhol = ql[1]
    rhor = qr[1]
    # 運動量
    mxl = ql[2]
    myl = ql[3]
    mzl = ql[4]
    mxr = qr[2]
    myr = qr[3]
    mzr = qr[4]
    # 磁場
    bxl = ql[5]
    byl = ql[6]
    bzl = ql[7]
    bxr = qr[5]
    byr = qr[6]
    bzr = qr[7]
    # エネルギー
    el = ql[8]
    er = qr[8]
    # 磁場発散抑制のための人工ポテンシャル (Dedner et al., 2002)
    psil = ql[9]
    psir = qr[9]
    # 速度
    ul = mxl/rhol
    vl = myl/rhol
    wl = mzl/rhol
    ur = mxr/rhor
    vr = myr/rhor
    wr = mzr/rhor
    # (熱的) 圧力
    pl = (gam-1.0)*(el-0.5*rhol*(ul*ul+vl*vl+wl*wl)-0.5*(bxl*bxl+byl*byl+bzl*bzl))
    pr = (gam-1.0)*(er-0.5*rhor*(ur*ur+vr*vr+wr*wr)-0.5*(bxr*bxr+byr*byr+bzr*bzr))
    ptl = pl+0.5*(bxl*bxl+byl*byl+bzl*bzl)
    # 総圧力 (熱的圧力 + 磁気圧)
    ptr = pr+0.5*(bxr*bxr+byr*byr+bzr*bzr)
    # 音速
    al = sqrt(gam*pl/rhol)
    ar = sqrt(gam*pr/rhor)
    # Alfvén速度
    cal = sqrt((bxl*bxl+byl*byl+bzl*bzl)/rhol)
    car = sqrt((bxr*bxr+byr*byr+bzr*bzr)/rhor)
    caxl = sqrt(bxl*bxl/rhol)
    caxr = sqrt(bxr*bxr/rhor)
    # 速進磁気音波速度
    cfl = sqrt(0.5*(al*al+cal*cal+sqrt((al*al+cal*cal)^2-4.0*(al*caxl)^2)))
    cfr = sqrt(0.5*(ar*ar+car*car+sqrt((ar*ar+car*car)^2-4.0*(ar*caxr)^2)))
    sfl = min(ul-cfl, ur-cfr)
    sfr = max(ul+cfl, ur+cfr)
    sm = ((sfr-ur)*rhor*ur-(sfl-ul)*rhol*ul-ptr+ptl)/((sfr-ur)*rhor-(sfl-ul)*rhol)
    um = sm
    bxm = bxl+0.5*(bxr-bxl)-0.5/ch*(psir-psil)
    psim = psil+0.5*(psir-psil)-0.5*ch*(bxr-bxl)
    # リーマンファン外側
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
    # リーマンファン内側
    vi = (srrhoml*vol+srrhomr*vor+(byor-byol)*sign(bxm))/(srrhoml+srrhomr)
    wi = (srrhoml*wol+srrhomr*wor+(bzor-bzol)*sign(bxm))/(srrhoml+srrhomr)
    byi = (srrhoml*byor+srrhomr*byol+srrhoml*srrhomr*(vor-vol)*sign(bxm))/(srrhoml+srrhomr)
    bzi = (srrhoml*bzor+srrhomr*bzol+srrhoml*srrhomr*(wor-wol)*sign(bxm))/(srrhoml+srrhomr)
    eil = eol-srrhoml*(vol*byol+wol*bzol-vi*byi-wi*bzi)*sign(bxm)
    eir = eor+srrhomr*(vor*byor+wor*bzor-vi*byi-wi*bzi)*sign(bxm)
    if sfl > 0 
        return mhdfx(rhol, ul, vl, wl, bxm, byl, bzl, el, psim, ptl, ch)
    elseif sal > 0
        return mhdfx(rhoml, um, vol, wol, bxm, byol, bzol, eol, psim, ptm, ch)
    elseif sm > 0
        return mhdfx(rhoml, um, vi, wi, bxm, byi, bzi, eil, psim, ptm, ch)
    elseif sar > 0
        return mhdfx(rhomr, um, vi, wi, bxm, byi, bzi, eir, psim, ptm, ch)
    elseif sfr > 0
        return mhdfx(rhomr, um, vor, wor, bxm, byor, bzor, eor, psim, ptm, ch)
    else
        return mhdfx(rhor, ur, vr, wr, bxm, byr, bzr, er, psim, ptr, ch)
    end
end
function hlldy(ql, qr, ch)
    # 密度
    rhol = ql[1]
    rhor = qr[1]
    # 運動量
    mxl = ql[2]
    myl = ql[3]
    mzl = ql[4]
    mxr = qr[2]
    myr = qr[3]
    mzr = qr[4]
    # 磁場
    bxl = ql[5]
    byl = ql[6]
    bzl = ql[7]
    bxr = qr[5]
    byr = qr[6]
    bzr = qr[7]
    # エネルギー
    el = ql[8]
    er = qr[8]
    # 磁場発散抑制のための人工ポテンシャル (Dedner et al., 2002)
    psil = ql[9]
    psir = qr[9]
    # 速度
    ul = mxl/rhol
    vl = myl/rhol
    wl = mzl/rhol
    ur = mxr/rhor
    vr = myr/rhor
    wr = mzr/rhor
    # (熱的) 圧力
    pl = (gam-1.0)*(el-0.5*rhol*(ul*ul+vl*vl+wl*wl)-0.5*(bxl*bxl+byl*byl+bzl*bzl))
    pr = (gam-1.0)*(er-0.5*rhor*(ur*ur+vr*vr+wr*wr)-0.5*(bxr*bxr+byr*byr+bzr*bzr))
    # 総圧力 (熱的圧力 + 磁気圧)
    ptl = pl+0.5*(bxl*bxl+byl*byl+bzl*bzl)
    ptr = pr+0.5*(bxr*bxr+byr*byr+bzr*bzr)
    # 音速
    al = sqrt(gam*pl/rhol)
    ar = sqrt(gam*pr/rhor)
    # Alfvén速度
    cal = sqrt((bxl*bxl+byl*byl+bzl*bzl)/rhol)
    car = sqrt((bxr*bxr+byr*byr+bzr*bzr)/rhor)
    cayl = sqrt(byl*byl/rhol)
    cayr = sqrt(byr*byr/rhor)
    # 速進磁気音波速度
    cfl = sqrt(0.5*(al*al+cal*cal+sqrt((al*al+cal*cal)^2-4.0*(al*cayl)^2)))
    cfr = sqrt(0.5*(ar*ar+car*car+sqrt((ar*ar+car*car)^2-4.0*(ar*cayr)^2)))
    sfl = min(vl-cfl, vr-cfr)
    sfr = max(vl+cfl, vr+cfr)
    sm = ((sfr-vr)*rhor*vr-(sfl-vl)*rhol*vl-ptr+ptl)/((sfr-vr)*rhor-(sfl-vl)*rhol)
    vm = sm
    bym = byl+0.5*(byr-byl)-0.5/ch*(psir-psil)
    psim = psil+0.5*(psir-psil)-0.5*ch*(byr-byl)
    # リーマンファン外側
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
    # リーマンファン内側
    wi = (srrhoml*wol+srrhomr*wor+(bzor-bzol)*sign(bym))/(srrhoml+srrhomr)
    ui = (srrhoml*uol+srrhomr*uor+(bxor-bxol)*sign(bym))/(srrhoml+srrhomr)
    bzi = (srrhoml*bzor+srrhomr*bzol+srrhoml*srrhomr*(wor-wol)*sign(bym))/(srrhoml+srrhomr)
    bxi = (srrhoml*bxor+srrhomr*bxol+srrhoml*srrhomr*(uor-uol)*sign(bym))/(srrhoml+srrhomr)
    eil = eol-srrhoml*(uol*bxol+wol*bzol-ui*bxi-wi*bzi)*sign(bym)
    eir = eor+srrhomr*(uor*bxor+wor*bzor-ui*bxi-wi*bzi)*sign(bym)
    if sfl > 0
        return mhdfy(rhol, ul, vl, wl, bxl, bym, bzl, el, psim, ptl, ch)
    elseif sal > 0
        return mhdfy(rhoml, uol, vm, wol, bxol, bym, bzol, eol, psim, ptm, ch)
    elseif sm > 0
        return mhdfy(rhoml, ui, vm, wi, bxi, bym, bzi, eil, psim, ptm, ch)
    elseif sar > 0
        return mhdfy(rhomr, ui, vm, wi, bxi, bym, bzi, eir, psim, ptm, ch)
    elseif sfr > 0
        return mhdfy(rhomr, uor, vm, wor, bxor, bym, bzor, eor, psim, ptm, ch)
    else
        return mhdfy(rhor, ur, vr, wr, bxr, bym, bzr, er, psim, ptr, ch)
    end
end
# 1次精度風上差分
function upwindx(q, i, j)
    ql = zeros(VN); qr = zeros(VN)
    il = (i-2+XN)%XN+1
    for k in 1:VN
        ql[k] = q[il, j, k]
        qr[k] = q[i, j, k]
    end
    return ql, qr
end
function upwindy(q, i, j)
    ql = zeros(VN); qr = zeros(VN)
    jl = (j-2+YN)%YN+1
    for k in 1:VN
        ql[k] = q[i, jl, k]
        qr[k] = q[i, j, k]
    end
    return ql, qr
end
# 2次精度MUSCL (minmod) (van Leer, 1979)
function musclx(q, i, j)
    ql = zeros(VN); qr = zeros(VN)
    ill = (i-3+XN)%XN+1
    il = (i-2+XN)%XN+1
    ir = i%XN+1
    for k in 1:VN
        ql[k] = q[il, j, k]+0.5*minmod(q[i, j, k]-q[il, j, k], q[il, j, k]-q[ill, j, k])
        qr[k] = q[i, j, k]-0.5*minmod(q[ir, j, k]-q[i, j, k], q[i, j, k]-q[il, j, k])
    end
    return ql, qr
end
function muscly(q, i, j)
    ql = zeros(VN); qr = zeros(VN)
    jll = (j-3+YN)%YN+1
    jl = (j-2+YN)%YN+1
    jr = j%YN+1
    for k in 1:VN
        ql[k] = q[i, jl, k]+0.5*minmod(q[i, j, k]-q[i, jl, k], q[i, jl, k]-q[i, jll, k])
        qr[k] = q[i, j, k]-0.5*minmod(q[i, jr, k]-q[i, j, k], q[i, j, k]-q[i, jl, k])
    end
    return ql, qr
end
# 5次精度MP5 (Suresh and Huynh, 1997)
function mp5x(q, i, j)
    ql = zeros(VN); qr = zeros(VN)
    illl = (i-4+XN)%XN+1
    ill = (i-3+XN)%XN+1
    il = (i-2+XN)%XN+1
    ir = i%XN+1
    irr = (i+1)%XN+1
    for k in 1:VN
        dll = q[illl, j, k]-2.0*q[ill, j, k]+q[il, j, k]
        dl = q[ill, j, k]-2.0*q[il, j, k]+q[i, j, k]
        d = q[il, j, k]-2.0*q[i, j, k]+q[ir, j, k]
        dr = q[i, j, k]-2.0*q[ir, j, k]+q[irr, j, k]
        dmml = minmod(dll, dl)
        dmm = minmod(dl, d)
        dmmr = minmod(d, dr)
        qull = q[il, j, k]+2.0*(q[il, j, k]-q[ill, j, k])
        qulr = q[i, j, k]+2.0*(q[i, j, k]-q[ir, j, k])
        #qull = q[il, j, k]+4.0*(q[il, j, k]-q[ill, j, k])
        #qulr = q[i, j, k]+4.0*(q[i, j, k]-q[ir, j, k])
        qav = 0.5*(q[il, j, k]+q[i, j, k])
        qmd = qav-0.5*dmm
        qlcl = q[il, j, k]+0.5*(q[il, j, k]-q[ill, j, k])+4.0/3.0*dmml
        qlcr = q[i, j, k]+0.5*(q[i, j, k]-q[ir, j, k])+4.0/3.0*dmmr
        qminl = max(min(q[il, j, k], q[i, j, k], qmd), min(q[il, j, k], qull, qlcl))
        qmaxl = min(max(q[il, j, k], q[i, j, k], qmd), max(q[il, j, k], qull, qlcl))
        qminr = max(min(q[i, j, k], q[il, j, k], qmd), min(q[i, j, k], qulr, qlcr))
        qmaxr = min(max(q[i, j, k], q[il, j, k], qmd), max(q[i, j, k], qulr, qlcr))
        q5l = (2.0*q[illl, j, k]-13.0*q[ill, j, k]+47.0*q[il, j, k]+27.0*q[i, j, k]-3.0*q[ir, j, k])/60.0
        q5r = (2.0*q[irr, j, k]-13.0*q[ir, j, k]+47.0*q[i, j, k]+27.0*q[il, j, k]-3.0*q[ill, j, k])/60.0
        ql[k] = median(q5l, qminl, qmaxl)
        qr[k] = median(q5r, qminr, qmaxr)
    end
    return ql, qr
end
function mp5y(q, i, j)
    ql = zeros(VN); qr = zeros(VN)
    jlll = (j-4+YN)%YN+1
    jll = (j-3+YN)%YN+1
    jl = (j-2+YN)%YN+1
    jr = j%YN+1
    jrr = (j+1)%YN+1
    for k in 1:VN
        dll = q[i, jlll, k]-2.0*q[i, jll, k]+q[i, jl, k]
        dl = q[i, jll, k]-2.0*q[i, jl, k]+q[i, j, k]
        d = q[i, jl, k]-2.0*q[i, j, k]+q[i, jr, k]
        dr = q[i, j, k]-2.0*q[i, jr, k]+q[i, jrr, k]
        dmml = minmod(dll, dl)
        dmm = minmod(dl, d)
        dmmr = minmod(d, dr)
        qull = q[i, jl, k]+2.0*(q[i, jl, k]-q[i, jll, k])
        qulr = q[i, j, k]+2.0*(q[i, j, k]-q[i, jr, k])
        #qull = q[i, jl, k]+4.0*(q[i, jl, k]-q[i, jll, k])
        #qulr = q[i, j, k]+4.0*(q[i, j, k]-q[i, jr, k])
        qav = 0.5*(q[i, jl, k]+q[i, j, k])
        qmd = qav-0.5*dmm
        qlcl = q[i, jl, k]+0.5*(q[i, jl, k]-q[i, jll, k])+4.0/3.0*dmml
        qlcr = q[i, j, k]+0.5*(q[i, j, k]-q[i, jr, k])+4.0/3.0*dmmr
        qminl = max(min(q[i, jl, k], q[i, j, k], qmd), min(q[i, jl, k], qull, qlcl))
        qmaxl = min(max(q[i, jl, k], q[i, j, k], qmd), max(q[i, jl, k], qull, qlcl))
        qminr = max(min(q[i, j, k], q[i, jl, k], qmd), min(q[i, j, k], qulr, qlcr))
        qmaxr = min(max(q[i, j, k], q[i, jl, k], qmd), max(q[i, j, k], qulr, qlcr))
        q5l = (2.0*q[i, jlll, k]-13.0*q[i, jll, k]+47.0*q[i, jl, k]+27.0*q[i, j, k]-3.0*q[i, jr, k])/60.0
        q5r = (2.0*q[i, jrr, k]-13.0*q[i, jr, k]+47.0*q[i, j, k]+27.0*q[i, jl, k]-3.0*q[i, jll, k])/60.0
        ql[k] = median(q5l, qminl, qmaxl)
        qr[k] = median(q5r, qminr, qmaxr)    
    end
    return ql, qr
end
# dq/dtの計算
function dqdt(q, ch)
    fx = zeros(XN, YN, VN); fy = zeros(XN, YN, VN)
    # x方向
    for j in 1:YN
        for i in 1:XN
            #ql, qr = upwindx(q, i, j)
            #ql, qr = musclx(q, i, j)
            ql, qr = mp5x(q, i, j)
            fx[i, j, :] = hlldx(ql, qr, ch)
        end
    end
    # y方向
    for j in 1:YN
        for i in 1:XN
            #ql, qr = upwindy(q, i, j)
            #ql, qr = muscly(q, i, j)
            ql, qr = mp5y(q, i, j)
            fy[i, j, :] = hlldy(ql, qr, ch)
        end
    end
    dqdt = zeros(XN, YN, VN)
    # 保存則の計算
    for j in 1:YN
        jr = j%YN+1
        for i in 1:XN
            ir = i%XN+1
            for k in 1:VN
                dqdt[i, j, k] = -(fx[ir, j, k]-fx[i, j, k])/dx-(fy[i, jr, k]-fy[i, j, k])/dy
            end
        end
    end
    return dqdt
end
# 1次精度Euler法
function euler(q, dt)
    ch = CFL*min(dx, dy)/dt
    q1 = q+dqdt(q, ch)*dt
    return q1
end
# 2次精度Runge-Kutta
function ssprk2(q, dt)
    ch = CFL*min(dx, dy)/dt
    q1 = q+dqdt(q, ch)*dt
    q2 = 0.5*q+0.5*(q1+dqdt(q1, ch)*dt)
    return q2
end
# 3次精度Runge-Kutta
function ssprk3(q, dt)
    ch = CFL*min(dx, dy)/dt
    q1 = q+dqdt(q, ch)*dt
    q2 = 0.75*q+0.25*(q1+dqdt(q1, ch)*dt)
    q3 = 1.0/3.0*q+2.0/3.0*(q2+dqdt(q2, ch)*dt)
    return q3
end
function main()
    # rho: 密度
    # p:   圧力
    # u:   x方向速度
    # v:   y方向速度
    # w:   z方向速度
    # mx:  x方向運動量
    # my:  y方向運動量
    # mz:  z方向運動量
    # bx:  x方向磁場
    # by:  y方向磁場
    # bz:  z方向磁場
    # e:   エネルギー
    # psi: 磁場発散抑制のための人工ポテンシャル (Dedner et al., 2002)
    # q: 保存変数ベクトル
    # -- 1: rho
    # -- 2: mx
    # -- 3: my
    # -- 4: mz
    # -- 5: bx
    # -- 6: by
    # -- 7: bz
    # -- 8: e
    # -- 9: psi
    x = 0.0:dx:dx*(XN-1)
    y = 0.0:dy:dy*(YN-1)
    # 初期値
    rho = ones(XN, YN)*gam^2
    p = ones(XN, YN)*gam
    u = repeat(-sin.(y)', XN, 1)
    v = repeat(sin.(x), 1, YN)
    w = zeros(XN, YN)
    bx = repeat(-sin.(y)', XN, 1)
    by = repeat(sin.(2.0*x), 1, YN)
    bz = zeros(XN, YN)
    psi = zeros(XN, YN)
    e = zeros(XN, YN)
    q = zeros(XN, YN, VN); q_n = zeros(XN, YN, VN)
    for j in 1:YN
        for i in 1:XN
            e[i, j] = p[i, j]/(gam-1.0)+0.5*rho[i, j]*(u[i, j]*u[i, j]+v[i, j]*v[i, j]+w[i, j]*w[i, j])+0.5*(bx[i, j]*bx[i, j]+by[i, j]*by[i, j]+bz[i, j]*bz[i, j])
            q[i, j, 1] = rho[i, j]
            q[i, j, 2] = rho[i, j]*u[i, j]
            q[i, j, 3] = rho[i, j]*v[i, j]
            q[i, j, 4] = rho[i, j]*w[i, j]
            q[i, j, 5] = bx[i, j]
            q[i, j, 6] = by[i, j]
            q[i, j, 7] = bz[i, j]
            q[i, j, 8] = e[i, j]
            q[i, j, 9] = psi[i, j]
        end
    end
    # 計算開始
    step = 0; n = 0
    t = 0.0
    @time while t <= TL
        lxmax = 0.0; lymax = 0.0
        for j in 1:YN
            for i in 1:XN
                # 密度
                rho[i, j] = q[i, j, 1]
                # 運動量
                mx = q[i, j, 2]
                my = q[i, j, 3]
                mz = q[i, j, 4]
                # 磁場
                bx[i, j] = q[i, j, 5]
                by[i, j] = q[i, j, 6]
                bz[i, j] = q[i, j, 7]
                # エネルギー
                e[i, j] = q[i, j, 8]
                # 磁場発散抑制のための人工ポテンシャル (Dedner et al., 2002)
                psi[i, j] = q[i, j, 9]
                # 速度
                u[i, j] = mx/rho[i, j]
                v[i, j] = my/rho[i, j]
                w[i, j] = mz/rho[i, j]
                # (熱的) 圧力
                p[i, j] = (gam-1.0)*(e[i, j]-0.5*rho[i, j]*(u[i, j]*u[i, j]+v[i, j]*v[i, j]+w[i, j]*w[i, j])-0.5*(bx[i, j]*bx[i, j]+by[i, j]*by[i, j]+bz[i, j]*bz[i, j]))
                # 音速
                a = sqrt(gam*p[i, j]/rho[i, j])
                # Alfvén速度
                ca = sqrt((bx[i, j]*bx[i, j]+by[i, j]*by[i, j]+bz[i, j]*bz[i, j])/rho[i, j])
                cax = sqrt(bx[i, j]*bx[i, j]/rho[i, j])
                cay = sqrt(by[i, j]*by[i, j]/rho[i, j])
                # 速進磁気音波速度
                cfx = sqrt(0.5*(a*a+ca*ca+sqrt((a*a+ca*ca)^2-4.0*(a*cax)^2)))
                cfy = sqrt(0.5*(a*a+ca*ca+sqrt((a*a+ca*ca)^2-4.0*(a*cay)^2)))
                lxmax = max(abs(u[i, j])+cfx, lxmax)
                lymax = max(abs(v[i, j])+cfy, lymax)
            end
        end
        # CFL数をもとにdtを設定
        dt = CFL*min(dx/lxmax, dy/lymax)
        # 磁場発散抑制のためのパラメータ (Dedner et al., 2002)
        ch = CFL*min(dx, dy)/dt
        cd = exp(-dt*ch/cr)
        # 時間発展
        #q_n = euler(q, dt)
        #q_n = ssprk2(q, dt)
        q_n = ssprk3(q, dt)
        q_n[:, :, 9] = cd*q_n[:, :, 9]
        q = q_n
        # ファイル出力
        if floor(t*PN/TL) != floor((t-dt)*PN/TL)
            @printf("n: %3d, t: %.2f\n", n, t)
            #varname = ["r", "p", "u", "v", "w", "bx", "by", "bz", "ps"]
            #vardata = [rho, p, u, v, w, bx, by, bz, psi]
            #for i in 1:VN
            #    DataFrame(vardata[i]') |> CSV.write(format("../data/{1}_{2:0>3}.csv", varname[i], n))
            #end
            n += 1
        end
        t += dt
        step += 1
    end
    # 計算終了
    # 計算失敗の検知
    if !isfinite(t)
        @printf("Computaion failed: step = %3d", step)
    end
    # y = πでの値
    println("   x    rho      p      u      v      w     bx     by     bz    psi")
    for i in 1:5:XN
        j = div(YN, 2)+1
        @printf("%4.2f  %5.2f  %5.2f  %5.2f  %5.2f  %5.2f  %5.2f  %5.2f  %5.2f  %5.2f\n",
                x[i], rho[i, j], p[i, j], u[i, j], v[i, j], w[i, j], bx[i, j], by[i, j], bz[i, j], psi[i, j])
    end
    # 結果をプロット
    fig = plt.figure()
    plt.axes().set_aspect("equal")
    plt.title("p")
    plt.xlabel("X")
    plt.ylabel("Y")
    #im = plt.pcolor(x, y, p', cmap = cm.jet, shading = "auto")
    im = plt.pcolor(x, y, p', cmap = cm.jet)
    im.set_clim(0.0, 6.5)
    fig.colorbar(im)
    plt.show()
end

main()
