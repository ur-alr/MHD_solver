import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# グリッド数
#XN = 200; YN = 200
XN = 100; YN = 100
# ファイル出力回数
PN = 100
# 保存変数の個数
VN = 9
# 計算領域のサイズ
XL = 2.0*np.pi; YL = 2.0*np.pi
TL = np.pi
# グリッドサイズ
dx = XL/XN; dy = YL/YN
# CFL数
CFL = 0.4
# 比熱比
gam = 5.0/3.0
# 磁場発散抑制のためのパラメータ (Dedner et al., 2002)
cr = 0.18

# MINMOD関数
def minmod(a, b):
    return 0.5*(np.sign(a)+np.sign(b))*np.minimum(np.abs(a), np.abs(b))
# 中央値
def median(a, b, c):
    return a+minmod(b-a, c-a)
# x方向フラックス
def mhdfx(rho, u, v, w, bx, by, bz, e, psi, pt, ch):
    fx = np.zeros((XN, YN, VN))
    fx[:, :, 0] = rho*u
    fx[:, :, 1] = rho*u*u+pt-bx*bx
    fx[:, :, 2] = rho*v*u-bx*by
    fx[:, :, 3] = rho*w*u-bx*bz
    fx[:, :, 4] = psi
    fx[:, :, 5] = by*u-bx*v
    fx[:, :, 6] = bz*u-bx*w
    fx[:, :, 7] = (e+pt)*u-bx*(u*bx+v*by+w*bz)
    fx[:, :, 8] = ch*ch*bx
    return fx
# y方向フラックス
def mhdfy(rho, u, v, w, bx, by, bz, e, psi, pt, ch):
    fy = np.zeros((XN, YN, VN))
    fy[:, :, 0] = rho*v
    fy[:, :, 1] = rho*u*v-by*bx
    fy[:, :, 2] = rho*v*v+pt-by*by
    fy[:, :, 3] = rho*w*v-by*bz
    fy[:, :, 4] = bx*v-by*u
    fy[:, :, 5] = psi
    fy[:, :, 6] = bz*v-by*w
    fy[:, :, 7] = (e+pt)*v-by*(u*bx+v*by+w*bz)
    fy[:, :, 8] = ch*ch*by
    return fy
# HLLDリーマンソルバ (Miyoshi and Kusano, 2005)
def hlldx(ql, qr, ch):
    # 密度
    rhol = ql[:, :, 0]
    rhor = qr[:, :, 0]
    # 運動量
    mxl = ql[:, :, 1]
    myl = ql[:, :, 2]
    mzl = ql[:, :, 3]
    mxr = qr[:, :, 1]
    myr = qr[:, :, 2]
    mzr = qr[:, :, 3]
    # 磁場
    bxl = ql[:, :, 4]
    byl = ql[:, :, 5]
    bzl = ql[:, :, 6]
    bxr = qr[:, :, 4]
    byr = qr[:, :, 5]
    bzr = qr[:, :, 6]
    # エネルギー
    el = ql[:, :, 7]
    er = qr[:, :, 7]
    # 磁場発散抑制のための人工ポテンシャル (Dedner et al., 2002)
    psil = ql[:, :, 8]
    psir = qr[:, :, 8]
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
    al = np.sqrt(gam*pl/rhol)
    ar = np.sqrt(gam*pr/rhor)
    # Alfvén速度
    cal = np.sqrt((bxl*bxl+byl*byl+bzl*bzl)/rhol)
    car = np.sqrt((bxr*bxr+byr*byr+bzr*bzr)/rhor)
    caxl = np.sqrt(bxl*bxl/rhol)
    caxr = np.sqrt(bxr*bxr/rhor)
    # 速進磁気音波速度
    cfl = np.sqrt(0.5*(al*al+cal*cal+np.sqrt((al*al+cal*cal)**2-4.0*(al*caxl)**2)))
    cfr = np.sqrt(0.5*(ar*ar+car*car+np.sqrt((ar*ar+car*car)**2-4.0*(ar*caxr)**2)))
    sfl = np.minimum(ul-cfl, ur-cfr)
    sfr = np.maximum(ul+cfl, ur+cfr)
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
    srrhoml = np.sqrt(rhoml)
    srrhomr = np.sqrt(rhomr)
    sal = sm-np.sqrt(bxm*bxm/rhoml)
    sar = sm+np.sqrt(bxm*bxm/rhomr)
    # リーマンファン内側
    vi = (srrhoml*vol+srrhomr*vor+(byor-byol)*np.sign(bxm))/(srrhoml+srrhomr)
    wi = (srrhoml*wol+srrhomr*wor+(bzor-bzol)*np.sign(bxm))/(srrhoml+srrhomr)
    byi = (srrhoml*byor+srrhomr*byol+srrhoml*srrhomr*(vor-vol)*np.sign(bxm))/(srrhoml+srrhomr)
    bzi = (srrhoml*bzor+srrhomr*bzol+srrhoml*srrhomr*(wor-wol)*np.sign(bxm))/(srrhoml+srrhomr)
    eil = eol-srrhoml*(vol*byol+wol*bzol-vi*byi-wi*bzi)*np.sign(bxm)
    eir = eor+srrhomr*(vor*byor+wor*bzor-vi*byi-wi*bzi)*np.sign(bxm)
    fxl = mhdfx(rhol, ul, vl, wl, bxm, byl, bzl, el, psim, ptl, ch)
    fxol = mhdfx(rhoml, um, vol, wol, bxm, byol, bzol, eol, psim, ptm, ch)
    fxil = mhdfx(rhoml, um, vi, wi, bxm, byi, bzi, eil, psim, ptm, ch)
    fxir = mhdfx(rhomr, um, vi, wi, bxm, byi, bzi, eir, psim, ptm, ch)
    fxor = mhdfx(rhomr, um, vor, wor, bxm, byor, bzor, eor, psim, ptm, ch)
    fxr = mhdfx(rhor, ur, vr, wr, bxm, byr, bzr, er, psim, ptr, ch)
    fx = np.zeros((XN, YN, VN))
    slindex = sfl > 0.0
    solindex = (0.0 >= sfl) & (sal > 0.0)
    silindex = (0.0 >= sal) & (sm > 0.0)
    sirindex = (0.0 >= sm) & (sar > 0.0)
    sorindex = (0.0 >= sar) & (sfr > 0.0)
    srindex = 0.0 >= sfr
    fx[slindex, :] = fxl[slindex, :]
    fx[solindex, :] = fxol[solindex, :]
    fx[silindex, :] = fxil[silindex, :]
    fx[sirindex, :] = fxir[sirindex, :]
    fx[sorindex, :] = fxor[sorindex, :]
    fx[srindex, :] = fxr[srindex, :]
    return fx
def hlldy(ql, qr, ch):
    # 密度
    rhol = ql[:, :, 0]
    rhor = qr[:, :, 0]
    # 運動量
    mxl = ql[:, :, 1]
    myl = ql[:, :, 2]
    mzl = ql[:, :, 3]
    mxr = qr[:, :, 1]
    myr = qr[:, :, 2]
    mzr = qr[:, :, 3]
    # 磁場
    bxl = ql[:, :, 4]
    byl = ql[:, :, 5]
    bzl = ql[:, :, 6]
    bxr = qr[:, :, 4]
    byr = qr[:, :, 5]
    bzr = qr[:, :, 6]
    # エネルギー
    el = ql[:, :, 7]
    er = qr[:, :, 7]
    # 磁場発散抑制のための人工ポテンシャル (Dedner et al., 2002)
    psil = ql[:, :, 8]
    psir = qr[:, :, 8]
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
    al = np.sqrt(gam*pl/rhol)
    ar = np.sqrt(gam*pr/rhor)
    # Alfvén速度
    cal = np.sqrt((bxl*bxl+byl*byl+bzl*bzl)/rhol)
    car = np.sqrt((bxr*bxr+byr*byr+bzr*bzr)/rhor)
    cayl = np.sqrt(byl*byl/rhol)
    cayr = np.sqrt(byr*byr/rhor)
    # 速進磁気音波速度
    cfl = np.sqrt(0.5*(al*al+cal*cal+np.sqrt((al*al+cal*cal)**2-4.0*(al*cayl)**2)))
    cfr = np.sqrt(0.5*(ar*ar+car*car+np.sqrt((ar*ar+car*car)**2-4.0*(ar*cayr)**2)))
    sfl = np.minimum(vl-cfl, vr-cfr)
    sfr = np.maximum(vl+cfl, vr+cfr)
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
    srrhoml = np.sqrt(rhoml)
    srrhomr = np.sqrt(rhomr)
    sal = sm-np.sqrt(bym*bym/rhoml)
    sar = sm+np.sqrt(bym*bym/rhomr)
    # リーマンファン内側
    wi = (srrhoml*wol+srrhomr*wor+(bzor-bzol)*np.sign(bym))/(srrhoml+srrhomr)
    ui = (srrhoml*uol+srrhomr*uor+(bxor-bxol)*np.sign(bym))/(srrhoml+srrhomr)
    bzi = (srrhoml*bzor+srrhomr*bzol+srrhoml*srrhomr*(wor-wol)*np.sign(bym))/(srrhoml+srrhomr)
    bxi = (srrhoml*bxor+srrhomr*bxol+srrhoml*srrhomr*(uor-uol)*np.sign(bym))/(srrhoml+srrhomr)
    eil = eol-srrhoml*(uol*bxol+wol*bzol-ui*bxi-wi*bzi)*np.sign(bym)
    eir = eor+srrhomr*(uor*bxor+wor*bzor-ui*bxi-wi*bzi)*np.sign(bym)
    fyl = mhdfy(rhol, ul, vl, wl, bxl, bym, bzl, el, psim, ptl, ch)
    fyol = mhdfy(rhoml, uol, vm, wol, bxol, bym, bzol, eol, psim, ptm, ch)
    fyil = mhdfy(rhoml, ui, vm, wi, bxi, bym, bzi, eil, psim, ptm, ch)
    fyir = mhdfy(rhomr, ui, vm, wi, bxi, bym, bzi, eir, psim, ptm, ch)
    fyor = mhdfy(rhomr, uor, vm, wor, bxor, bym, bzor, eor, psim, ptm, ch)
    fyr = mhdfy(rhor, ur, vr, wr, bxr, bym, bzr, er, psim, ptr, ch)
    fy = np.zeros((XN, YN, VN))
    slindex = sfl > 0.0
    solindex = (0.0 >= sfl) & (sal > 0.0)
    silindex = (0.0 >= sal) & (sm > 0.0)
    sirindex = (0.0 >= sm) & (sar > 0.0)
    sorindex = (0.0 >= sar) & (sfr > 0.0)
    srindex = 0.0 >= sfr
    fy[slindex, :] = fyl[slindex, :]
    fy[solindex, :] = fyol[solindex, :]
    fy[silindex, :] = fyil[silindex, :]
    fy[sirindex, :] = fyir[sirindex, :]
    fy[sorindex, :] = fyor[sorindex, :]
    fy[srindex, :] = fyr[srindex, :]
    return fy
# 1次精度風上差分
def upwindx(q):
    qil = np.roll(q, 1, axis = 0)
    qi = q
    ql = qil
    qr = qi
    return ql, qr
def upwindy(q):
    qjl = np.roll(q, 1, axis = 1)
    qj = q
    ql = qjl
    qr = qj
    return ql, qr
# 2次精度MUSCL (minmod) (van Leer, 1979)
def musclx(q):
    qill = np.roll(q, 2, axis = 0)
    qil = np.roll(q, 1, axis = 0)
    qi = q
    qir = np.roll(q, -1, axis = 0)
    ql = qil+0.5*minmod(qi-qil, qil-qill)
    qr = qi-0.5*minmod(qir-qi, qi-qil)
    return ql, qr
def muscly(q):
    qjll = np.roll(q, 2, axis = 1)
    qjl = np.roll(q, 1, axis = 1)
    qj = q
    qjr = np.roll(q, -1, axis = 1)
    ql = qjl+0.5*minmod(qj-qjl, qjl-qjll)
    qr = qj-0.5*minmod(qjr-qj, qj-qjl)
    return ql, qr
# 5次精度MP5 (Suresh and Huynh, 1997)
def mp5x(q):
    qilll = np.roll(q, 3, axis = 0)
    qill = np.roll(q, 2, axis = 0)
    qil = np.roll(q, 1, axis = 0)
    qi = q
    qir = np.roll(q, -1, axis = 0)
    qirr = np.roll(q, -2, axis = 0)
    dll = qilll-2.0*qill+qil
    dl = qill-2.0*qil+qi
    d = qil-2.0*qi+qir
    dr = qi-2.0*qir+qirr
    dmml = minmod(dll, dl)
    dmm = minmod(dl, d)
    dmmr = minmod(d, dr)
    qull = qil+2.0*(qil-qill)
    qulr = qi+2.0*(qi-qir)
    #qull = qil+4.0*(qil-qill)
    #qulr = qi+4.0*(qi-qir)
    qav = 0.5*(qil+qi)
    qmd = qav-0.5*dmm
    qlcl = qil+0.5*(qil-qill)+4.0/3.0*dmml
    qlcr = qi+0.5*(qi-qir)+4.0/3.0*dmmr
    qminl = np.maximum(np.minimum(qil, qi, qmd), np.minimum(qil, qull, qlcl))
    qmaxl = np.minimum(np.maximum(qil, qi, qmd), np.maximum(qil, qull, qlcl))
    qminr = np.maximum(np.minimum(qi, qil, qmd), np.minimum(qi, qulr, qlcr))
    qmaxr = np.minimum(np.maximum(qi, qil, qmd), np.maximum(qi, qulr, qlcr))
    q5l = (2.0*qilll-13.0*qill+47.0*qil+27.0*qi-3.0*qir)/60.0
    q5r = (2.0*qirr-13.0*qir+47.0*qi+27.0*qil-3.0*qill)/60.0
    ql = median(q5l, qminl, qmaxl)
    qr = median(q5r, qminr, qmaxr)
    return ql, qr
def mp5y(q):
    qjlll = np.roll(q, 3, axis = 1)
    qjll = np.roll(q, 2, axis = 1)
    qjl = np.roll(q, 1, axis = 1)
    qj = q
    qjr = np.roll(q, -1, axis = 1)
    qjrr = np.roll(q, -2, axis = 1)
    dll = qjlll-2.0*qjll+qjl
    dl = qjll-2.0*qjl+qj
    d = qjl-2.0*qj+qjr
    dr = qj-2.0*qjr+qjrr
    dmml = minmod(dll, dl)
    dmm = minmod(dl, d)
    dmmr = minmod(d, dr)
    qull = qjl+2.0*(qjl-qjll)
    qulr = qj+2.0*(qj-qjr)
    #qull = qjl+4.0*(qjl-qjll)
    #qulr = qj+4.0*(qj-qjr)
    qav = 0.5*(qjl+qj)
    qmd = qav-0.5*dmm
    qlcl = qjl+0.5*(qjl-qjll)+4.0/3.0*dmml
    qlcr = qj+0.5*(qj-qjr)+4.0/3.0*dmmr
    qminl = np.maximum(np.minimum(qjl, qj, qmd), np.minimum(qjl, qull, qlcl))
    qmaxl = np.minimum(np.maximum(qjl, qj, qmd), np.maximum(qjl, qull, qlcl))
    qminr = np.maximum(np.minimum(qj, qjl, qmd), np.minimum(qj, qulr, qlcr))
    qmaxr = np.minimum(np.maximum(qj, qjl, qmd), np.maximum(qj, qulr, qlcr))
    q5l = (2.0*qjlll-13.0*qjll+47.0*qjl+27.0*qj-3.0*qjr)/60.0
    q5r = (2.0*qjrr-13.0*qjr+47.0*qj+27.0*qjl-3.0*qjll)/60.0
    ql = median(q5l, qminl, qmaxl)
    qr = median(q5r, qminr, qmaxr)    
    return ql, qr
# dq/dtの計算
def dqdt(q, ch):
    # x方向
    #ql, qr = upwindx(q)
    #ql, qr = musclx(q)
    ql, qr = mp5x(q)
    fx = hlldx(ql, qr, ch)
    # y方向
    #ql, qr = upwindy(q)
    #ql, qr = muscly(q)
    ql, qr = mp5y(q)
    fy = hlldy(ql, qr, ch)
    # 保存則の計算
    fxi = fx
    fyj = fy
    fxir = np.roll(fx, -1, axis = 0)
    fyjr = np.roll(fy, -1, axis = 1)
    dqdt = -(fxir-fxi)/dx-(fyjr-fyj)/dy
    return dqdt
# 1次精度Euler法
def euler(q, dt):
    ch = CFL*np.minimum(dx, dy)/dt
    q1 = q+dqdt(q, ch)*dt
    return q1
# 2次精度Runge-Kutta
def ssprk2(q, dt):
    ch = CFL*np.minimum(dx, dy)/dt
    q1 = q+dqdt(q, ch)*dt
    q2 = 0.5*q+0.5*(q1+dqdt(q1, ch)*dt)
    return q2
# 3次精度Runge-Kutta
def ssprk3(q, dt):
    ch = CFL*np.minimum(dx, dy)/dt
    q1 = q+dqdt(q, ch)*dt
    q2 = 0.75*q+0.25*(q1+dqdt(q1, ch)*dt)
    q3 = 1.0/3.0*q+2.0/3.0*(q2+dqdt(q2, ch)*dt)
    return q3

def main():
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
    # -- 0: rho
    # -- 1: mx
    # -- 2: my
    # -- 3: mz
    # -- 4: bx
    # -- 5: by
    # -- 6: bz
    # -- 7: e
    # -- 8: psi
    x = np.arange(0.0, XL, dx)
    y = np.arange(0.0, YL, dy)
    # 初期値
    rho = np.ones((XN, YN))*gam**2
    p = np.ones((XN, YN))*gam
    u = np.tile(-np.sin(y), (XN, 1))
    v = np.tile(np.sin(x), (YN, 1)).T
    w = np.zeros((XN, YN))
    bx = np.tile(-np.sin(y), (XN, 1))
    by = np.tile(np.sin(2.0*x), (YN, 1)).T
    bz = np.zeros((XN, YN))
    e = p/(gam-1.0)+0.5*rho*(u*u+v*v+w*w)+0.5*(bx*bx+by*by+bz*bz)
    psi = np.zeros((XN, YN))
    q = np.zeros((XN, YN, VN)); q_n = np.zeros((XN, YN, VN))
    q[:, :, 0] = rho
    q[:, :, 1] = rho*u
    q[:, :, 2] = rho*v
    q[:, :, 3] = rho*w
    q[:, :, 4] = bx
    q[:, :, 5] = by
    q[:, :, 6] = bz
    q[:, :, 7] = e
    q[:, :, 8] = psi
    # 計算開始
    n = 0
    t = 0.0
    step = 0
    start = time.time()
    while (t <= TL):
        # 密度
        rho = q[:, :, 0]
        # 運動量
        mx = q[:, :, 1]
        my = q[:, :, 2]
        mz = q[:, :, 3]
        # 磁場
        bx = q[:, :, 4]
        by = q[:, :, 5]
        bz = q[:, :, 6]
        # エネルギー
        e = q[:, :, 7]
        # 磁場発散抑制のための人工ポテンシャル (Dedner et al., 2002)
        psi = q[:, :, 8]
        # 速度
        u = mx/rho
        v = my/rho
        w = mz/rho
        # (熱的) 圧力
        p = (gam-1.0)*(e-0.5*rho*(u*u+v*v+w*w)-0.5*(bx*bx+by*by+bz*bz))
        # 音速
        a = np.sqrt(gam*p/rho)
        # Alfvén速度
        ca = np.sqrt((bx*bx+by*by+bz*bz)/rho)
        cax = np.sqrt(bx*bx/rho)
        cay = np.sqrt(by*by/rho)
        # 速進磁気音波速度
        cfx = np.sqrt(0.5*(a*a+ca*ca+np.sqrt((a*a+ca*ca)**2-4.0*(a*cax)**2)))
        cfy = np.sqrt(0.5*(a*a+ca*ca+np.sqrt((a*a+ca*ca)**2-4.0*(a*cay)**2)))
        # CFL数をもとにdtを設定
        dt = CFL*np.minimum(dx/np.max(np.abs(u)+cfx), dy/np.max(np.abs(v)+cfy))
        # 磁場発散抑制のためのパラメータ (Dedner et al., 2002)
        ch = CFL*np.minimum(dx, dy)/dt
        cd = np.exp(-dt*ch/cr)
        # 時間発展
        #q_n = euler(q, dt)
        #q_n = ssprk2(q, dt)
        q_n = ssprk3(q, dt)
        q_n[:, :, 8] = cd*q_n[:, :, 8]
        q = q_n
        # ファイル出力
        if (np.floor(t*PN/TL) != np.floor((t-dt)*PN/TL)):
            print("n: {0:3d}, t: {1:.2f}".format(n, t))
            #varname = ["r", "p", "u", "v", "w", "bx", "by", "bz", "ps"]
            #vardata = [rho, p, u, v, w, bx, by, bz, psi]
            #for i in range(VN):
            #    np.savetxt("../data/"+varname[i]+"_{0:0>3}.csv".format(n), vardata[i].T, delimiter = ",")
            n += 1
        t += dt
        step += 1
    elapsed_time = time.time()-start
    print("Elapsed time: {0} s".format(elapsed_time))
    # 計算終了
    # 計算失敗の検知
    if (not np.isfinite(t)):
        print("Computaion failed: step = {0:3d}".format(step))
    # y = πでの値
    print("   x    rho      p      u      v      w     bx     by     bz    psi")
    for i in range(0, XN, 5):
        j = YN//2
        print("{0:4.2f}  {1:5.2f}  {2:5.2f}  {3:5.2f}  {4:5.2f}  {5:5.2f}  {6:5.2f}  {7:5.2f}  {8:5.2f}  {9:5.2f}"
              .format(dx*i, rho[i, j], p[i, j], u[i, j], v[i, j], w[i, j], bx[i, j], by[i, j], bz[i, j], psi[i, j]))
    # 結果をプロット
    fig = plt.figure()
    plt.axes().set_aspect("equal")
    plt.title("p (t = {0:.2f})".format(t))
    plt.xlabel("X")
    plt.ylabel("Y")
    im = plt.pcolor(x, y, p.T, cmap = cm.jet, shading = "auto")
    im.set_clim(0.0, 6.5)
    fig.colorbar(im)
    plt.show()

if __name__ == "__main__":
    main()