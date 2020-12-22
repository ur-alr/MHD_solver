#include "mhd.hh"

double minmod(double a, double b) {
    return 0.5*(sign(a)+sign(b))*min(std::abs(a), std::abs(b));
}
double median(double a, double b, double c) {
    return a+minmod(b-a, c-a);
}
auto hlldx(const ndarray_t<double, VN> &ql, const ndarray_t<double, VN> &qr, double ch) {
    // 密度
    double rhol = ql[0];
    double rhor = qr[0];
    // 運動量
    double mxl = ql[1];
    double myl = ql[2];
    double mzl = ql[3];
    double mxr = qr[1];
    double myr = qr[2];
    double mzr = qr[3];
    // 磁場
    double bxl = ql[4];
    double byl = ql[5];
    double bzl = ql[6];
    double bxr = qr[4];
    double byr = qr[5];
    double bzr = qr[6];
    // エネルギー
    double el = ql[7];
    double er = qr[7];
    // 磁場発散抑制のための人工ポテンシャル (Dedner et al., 2002)
    double psil = ql[8];
    double psir = qr[8];
    // 速度
    double ul = mxl/rhol;
    double vl = myl/rhol;
    double wl = mzl/rhol;
    double ur = mxr/rhor;
    double vr = myr/rhor;
    double wr = mzr/rhor;
    // (熱的) 圧力
    double pl = (gam-1.0)*(el-0.5*rhol*(ul*ul+vl*vl+wl*wl)-0.5*(bxl*bxl+byl*byl+bzl*bzl));
    double pr = (gam-1.0)*(er-0.5*rhor*(ur*ur+vr*vr+wr*wr)-0.5*(bxr*bxr+byr*byr+bzr*bzr));
    double ptl = pl+0.5*(bxl*bxl+byl*byl+bzl*bzl);
    // 総圧力 (熱的圧力 + 磁気圧)
    double ptr = pr+0.5*(bxr*bxr+byr*byr+bzr*bzr);
    // 音速
    double al = std::sqrt(gam*pl/rhol);
    double ar = std::sqrt(gam*pr/rhor);
    // Alfvén速度
    double cal = std::sqrt((bxl*bxl+byl*byl+bzl*bzl)/rhol);
    double car = std::sqrt((bxr*bxr+byr*byr+bzr*bzr)/rhor);
    double caxl = std::sqrt(bxl*bxl/rhol);
    double caxr = std::sqrt(bxr*bxr/rhor);
    // 速進磁気音波速度
    double cfl = std::sqrt(0.5*(al*al+cal*cal+std::sqrt((al*al+cal*cal)*(al*al+cal*cal)-4.0*al*al*caxl*caxl)));
    double cfr = std::sqrt(0.5*(ar*ar+car*car+std::sqrt((ar*ar+car*car)*(ar*ar+car*car)-4.0*ar*ar*caxr*caxr)));
    double sfl = min(ul-cfl, ur-cfr);
    double sfr = max(ul+cfl, ur+cfr);
    double sm = ((sfr-ur)*rhor*ur-(sfl-ul)*rhol*ul-ptr+ptl)/((sfr-ur)*rhor-(sfl-ul)*rhol);
    double um = sm;
    double bxm = bxl+0.5*(bxr-bxl)-0.5/ch*(psir-psil);
    double psim = psil+0.5*(psir-psil)-0.5*ch*(bxr-bxl);
    // リーマンファン外側
    double ptm = ((sfr-ur)*rhor*ptl-(sfl-ul)*rhol*ptr+rhol*rhor*(sfr-ur)*(sfl-ul)*(ur-ul))/((sfr-ur)*rhor-(sfl-ul)*rhol);
    double rhoml = rhol*(sfl-ul)/(sfl-sm);
    double rhomr = rhor*(sfr-ur)/(sfr-sm);
    double vol = vl-bxm*byl*(sm-ul)/(rhol*(sfl-ul)*(sfl-sm)-bxm*bxm);
    double vor = vr-bxm*byr*(sm-ur)/(rhor*(sfr-ur)*(sfr-sm)-bxm*bxm);
    double wol = wl-bxm*bzl*(sm-ul)/(rhol*(sfl-ul)*(sfl-sm)-bxm*bxm);
    double wor = wr-bxm*bzr*(sm-ur)/(rhor*(sfr-ur)*(sfr-sm)-bxm*bxm);
    double byol = byl*(rhol*(sfl-ul)*(sfl-ul)-bxm*bxm)/(rhol*(sfl-ul)*(sfl-sm)-bxm*bxm);
    double byor = byr*(rhor*(sfr-ur)*(sfr-ur)-bxm*bxm)/(rhor*(sfr-ur)*(sfr-sm)-bxm*bxm);
    double bzol = bzl*(rhol*(sfl-ul)*(sfl-ul)-bxm*bxm)/(rhol*(sfl-ul)*(sfl-sm)-bxm*bxm);
    double bzor = bzr*(rhor*(sfr-ur)*(sfr-ur)-bxm*bxm)/(rhor*(sfr-ur)*(sfr-sm)-bxm*bxm);
    double eol = ((sfl-ul)*el-ptl*ul+ptm*sm+bxm*(ul*bxl+vl*byl+wl*bzl-um*bxm-vol*byol-wol*bzol))/(sfl-sm);
    double eor = ((sfr-ur)*er-ptr*ur+ptm*sm+bxm*(ur*bxr+vr*byr+wr*bzr-um*bxm-vor*byor-wor*bzor))/(sfr-sm);
    double srrhoml = std::sqrt(rhoml);
    double srrhomr = std::sqrt(rhomr);
    double sal = sm-std::sqrt(bxm*bxm/rhoml);
    double sar = sm+std::sqrt(bxm*bxm/rhomr);
    // リーマンファン内側
    double vi = (srrhoml*vol+srrhomr*vor+(byor-byol)*sign(bxm))/(srrhoml+srrhomr);
    double wi = (srrhoml*wol+srrhomr*wor+(bzor-bzol)*sign(bxm))/(srrhoml+srrhomr);
    double byi = (srrhoml*byor+srrhomr*byol+srrhoml*srrhomr*(vor-vol)*sign(bxm))/(srrhoml+srrhomr);
    double bzi = (srrhoml*bzor+srrhomr*bzol+srrhoml*srrhomr*(wor-wol)*sign(bxm))/(srrhoml+srrhomr);
    double eil = eol-srrhoml*(vol*byol+wol*bzol-vi*byi-wi*bzi)*sign(bxm);
    double eir = eor+srrhomr*(vor*byor+wor*bzor-vi*byi-wi*bzi)*sign(bxm);
    if (sfl > 0) {
        return ndarray_t<double, VN>{
            rhol*ul,
            rhol*ul*ul+ptl-bxm*bxm,
            rhol*vl*ul-bxm*byl,
            rhol*wl*ul-bxm*bzl,
            psim,
            byl*ul-bxm*vl,
            bzl*ul-bxm*wl,
            (el+ptl)*ul-bxm*(ul*bxm+vl*byl+wl*bzl),
            ch*ch*bxm
        };
    } else if (sal > 0) {
        return ndarray_t<double, VN>{
            rhoml*um,
            rhoml*um*um+ptm-bxm*bxm,
            rhoml*vol*um-bxm*byol,
            rhoml*wol*um-bxm*bzol,
            psim,
            byol*um-bxm*vol,
            bzol*um-bxm*wol,
            (eol+ptm)*um-bxm*(um*bxm+vol*byol+wol*bzol),
            ch*ch*bxm
        };
    } else if (sm > 0) {
        return ndarray_t<double, VN>{
            rhoml*um,
            rhoml*um*um+ptm-bxm*bxm,
            rhoml*vi*um-bxm*byi,
            rhoml*wi*um-bxm*bzi,
            psim,
            byi*um-bxm*vi,
            bzi*um-bxm*wi,
            (eil+ptm)*um-bxm*(um*bxm+vi*byi+wi*bzi),
            ch*ch*bxm
        };
    } else if (sar > 0) {
        return ndarray_t<double, VN>{
            rhomr*um,
            rhomr*um*um+ptm-bxm*bxm,
            rhomr*vi*um-bxm*byi,
            rhomr*wi*um-bxm*bzi,
            psim,
            byi*um-bxm*vi,
            bzi*um-bxm*wi,
            (eir+ptm)*um-bxm*(um*bxm+vi*byi+wi*bzi),
            ch*ch*bxm
        };
    } else if (sfr > 0) {
        return ndarray_t<double, VN>{
            rhomr*um,
            rhomr*um*um+ptm-bxm*bxm,
            rhomr*vor*um-bxm*byor,
            rhomr*wor*um-bxm*bzor,
            psim,
            byor*um-bxm*vor,
            bzor*um-bxm*wor,
            (eor+ptm)*um-bxm*(um*bxm+vor*byor+wor*bzor),
            ch*ch*bxm
        };
    } else {
        return ndarray_t<double, VN>{
            rhor*ur,
            rhor*ur*ur+ptr-bxm*bxm,
            rhor*vr*ur-bxm*byr,
            rhor*wr*ur-bxm*bzr,
            psim,
            byr*ur-bxm*vr,
            bzr*ur-bxm*wr,
            (er+ptr)*ur-bxm*(ur*bxm+vr*byr+wr*bzr),
            ch*ch*bxm
        };
    }
}
auto hlldy(const ndarray_t<double, VN> &ql, const ndarray_t<double, VN> &qr, double ch) {
    // 密度
    double rhol = ql[0];
    double rhor = qr[0];
    // 運動量
    double mxl = ql[1];
    double myl = ql[2];
    double mzl = ql[3];
    double mxr = qr[1];
    double myr = qr[2];
    double mzr = qr[3];
    // 磁場
    double bxl = ql[4];
    double byl = ql[5];
    double bzl = ql[6];
    double bxr = qr[4];
    double byr = qr[5];
    double bzr = qr[6];
    // エネルギー
    double el = ql[7];
    double er = qr[7];
    // 磁場発散抑制のための人工ポテンシャル (Dedner et al., 2002)
    double psil = ql[8];
    double psir = qr[8];
    // 速度
    double ul = mxl/rhol;
    double vl = myl/rhol;
    double wl = mzl/rhol;
    double ur = mxr/rhor;
    double vr = myr/rhor;
    double wr = mzr/rhor;
    // (熱的) 圧力
    double pl = (gam-1.0)*(el-0.5*rhol*(ul*ul+vl*vl+wl*wl)-0.5*(bxl*bxl+byl*byl+bzl*bzl));
    double pr = (gam-1.0)*(er-0.5*rhor*(ur*ur+vr*vr+wr*wr)-0.5*(bxr*bxr+byr*byr+bzr*bzr));
    // 総圧力 (熱的圧力 + 磁気圧)
    double ptl = pl+0.5*(bxl*bxl+byl*byl+bzl*bzl);
    double ptr = pr+0.5*(bxr*bxr+byr*byr+bzr*bzr);
    // 音速
    double al = std::sqrt(gam*pl/rhol);
    double ar = std::sqrt(gam*pr/rhor);
    // Alfvén速度
    double cal = std::sqrt((bxl*bxl+byl*byl+bzl*bzl)/rhol);
    double car = std::sqrt((bxr*bxr+byr*byr+bzr*bzr)/rhor);
    double cayl = std::sqrt(byl*byl/rhol);
    double cayr = std::sqrt(byr*byr/rhor);
    // 速進磁気音波速度
    double cfl = std::sqrt(0.5*(al*al+cal*cal+std::sqrt((al*al+cal*cal)*(al*al+cal*cal)-4.0*al*al*cayl*cayl)));
    double cfr = std::sqrt(0.5*(ar*ar+car*car+std::sqrt((ar*ar+car*car)*(ar*ar+car*car)-4.0*ar*ar*cayr*cayr)));
    double sfl = min(vl-cfl, vr-cfr);
    double sfr = max(vl+cfl, vr+cfr);
    double sm = ((sfr-vr)*rhor*vr-(sfl-vl)*rhol*vl-ptr+ptl)/((sfr-vr)*rhor-(sfl-vl)*rhol);
    double vm = sm;
    double bym = byl+0.5*(byr-byl)-0.5/ch*(psir-psil);
    double psim = psil+0.5*(psir-psil)-0.5*ch*(byr-byl);
    // リーマンファン外側
    double ptm = ((sfr-vr)*rhor*ptl-(sfl-vl)*rhol*ptr+rhol*rhor*(sfr-vr)*(sfl-vl)*(vr-vl))/((sfr-vr)*rhor-(sfl-vl)*rhol);
    double rhoml = rhol*(sfl-vl)/(sfl-sm);
    double rhomr = rhor*(sfr-vr)/(sfr-sm);
    double wol = wl-bym*bzl*(sm-vl)/(rhol*(sfl-vl)*(sfl-sm)-bym*bym);
    double wor = wr-bym*bzr*(sm-vr)/(rhor*(sfr-vr)*(sfr-sm)-bym*bym);
    double uol = ul-bym*bxl*(sm-vl)/(rhol*(sfl-vl)*(sfl-sm)-bym*bym);
    double uor = ur-bym*bxr*(sm-vr)/(rhor*(sfr-vr)*(sfr-sm)-bym*bym);
    double bzol = bzl*(rhol*(sfl-vl)*(sfl-vl)-bym*bym)/(rhol*(sfl-vl)*(sfl-sm)-bym*bym);
    double bzor = bzr*(rhor*(sfr-vr)*(sfr-vr)-bym*bym)/(rhor*(sfr-vr)*(sfr-sm)-bym*bym);
    double bxol = bxl*(rhol*(sfl-vl)*(sfl-vl)-bym*bym)/(rhol*(sfl-vl)*(sfl-sm)-bym*bym);
    double bxor = bxr*(rhor*(sfr-vr)*(sfr-vr)-bym*bym)/(rhor*(sfr-vr)*(sfr-sm)-bym*bym);
    double eol = ((sfl-vl)*el-ptl*vl+ptm*sm+bym*(ul*bxl+vl*byl+wl*bzl-uol*bxol-vm*bym-wol*bzol))/(sfl-sm);
    double eor = ((sfr-vr)*er-ptr*vr+ptm*sm+bym*(ur*bxr+vr*byr+wr*bzr-uor*bxor-vm*bym-wor*bzor))/(sfr-sm);
    double srrhoml = std::sqrt(rhoml);
    double srrhomr = std::sqrt(rhomr);
    double sal = sm-std::sqrt(bym*bym/rhoml);
    double sar = sm+std::sqrt(bym*bym/rhomr);
    // リーマンファン内側
    double wi = (srrhoml*wol+srrhomr*wor+(bzor-bzol)*sign(bym))/(srrhoml+srrhomr);
    double ui = (srrhoml*uol+srrhomr*uor+(bxor-bxol)*sign(bym))/(srrhoml+srrhomr);
    double bzi = (srrhoml*bzor+srrhomr*bzol+srrhoml*srrhomr*(wor-wol)*sign(bym))/(srrhoml+srrhomr);
    double bxi = (srrhoml*bxor+srrhomr*bxol+srrhoml*srrhomr*(uor-uol)*sign(bym))/(srrhoml+srrhomr);
    double eil = eol-srrhoml*(uol*bxol+wol*bzol-ui*bxi-wi*bzi)*sign(bym);
    double eir = eor+srrhomr*(uor*bxor+wor*bzor-ui*bxi-wi*bzi)*sign(bym);
    if (sfl > 0) {
        return ndarray_t<double, VN>{
            rhol*vl,
            rhol*ul*vl-bym*bxl,
            rhol*vl*vl+ptl-bym*bym,
            rhol*wl*vl-bym*bzl,
            bxl*vl-bym*ul,
            psim,
            bzl*vl-bym*wl,
            (el+ptl)*vl-bym*(ul*bxl+vl*bym+wl*bzl),
            ch*ch*bym
        };
    } else if (sal > 0) { 
        return ndarray_t<double, VN>{
            rhoml*vm,
            rhoml*uol*vm-bym*bxol,
            rhoml*vm*vm+ptm-bym*bym,
            rhoml*wol*vm-bym*bzol,
            bxol*vm-bym*uol,
            psim,
            bzol*vm-bym*wol,
            (eol+ptm)*vm-bym*(uol*bxol+vm*bym+wol*bzol),
            ch*ch*bym
        };
    } else if (sm > 0) {
        return ndarray_t<double, VN>{
            rhoml*vm,
            rhoml*ui*vm-bym*bxi,
            rhoml*vm*vm+ptm-bym*bym,
            rhoml*wi*vm-bym*bzi,
            bxi*vm-bym*ui,
            psim,
            bzi*vm-bym*wi,
            (eil+ptm)*vm-bym*(ui*bxi+vm*bym+wi*bzi),
            ch*ch*bym
        };
    } else if (sar > 0) {
        return ndarray_t<double, VN>{
            rhomr*vm,
            rhomr*ui*vm-bym*bxi,
            rhomr*vm*vm+ptm-bym*bym,
            rhomr*wi*vm-bym*bzi,
            bxi*vm-bym*ui,
            psim,
            bzi*vm-bym*wi,
            (eir+ptm)*vm-bym*(ui*bxi+vm*bym+wi*bzi),
            ch*ch*bym
        };
    } else if (sfr > 0) {
        return ndarray_t<double, VN>{
            rhomr*vm,
            rhomr*uor*vm-bym*bxor,
            rhomr*vm*vm+ptm-bym*bym,
            rhomr*wor*vm-bym*bzor,
            bxor*vm-bym*uor,
            psim,
            bzor*vm-bym*wor,
            (eor+ptm)*vm-bym*(uor*bxor+vm*bym+wor*bzor),
            ch*ch*bym
        };
    } else {
        return ndarray_t<double, VN>{
            rhor*vr,
            rhor*ur*vr-bym*bxr,
            rhor*vr*vr+ptr-bym*bym,
            rhor*wr*vr-bym*bzr,
            bxr*vr-bym*ur,
            psim,
            bzr*vr-bym*wr,
            (er+ptr)*vr-bym*(ur*bxr+vr*bym+wr*bzr),
            ch*ch*bym
        };
    }
}
auto upwindx(const ndarray_t<double, XN, YN, VN> &q, int i, int j) {
    ndarray_t<double, VN> ql, qr;
    int il = (i-1+XN)%XN;
    for (int k = 0; k < VN; k++) {
        ql[k] = q[il][j][k];
        qr[k] = q[i][j][k];
    }
    return std::make_tuple(ql, qr);
}
auto upwindy(const ndarray_t<double, XN, YN, VN> &q, int i, int j) {
    ndarray_t<double, VN> ql, qr;
    int jl = (j-1+YN)%YN;
    for (int k = 0; k < VN; k++) {
        ql[k] = q[i][jl][k];
        qr[k] = q[i][j][k];
    }
    return std::make_tuple(ql, qr);
}
auto musclx(const ndarray_t<double, XN, YN, VN> &q, int i, int j) {
    ndarray_t<double, VN> ql, qr;
    int ill = (i-2+XN)%XN;
    int il = (i-1+XN)%XN;
    int ir = (i+1)%XN;
    for (int k = 0; k < VN; k++) {
        ql[k] = q[il][j][k]+0.5*minmod(q[i][j][k]-q[il][j][k], q[il][j][k]-q[ill][j][k]);
        qr[k] = q[i][j][k]-0.5*minmod(q[ir][j][k]-q[i][j][k], q[i][j][k]-q[il][j][k]);
    }
    return std::make_tuple(ql, qr);
}
auto muscly(const ndarray_t<double, XN, YN, VN> &q, int i, int j) {
    ndarray_t<double, VN> ql, qr;
    int jll = (j-2+YN)%YN;
    int jl = (j-1+YN)%YN;
    int jr = (j+1)%YN;
    for (int k = 0; k < VN; k++) {
        ql[k] = q[i][jl][k]+0.5*minmod(q[i][j][k]-q[i][jl][k], q[i][jl][k]-q[i][jll][k]);
        qr[k] = q[i][j][k]-0.5*minmod(q[i][jr][k]-q[i][j][k], q[i][j][k]-q[i][jl][k]);
    }
    return std::make_tuple(ql, qr);
}
auto mp5x(const ndarray_t<double, XN, YN, VN> &q, int i, int j) {
    ndarray_t<double, VN> ql, qr;
    int illl = (i-3+XN)%XN;
    int ill = (i-2+XN)%XN;
    int il = (i-1+XN)%XN;
    int ir = (i+1)%XN;
    int irr = (i+2)%XN;
    for (int k = 0; k < VN; k++) {
        double dll = q[illl][j][k]-2.0*q[ill][j][k]+q[il][j][k];
        double dl = q[ill][j][k]-2.0*q[il][j][k]+q[i][j][k];
        double d = q[il][j][k]-2.0*q[i][j][k]+q[ir][j][k];
        double dr = q[i][j][k]-2.0*q[ir][j][k]+q[irr][j][k];
        double dmml = minmod(dll, dl);
        double dmm = minmod(dl, d);
        double dmmr = minmod(d, dr);
        double qull = q[il][j][k]+2.0*(q[il][j][k]-q[ill][j][k]);
        double qulr = q[i][j][k]+2.0*(q[i][j][k]-q[ir][j][k]);
        //double qull = q[il][j][k]+4.0*(q[il][j][k]-q[ill][j][k]);
        //double qulr = q[i][j][k]+4.0*(q[i][j][k]-q[ir][j][k]);
        double qav = 0.5*(q[il][j][k]+q[i][j][k]);
        double qmd = qav-0.5*dmm;
        double qlcl = q[il][j][k]+0.5*(q[il][j][k]-q[ill][j][k])+4.0/3.0*dmml;
        double qlcr = q[i][j][k]+0.5*(q[i][j][k]-q[ir][j][k])+4.0/3.0*dmmr;
        double qminl = max(min(q[il][j][k], q[i][j][k], qmd), min(q[il][j][k], qull, qlcl));
        double qmaxl = min(max(q[il][j][k], q[i][j][k], qmd), max(q[il][j][k], qull, qlcl));
        double qminr = max(min(q[i][j][k], q[il][j][k], qmd), min(q[i][j][k], qulr, qlcr));
        double qmaxr = min(max(q[i][j][k], q[il][j][k], qmd), max(q[i][j][k], qulr, qlcr));
        double q5l = (2.0*q[illl][j][k]-13.0*q[ill][j][k]+47.0*q[il][j][k]+27.0*q[i][j][k]-3.0*q[ir][j][k])/60.0;
        double q5r = (2.0*q[irr][j][k]-13.0*q[ir][j][k]+47.0*q[i][j][k]+27.0*q[il][j][k]-3.0*q[ill][j][k])/60.0;
        ql[k] = median(q5l, qminl, qmaxl);
        qr[k] = median(q5r, qminr, qmaxr);
    }
    return std::make_tuple(ql, qr);
}
auto mp5y(const ndarray_t<double, XN, YN, VN> &q, int i, int j) {
    ndarray_t<double, VN> ql, qr;
    int jlll = (j-3+YN)%YN;
    int jll = (j-2+YN)%YN;
    int jl = (j-1+YN)%YN;
    int jr = (j+1)%YN;
    int jrr = (j+2)%YN;
    for (int k = 0; k < VN; k++) {
        double dll = q[i][jlll][k]-2.0*q[i][jll][k]+q[i][jl][k];
        double dl = q[i][jll][k]-2.0*q[i][jl][k]+q[i][j][k];
        double d = q[i][jl][k]-2.0*q[i][j][k]+q[i][jr][k];
        double dr = q[i][j][k]-2.0*q[i][jr][k]+q[i][jrr][k];
        double dmml = minmod(dll, dl);
        double dmm = minmod(dl, d);
        double dmmr = minmod(d, dr);
        double qull = q[i][jl][k]+2.0*(q[i][jl][k]-q[i][jll][k]);
        double qulr = q[i][j][k]+2.0*(q[i][j][k]-q[i][jr][k]);
        //double qull = q[i][jl][k]+4.0*(q[i][jl][k]-q[i][jll][k]);
        //double qulr = q[i][j][k]+4.0*(q[i][j][k]-q[i][jr][k]);
        double qav = 0.5*(q[i][jl][k]+q[i][j][k]);
        double qmd = qav-0.5*dmm;
        double qlcl = q[i][jl][k]+0.5*(q[i][jl][k]-q[i][jll][k])+4.0/3.0*dmml;
        double qlcr = q[i][j][k]+0.5*(q[i][j][k]-q[i][jr][k])+4.0/3.0*dmmr;
        double qminl = max(min(q[i][jl][k], q[i][j][k], qmd), min(q[i][jl][k], qull, qlcl));
        double qmaxl = min(max(q[i][jl][k], q[i][j][k], qmd), max(q[i][jl][k], qull, qlcl));
        double qminr = max(min(q[i][j][k], q[i][jl][k], qmd), min(q[i][j][k], qulr, qlcr));
        double qmaxr = min(max(q[i][j][k], q[i][jl][k], qmd), max(q[i][j][k], qulr, qlcr));
        double q5l = (2.0*q[i][jlll][k]-13.0*q[i][jll][k]+47.0*q[i][jl][k]+27.0*q[i][j][k]-3.0*q[i][jr][k])/60.0;
        double q5r = (2.0*q[i][jrr][k]-13.0*q[i][jr][k]+47.0*q[i][j][k]+27.0*q[i][jl][k]-3.0*q[i][jll][k])/60.0;
        ql[k] = median(q5l, qminl, qmaxl);
        qr[k] = median(q5r, qminr, qmaxr);
    }
    return std::make_tuple(ql, qr);
}
auto dqdt(const ndarray_t<double, XN, YN, VN> &q, double ch) {
    ndarray_t<double, XN, YN, VN> fx, fy;
    // x方向
    #pragma omp parallel for
    for (int i = 0; i < XN; i++) {
        for (int j = 0; j < YN; j++) {
            //auto [ql, qr] = upwindx(q, i, j);
            //auto [ql, qr] = musclx(q, i, j);
            auto [ql, qr] = mp5x(q, i, j);
            auto fxi = hlldx(ql, qr, ch);
            for (int k = 0; k < VN; k++) {
                fx[i][j][k] = fxi[k];
            }
        }
    }
    // y方向
    #pragma omp parallel for
    for (int j = 0; j < YN; j++) {
        for (int i = 0; i < XN; i++) {
            //auto [ql, qr] = upwindy(q, i, j);
            //auto [ql, qr] = muscly(q, i, j);
            auto [ql, qr] = mp5y(q, i, j);
            auto fyj = hlldy(ql, qr, ch);
            for (int k = 0; k < VN; k++) {
                fy[i][j][k] = fyj[k];
            }
        }
    }
    ndarray_t<double, XN, YN, VN> dqdt;
    // 保存則の計算
    #pragma omp parallel for
    for (int i = 0; i < XN; i++) {
        int ir = (i+1)%XN;
        for (int j = 0; j < YN; j++) {
            int jr = (j+1)%YN;
            for (int k = 0; k < VN; k++) {
                dqdt[i][j][k] = -(fx[ir][j][k]-fx[i][j][k])/dx-(fy[i][jr][k]-fy[i][j][k])/dy;
            }
        }
    }
    return dqdt;
}
ndarray_t<double, XN, YN, VN> euler(const ndarray_t<double, XN, YN, VN> &q, double dt) {
    ndarray_t<double, XN, YN, VN> q1;
    double ch = CFL*min(dx, dy)/dt;
    auto dqdtij = dqdt(q, ch);
    #pragma omp parallel for
    for (int i = 0; i < XN; i++) {
        for (int j = 0; j < YN; j++) {
            for (int k = 0; k < VN; k++) {
                q1[i][j][k] = q[i][j][k]+dqdtij[i][j][k]*dt;
            }
        }
    }
    return q1;
}
ndarray_t<double, XN, YN, VN> ssprk2(const ndarray_t<double, XN, YN, VN> &q, double dt) {
    ndarray_t<double, XN, YN, VN> q1, q2;
    double ch = CFL*min(dx, dy)/dt;
    auto dqdtij = dqdt(q, ch);
    #pragma omp parallel for
    for (int i = 0; i < XN; i++) {
        for (int j = 0; j < YN; j++) {
            for (int k = 0; k < VN; k++) {
                q1[i][j][k] = q[i][j][k]+dqdtij[i][j][k]*dt;
            }
        }
    }
    dqdtij = dqdt(q1, ch);
    #pragma omp parallel for
    for (int i = 0; i < XN; i++) {
        for (int j = 0; j < YN; j++) {
            for (int k = 0; k < VN; k++) {
                q2[i][j][k] = 0.5*q[i][j][k]+0.5*(q1[i][j][k]+dqdtij[i][j][k]*dt);
            }
        }
    }
    return q2;
}
ndarray_t<double, XN, YN, VN> ssprk3(const ndarray_t<double, XN, YN, VN> &q, double dt) {
    ndarray_t<double, XN, YN, VN> q1, q2, q3;
    double ch = CFL*min(dx, dy)/dt;
    auto dqdtij = dqdt(q, ch);
    #pragma omp parallel for
    for (int i = 0; i < XN; i++) {
        for (int j = 0; j < YN; j++) {
            for (int k = 0; k < VN; k++) {
                q1[i][j][k] = q[i][j][k]+dqdtij[i][j][k]*dt;
            }
        }
    }
    dqdtij = dqdt(q1, ch);
    #pragma omp parallel for
    for (int i = 0; i < XN; i++) {
        for (int j = 0; j < YN; j++) {
            for (int k = 0; k < VN; k++) {
                q2[i][j][k] = 0.75*q[i][j][k]+0.25*(q1[i][j][k]+dqdtij[i][j][k]*dt);
            }
        }
    }
    dqdtij = dqdt(q2, ch);
    #pragma omp parallel for
    for (int i = 0; i < XN; i++) {
        for (int j = 0; j < YN; j++) {
            for (int k = 0; k < VN; k++) {
                q3[i][j][k] = 1.0/3.0*q[i][j][k]+2.0/3.0*(q2[i][j][k]+dqdtij[i][j][k]*dt);
            }
        }
    }
    return q3;
}
