#include "mhd.hh"
#include <chrono>

int main() {
    // rho: 密度
    // p:   圧力
    // u:   x方向速度
    // v:   y方向速度
    // w:   z方向速度
    // mx:  x方向運動量
    // my:  y方向運動量
    // mz:  z方向運動量
    // bx:  x方向磁場
    // by:  y方向磁場
    // bz:  z方向磁場
    // e:   エネルギー
    // psi: 磁場発散抑制のための人工ポテンシャル (Dedner et al., 2002)
    // q: 保存変数ベクトル
    // -- 0: rho
    // -- 1: mx
    // -- 2: my
    // -- 3: mz
    // -- 4: bx
    // -- 5: by
    // -- 6: bz
    // -- 7: e
    // -- 8: psi
    ndarray_t<double, XN, YN> rho, p, u, v, w, mx, my, mz, e, bx, by, bz, psi;
    ndarray_t<double, XN, YN, VN> q, q_n;
    // 初期値
    #pragma omp parallel for
    for (int i = 0; i < XN; i++) {
        double x = i*dx;
        for (int j = 0; j < YN; j++) {
            double y = j*dy;
            rho[i][j] = gam*gam;
            p[i][j] = gam;
            u[i][j] = -std::sin(y);
            v[i][j] = std::sin(x);
            w[i][j] = 0.0;
            bx[i][j] = -std::sin(y);
            by[i][j] = std::sin(2.0*x);
            bz[i][j] = 0.0;
            e[i][j] = p[i][j]/(gam-1.0)
                     +0.5*rho[i][j]*(u[i][j]*u[i][j]+v[i][j]*v[i][j]+w[i][j]*w[i][j])
                     +0.5*(bx[i][j]*bx[i][j]+by[i][j]*by[i][j]+bz[i][j]*bz[i][j]);
            psi[i][j] = 0.0;
            q[i][j][0] = rho[i][j];
            q[i][j][1] = rho[i][j]*u[i][j];
            q[i][j][2] = rho[i][j]*v[i][j];
            q[i][j][3] = rho[i][j]*w[i][j];
            q[i][j][4] = bx[i][j];
            q[i][j][5] = by[i][j];
            q[i][j][6] = bz[i][j];
            q[i][j][7] = e[i][j];
            q[i][j][8] = psi[i][j];
        }
    }
    // 計算開始
    size_t step = 0, n = 0;
    double t = 0.0;
    auto start = std::chrono::system_clock::now();
    while (t <= TL) {
        double lxmax = 0.0, lymax = 0.0;
        #pragma omp parallel for
        for (int i = 0; i < XN; i++) {
            for (int j = 0; j < YN; j++) {
                // 密度
                rho[i][j] = q[i][j][0];
                // 運動量
                mx[i][j] = q[i][j][1];
                my[i][j] = q[i][j][2];
                mz[i][j] = q[i][j][3];
                // 磁場
                bx[i][j] = q[i][j][4];
                by[i][j] = q[i][j][5];
                bz[i][j] = q[i][j][6];
                // エネルギー
                e[i][j] = q[i][j][7];
                // 磁場発散抑制のための人工ポテンシャル (Dedner et al., 2002)
                psi[i][j] = q[i][j][8];
                // 速度
                u[i][j] = mx[i][j]/rho[i][j];
                v[i][j] = my[i][j]/rho[i][j];
                w[i][j] = mz[i][j]/rho[i][j];
                // (熱的) 圧力
                p[i][j] = (gam-1.0)*(e[i][j]
                         -0.5*rho[i][j]*(u[i][j]*u[i][j]+v[i][j]*v[i][j]+w[i][j]*w[i][j])
                         -0.5*(bx[i][j]*bx[i][j]+by[i][j]*by[i][j]+bz[i][j]*bz[i][j]));
                // 音速
                double a = std::sqrt(gam*p[i][j]/rho[i][j]);
                // Alfvén速度
                double ca = std::sqrt((bx[i][j]*bx[i][j]+by[i][j]*by[i][j]+bz[i][j]*bz[i][j])/rho[i][j]);
                double cax = std::sqrt(bx[i][j]*bx[i][j]/rho[i][j]);
                double cay = std::sqrt(by[i][j]*by[i][j]/rho[i][j]);
                // 速進磁気音波速度
                double cfx = std::sqrt(0.5*(a*a+ca*ca+std::sqrt((a*a+ca*ca)*(a*a+ca*ca)-4.0*a*a*cax*cax)));
                double cfy = std::sqrt(0.5*(a*a+ca*ca+std::sqrt((a*a+ca*ca)*(a*a+ca*ca)-4.0*a*a*cay*cay)));
                // 最大固有値
                #pragma omp critical
                {
                    lxmax = max(std::abs(u[i][j])+cfx, lxmax);
                    lymax = max(std::abs(v[i][j])+cfy, lymax);
                }
            }
        }
        // CFL数をもとにdtを設定
        double dt = CFL*min(dx/lxmax, dy/lymax);
        // 磁場発散抑制のためのパラメータ (Dedner et al., 2002)
        double ch = CFL*min(dx, dy)/dt;
        double cd = std::exp(-dt*ch/cr);
        // 時間発展
        //q_n = euler(q, dt);
        //q_n = ssprk2(q, dt);
        q_n = ssprk3(q, dt);
        #pragma omp parallel for
        for (int i = 0; i < XN; i++) {
            for (int j = 0; j < YN; j++) {
                q_n[i][j][8] = cd*q_n[i][j][8];
                for (int k = 0; k < VN; k++) {
                    q[i][j][k] = q_n[i][j][k];
                }
            }
        }
        // ファイル出力
        if (static_cast<int>(t*PN/TL+1) != static_cast<int>((t-dt)*PN/TL+1)) {
            printf("n: %3lu, t: %.2f\n", n, t);
            const char *var[VN] = {"r", "p", "u", "v", "w", "bx", "by", "bz", "ps"};
            const char *dir = "../data/";
            char name[64][VN];
            char footer[64];
            std::ofstream os[VN];
            sprintf(footer, "_%03lu.csv", n);
            for (int k = 0; k < VN; k++) {
                sprintf(name[k], "%s%s%s", dir, var[k], footer);
                os[k].open(name[k]);
            }
            for (int j = 0; j < YN; j++) {
                for (int i = 0; i < XN-1; i++) {
                    os[0] << rho[i][j] << ",";
                    os[1] << p[i][j] << ",";
                    os[2] << u[i][j] << ",";
                    os[3] << v[i][j] << ",";
                    os[4] << w[i][j] << ",";
                    os[5] << bx[i][j] << ",";
                    os[6] << by[i][j] << ",";
                    os[7] << bz[i][j] << ",";
                    os[8] << psi[i][j] << ",";
                }
                os[0] << rho[XN-1][j] << "\n";
                os[1] << p[XN-1][j] << "\n";
                os[2] << u[XN-1][j] << "\n";
                os[3] << v[XN-1][j] << "\n";
                os[4] << w[XN-1][j] << "\n";
                os[5] << bx[XN-1][j] << "\n";
                os[6] << by[XN-1][j] << "\n";
                os[7] << bz[XN-1][j] << "\n";
                os[8] << psi[XN-1][j] << "\n";
            }
            n++;
        }
        t += dt;
        step++;
    }
    auto elapsed_time = std::chrono::system_clock::now()-start;
    std::cout << "Elapsed time: " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed_time).count()/1000.0 << " s\n";
    // 計算終了
    // 計算失敗の検知
    if (!std::isfinite(t)) {
        std::cout << "Computaion failed: step = " << step << "\n";
    }
    // y = πでの値
    std::cout << "   x    rho      p      u      v      w     bx     by     bz    psi\n";
    for (int i = 0; i < XN; i += 5) {
        int j = YN/2;
        std::printf("%4.2f  %5.2f  %5.2f  %5.2f  %5.2f  %5.2f  %5.2f  %5.2f  %5.2f  %5.2f\n",
                    dx*i, rho[i][j], p[i][j], u[i][j], v[i][j], w[i][j], bx[i][j], by[i][j], bz[i][j], psi[i][j]);
    }
    return 0;
}
