#ifndef MHD_HH
#define MHD_HH

#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <array>
#include <tuple>
#include <cmath>

// グリッド数
constexpr size_t XN = 200, YN = 200;
// ファイル出力回数
constexpr size_t PN = 100;
// 保存変数の個数
constexpr size_t VN = 9;
// 計算領域のサイズ
constexpr double XL = 2.0*M_PI, YL = 2.0*M_PI;
constexpr double TL = M_PI;
// グリッドサイズ
constexpr double dx = XL/XN, dy = YL/YN;
// CFL数
constexpr double CFL = 0.4;
// 比熱比
constexpr double gam = 5.0/3.0;
// 磁場発散抑制のためのパラメータ (Dedner et al., 2002)
constexpr double cr = 0.18;

// 多次元配列
template<typename T, size_t... N>
struct ndarray;
template<typename T>
struct ndarray<T> {
    using type = T;
};
template<typename T, size_t N, size_t... Rest>
struct ndarray<T, N, Rest...> {
    using type = std::array<typename ndarray<T, Rest...>::type, N>;
};
template<typename T, size_t... N>
using ndarray_t = typename ndarray<T, N...>::type;

// 可変長max
template<typename T, typename... Rest>
constexpr T max(T a, Rest... rest) {
    if constexpr (sizeof...(rest) == 0) {
        return a;
    } else {
        auto b = max(rest...);
        return a > b ? a : b;
    }
}
// 可変長min
template<typename T, typename... Rest>
constexpr T min(T a, Rest... rest) {
    if constexpr (sizeof...(rest) == 0) {
        return a;
    } else {
        auto b = min(rest...);
        return a < b ? a : b;
    }
}
// 符号関数
template<typename T>
constexpr T sign(T a) {
    return a > 0 ? 1 : a < 0 ? -1 : 0;
}
// MINMOD関数
double minmod(double a, double b);
// 中央値
double median(double a, double b, double c);
// HLLDリーマンソルバ (Miyoshi and Kusano, 2005)
auto hlldx(const ndarray_t<double, VN> &ql, const ndarray_t<double, VN> &qr, double ch);
auto hlldy(const ndarray_t<double, VN> &ql, const ndarray_t<double, VN> &qr, double ch);
// 1次精度風上差分
auto upwindx(const ndarray_t<double, XN, YN, VN> &q, int i, int j);
auto upwindy(const ndarray_t<double, XN, YN, VN> &q, int i, int j);
// 2次精度MUSCL (minmod) (van Leer, 1979)
auto musclx(const ndarray_t<double, XN, YN, VN> &q, int i, int j);
auto muscly(const ndarray_t<double, XN, YN, VN> &q, int i, int j);
// 5次精度MP5 (Suresh and Huynh, 1997)
auto mp5x(const ndarray_t<double, XN, YN, VN> &q, int i, int j);
auto mp5y(const ndarray_t<double, XN, YN, VN> &q, int i, int j);
// dq/dtの計算
auto dqdt(const ndarray_t<double, XN, YN, VN> &q, double ch);
// 1次精度Euler法
ndarray_t<double, XN, YN, VN> euler(const ndarray_t<double, XN, YN, VN> &q, double dt);
// 2次精度Runge-Kutta
ndarray_t<double, XN, YN, VN> ssprk2(const ndarray_t<double, XN, YN, VN> &q, double dt);
// 3次精度Runge-Kutta
ndarray_t<double, XN, YN, VN> ssprk3(const ndarray_t<double, XN, YN, VN> &q, double dt);

#endif // MHD_HH
