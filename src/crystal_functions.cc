#ifdef _MSC_VER
#pragma warning(disable : 4996) // assignment operator could not be generated.
#endif


#include <iostream>
#include <iomanip>
#include <limits>
#include <type_traits>
#include <cmath>

#include <boost/math/tools/assert.hpp>
#include <boost/math/special_functions/next.hpp> // for float_distance

#include <boost/multiprecision/cpp_dec_float.hpp>
using boost::multiprecision::number;
using boost::multiprecision::cpp_dec_float;

// Re-compute using 5 extra decimal digits precision (22) than double (17).
#define MP_DIGITS10 unsigned(std::numeric_limits<double>::max_digits10 + 5)

typedef cpp_dec_float<MP_DIGITS10> mp_backend;
typedef number<mp_backend> mp_type;

#include "crystal_functions.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

/**
* Функция вычисления производной потенциала кристалла в конкретной точке, домноженной на коэффициент
*/
ld1 Crystal_Functions::ducr(ld1 t, ld1 x, array<ld1, 2*nmax+1> cpot, ld1 dp, string cr_t)
{
 	ld1 u0 = 0;
    ld1 x0;
    if (cr_t == "w")
        x0 = 0.0;
    else
        x0 = dp / 2;
    for (int m = 2; m <= 21; m++)
    {
        const mp_type mp =
            derivative(mp_type(mp_type(x)),              // x = 3/2
                       mp_type(mp_type(1) / 100000000U), // Step size 10^-8.
                       [&](const mp_type &x) -> mp_type
                       {
                           mp_type pf = 2 * pi * (m - 1) * (x + x0) / dp;
                           mp_type cs = cos(pf);
                           return 2 * cpot[19 + m] * cs; // Function
                       });
        ld1 dh = mp.convert_to<ld1>(); // Convert to closest double.
        u0 += dh;
    }
   return -u0*coeff;
}

/**
* Функция вычисления производной потенциала кристалла в конкретной точке
*/
ld1 Crystal_Functions::ducr_2(ld1 t, ld1 x, array<ld1, 2*nmax+1> cpot, ld1 dp, string cr_t)
{
 	ld1 u0 = 0;
    ld1 x0;
    if (cr_t == "w")
        x0 = 0.0;
    else
        x0 = dp / 2;
    for (int m = 2; m <= 21; m++)
    {
        const mp_type mp =
            derivative(mp_type(mp_type(x)),              // x = 3/2
                       mp_type(mp_type(1) / 100000000U), // Step size 10^-8.
                       [&](const mp_type &x) -> mp_type
                       {
                           mp_type pf = 2 * pi * (m - 1) * (x + x0) / dp;
                           mp_type cs = cos(pf);
                           return 2 * cpot[19 + m] * cs; // Function
                       });
        ld1 dh = mp.convert_to<ld1>(); // Convert to closest double.
        u0 += dh;
    }
   return u0;
}

/**
* Функция вычисления потенциала кристалла в конкретной точке
*/
ld1 Crystal_Functions::ucr(ld1 x, array<ld1, 2*nmax+1> cpot, ld1 dp, string cr_t)
{
    ld1 x0;
    if (cr_t == "w")
        x0 = 0.0;
    else
        x0 = dp / 2;
    ld1 u = cpot[20];
    for (int m = 2; m <= 21; m++)
        u += 2 * cos(2 * pi * (m - 1) * (x + x0) / dp) * cpot[m + 19];
    return u;
}

ld1 Crystal_Functions::eve(ld1 i, array<ld1, 2*nmax+1> cpot, ld1 dp, ld1 tetta1, string cr_t)
{
    return ucr(1.0 / (nbeam + 1) * i * dp / 2, cpot, dp, cr_t) + pow(mom, 2) * pow(10, 12) * pow(sin(tetta1), 2) / (2 * gama * m0);
}

/**
* Функция правой части систесы дифференциальных уравнений
*/
ld1 Crystal_Functions::functz(ld1 t1z, ld1 z1, ld1 z2, ld1 x2, ld1 k0) // the right part of the eq: z''=...
{
    ld1 numirator = -z2 * x2 * k0;
    ld1 square_root = 1 - (pow(x2, 2) + pow(z2, 2)) / (pow(c0, 2) * pow(10, 20));
    ld1 denominator = (pow(10, 20) * pow(c0, 2) * square_root) + pow(z2, 2);
    return numirator / denominator;
}

/**
*Тета-функция
*/
ld1 Crystal_Functions::theta(ld1 x)
{
    if (x < 0)
        return 0;
    else
        return 1;
}

/**
* Функция вычисления частоты
*/
ld1 Crystal_Functions::omega(ld1 n, ld1 TT)
{
    return n * 2 * pi / TT;
}

void Crystal_Functions::print(vector<ld1> &vec)
{
    for (auto j = vec.begin(); j != vec.end(); ++j) //vector<ld1>::iterator
        cout << *j << ", ";
}

/**
* Функция расчета излучения
*/
ld1 Crystal_Functions::dWdE(ld1 Nt, ld1 w, ld1 TT, ld1 integral[nbeam])
{
    ld1 multiplier = pow(10, -2) * pow(echarge, 2) / (hbar * pow(c0, 2)) * Nt;
    ld1 sum = 0;
    for (int n = 1; n <= 5; n++)    //число гармоник
    {
        sum += theta(2 * pow(gama, 2) * omega(n, TT) - w) * w * (1 - 2 * w / (2 * pow(gama, 2) * omega(n, TT)) + 2 * pow(w / (2 * pow(gama, 2) * omega(n, TT)), 2)) / pow(TT, 2) * integral[n - 1];
    }
    sum *= multiplier;
    return sum;
}

/**
* Функция которая впоследствии интегрируется
*/
ld1 Crystal_Functions::integr_f(ld1 x, ld1 eve2, array<ld1, 2*nmax+1> cpot, ld1 dp, string cr_t)
{
    return 1 / sqrt(eve2 - ucr(x, cpot, dp, cr_t));
}

/**
* Функция интегрирования методом Симпсона
*/
ld1 Crystal_Functions::simpson(ld1 a, ld1 b, int n0, ld1 eve2, array<ld1, 2*nmax+1> cpot, ld1 dp, string cr_t)
{
    ld1 tmp1 = 0;
    for (int i = 1; i <= n0 - 1; ++i)
        tmp1 += integr_f(a + (b - a) / n0 * i, eve2, cpot, dp, cr_t);
    ld1 tmp2 = 0;
    for (int i = 1; i <= n0; ++i)
        tmp2 += integr_f(a + (b - a) / 2 / n0 * (2 * i - 1), eve2, cpot, dp, cr_t);

    return (b - a) / 6 / n0 * (integr_f(a + 0.00001, eve2, cpot, dp, cr_t) + integr_f(b - 0.00001, eve2, cpot, dp, cr_t) + 2 * tmp1 + 4 * tmp2);
}

/**
* Функция нахождения корня нелинейного уравнения методом Ньютона
*/
ld1 Crystal_Functions::find(ld1 x, ld1 eps, ld1 eve1, array<ld1, 2*nmax+1> cpot, ld1 dp, string cr_t)
{
    ld1 f, df;
    int iter = 0;
    cout << "x0= " << x << " ";
    do
    {
        f = ucr(x, cpot, dp, cr_t) - eve1;
        df = ducr_2(0,x, cpot, dp, cr_t);
        x = x - f / df;
        iter++;
    } while (fabs(f) > eps && iter < 1000); // 20000
    cout << iter << " iterations" << endl;
    return x;
}

/**
* Функция интегрирования методом трапеций
*/
ld1 Crystal_Functions::trapezoid(vector<ld1> times, vector<ld1> velocities, int n, ld1 TT)
{
    ld1 omega1 = omega(n, TT); //число гармоник = 5 n=[0,5]
    ld1 tmp1 = 0;
    ld1 tmp2 = 0, i = 1;
    while (times.at(i) <= TT)
    {
        tmp1 += (velocities[i - 1] * cos(omega1 * times[i - 1]) + velocities[i] * cos(omega1 * times[i])) / 2 * (times[i] - times[i - 1]) / (c0 * pow(10, 10));
        tmp2 += (velocities[i - 1] * sin(omega1 * times[i - 1]) + velocities[i] * sin(omega1 * times[i])) / 2 * (times[i] - times[i - 1]) / (c0 * pow(10, 10));
        i++;
    }
    return (pow(tmp1, 2) + pow(tmp2, 2));
}
