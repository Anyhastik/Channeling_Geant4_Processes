#ifndef CRYSTAL_FUNCTIONS_HPP
#define CRYSTAL_FUNCTIONS_HPP

#include "phys_const.hh"

/**
* Класс содержит необходимые для работы с кристаллом функции
*/
class Crystal_Functions
{
public:
    Crystal_Functions() {}
    virtual ~Crystal_Functions() {}

    template <typename value_type, typename function_type>
    value_type derivative(const value_type x, const value_type dx, function_type function)
    {
        /*! \brief Compute the derivative of function using a 3-point central difference rule of O(dx^6).
          \tparam value_type, floating-point type, for example: `double` or `cpp_dec_float_50`
          \tparam function_type

          \param x Value at which to evaluate derivative.
          \param dx Incremental step-size.
          \param function Function whose derivative is to computed.

          \return derivative at x.
        */

        static_assert(false == std::numeric_limits<value_type>::is_integer, "value_type must be a floating-point type!");

        const value_type dx2(dx * 2U);
        const value_type dx3(dx * 3U);
        // Difference terms.
        const value_type m1((function(x + dx) - function(x - dx)) / 2U);
        const value_type m2((function(x + dx2) - function(x - dx2)) / 4U);
        const value_type m3((function(x + dx3) - function(x - dx3)) / 6U);
        const value_type fifteen_m1(m1 * 15U);
        const value_type six_m2(m2 * 6U);
        const value_type ten_dx(dx * 10U);
        return ((fifteen_m1 - six_m2) + m3) / ten_dx; // Derivative.
    }

    ld1 ducr(ld1 t, ld1 x, array<ld1, 2*nmax+1> cpot, ld1 dp, string cr_t);

    ld1 ducr_2(ld1 t, ld1 x, array<ld1, 2*nmax+1> cpot, ld1 dp, string cr_t);

    ld1 ucr(ld1 x, array<ld1, 2*nmax+1> cpot, ld1 dp, string cr_t);

    ld1 eve(ld1 i, array<ld1, 2*nmax+1> cpot, ld1 dp, ld1 tetta1, string cr_t);

    ld1 functz(ld1 t1z, ld1 z1, ld1 z2, ld1 x2, ld1 k0);

    static void print(vector<ld1> &vec);

    ld1 dWdE(ld1 Nt, ld1 w, ld1 TT, ld1 integral[nbeam]);

    ld1 integr_f(ld1 x, ld1 eve2, array<ld1, 2*nmax+1> cpot, ld1 dp, string cr_t);

    ld1 simpson(ld1 a, ld1 b, int n0, ld1 eve2, array<ld1, 2*nmax+1> cpot,ld1 dp, string cr_t);

    ld1 find(ld1 x, ld1 eps, ld1 eve1, array<ld1, 2*nmax+1> cpot, ld1 dp, string cr_t);

    ld1 theta(ld1 x);

    ld1 omega(ld1 n, ld1 TT);

    ld1 trapezoid(vector<ld1> times, vector<ld1> velocities, int n, ld1 TT);
};

#endif /*CRYSTAL_FUNCTIONS_HPP*/