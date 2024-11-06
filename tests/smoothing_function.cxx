#include "sys.h"
#include "mpreal/mpreal.h"
#include "cairowindow/QuickGraph.h"
#include "cairowindow/Line.h"
#include "math/Polynomial.h"
#include "utils/ColorPool.h"
#include <glpk.h>
#include <boost/math/tools/polynomial.hpp>
#include <Eigen/Dense>
#include <stdexcept>
#include <iostream>
#include <cassert>
#include "debug.h"

#define SHOW_POLYNOMIAL_BASIS 1

using namespace mpfr;

void print(mpreal const& val)
{
  std::cout << val << std::endl;
}

// Define a few constants.
constexpr mp_prec_t precision_in_bits = 128;
constexpr double tolerance_dbl = 1e-30;
constexpr int polynomial_degree = 3;

// Make Eigen work with mpreal.
namespace Eigen {

template<>
struct NumTraits<mpreal> : NumTraits<double>      // Use double as the base.
{
  using Real = mpreal;
  using NonInteger = mpreal;
  using Nested = mpreal;
  enum
  {
    IsComplex = 0,
    IsInteger = 0,
    IsSigned = 1,
    RequireInitialization = 1,
    ReadCost = 1,
    AddCost = 3,
    MulCost = 3
  };
  static Real epsilon() { return machine_epsilon(); }
  static Real dummy_precision() { return tolerance_dbl; }
  static Real highest() { return pow(mpreal(10), mpreal::get_emax_max()); }
  static Real lowest() { return -pow(mpreal(10), mpreal::get_emax_max());; }
  static int digits10() { return mpreal::get_default_prec(); }
};

} // namespace Eigen

// Determining the value of S(C0) where C0 >= 0.
//
// Let r be the smallest root of C0 - 3u + u^3,
// then S(C0) is defined such that r = (1 - S) * (-sqrt(3)) + S * (-(cbrt(C0) + 1/cbrt(C0))).
//
// Note that if C0 = 0, the smallest root of -3u + u^3 is -sqrt(3), and thus S will be 0.
// While if C0 goes to infinity, then the smallest root of C0 - 3u + u^3 is -(cbrt(C0) + 1/cbrt(C0)), and this S will be 1.
//
// S(C0) smoothly and monotonically goes from 0 to 1, for C0 going from 0 to +infinity,
// where S(0.90375741845959156233304814223072905692...) ≈ 0.5.
//
// This function also returns the root r.
//
mpfr::mpreal S(mpfr::mpreal const& C0, mpfr::mpreal& root_out)
{
  static mpreal const sqrt_3(sqrt(3));
  static mpreal const C1(-3);

  // Calculate a reasonable guess for the root of C0 - 3u + u³.
  mpreal u = -sqrt_3;    // Just use -√3 because C0 < 1.

  // Use Halley's method to find the root.
  mpreal prev_u = 0;
  mpreal prev_prev_u;
  mpreal step = 0;
  mpreal prev_step;
  do
  {
    prev_prev_u = prev_u;
    prev_u = u;
    prev_step = step;
    mpreal Q_u = C0 + u * (u * u + C1);
    mpreal half_Qpp_u = 3 * u;
    mpreal Qp_u = half_Qpp_u * u + C1;
    step = -Q_u * Qp_u / (Qp_u * Qp_u - Q_u * half_Qpp_u);
    u += step;
  }
  while (step != 0.0 && u != prev_u && u != prev_prev_u);

  // Calculate the smoothing value S(C0).
  mpreal c = cbrt(C0);
  mpreal SC0 = (sqrt_3 + u) / (sqrt_3 - (c + 1.0 / c));

  return SC0;
}

// Find and return the value x in the range [x_min, x_max] such that S(x) = target +/- tolerance.
mpreal bracket_S(mpreal x_min, mpreal x_max, mpreal const& target, mpreal const& tolerance)
{
  mpreal root_dummy;
  mpreal y_min = S(x_min, root_dummy);
  mpreal y_max = S(x_max, root_dummy);

  // Check if the solution is bracketed.
  if (y_min > target || y_max < target)
    throw std::runtime_error("target is not bracketed within x_min and x_max.");

  mpreal y_mid = 0;
  while (x_max - x_min > tolerance)
  {
    mpreal x_mid = 0.5 * (x_min + x_max);
    mpreal prev_y_mid = y_mid;
    y_mid = S(x_mid, root_dummy);

    if (y_mid < target)
      x_min = x_mid;
    else
      x_max = x_mid;
  }

  // This appears to return a value of target +/- 0.5 * tolerance, but due to floating point
  // round off error, it doesn't guarantee that (if tolerance is of the order the current
  // resolution at the current value).
  return 0.5 * (x_min + x_max);
}

// Function to compute the Jacobi polynomial P_N^{(alpha, beta)}(t)
// Jacobi polynomials are a family of orthogonal polynomials on the interval
// [−1,1] with respect to the weight function: w(t) = (1 - t)^alpha (1 + t)^beta
// where alpha, beta > -1.
// In our case we need an orthogonal basis that is zero in t = -1, which maps
// to our C0 = 0. Therefore we need beta > 0.
mpreal jacobiP(int N, mpreal alpha, mpreal beta, mpreal t)
{
  mpreal const P_0{1};
  mpreal const P_1{alpha + 1.0 + (alpha + beta + 2.0) * (t - 1.0) / 2.0};

  // Special cases.
  if (N == 0)
    return P_0;
  if (N == 1)
    return P_1;

  // Initialize previous two values (where n = 2).
  mpreal Pnm2 = P_0;    // P_{n-2}
  mpreal Pnm1 = P_1;    // P_{n-1}
  mpreal Pn = 0;

  for (int n = 2; n <= N; ++n)
  {
                                                        // Legendre (alpha = beta = 0):
    mpreal const a = n + alpha;                         // a = n
    mpreal const b = n + beta;                          // b = n
    mpreal const c = 2 * n + alpha + beta;              // c = 2n

    mpreal a1 = (c - 1) * c * (c - 2);                  // 2n (2n - 1) (2n - 2)
    mpreal a2 = (c - 1) * (a - b) * (c - 2 * n);
    mpreal a3 = 2 * n * (c - n) * (c - 2);              // 2n * n * (2n - 2)
    mpreal a4 = 2 * (a - 1) * (b - 1) * c;

    a1 = a1 / a3;
    a2 = a2 / a3;
    a4 = a4 / a3;

    Pn = (a1 * t + a2) * Pnm1 - a4 * Pnm2;

    // Update Pnm2 and Pnm1 for next iteration.
    Pnm2 = Pnm1;
    Pnm1 = Pn;
  }

  return Pn;
}

math::Polynomial jacobiP_sym(int N, double alpha, double beta)
{
  // The constructor takes the number of coefficients, thus - the degree plus one.
  math::Polynomial P_N{N + 1 COMMA_CWDEBUG_ONLY("t")};                  // The result.

  // Special case.
  if (N == 0)
    P_N[0] = 1.0;
  else
  {
    // Initialize as-if N=1, thus as P_1.
    P_N[0] = 0.5 * (alpha - beta);
    P_N[1] = 0.5 * (alpha + beta) + 1.0;

    math::Polynomial P_N_minus_one{N + 1 COMMA_CWDEBUG_ONLY("t")};      // Initialize as-if N=1, thus as P_0.
    P_N_minus_one[0] = 1.0;

    math::Polynomial P_N_minus_two{N + 1 COMMA_CWDEBUG_ONLY("t")};

    for (int n = 2; n <= N; ++n)
    {
      // Update P_N_minus_two and P_N_minus_one.
      P_N_minus_two = std::move(P_N_minus_one);
      P_N_minus_one = std::move(P_N);
      P_N = math::Polynomial{N + 1 COMMA_CWDEBUG_ONLY("t")};

      double const a = n + alpha;
      double const b = n + beta;
      double const c = 2 * n + alpha + beta;

      double a1 = (c - 1) * c * (c - 2);
      double a2 = (c - 1) * (a - b) * (c - 2 * n);
      double a3 = 2 * n * (c - n) * (c - 2);
      double a4 = 2 * (a - 1) * (b - 1) * c;

      a1 = a1 / a3;
      a2 = a2 / a3;
      a4 = a4 / a3;

      // Assign P_N = a1 * P_N_minus_one * t.
      for (int i = 0; i < n; ++i)
        P_N[i + 1] = P_N_minus_one[i] * a1;

      P_N += a2 * P_N_minus_one - a4 * P_N_minus_two;
    }
  }

  return P_N;
}

#if 0
// Compute Chebyshev zeroes in the interval [a, b].
std::vector<mpreal> compute_chebyshev_zeros(int n, mpreal a, mpreal b)
{
  static mpreal const pi = const_pi();

  // First calculate the zeroes in the interval [-1, 1].
  std::vector<mpreal> zeroes;
  for (int k = 1; k <= n; ++k)
  {
    mpreal two_k_minus_one = 2 * k - 1;
    mpreal x_k = cos(two_k_minus_one / (2 * n) * pi);
    zeroes.push_back(x_k);
  }
  std::sort(zeroes.begin(), zeroes.end());

  // Next map the zeroes so that -1 --> a and 1 --> b.
  for (mpreal& zero : zeroes)
    zero = a + (zero + 1) / 2 * (b - a);

  return zeroes;
}

// Compute Newton's divided differences coefficients.
std::vector<mpreal> compute_newton_coefficients(std::vector<mpreal> const& x, std::vector<mpreal> const& y)
{
  int n = x.size();
  std::vector<mpreal> coef = y; // Divided difference coefficients.

  for (int j = 1; j < n; ++j)
    for (int i = n - 1; i >= j; --i)
      coef[i] = (coef[i] - coef[i - 1]) / (x[i] - x[i - j]);

  return coef;
}
#endif

using Polynomial = boost::math::tools::polynomial<mpreal>;

// Expand Newton's polynomial into monomial form.
Polynomial compute_newton_polynomial(std::vector<mpreal> const& x, std::vector<mpreal> const& coef)
{
  int n = coef.size();
  Polynomial P(coef[0]);        // Start with constant term c0.

  Polynomial term(1.0);         // Initialize term to 1
  for (int k = 1; k < n; ++k)
  {
    // Multiply term by (x - x_{k-1}).
    Polynomial factor({-x[k - 1], 1});  // Represents (x - x_{k-1}).
    term *= factor;

    // Add c_k * term to polynomial P.
    P += term * coef[k];
  }

  return P;
}

// Access coefficients from the polynomial
void print_polynomial_coefficients(Polynomial const& P)
{
  std::cout << "Polynomial coefficients (from lowest degree to highest):\n";
  for (size_t i = 0; i < P.size(); ++i)
    std::cout << "Coefficient of x^" << i << " is " << P[i] << '\n';
}

struct Extremum
{
  mpreal x_;
  mpreal y_;
  int sign_;            // -1: minimum; +1: maximum.
};

// Find and return the value x in the range [x_min, x_max] such that g(x) = is minimized/maximized within tolerance.
Extremum find_extremum(std::function<mpreal(mpfr::mpreal const&)> const& g, mpreal x_min, mpreal x_max, mpreal const& tolerance, int sign = 0)
{
  DoutEntering(dc::notice|continued_cf, "find_extremum(g, " << x_min << ", " << x_max << ", " << tolerance << ") --> ");

  Extremum result;

  static mpreal const phi = (sqrt(mpreal{5}) + 1) / 2;  // φ
  static mpreal const r = phi - 1;                      // φ − 1 = 0.61803398875
  static mpreal const c = 1 - r;                        // 2 - φ = 0.38196601125

  mpreal const tau = sqrt(tolerance);
  mpreal const x_epsilon = tau * (x_max - x_min);

  // Lambda for obtaining the x value of an interior point in the interval [x_min, x_max].
  auto x_at = [&x_min, &x_max](mpreal const& ratio){ return x_min + ratio * (x_max - x_min); };

  // Determine the g(x) value at c.
  mpreal const x_at_c = x_at(c);
  mpreal const gc = g(x_at_c);

  // Since this value is inbetween two zeroes (or one being an edge, but then it works too),
  // the sign tells us if we have to look for a maximum or a minimum.
  result.sign_ = sign == 0 ? (gc < 0 ? -1 : 1) : sign;
  Dout(dc::continued, "{sign_ = " << std::boolalpha << result.sign_ << "} = ");

  // Implement the Golden-section search algorithm, see https://en.wikipedia.org/wiki/Golden-section_search
  // This algorithm searches for a maximum, therefore, if sign_ is -1, we use -g(x) instead of g(x).

  // Specify the function to be minimized, f(x):
  auto f = [sign_ = result.sign_, &g](mpreal const& x){ return sign_ * g(x); };

  // The interval to be searched is {X0, X3}, X1 and X2 are interior points such that the three intervals have the ratios c : cr : c.
  std::array<mpreal, 4> X = {{x_min, x_at_c, x_at(c * phi), x_max}};

  // Calculate the functional values of F[i] and keep track of the index of the smallest value.
  std::array<mpreal, 4> F;
  // Initialize index 1 first because we already calculated that anyway.
  F[1] = result.sign_ * gc;
  int largest_index = 1;
  for (int i = 0; i < F.size(); ++i)
  {
    if (i == 1)
      continue;
    F[i] = f(X[i]);
    if (F[i] > F[largest_index])
      largest_index = i;
  }

  do
  {
    if (largest_index <= 1)
    {
      X[3] = std::move(X[2]);
      F[3] = std::move(F[2]);
      X[2] = std::move(X[1]);
      F[2] = std::move(F[1]);
      if (largest_index == 1)
        largest_index = 2;
      X[1] = r * X[2] + c * X[0];
      F[1] = f(X[1]);
      if (F[1] > F[largest_index])
        largest_index = 1;
    }
    else
    {
      X[0] = std::move(X[1]);
      F[0] = std::move(F[1]);
      X[1] = std::move(X[2]);
      F[1] = std::move(F[2]);
      if (largest_index == 2)
        largest_index = 1;
      X[2] = r * X[1] + c * X[3];
      F[2] = f(X[2]);
      if (F[2] > F[largest_index])
        largest_index = 2;
    }
  }
  while (X[3] - X[0] > x_epsilon || abs(F[0] - F[3]) > tolerance);

  result.x_ = (X[0] + X[3]) / 2;
  result.y_ = result.sign_ * (F[0] + F[3]) / 2;
  Dout(dc::finish, result.x_ << ", " << result.y_);
  return result;
}

double find_zero(std::function<mpreal(mpreal)> const& f, double x_min, double x_max, double tolerance)
{
  // Return the x-coordinate, in the range [x_min, x_max] at which f(x) ≈ 0.
  double f_min = f(x_min).toDouble();
  double f_max = f(x_max).toDouble();

  double x_mid, f_mid;

  while ((x_max - x_min) > tolerance)
  {
    x_mid = (x_min + x_max) / 2.0;
    f_mid = f(x_mid).toDouble();

    if (f_mid * f_min < 0)
    {
      x_max = x_mid;
      f_max = f_mid;
    }
    else
    {
      x_min = x_mid;
      f_min = f_mid;
    }
  }

  // Return the midpoint of the final interval.
  return (x_min + x_max) / 2.0;
}

void find_extrema(std::function<mpreal(mpreal)> const& E, mpreal const& zero, mpreal const& C0_half,
    std::vector<mpreal> const& zeroes, mpreal const& tolerance, std::vector<Extremum>& extrema,
    cairowindow::QuickGraph& graph2, utils::ColorPool<32>& color_pool)
{
  cairowindow::draw::PointStyle min_style({.color_index = color_pool.get_and_use_color(), .filled_shape = 1});
  cairowindow::draw::PointStyle max_style({.color_index = color_pool.get_and_use_color(), .filled_shape = 7});

  // The first extrema depends on the sign of the value at the start of the interval.
  int expected_sign = E(zero) > 0 ? 1 : -1;

  // For each interval between zeroes[i] and zeroes[i+1].
  for (size_t i = 0; i <= zeroes.size(); ++i)
  {
    mpreal c_min = i == 0 ? zero : zeroes[i - 1];
    mpreal c_max = i == zeroes.size() ? C0_half : zeroes[i];

    // Find the extremum in the interval.
    graph2.add_point({c_min.toDouble(), E(c_min).toDouble()}, min_style);
    graph2.add_point({c_max.toDouble(), E(c_max).toDouble()}, max_style);
    Extremum extremum = find_extremum(E, c_min, c_max, tolerance);

    if (extremum.sign_ != expected_sign)
    {
      // This should't fail if the initial expected_sign is correct.
      ASSERT(!extrema.empty());

      // We skipped over an extreme because the error function has more zeroes than N+1 zeroes (which is just the minimum)!
      // Look for the missed extreme between the last two extrema found. It is assumed we can only miss one extreme at a time
      // (two wouldn't be detected and three or more aint gonna happen).
      c_min = extrema.back().x_;
      c_max = extremum.x_;

      Dout(dc::notice, "Missed an extreme! Looking for a " << (expected_sign == 1 ? "maximum" : "minimum") <<
          " between " << c_min << " and " << c_max << ".");

      Extremum missed_extremum = find_extremum(E, c_min, c_max, tolerance, expected_sign);

      extrema.push_back(missed_extremum);
      expected_sign = -expected_sign;
    }

    extrema.push_back(extremum);

    // We expect minima and maxima to alternate.
    expected_sign = -expected_sign;
  }

  // Show the extrema. If found Z extra zeroes, we now have N+2+Z extrema.
  for (auto const& extremum : extrema)
    graph2.add_point({extremum.x_.toDouble(), extremum.y_.toDouble()});
}

int main()
{
  Debug(NAMESPACE_DEBUG::init());

  // Set the default floating-point precision.
  mpfr_set_default_prec(precision_in_bits);

  // Print everything with 38 digits precision.
  std::streamsize const old_precision = std::cout.precision(38);
  std::ios::fmtflags const old_flags = std::cout.setf(std::ios::fixed);

#if SHOW_POLYNOMIAL_BASIS
  // Jacobi polynomials test.
  cairowindow::QuickGraph graph_jacobi("Jacobi polynomials", "t", "J_N^{0,1}", {-1.0, 1.0}, {-1.0, 2.0});
  cairowindow::draw::LineStyle jacobi_style({.line_width = 1.0});
  namespace color = cairowindow::color;

  std::array<cairowindow::Color, 6> colors = {{ color::darkcyan, color::orange, color::green, color::red, color::purple, color::brown }};
#endif

  // Define the weight function.
  math::Polynomial sqrt_w(2 COMMA_CWDEBUG_ONLY("t"));   // 1 + t.
  sqrt_w[0] = 1.0;
  sqrt_w[1] = 1.0;
  math::Polynomial w = utils::square(sqrt_w);           // (1 - t)^alpha * (1 + t)^beta = (1 + t)^2.

  std::vector<math::Polynomial> polynomial_basis;
  for (int N = 0; N <= 5; ++N)
  {
    polynomial_basis.push_back(sqrt_w * jacobiP_sym(N, 0.0, 2.0));
#if SHOW_POLYNOMIAL_BASIS
    graph_jacobi.add_function([&, N](double t) -> double { return polynomial_basis[N](t); }, jacobi_style({.line_color = colors[N]}));
#endif
  }

  // Make sure the basis is orthogonal.
  for (int N1 = 0; N1 < 5; ++N1)
    for (int N2 = N1 + 1; N2 <= 5; ++N2)
    {
      auto product = polynomial_basis[N1] * polynomial_basis[N2];       // J_{N1}^{0,2}(t) J_{N2}^{0,2}(t) w(t)
      auto indefinite_integral = product.integrate();
      double definite_integral = indefinite_integral(1.0) - indefinite_integral(-1.0);
      ASSERT(std::abs(definite_integral) < 1e-15);
    }

#if SHOW_POLYNOMIAL_BASIS
  {
    Dout(dc::notice, "Orthogonal polynomial basis:");
    debug::Indent indent(2);
    for (math::Polynomial const& polynomial : polynomial_basis)
      Dout(dc::notice, "1/8 * (" << (8.0 * polynomial) << ")");
  }
  graph_jacobi.wait_for_keypress();
#endif

  // Loop over all values of C0 between 0 and 1 in steps of 1/128th.
  mpreal C0(0);
  mpreal const step(1.0 / 128.0);

  // Find the value of C0 for which S(C0) returns 0.5.
  long double const C0_half_approximation = 0.90375741845959156233304814223072905692L;
  mpreal const C0_half = bracket_S(std::nextafter(C0_half_approximation, 0.0L), std::nextafter(C0_half_approximation, 1.0L), 0.5, tolerance_dbl);

  std::cout << "C0_half = " << C0_half << std::endl;

  // Draw S(C0) on the interval [0, 0.9038].
  mpreal root_dummy;
  cairowindow::QuickGraph graph1("The smoothing function S(C0)", "C0", "S",
      {0.0, C0_half.toDouble()}, [&root_dummy](double C0) -> double { return S(C0, root_dummy).toDouble(); });

#if 0
  // Define zero.
  mpreal const zero(0);
  mpreal const sqrt_3(sqrt(3));
  mpreal const eps("1e-38");

  // Set the tolerance for the search.
  mpreal const tolerance{tolerance_dbl};

  // Initialize a vector to store extrema of E.
  std::vector<Extremum> extrema;

  // Draw E(C0) on the interval [0, 0.9038].
  cairowindow::QuickGraph graph2("The error function E(C0)", "C0", "E", {0.0, C0_half.toDouble()});

  utils::ColorPool<32> color_pool;
  namespace color = cairowindow::color;

  {
    // We will use an orthogal polynomial basis that is zero in x = 0; therefore we
    // don't need to fit through the point (0, 0) as that is already guaranteed.
    // However, we do want more density near zero. Lets use a distribution like
    // Chebyshev–Gauss–Lobatto nodes on the interval [0, 2 * C0_half] and then
    // only use the nodes on the interval [0, C0_half]. To guarantee that a node
    // at C0_half exists, there should be an odd number of Chebyshev–Gauss–Lobatto
    // nodes on the interval [0, 2 * C0_half]. Moreover, we need N+1 nodes to fit
    // an N-degree polynomial and we'll be ignoring the node at 0.
    // Therefore, the number of nodes on the [0, 2 * C0_half] must be 2(N+2)+1 = 2N+5.

    // Compute N+1 Chebyshev zeroes, where N is the degree of the polynomial that we want to fit to S.
    std::vector<mpreal> zeroes = compute_chebyshev_zeros(polynomial_degree + 1, 0, C0_half);
    zeroes[0] = 0.0;

    // Evaluate S(c) at the mapped nodes.
    std::vector<mpreal> S_values;
    for (mpreal c : zeroes)
      S_values.push_back(S(c, root_dummy));

    // Define our initial polynomial approximation.
    auto coefficients = compute_newton_coefficients(zeroes, S_values);
    Polynomial P = compute_newton_polynomial(zeroes, coefficients);

    graph1.add_function([&P](double C0) -> double { return (P.evaluate(C0) * C0).toDouble(); }, color::red);

    cairowindow::QuickGraph graph1_zoom("The smoothing function S(C0)", "C0", "S",
        {0.0, C0_half.toDouble()},
        [&root_dummy](double C0) -> double { return S(C0, root_dummy).toDouble(); });
    graph1_zoom.add_function([&P, &root_dummy](double C0) -> double
        {
          mpreal SC0 = S(C0, root_dummy);
          mpreal diff = P.evaluate(C0) - SC0;
          return (SC0 + 100 * diff).toDouble();
        }, color::red);

    {
      // Define the error function E(c) = S(c) - c * P(c).
      mpreal root_value;
      auto E = [&P, &root_value, &sqrt_3, &eps](mpreal const& c)
        {
          mpreal Sc = S(c, root_value);
          mpreal Pc = P.evaluate(c);
          mpreal cbrtc = cbrt(c);
          if (cbrtc < 1e-6)
            cbrtc = 1e-6;
          if (abs(Pc) < 1e-30)
          {
            if (Pc < 0.0)
              Pc = -1e-30;
            else
              Pc = 1e-30;
          }
          mpreal result = (Sc - Pc) * ((sqrt_3 - cbrtc) * cbrtc - 1) / abs(Pc * cbrtc);
          ASSERT(std::isfinite(result.toDouble()));
          return result;
        };

      // Draw E(C0) on the interval [0, 0.9038].
      graph2.add_function([&E](double C0) -> double { return E(C0).toDouble(); });

      find_extrema(E, zero, C0_half, zeroes, tolerance, extrema, graph2, color_pool);
    }
  }

  do
  {
    // Lets define M = N+Z, so that M=N if we didn't find any extra zeroes.
    int const M = extrema.size() - 2;

    // Calculate the initial M+2, or more, Chebyshev control points.
    std::vector<mpreal> extrema_S_values;
    std::vector<mpreal> extrema_root_values;
    mpreal extrema_root;
    for (auto const& extremum : extrema)
    {
      extrema_S_values.push_back(S(extremum.x_, extrema_root));
      extrema_root_values.push_back(extrema_root);
    }

    // Find the next polynomial P(x) = c0 + c1 x + c2 x^2 + ... + c_M x^M
    // such that P(extrema[i].x_) + extrema[i].sign_ * E = extrema[i].y_,
    // with unknowns c_i and E.
    int const num_equations = M + 2;

    // Set up the system using Eigen matrices and vectors.
    using MatrixXmp = Eigen::Matrix<mpreal, Eigen::Dynamic, Eigen::Dynamic>;
    using VectorXmp = Eigen::Matrix<mpreal, Eigen::Dynamic, 1>;
    MatrixXmp A(num_equations, num_equations);
    VectorXmp b(num_equations);

    // Fill in the matrix A and vector b.
    for (int i = 0; i < num_equations; ++i)
    {
      mpreal x_i = extrema[i].x_;
      mpreal sign = extrema[i].sign_;

      // The value of C0 at point i is x_i.
      mpreal cbrt_C0 = cbrt(x_i);
      mpreal scaling = abs(extrema_root_values[i] / (sqrt_3 - cbrt_C0 - 1 / cbrt_C0));

      // Fill in the polynomial terms.
      for (int j = 0; j <= M; ++j)
        A(i, j) = pow(x_i, j);

      // Fill in the term for E.
      A(i, M + 1) = sign * scaling;

      // Fill in b_i.
      b(i) = extrema_S_values[i];
    }

    // Solve the equation A⋅x = b.
    VectorXmp x = A.colPivHouseholderQr().solve(b);

    // Extract coefficients c0 to cM and E.
    std::vector<mpreal> coefficients(M + 1);
    for (int i = 0; i <= M; ++i)
      coefficients[i] = x(i);
    mpreal Err = x(M + 1);

    Polynomial P{coefficients};
    std::cout << "Err = " << Err << std::endl;
    std::cout << "Polynomial:\n";
    for (int i = 0; i <= M; ++i)
      std::cout << "  " << coefficients[i] << " * x^" << i << std::endl;

    graph1.add_function([&P](double C0) -> double { return P.evaluate(C0).toDouble(); }, color::green);

    // Define the error function E(c) = S(c) - P(c).
    mpreal root_value;
    auto E = [&P, &root_value, &sqrt_3, &eps](mpreal const& c)
      {
        mpreal Sc = S(c, root_value);
        mpreal Pc = P.evaluate(c);
        mpreal cbrtc = cbrt(c);
        return (Sc - Pc) * ((sqrt_3 - cbrtc) * cbrtc - 1) / max(abs(Pc * cbrtc), eps);
      };

    // Draw E(C0) on the interval [0, 0.9038].
    graph2.add_function([&E](double C0) -> double { return E(C0).toDouble(); }, color::blue);

    // Find the roots of E.
    std::vector<mpreal> roots;

    // Run over the M+1 intervals between the previously found extrema and store the value of the new E.
    std::vector<math::Point> points;
    for (int i  = 0; i < M + 1; ++i)
    {
      // The begin and end of interval i.
      double x_min = extrema[i].x_.toDouble();
      double x_max = extrema[i + 1].x_.toDouble();
      // Evaluate the function S at three intermediate points, as well the boundaries of the intervals.
      for (int j = 0; j <= 3; ++j)
      {
        double x = x_min + 0.25 * j * (x_max - x_min);
        points.emplace_back(x, E(x).toDouble());
      }
    }
    // add the last extremum too.
    {
      double x = extrema[M + 1].x_.toDouble();
      points.emplace_back(x, E(x).toDouble());
    }

    // Find all intervals that contain a zero of E.
    cairowindow::draw::PointStyle zero_style({.color_index = color_pool.get_and_use_color(), .filled_shape = 5});
    std::vector<mpreal> zeroes;
    for (int i = 0; i < points.size() - 1; ++i)
    {
      if (std::copysign(1.0, points[i].y()) != std::copysign(1.0, points[i + 1].y()))
      {
        zeroes.emplace_back(find_zero(E, points[i].x(), points[i + 1].x(), 1e-6));
        graph2.add_point({zeroes.back().toDouble(), 0.0}, zero_style);
      }
    }

    // Find and store all extrema.
    extrema.clear();
    find_extrema(E, zero, C0_half, zeroes, tolerance, extrema, graph2, color_pool);

    // Give a chance to view the graphs before clearing.
    std::cin.get();
    graph2.clear();
  }
  while (true);
#endif

  graph1.wait_for_keypress();
}
