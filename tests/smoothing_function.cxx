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

#define SHOW_POLYNOMIAL_BASIS 0
#define SHOW_SMOOTHING_FUNCTION 1
#define SHOW_SMOOTHING_FUNCTION_ZOOM 0
#define SHOW_ERROR_FUNCTION 1
#define SHOW_WEIGH_FUNCTION 0
#define SHOW_REL_ERROR_ROOT 1
#define TILL_C0_half 1

using namespace mpfr;

void print(mpreal const& val)
{
  std::cout << val << std::endl;
}

// Define a few constants.
constexpr mp_prec_t precision_in_bits = 128;
constexpr double tolerance_dbl = 1e-30;
#if TILL_C0_half
constexpr int polynomial_degree = 2;
#else
constexpr int polynomial_degree = 3;
#endif

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
  mpreal SC0 = (sqrt_3 + u) * c / (c * (sqrt_3 - c) - 1.0);
  root_out = u;

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

using Polynomial = boost::math::tools::polynomial<mpreal>;

// Function to compute the Chebyshev polynomial P_N(t)
//
// Chebyshev polynomials are orthogonal polynomials on the interval [-1, 1] under
// the weight function 1/sqrt(1 - t^2), but are not themselves orthogonal.
//
Polynomial chebyshevP(int N)
{
  Polynomial P_N;                               // The result.

  // Special case.
  if (N == 0)
    P_N = Polynomial({1});
  else
  {
    // Initialize as-if N=1, thus as P_1.
    P_N = Polynomial({0, 1});

    Polynomial P_N_minus_one({1});              // Initialize as-if N=1, thus as P_0.
    Polynomial P_N_minus_two;

    for (int n = 2; n <= N; ++n)
    {
      // Update P_N_minus_two and P_N_minus_one.
      P_N_minus_two = std::move(P_N_minus_one);
      P_N_minus_one = std::move(P_N);

      P_N = P_N_minus_one * Polynomial({0, 2}) - P_N_minus_two;
    }
  }

  return P_N;
}

// Function to compute the Jacobi polynomial P_N^{(alpha, beta)}(t)
//
// Jacobi polynomials are a family of orthogonal polynomials on the interval
// [−1,1] with respect to the weight function: w(t) = (1 - t)^alpha (1 + t)^beta
// where alpha, beta > -1.
//
// In our case we need an orthogonal basis that is zero in t = -1, which maps
// to our C0 = 0. Therefore we need beta > 0.
Polynomial jacobiP(int N, double alpha, double beta)
{
  Polynomial P_N;                               // The result.

  // Special case.
  if (N == 0)
    P_N = Polynomial({1});
  else
  {
    // Initialize as-if N=1, thus as P_1.
    P_N = Polynomial({0.5 * (alpha - beta), 0.5 * (alpha + beta) + 1.0});

    Polynomial P_N_minus_one({1});              // Initialize as-if N=1, thus as P_0.
    Polynomial P_N_minus_two;

    for (int n = 2; n <= N; ++n)
    {
      // Update P_N_minus_two and P_N_minus_one.
      P_N_minus_two = std::move(P_N_minus_one);
      P_N_minus_one = std::move(P_N);
      P_N = Polynomial(0);

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

      P_N = P_N_minus_one * Polynomial({0, a1}) + a2 * P_N_minus_one - a4 * P_N_minus_two;
    }
  }

  return P_N;
}

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

#if 0
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
#endif

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

void find_extrema(mpreal const& begin, mpreal const& end, std::function<mpreal(mpreal)> const& E,
    std::vector<mpreal> const& zeroes, mpreal const& tolerance, std::vector<Extremum>& extrema
#if SHOW_ERROR_FUNCTION
    , cairowindow::QuickGraph& error_function_graph, utils::ColorPool<32>& color_pool
#endif
    )
{
  DoutEntering(dc::notice, "find_extrema(" << begin << ", " << end << ", E(c), " << zeroes << ", " << tolerance << ", ...)");

  cairowindow::draw::PointStyle min_style({.color_index = color_pool.get_and_use_color(), .filled_shape = 1});
  cairowindow::draw::PointStyle max_style({.color_index = color_pool.get_and_use_color(), .filled_shape = 7});

  // The first extrema depends on the sign of the value at the start of the interval.
  int expected_sign = E(begin) > 0 ? 1 : -1;

  // For each interval between zeroes[i] and zeroes[i+1].
  for (size_t i = 0; i <= zeroes.size(); ++i)
  {
    mpreal c_min = i == 0 ? begin : zeroes[i - 1];
    mpreal c_max = i == zeroes.size() ? end : zeroes[i];

    // Find the extremum in the interval.
#if SHOW_ERROR_FUNCTION
    error_function_graph.add_point({c_min.toDouble(), E(c_min).toDouble()}, min_style);
    error_function_graph.add_point({c_max.toDouble(), E(c_max).toDouble()}, max_style);
#endif
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

#if SHOW_ERROR_FUNCTION
  // Show the extrema. If found Z extra zeroes, we now have N+2+Z extrema.
  for (auto const& extremum : extrema)
    error_function_graph.add_point({extremum.x_.toDouble(), extremum.y_.toDouble()});
#endif
}

// Function to compose two polynomials: Q(x) = P(t(x))
Polynomial compose(Polynomial const& P, Polynomial const& t)
{
  Polynomial Q(0);              // Initialize Q(x) = 0.
  Polynomial t_power(1);        // Start with t(x)^0 = 1.

  for (std::size_t k = 0; k <= P.degree(); ++k)
  {
    if (k > 0)
      t_power = t_power * t;    // Compute t(x)^k.
    Q += P[k] * t_power;        // Accumulate Q(x) += p_k * t(x)^k.
  }

  return Q;
}

// Function to transform P(t) to Q(x) where t = 2 * (x - begin) / (end - begin) - 1.
Polynomial transform_Pt_to_Qx(mpreal begin, mpreal end, Polynomial const& Pt)
{
  // Define the coefficients a and b for t(x).
  mpreal const a = 2 / (end - begin);
  mpreal const b = -(a * begin + 1);

  // Create t(x) = b + a x.
  Polynomial t({b, a});

  // Compose P(t(x)) to get Q(x).
  Polynomial Qx = compose(Pt, t);

  return Qx;
}

int main()
{
  Debug(NAMESPACE_DEBUG::init());

  // Set the default floating-point precision.
  mpfr_set_default_prec(precision_in_bits);

  // Print everything with 38 digits precision.
  std::streamsize const old_precision = std::cout.precision(38);
  std::ios::fmtflags const old_flags = std::cout.setf(std::ios::fixed);

  // Find the value of C0 for which S(C0) returns 0.5.
  long double const C0_half_approximation = 0.90375741845959156233304814223072905692L;
  mpreal const C0_half = bracket_S(std::nextafter(C0_half_approximation, 0.0L), std::nextafter(C0_half_approximation, 1.0L), 0.5, tolerance_dbl);

  std::cout << "0.5 = S(" << C0_half << ")" << std::endl;

  // Define the interval over which we need to fit a polynomial.
#if TILL_C0_half
  mpreal const begin = 0;
  mpreal const end = C0_half;
#else
  mpreal const begin = C0_half;
  mpreal const end = 6.82;
#endif

#if SHOW_POLYNOMIAL_BASIS
  // Jacobi polynomials test.
  cairowindow::QuickGraph orthogonal_basis_graph("Orthogonal basis", "t",
#if TILL_C0_half
      "JacobiP_n^{0,2}(t) (1 + t)",
#else
      "ChebyshevP_n(t)",
#endif
      {begin.toDouble(), end.toDouble()}, {-1.0, 2.0});
  cairowindow::draw::LineStyle orthogonal_basis_style({.line_width = 1.0});
  namespace color = cairowindow::color;

  std::array<cairowindow::Color, 6> colors = {{ color::darkcyan, color::orange, color::green, color::red, color::purple, color::brown }};
#endif

  std::vector<Polynomial> polynomial_basis;

#if TILL_C0_half
  // Define the weight function.
  Polynomial sqrt_w({1, 1});                            // 1 + t.
  Polynomial w = utils::square(sqrt_w);                 // w(t) = (1 - t)^alpha * (1 + t)^beta = (1 + t)^2.

  // N is the degree of the jacobi polynomials, which are multiplied with 1+t to get
  // the orthogonal basis, which is thus one degree higher.
  // Therefore here N runs from 0 till (but not including) polynomial_degree.
  for (int N = 0; N < polynomial_degree; ++N)
  {
    polynomial_basis.push_back(transform_Pt_to_Qx(begin, end, sqrt_w * jacobiP(N, 0.0, 2.0)));
#else // TILL_C0_half
  for (int N = 0; N <= polynomial_degree; ++N)
  {
    polynomial_basis.push_back(transform_Pt_to_Qx(begin, end, chebyshevP(N)));
#endif // TILL_C0_half
#if SHOW_POLYNOMIAL_BASIS
    orthogonal_basis_graph.add_function(
        [&, N](double t) -> double { return polynomial_basis[N](t).toDouble(); },
        orthogonal_basis_style({.line_color = colors[N]}));
#endif
  }
//}

// The Chebychev polynomials of the first kind are not orthogonal.
#if TILL_C0_half
  // Make sure the basis is orthogonal.
  // N1 and N2 are the indices into polynomial_basis and therefore the degree of the Jacobi polynomial.
  // Hence, again, N2 only runs till (but not including) polynomial_degree.
  for (int N1 = 0; N1 < polynomial_basis.size() - 1; ++N1)
    for (int N2 = N1 + 1; N2 < polynomial_basis.size(); ++N2)
    {
      auto product = polynomial_basis[N1] * polynomial_basis[N2];       // J_{N1}^{0,2}(t) J_{N2}^{0,2}(t) w(t)
      auto indefinite_integral = product.integrate();
      mpreal definite_integral = indefinite_integral(end) - indefinite_integral(begin);
      ASSERT(std::abs(definite_integral.toDouble()) < 1e-15);
    }
#endif

#if SHOW_POLYNOMIAL_BASIS
  {
    Dout(dc::notice, "Polynomial basis:");
    debug::Indent indent(2);
    for (Polynomial const& polynomial : polynomial_basis)
      Dout(dc::notice, polynomial);
  }
  orthogonal_basis_graph.wait_for_keypress();
#endif

  // Draw S(C0) on the interval [0, 0.9038] or [0, 4].
  mpreal root_dummy;
#if SHOW_SMOOTHING_FUNCTION
  cairowindow::QuickGraph smoothing_function_graph("The smoothing function S(C0)", "C0", "S",
#if TILL_C0_half
      {0.0, C0_half.toDouble()},
#else
      {begin.toDouble(), end.toDouble()},
#endif
      [&root_dummy](double C0) -> double { return S(C0, root_dummy).toDouble(); });
#endif

  // Define zero.
  mpreal const zero(0);
  mpreal const sqrt_3(sqrt(3));
  mpreal const eps("1e-38");

  // Set the tolerance for the search.
  mpreal const tolerance{tolerance_dbl};

  // Initialize a vector to store extrema of E.
  std::vector<Extremum> extrema;

#if SHOW_ERROR_FUNCTION
  // Draw E(C0) on the interval [0, 0.9038] or [0, 10].
  // If the absolute value stays under 0.012 then we can find the actual root to the full resolution of a double in only two Halley iterations.
  cairowindow::QuickGraph error_function_graph("Relative error in the root", "C0", "E",
#if TILL_C0_half
      {0.0, C0_half.toDouble()},
#else
      {0.0, 10.0},
#endif
      {-0.02, 0.02});

  utils::ColorPool<32> color_pool;
  namespace color = cairowindow::color;
#if SHOW_REL_ERROR_ROOT
//  cairowindow::QuickGraph rel_error_root_graph("Relative error of cbrt(C0) + 1/cbrt(C0) vs real root", "C0", "rel err",
//      {0.0, 10.0}, {-0.02, 0.005});
  error_function_graph.add_function([&](double c) -> double {
      mpreal root_value;
      mpreal Sc = S(c, root_value);
      mpreal cbrtc = cbrt(c);
      return ((-(cbrtc + 1/cbrtc) - root_value) / abs(root_value)).toDouble();
    }, color::blue);
//  rel_error_root_graph.wait_for_keypress();
#endif
#endif

  {
    int const N = polynomial_degree;
#if TILL_C0_half
    // Compute N nodes, where N is the degree of the polynomial that we want to fit to S.
    int const number_of_nodes = N;
#else
    // Compute N+1 nodes, where N is the degree of the polynomial that we want to fit to S.
    int const number_of_nodes = N + 1;
#endif
    // The x-coordinates of the nodes, in the interval (begin, end].
    std::vector<mpreal> nodes = compute_chebyshev_zeros(number_of_nodes, begin, end);

    // Evaluate S(c) at the mapped nodes.
    std::vector<mpreal> S_values;
    for (mpreal c : nodes)
      S_values.push_back(S(c, root_dummy));

    // Set up the system using Eigen matrices and vectors.
    using MatrixXmp = Eigen::Matrix<mpreal, Eigen::Dynamic, Eigen::Dynamic>;
    using VectorXmp = Eigen::Matrix<mpreal, Eigen::Dynamic, 1>;
    int const number_of_equations = number_of_nodes;
    MatrixXmp A(number_of_equations, number_of_equations);
    VectorXmp b(number_of_equations);

    // Fill vector b with S_values.
    for (int i = 0; i < number_of_equations; ++i)
      b(i) = S_values[i];

    // Fill matrix A with P_j(x_i).
    for (int i = 0; i < number_of_equations; ++i)
      for (int j = 0; j < number_of_equations; ++j)
        A(i, j) = polynomial_basis[j](nodes[i]);

    // Solve the Linear System.
    VectorXmp c = A.colPivHouseholderQr().solve(b);

    // Define our initial polynomial approximation.
    Polynomial P_initial({0});
    for (int i = 0; i < number_of_equations; ++i)
      P_initial += c(i) * polynomial_basis[i];

#if SHOW_SMOOTHING_FUNCTION
    // Add this initial approximation to the graph.
    smoothing_function_graph.add_function([&](double C0) -> double { return P_initial.evaluate(C0).toDouble(); }, color::red);
#endif

#if SHOW_SMOOTHING_FUNCTION_ZOOM
    // The initial approximation is already so accurate that we also show a graph where the difference is magnified.
    cairowindow::QuickGraph smoothing_function_graph_zoom("S(C0) + zoomed in approximation", "C0", "S and S + 100 E",
#if TILL_C0_half
        {0.0, C0_half.toDouble()},
        {-0.05, 0.5},
#else
        {0.0, end},
        {-0.05, 1.05},
#endif
        [&](double C0) -> double { return S(C0, root_dummy).toDouble(); });
    smoothing_function_graph_zoom.add_function([&](double C0) -> double
        {
          mpreal SC0 = S(C0, root_dummy);
          mpreal diff = P_initial.evaluate(C0) - SC0;
          return (SC0 + 100 * diff).toDouble();
        }, color::red);
#if !TILL_C0_half
    math::Line line(begin.toDouble(), 0.0}, math::Direction::up);
    smoothing_function_graph_zoom.add_line(line);
#endif
    for (mpreal const& node : nodes)
    {
      math::Line line({node.toDouble(), 0.0}, math::Direction::up);
      smoothing_function_graph_zoom.add_line(line);
    }
    smoothing_function_graph_zoom.wait_for_keypress();
#endif // SHOW_SMOOTHING_FUNCTION_ZOOM

    {
      // Define the error function E(c) = P(c) - S(c).
      mpreal root_value;
      auto E = [&](mpreal const& c)
        {
          mpreal Sc = S(c, root_value);
          mpreal Pc = P_initial.evaluate(c);
          mpreal cbrtc = cbrt(c);
          if (cbrtc < 1e-6)
            cbrtc = 1e-6;
          mpreal root_S = (Sc - 1) * sqrt_3 - Sc * (cbrtc + 1 / cbrtc);
          mpreal root_P = (Pc - 1) * sqrt_3 - Pc * (cbrtc + 1 / cbrtc);
          mpreal result = (root_P - root_S) / abs(root_S);
          ASSERT(std::isfinite(result.toDouble()));
          return result;
        };

#if SHOW_ERROR_FUNCTION
      // Draw E(C0) on the interval [0, 0.9038].
      error_function_graph.add_function([&E](double C0) -> double { return E(C0).toDouble(); });
#endif

      find_extrema(begin, end, E, nodes, tolerance, extrema
#if SHOW_ERROR_FUNCTION
          , error_function_graph, color_pool
#endif
          );
    }
  }

  do
  {
    std::cout << "Extrema found:\n";
    for (Extremum const& extremum : extrema)
      std::cout << "  " << extremum.x_ << std::endl;

    // Lets define M = N+Z, so that M=N if we didn't find any extra zeroes.
#if TILL_C0_half
    int const M = extrema.size() - 1;
#else
    int const M = extrema.size() - 2;
#endif

    // Calculate the initial M+1, Chebyshev control points.
    std::vector<mpreal> extrema_S_values;
    std::vector<mpreal> extrema_root_values;
    mpreal extrema_root;
    for (auto const& extremum : extrema)
    {
      extrema_S_values.push_back(S(extremum.x_, extrema_root));
      extrema_root_values.push_back(extrema_root);
    }

    // Implement the Remez algorithm.

    // Find the next polynomial P(x) = c1 x + c2 x^2 + ... + c_M x^M
    // such that P(extrema[i].x_) + extrema[i].sign_ * abs(w(extrema[i].x_)) * E = extrema[i].y_,
    // with unknowns c_i and E (a total of M+1 unknowns).
    //
    // Here w(x) is a weight that adjusts the height of the error term in order to make
    // the relative error in the corresponding root the same.
    //
    // Since P(c) approximates S(c), the approximation of the root corresponding to P(c),
    // for c > 0, is:
    //
    //     root_P(c) = (P(c) - 1) * sqrt_3 - P(c) * (cbrt(c) + 1 / cbrt(c))
    //
    // thus the sensitivity of root_P(c) to changes in P(c) is
    //
    //     ∂root_P(c)/∂P = sqrt_3 - (cbrt(c) + 1 / cbrt(c))
    //
    // To ensure that Δroot_P(c)/root(c) is constant across all extrema (for some ΔP(c)),
    // we must set:
    //
    //                root(c)               root(c) cbrt(c)
    //     w(c) = --------------- = ---------------------------------
    //             ∂root_P(c)/∂P     cbrt(c) (sqrt_3 - cbrt(c)) - 1
    //
    // Note that the denominator never reaches zero (its maximum value is -0.25 at c ≈ 0.64951891).

    mpreal root_value;
    auto w = [&](double c) -> double {
      mpreal Sc = S(c, root_value);
      mpreal cbrtc = cbrt(c);
      mpreal root = (c > 0.0) ? (Sc - 1) * sqrt_3 - Sc * (cbrtc + 1 / cbrtc) : -sqrt_3;
      mpreal w = root * cbrtc / (cbrtc * (sqrt_3 - cbrtc) - 1);
      return w.toDouble();
    };
#if SHOW_WEIGH_FUNCTION
    static bool weigh_function_shown = false;
    if (!weigh_function_shown)
    {
      weigh_function_shown = true;
      cairowindow::QuickGraph weigh_function_graph("Weight function", "c", "w(c)", {0.0, 1.0}, w);
      weigh_function_graph.wait_for_keypress();
    }
#endif

    int const number_of_equations = extrema.size();
#if TILL_C0_half
    ASSERT(number_of_equations == M + 1);
#else
    ASSERT(number_of_equations == M + 2);
#endif

    // Set up the system using Eigen matrices and vectors.
    using MatrixXmp = Eigen::Matrix<mpreal, Eigen::Dynamic, Eigen::Dynamic>;
    using VectorXmp = Eigen::Matrix<mpreal, Eigen::Dynamic, 1>;
    MatrixXmp A(number_of_equations, number_of_equations);
    VectorXmp b(number_of_equations);

    // Fill in the matrix A and vector b.
    for (int i = 0; i < number_of_equations; ++i)
    {
      mpreal c = extrema[i].x_;
      mpreal sign = extrema[i].sign_;

      mpreal cbrt_c = cbrt(c);
      mpreal w = extrema_root_values[i] * cbrt_c / (cbrt_c * (sqrt_3 - cbrt_c) - 1);

      // Fill in the polynomial terms.
      for (int j = 0; j < number_of_equations - 1; ++j)
#if TILL_C0_half
        A(i, j) = pow(c, j + 1);
#else
        A(i, j) = pow(c, j);
#endif

      // Fill in the term for E.
      A(i, number_of_equations - 1) = sign * abs(w);

      // Fill in b_i.
      b(i) = extrema_S_values[i];
    }

    // Solve the equation A⋅x = b.
    VectorXmp x = A.colPivHouseholderQr().solve(b);

    std::vector<mpreal> coefficients(M + 1);
#if TILL_C0_half
    // Extract coefficients c_1 to c_{M}.
    for (int i = 1; i <= M; ++i)
      coefficients[i] = x(i - 1);
#else
    // Extract coefficients c_0 to c_{M}.
    for (int i = 0; i <= M; ++i)
      coefficients[i] = x(i);
#endif
    // Extract E.
    mpreal Err = x(number_of_equations - 1);

    Polynomial P{coefficients};
    std::cout << "Err = " << Err << std::endl;
    std::cout << "Polynomial:\n";
    for (int i = 0; i <= M; ++i)
      std::cout << "  " << coefficients[i] << " * x^" << i << std::endl;

#if SHOW_ERROR_FUNCTION
    error_function_graph.add_function(
        [&](double c) -> double
        {
          mpreal Sc = S(c, root_value);
          mpreal cbrtc = cbrt(c);
          mpreal droot_dP = sqrt_3 - (cbrtc + 1 / cbrtc);       // ∂root_P(c)/∂P = sqrt_3 - (cbrt(c) + 1 / cbrt(c))
          mpreal delta_P = abs(Err * w(c));
          mpreal delta_root = droot_dP * delta_P;
          mpreal rel_err = abs(delta_root / root_value);
          return rel_err.toDouble();
        }, color::red);
    error_function_graph.add_function(
        [&](double c) -> double
        {
          mpreal Sc = S(c, root_value);
          mpreal cbrtc = cbrt(c);
          mpreal droot_dP = sqrt_3 - (cbrtc + 1 / cbrtc);       // ∂root_P(c)/∂P = sqrt_3 - (cbrt(c) + 1 / cbrt(c))
          mpreal delta_P = abs(Err * w(c));
          mpreal delta_root = droot_dP * delta_P;
          mpreal rel_err = -abs(delta_root / root_value);
          return rel_err.toDouble();
        }, color::red);

    error_function_graph.add_function([&](double c) -> double {
        mpreal root_value;
        mpreal Sc = S(c, root_value);
        mpreal cbrtc = cbrt(c);
        return ((-(cbrtc + 1/cbrtc) - root_value) / abs(root_value)).toDouble();
      }, color::blue);
#endif

#if SHOW_SMOOTHING_FUNCTION
    smoothing_function_graph.add_function([&P](double C0) -> double { return P.evaluate(C0).toDouble(); }, color::green);
#endif

    // Define the error function E(c) = S(c) - P(c).
    auto E = [&](mpreal const& c)
      {
        mpreal Sc = S(c, root_value);
        mpreal Pc = P.evaluate(c);
        mpreal cbrtc = cbrt(c);
        if (cbrtc < 1e-6)
          cbrtc = 1e-6;
        mpreal root_S = (Sc - 1) * sqrt_3 - Sc * (cbrtc + 1 / cbrtc);
        mpreal root_P = (Pc - 1) * sqrt_3 - Pc * (cbrtc + 1 / cbrtc);
        mpreal result = (root_P - root_S) / abs(root_S);
        ASSERT(std::isfinite(result.toDouble()));
        return result;
      };

#if SHOW_ERROR_FUNCTION
    // Draw E(C0) on the interval [0, 0.9038].
    error_function_graph.add_function([&E](double C0) -> double { return E(C0).toDouble(); }, color::blue);
#endif

    // Find the roots of E.
    std::vector<mpreal> roots;

    // Run over the M intervals between the previously found extrema and store the value of the new E.
    std::vector<math::Point> points;
    for (int i  = 0; i < extrema.size() - 1; ++i)
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
      double x = extrema[extrema.size() - 1].x_.toDouble();
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
#if SHOW_ERROR_FUNCTION
        error_function_graph.add_point({zeroes.back().toDouble(), 0.0}, zero_style);
#endif
      }
    }

    // Find and store all extrema.
    extrema.clear();
    find_extrema(begin, end, E, zeroes, tolerance, extrema
#if SHOW_ERROR_FUNCTION
        , error_function_graph, color_pool
#endif
        );

    // Give a chance to view the graphs before clearing.
    std::cin.get();
#if SHOW_ERROR_FUNCTION
    error_function_graph.clear();
#endif
  }
  while (true);

#if SHOW_ERROR_FUNCTION
  error_function_graph.wait_for_keypress();
#elif SHOW_SMOOTHING_FUNCTION
  smoothing_function_graph.wait_for_keypress();
#endif
}
