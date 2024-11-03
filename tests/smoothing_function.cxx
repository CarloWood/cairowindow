#include "sys.h"
#include "mpreal/mpreal.h"
#include "cairowindow/QuickGraph.h"
#include "utils/ColorPool.h"
#include <glpk.h>
#include <boost/math/tools/polynomial.hpp>
#include <Eigen/Dense>
#include <stdexcept>
#include <iostream>
#include <cassert>
#include "debug.h"

using namespace mpfr;

void print(mpreal const& val)
{
  std::cout << val << std::endl;
}

// Define a few constants.
constexpr mp_prec_t precision_in_bits = 128;
constexpr double tolerance_dbl = 5e-39;
constexpr int polynomial_degree = 4;

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
mpfr::mpreal S(mpfr::mpreal const& C0)
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
  mpreal y_min = S(x_min);
  mpreal y_max = S(x_max);

  // Check if the solution is bracketed.
  if (y_min > target || y_max < target)
    throw std::runtime_error("target is not bracketed within x_min and x_max.");

  mpreal y_mid = 0;
  while (x_max - x_min > tolerance)
  {
    mpreal x_mid = 0.5 * (x_min + x_max);
    mpreal prev_y_mid = y_mid;
    y_mid = S(x_mid);

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
Extremum find_extremum(std::function<mpreal(mpfr::mpreal const&)> const& g, mpreal x_min, mpreal x_max, mpreal const& tolerance)
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
  result.sign_ = gc < 0 ? -1 : 1;
  Dout(dc::continued, "{sign_ = " << std::boolalpha << result.sign_ << "} = ");

  // Implement the Golden-section search algorithm, see https://en.wikipedia.org/wiki/Golden-section_search
  // This algorithm searches for a maximum, therefore, if sign_ is -1, we use -g(x) instead of g(x).

  // Specify the function to be minimized, f(x):
  auto f = [sign = result.sign_, &g](mpreal const& x){ return sign * g(x); };

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

  // If this fails then sign_ has the wrong value!
  // The reason for this is that we are always using at least one zero for points 0 and 3
  // and if the other point isn't a zero then the function is monotonically rising or descending.
  ASSERT(F[1] > 0 && F[2] > 0);

  std::cout << "Entering loop.\n";
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

    for (int i = 0; i < X.size(); ++i)
      std::cout << "X[" << i << "] = " << X[i] << std::endl;
  }
  while (X[3] - X[0] > x_epsilon || abs(F[0] - F[3]) > tolerance);

  result.x_ = (X[0] + X[3]) / 2;
  result.y_ = result.sign_ * (F[0] + F[3]) / 2;
  Dout(dc::finish, result.x_ << ", " << result.y_);
  return result;
}

int main()
{
  Debug(NAMESPACE_DEBUG::init());

  // Set the default floating-point precision.
  mpfr_set_default_prec(precision_in_bits);

  // Print everything with 30 digits precision.
  std::streamsize const old_precision = std::cout.precision(38);
  std::ios::fmtflags const old_flags = std::cout.setf(std::ios::fixed);

  // Loop over all values of C0 between 0 and 1 in steps of 1/128th.
  mpreal C0(0);
  mpreal const step(1.0 / 128.0);

  // Find the value of C0 for which S(C0) returns 0.5.
  long double const C0_half_approximation = 0.90375741845959156233304814223072905692L;
  mpreal const C0_half = bracket_S(std::nextafter(C0_half_approximation, 0.0L), std::nextafter(C0_half_approximation, 1.0L), 0.5, tolerance_dbl);

  std::cout << "C0_half = " << C0_half << std::endl;

  // Draw S(C0) on the interval [0, 0.9038].
  cairowindow::QuickGraph graph1("The smoothing function S(C0)", "C0", "S",
      {0.0, C0_half.toDouble()}, [](double C0) -> double { return S(C0).toDouble(); });

  // Compute N+1 Chebyshev zeroes, where N is the degree of the polynomial that we want to fit to S.
  std::vector<mpreal> zeroes = compute_chebyshev_zeros(polynomial_degree + 1, 0, C0_half);

  // Evaluate S(c) at the mapped nodes.
  std::vector<mpreal> S_values;
  for (mpreal c : zeroes)
    S_values.push_back(S(c));

  auto coefficients = compute_newton_coefficients(zeroes, S_values);
  auto P = compute_newton_polynomial(zeroes, coefficients);

  namespace color = cairowindow::color;
  graph1.add_function([&P](double C0) -> double { return P.evaluate(C0).toDouble(); }, color::red);

  // Define the error function E(c) = S(c) - P(c).
  auto E = [&](mpreal const& c) { return S(c) - P.evaluate(c); };

  // Draw E(C0) on the interval [0, 0.9038].
  cairowindow::QuickGraph graph2("The error function E(C0)", "C0", "E",
      {0.0, C0_half.toDouble()}, [&E](double C0) -> double { return E(C0).toDouble(); });

  // Initialize vectors to store extrema.
  std::vector<Extremum> extrema;

  // Set the tolerance for the search.
  mpreal const tolerance{tolerance_dbl};

  // Define a zero.
  mpreal const zero(0);

  utils::ColorPool<32> color_pool;
  cairowindow::draw::PointStyle min_style({.color_index = color_pool.get_and_use_color(), .filled_shape = 1});
  cairowindow::draw::PointStyle max_style({.color_index = color_pool.get_and_use_color(), .filled_shape = 7});

  // For each interval between zeroes[i] and zeroes[i+1].
  for (size_t i = 0; i <= zeroes.size(); ++i)
  {
    mpreal c_min = i == 0 ? zero : zeroes[i - 1];
    mpreal c_max = i == zeroes.size() ? C0_half : zeroes[i];

    // Find the extremum in the interval.
    graph2.add_point({c_min.toDouble(), E(c_min).toDouble()}, min_style);
    graph2.add_point({c_max.toDouble(), E(c_max).toDouble()}, max_style);
    Extremum extremum = find_extremum(E, c_min, c_max, tolerance);

    extrema.push_back(extremum);
  }

  // Show the extrema.
  for (auto const& extremum : extrema)
    graph2.add_point({extremum.x_.toDouble(), extremum.y_.toDouble()});

  // Give a chance to view the graphs before exiting.
  graph2.wait_for_keypress();

  // Calculate the initial N+2 Chebyshev control points.
  std::vector<mpreal> extrema_S_values;
  for (auto const& extremum : extrema)
    extrema_S_values.push_back(S(extremum.x_));

  // Find the next polynomial P(x) = c0 + c1 x + c2 x^2 + ... + c_N x^N
  // such that P(extrema[i].x_) + extrema[i].sign_ * E = extrema[i].y_,
  // with unknowns c_i and E.
  int const num_equations = polynomial_degree + 2;      // Number of equations and unknowns.

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

    // Fill in the polynomial terms.
    for (int j = 0; j <= polynomial_degree; ++j)
      A(i, j) = pow(x_i, j);

    // Fill in the term for E.
    A(i, polynomial_degree + 1) = sign;

    // Fill in b_i.
    b(i) = extrema_S_values[i];
  }

  // Solve the equation A⋅x = b.
  VectorXmp x = A.colPivHouseholderQr().solve(b);

  // Extract coefficients c0 to cN and E.
  std::vector<mpreal> c_coefficients(polynomial_degree + 1);
  for (int i = 0; i <= polynomial_degree; ++i)
    c_coefficients[i] = x(i);
  mpreal Err = x(polynomial_degree + 1);

  std::cout << "Err = " << Err << std::endl;

  // Temporarily convert the coefficients of the polynomial to double.
  std::vector<double> c_coefficients_double;
  for (mpreal const& coeff : c_coefficients)
    c_coefficients_double.push_back(coeff.toDouble());

  //  Eigen::Matrix<double, 5, 5> C;
  //  C <<  0, 0, 0, 0, -a[0] / a[5],
  //        1, 0, 0, 0, -a[1] / a[5],
  //        0, 1, 0, 0, -a[2] / a[5],
  //        0, 0, 1, 0, -a[3] / a[5],
  //        0, 0, 0, 1, -a[4] / a[5];

  // Construct the companion matrix for this polynomial.
  int degree = c_coefficients_double.size() - 1;
  Eigen::MatrixXd companion_matrix = Eigen::MatrixXd::Zero(degree, degree);

  // Fill the last column with the negative coefficients.
  double const leading_coefficient = c_coefficients_double.back();
  for (int i = 0; i < degree; ++i)
    companion_matrix(i, degree - 1) = -c_coefficients_double[i] / leading_coefficient;

  // Fill the sub-diagonal with ones.
  for (int i = 1; i < degree; ++i)
    companion_matrix(i, i - 1) = 1.0;

  // Calculate the roots.
  Eigen::EigenSolver<Eigen::MatrixXd> solver(companion_matrix);
  Eigen::VectorXcd eigenvalues = solver.eigenvalues();
  std::vector<double> roots;

  for (int i = 0; i < eigenvalues.size(); ++i)
  {
    std::complex<double> root = eigenvalues[i];
    assert(std::abs(root.imag()) < 1e-8);
    roots.push_back(root.real());
//    assert(0.0 <= roots.back() && roots.back() <= 0.903757418459591562);
  }

  // Output the roots.
  std::cout << "Roots within the interval [0, C0_half]:\n";
  for (double root : roots)
    std::cout << root << std::endl;
}
