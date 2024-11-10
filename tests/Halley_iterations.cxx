#include "sys.h"
#include "cairowindow/QuickGraph.h"
#include "math/CubicPolynomial.h"
#include "debug.h"

// Let P(x) be the cubic:
//
//   P(x) = C0 - 3x + x^3
//
// Let r(C0) be the smallest (negative) root of P, where C0 >= 0.
//
// Let
//
//             |r - r₀|
//   relerr₀ = --------
//                r
//
// be the relative error of r₀ with respect to r.
// Hence r₀ = r * (1 ± relerr₀)
//
// Let N be the number of Halley iterations required to
// find r, from a starting point of r₀, with a precision
// of a double.
//
// Then we are interested in the intervals of relerr₀
// within which N has the given value 1, 2, 3 or 4, as
// function of C0.

class MaxRelerr
{
 private:
  int N_;
  double tolerance_;

 public:
  MaxRelerr(int N, double tolerance) : N_(N), tolerance_(tolerance) { }

  double operator()(double const C0) const;
  int iterations(double relerr, double C0, double r) const;
};

int MaxRelerr::iterations(double relerr, double C0, double r) const
{
  DoutEntering(dc::notice, "MaxRelerr::iterations(" << relerr << ", " << C0 << ", " << r << ") [with N_ = " << N_ << "]");

  std::array<double, 2> r0 = {{ r * (1 - relerr), r * (1 + relerr) }};
  std::array<int, 2> N;

  for (int minus_plus = 0; minus_plus <= 1; ++minus_plus)
  {
    double const C1 = -3.0;
    double u = r0[minus_plus];
    int iterations;
    ASSERT(u < -1.0);
    Dout(dc::notice, "Initial u = " << std::setprecision(18) << u);

#define HALLEY_ITERATIONS_TEST
#   include "math/CubicPolynomial_get_roots.cpp"
#undef HALLEY_ITERATIONS_TEST

    N[minus_plus] = iterations;
  }

  // Return the worst cases.
  int max_iterations = std::max(N[0], N[1]);
  Dout(dc::notice, "Number of iterations: " << max_iterations);

  return max_iterations;
}

double MaxRelerr::operator()(double const C0) const
{
  DoutEntering(dc::notice, "MaxRelerr::operator()(" << C0 << ") [with N_ = " << N_ << ", tolerance_ = " << tolerance_ << "]");

  math::CubicPolynomial cubic(C0, -3.0, 0.0, 1.0);

  std::array<double, 3> roots;
//  Debug(dc::notice.off());
  cubic.get_roots(roots);
//  Debug(dc::notice.on());

  double r = roots[0];

  if (C0 == 0.0)                // In this case roots[0] is zero, but we need r to be the most negative root.
    r = -std::sqrt(3);

  ASSERT(r <= -std::sqrt(3.0));

  double relerr_min = 0.0;
  // r0 = r * (1 ± relerr).
  // We want to set relerr_max such that one of those values is slightly less than -1 (the maximum of the cubic).
  // Thus -1.00001 = r * (1 ± relerr_max) --> relerr_max = ±(-1.00001 / r - 1), picking the sign that makes it positive of course.
  double relerr_max = std::abs(-1.00001 / r - 1.0);

  int N_min = iterations(relerr_min, C0, r);    // We're looking for the largest value such that N_min == N_.
  int N_max = iterations(relerr_max, C0, r);    // We're looking for the smallest value, larger than N_min, such that N_max == N_ + 1.

  // Check if the solution is bracketed.
  if (!(N_min <= N_ && N_max > N_))
    throw std::runtime_error("target is not bracketed within x_min and x_max.");

  while (relerr_max - relerr_min > tolerance_)
  {
    double relerr_mid = 0.5 * (relerr_min + relerr_max);
    int N_mid = iterations(relerr_mid, C0, r);

    if (N_mid <= N_)
      relerr_min = relerr_mid;
    else
      relerr_max = relerr_mid;
  }

  // Return the largest value for which the number of iterations is still (just) N_.
  return relerr_min;
}

int main()
{
  Debug(NAMESPACE_DEBUG::init());
  Dout(dc::notice, "Entering main()");

  using namespace cairowindow;

  QuickGraph graph("Maximum relative error in initial guess to get 2 iterations",
      "C0", "relerr₀", {0.0, 6.82});

//  Debug(dc::notice.off());
  {
    MaxRelerr max_relerr2(2, 1e-9);
//    MaxRelerr max_relerr3(3, 1e-9);
//    graph.add_function(max_relerr3, color::red);
    graph.add_function(max_relerr2, color::blue);
  }
//  Debug(dc::notice.on());

  graph.wait_for_keypress();
}
