#include "sys.h"
#include "utils/square.h"
#include "cairowindow/QuickGraph.h"
#include "mpreal/mpreal.h"
#include <array>
#include <vector>
#include <cmath>
#include <iostream>
#include "debug.h"

#define MINUS_THREE 0

// Let P(x) be the cubic:
//
//   P(x) = C0 + C1 x + x^3
//
// Let r(C0, C1) be the smallest (negative) root of P, where C0 >= 0.
//
// Let
//
//              r₀ - r   r - r₀
//   relerr₀ = ------- = ------
//               |r|       r
//
// be the relative error of r₀ with respect to r.
// Note that |r| = -r because r <= 0.
//
// Hence r₀ = r (1 - relerr₀)
//
// Let N be the number of Halley iterations required to
// find r, from a starting point of r₀, with a precision
// of a double.
//
// Then we are interested in the intervals of relerr₀
// within which N has the given value 1, 2, 3 or 4, as
// function of C0.

using namespace mpfr;

#ifdef CWDEBUG
// This can be used from within gdb to print an mpreal value.
void print(mpreal const& val)
{
  std::cout << val << std::endl;
}
#endif

constexpr mp_prec_t precision_in_bits = 128;

// Function to calculate the real root of P(u) = C0 + C1 u + u^3 for given C0 and C1.
mpreal calculate_root(mpreal const& C0, mpreal const& C1)
{
  mpreal u = C1 == mpreal{3} ? mpreal{0} : -sqrt(3);
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

  return u;
}

// Function that returns a double with the root of P(u) = C0 + C1 u + u^3 with maximum possible precision.
double exact_root(double C0, double C1)
{
  return calculate_root(C0, C1).toDouble();
}

// Returns the minimum number of times a Halley iteration must be performed on
// initial_guess in order to reach the root of C0 + C1 x + x^3 as close as possible
// (meaning that additional iterations do not improve the accuracy).
int iterations(double C0, double C1, double initial_guess)
{
  double correct_root = exact_root(C0, C1);

  constexpr int large_count = 5;        // Five means that we definitely need more iterations.
  double u = initial_guess;
  int number_of_iterations = -2;        // The number of iterations of the previous count.
  int prev_count = large_count;         // The number of std::nextafter steps that u is away from correct_root.
  do
  {
    // Count how many iterations have been performed.
    ++number_of_iterations;

    Dout(dc::notice, "After " << (number_of_iterations + 1) << " iterations, u = " << std::setprecision(18) << u);

    // Calculate the number of times `nextafter` must be called to reach `correct_root`.
    int count;
    if (u == correct_root)
      count = 0;
    else
    {
      double nup = u, num = u;
      count = 1;
      for (; count <= 4; ++count)
      {
        nup = std::nextafter(nup, std::numeric_limits<double>::infinity());
        num = std::nextafter(num, -std::numeric_limits<double>::infinity());
        if (nup == correct_root || num == correct_root)
          break;
      }
    }

    // Stop if the count is greater than or equal to the previous count (but less than large_count).
    if (prev_count != large_count && count >= prev_count - 1)
      break;

    // Remember the distance of the previous iteration.
    prev_count = count;

    // Perform one Halley iteration.
    double Q_u = C0 + u * (utils::square(u) + C1);
    double half_Qpp_u = 3.0 * u;
    double Qp_u = half_Qpp_u * u + C1;
    u += -Q_u * Qp_u / (utils::square(Qp_u) - Q_u * half_Qpp_u);
  }
  while (true);

  return number_of_iterations;
}

class MaxRelerr
{
 private:
  int const N_;
  double const tolerance_;

 public:
  MaxRelerr(int N, double tolerance) : N_(N), tolerance_(tolerance) { }

  int iterations(double relerr0, double C0, double C1, double r) const;
  std::array<double, 2> max_relerr(double C0, double C1) const;
};

int MaxRelerr::iterations(double relerr0, double C0, double C1, double r) const
{
  DoutEntering(dc::notice, "MaxRelerr::iterations(" << relerr0 << ", " << C0 << ", " << C1 << ", " << r << ") [with N_ = " << N_ << "]");
  double r0 = r * (1.0 - relerr0);
  return ::iterations(C0, C1, r0);
}

std::array<double, 2> MaxRelerr::max_relerr(double C0, double C1) const
{
  DoutEntering(dc::notice, "MaxRelerr::max_relerr(" << C0 << ", " << C1 << ") [with N_ = " << N_ << ", tolerance_ = " << tolerance_ << "]");

  std::array<double, 2> result;

  double r = exact_root(C0, C1);
  ASSERT(r < 0.0);

  auto r0_to_relerr = [=](double r0) -> double { return (r - r0) / r; };

  // Index 0: negative relerr, 1: positive relerr.
  std::array<double, 2> relerr_min = {{ 0.0, 0.0 }};
  std::array<double, 2> relerr_max = {{ r0_to_relerr(-1e10), r0_to_relerr(1e10) }};
  if (C1 == -3.0)
  {
    // We want to set relerr_max such that one of those values is slightly less than -1 (the maximum of the cubic).
    relerr_max[1] = r0_to_relerr(-1.00001);
  }
  ASSERT(relerr_max[0] < 0.0);
  ASSERT(relerr_max[1] > 0.0);

  for (int neg_pos = 0; neg_pos <= 1; ++neg_pos)
  {
    // We're looking for the largest value such that N_min == N_.
    int N_min = iterations(relerr_min[neg_pos], C0, C1, r);
    // We're looking for the smallest value, larger than N_min, such that N_max == N_ + 1.
    int N_max = iterations(relerr_max[neg_pos], C0, C1, r);

    // Check if the solution is bracketed.
    if (!(N_min <= N_ && N_ < N_max))
      throw std::runtime_error("Target N_ is not bracketed within N_min and N_max.");

    while (std::abs(relerr_max[neg_pos] - relerr_min[neg_pos]) > tolerance_)
    {
      double relerr_mid = 0.5 * (relerr_min[neg_pos] + relerr_max[neg_pos]);
      int N_mid = iterations(relerr_mid, C0, C1, r);

      if (N_mid <= N_)
        relerr_min[neg_pos] = relerr_mid;
      else
        relerr_max[neg_pos] = relerr_mid;
    }

    // Return the largest value for which the number of iterations is still (just) N_.
    Dout(dc::notice, "C0 = " << C0 << "; relerr_min[" << neg_pos << "] = " << relerr_min[neg_pos]);
    result[neg_pos] = relerr_min[neg_pos];
  }

  return result;
}

// Function to calculate the approximate root0 = -C0/3 + C0^3/81.
double approximate_root(double C0)
{
  return -C0 / 3.0 + C0 * C0 * C0 / 81.0;
}

int main()
{
  Debug(NAMESPACE_DEBUG::init());
  Dout(dc::notice, "Entering main()");

  // Set the default floating-point precision.
  mpfr_set_default_prec(precision_in_bits);

  // Print everything with 38 digits precision.
  std::streamsize const old_precision = std::cout.precision(38);
  std::ios::fmtflags const old_flags = std::cout.setf(std::ios::fixed);

  using namespace cairowindow;

#if MINUS_THREE
  Range xrange{0.0, 6.82};      // Larger than 6.82 is approached with -(cbrt(C0) + 1 / cbrt(C0)).
  Range yrange{-0.02, 0.02};
#else
  Range xrange{0.0, 6.82};
  Range yrange{-0.2, 0.2};
#endif
  QuickGraph graph("Maximum relative error in initial guess still resulting in 2 iterations", "C0", "relerr₀", xrange, yrange);
#if MINUS_THREE
  graph.add_function([](double C0) -> double {
      if (C0 <= 2.0)
        return 0.0125 + 0.00078 * C0;
      //return (0.0122 + 0.00078 * 5.0) + 0.0003 * (C0 - 5.0);
      return 0.0145 + 0.0003 * C0;
    });
  graph.add_function([](double C0) -> double {
      if (C0 <= 2.0)
        return -(0.0128 + 0.00078 * C0);
      return -0.0149 - 0.0003 * C0;
    });
#else
  graph.add_function([](double C0) -> double { return C0 < 1e-9 ? 1e6 : 0.075 * std::pow(C0, -0.75); });
  graph.add_function([](double C0) -> double { return C0 < 1e-9 ? -1e6 : -0.075 * std::pow(C0, -0.75); });

  // Add function to plot the relative error for the approximation of root0 near zero.
  graph.add_function(
    [](double C0) -> double {
      double root = exact_root(C0, 3.0);
      double root0 = approximate_root(C0);
      return (root0 - root) / std::abs(root);
    }
  );
#endif

  draw::PointStyle point_style0({.color_index = 0, .filled_shape = 15});
  draw::PointStyle point_style1({.color_index = 1, .filled_shape = 15});

  {
    Debug(dc::notice.off());
    MaxRelerr max_relerr(2, 1e-9);
    double worst_relerr = 100.0, worst_C0;
    double const step = xrange.size() / (1 << 16);
    for (double C0 = step; C0 < xrange.max(); C0 += step)
    {
      std::array<double, 2> relerr =
#if MINUS_THREE
        max_relerr.max_relerr(C0, -3.0);
#else
        max_relerr.max_relerr(C0, 3.0);
#endif
      graph.add_point({C0, relerr[0]}, point_style0);
      graph.add_point({C0, relerr[1]}, point_style1);
      if (std::abs(relerr[0]) < std::abs(worst_relerr))
      {
        worst_relerr = relerr[0];
        worst_C0 = C0;
      }
      if (std::abs(relerr[1]) < std::abs(worst_relerr))
      {
        worst_relerr = relerr[1];
        worst_C0 = C0;
      }
    }
    Debug(dc::notice.on());
    Dout(dc::notice, "Worst C0 = " << std::setprecision(18) << worst_C0 << " with a min. relerr of " << worst_relerr << ".");
  }

  graph.wait_for_keypress();
}
