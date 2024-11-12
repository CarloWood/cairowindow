#include "sys.h"
#include "mpreal/mpreal.h"
#include "cairowindow/QuickGraph.h"
#include "cairowindow/Line.h"
#include "debug.h"
#include <cmath>

// This file analyses P(u) = C0 + 3u + u^3 and its roots.

using namespace mpfr;

void print(mpreal const& val)
{
  std::cout << val << std::endl;
}

// Define a few constants.
constexpr mp_prec_t precision_in_bits = 128;

// Function to calculate the real root of P(u) = C0 + 3u + u^3 for a given C0.
mpreal calculate_root(mpreal const& C0)
{
  static mpreal const C1(3);

  mpreal u = 0;
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

// Function to calculate the approximate root0 = -C0/3 + C0^3/81.
mpreal approximate_root(mpreal const& C0)
{
  return -C0 / 3 + C0 * C0 * C0 / 81;
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

  // Set up graph for plotting the root as a function of C0 in the range [0, 10].
  cairowindow::QuickGraph root_graph("Root of P(u) = C0 + 3u + u^3", "C0", "Root u", {0.0, 10.0});

  // Add function to plot the root.
  root_graph.add_function(
    [](double C0) -> double {
      mpreal C0_mpreal(C0);
      return calculate_root(C0_mpreal).toDouble();
    }
  );

  // Set up graph for plotting the relative error of the approximation.
  cairowindow::QuickGraph error_graph("Relative error of approximation", "C0", "(root0 - root) / |root|", {0.1, 1.0});

  // Add function to plot the relative error.
  error_graph.add_function(
    [](double C0) -> double {
      mpreal C0_mpreal(C0);
      mpreal root = calculate_root(C0_mpreal);
      mpreal root0 = approximate_root(C0_mpreal);
      mpreal relative_error = (root0 - root) / abs(root);
      Dout(dc::notice, "lambda(" << C0 << ") = " << relative_error);
      return relative_error.toDouble();
    }
  );

  Dout(dc::notice, "Press any key to close the graphs.");
  root_graph.wait_for_keypress();
  error_graph.wait_for_keypress();
}
