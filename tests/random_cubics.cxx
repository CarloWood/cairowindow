#include "sys.h"
#include "math/CubicPolynomial.h"
#include "math/AnalyzedCubic.h"
#include "utils/almost_equal.h"
#include <random>
#include <chrono>
#include <sstream>
#include <array>
#include "debug.h"

constexpr int number_of_cubics = 10000;

void transform(std::array<double, 4>&c, double x_scale, double y_scale, double x_shift, double y_shift)
{
  // First apply the scaling:
  for (int i = 0; i < c.size(); ++i)
    c[i] *= y_scale * std::pow(x_scale, -i);
  // Then the offset.
  c[0] += y_shift + (-c[1] + (c[2] - c[3] * x_shift) * x_shift) * x_shift;
  c[1] += (-2.0 * c[2] + 3.0 * c[3] * x_shift) * x_shift;
  c[2] += -3.0 * c[3] * x_shift;
}

int main(int argc, char* argv[])
{
  Debug(NAMESPACE_DEBUG::init());

#if 0
  // 17.296348 x^2 + 99.720970 x + 10.350319
  math::CubicPolynomial p(10.350319, 99.720970, 17.296348, 1.0);
  std::array<double, 3> roots2;
  int n = p.get_roots(roots2);
  Dout(dc::notice, "The roots of " << p << " are:");
  for (int i = 0; i < n; ++i)
    Dout(dc::notice, std::setprecision(12) << roots2[i]);
  return 0;
#endif

  // Handle random seed.
  unsigned int seed;
  if (argc > 1)
  {
    std::istringstream ss(argv[1]);
    ss >> seed;
  }
  else
    seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::cout << "Seed: " << seed << '\n';

  // Use the mersenne twister engine.
  std::mt19937 engine(seed);

  // Keep all numbers such that their cube still fits in a double.
  double const max_value = std::cbrt(std::numeric_limits<double>::max());
  double const min_value = std::cbrt(std::numeric_limits<double>::min());
  Dout(dc::notice, "min_value = " << min_value << "; max_value = " << max_value);

  // We'll deal with all numbers logarithmicly, because of their huge spread.
  double const log_max_value = std::log(max_value);
  double const log_min_value = std::log(min_value);
  Dout(dc::notice, "log(min_value) = " << log_min_value << "; log(max_value) = " << log_max_value);

  // Truncate these exponents towards zero, with a tiny little bit more leeway.
  double const log_min_value_trunc = std::trunc(log_min_value) + 1;
  double const log_max_value_trunc = std::trunc(log_max_value) - 1;
  Dout(dc::notice, "log(min_value) = " << log_min_value_trunc << "; log(max_value) = " << log_max_value_trunc);

  std::normal_distribution normal_dist(0.0, log_max_value_trunc / 5.0); // 1 in 3000000 will have the most extreme value.
  auto random_value = [&](){
    double exponent;
    do { exponent = normal_dist(engine); } while (exponent < log_min_value_trunc || log_max_value_trunc > log_max_value_trunc);
    return std::exp(exponent);
  };

  // Generate random cubics.
  std::cout << "Generating random cubic polynomials..." << std::endl;
  std::vector<math::CubicPolynomial> cubics;
  for (int i = 0; i < number_of_cubics; ++i)
  {
    // Generate the position of the inflection point.
    double Ix = random_value();
    double Iy = random_value();
    // Generate the position of the local minimum.
    double Ex = random_value();
    double Ey = random_value();
    // Create an array with coefficients for the cubic x * (x^2 - 3) = -3 x + x^3 if Ey < Iy, or +3 x + x^3 if Ey >= Iy.
    std::array<double, 4> c = { 0, std::copysign(3.0, Ey - Iy), 0, 1 };
    // Scale the cubic such that the minimum of -3 x + x^3 (at x = 1) goes to Ex - Ix:
    double x_scale = Ex - Ix;
    // Scale the cubic such that the minimum, at y = -2, goes to Ey - Iy.
    double y_scale = 0.5 * (Iy - Ey);
    // Shift the cubic such that the inflection point (still at the origin) moves to (Ix, Iy).
    double x_shift = Ix;
    double y_shift = Iy;

    // Apply the linear transformation.
    transform(c, x_scale, y_scale, x_shift, y_shift);

    cubics.emplace_back(c[0], c[1], c[2], c[3]);
  }

  std::array<double, 3> roots;
  std::cout << "Calculating roots..." << std::endl;
  int max_iterations = 0;
//  Debug(libcw_do.off());
  unsigned long total_iterations = 0;
  for (math::CubicPolynomial const& cubic : cubics)
  {
    int iterations;
    int n = cubic.get_roots(roots, iterations);
    total_iterations += iterations;

    if (iterations == 12 || (iterations < 100 && iterations > max_iterations))
    {
      max_iterations = iterations;
//      Debug(libcw_do.on());

      bool inflection_point_y_larger_than_zero =
        13.5 * cubic[0] + cubic[2] * (utils::square(cubic[2] / cubic[3]) - 4.5 * (cubic[1] / cubic[3])) > 0.0;

      math::AnalyzedCubic acubic;
      // Calculate the minimum (-1) if the inflection point lays above the x-axis, and the maximum otherwise.
      // This way acubic.get_extreme() (E above) becomes the extreme that is the closest to the x-axis.
      acubic.initialize(cubic, inflection_point_y_larger_than_zero ? -1 : 1);

      // Obtain the calculated inflection point.
      double const inflection_point_x = acubic.inflection_point();

      // Remember if we have local extrema or not.
      bool const cubic_has_local_extrema = acubic.has_extrema();

      // Avoid the local extrema and the inflection point because the derivative might be zero there too.
      double x = cubic_has_local_extrema ?
          3 * inflection_point_x - 2 * acubic.get_extreme() :
          inflection_point_x + (((cubic[3] > 0.0) == inflection_point_y_larger_than_zero) ? -1.0 : 1.0);

      Dout(dc::notice|continued_cf, cubic << " roots: ");
      for (int r = 0; r < n; ++r)
        Dout(dc::continued, roots[r] << ", ");
      Dout(dc::continued, "inflection point: " << inflection_point_x);
      if (cubic_has_local_extrema)
        Dout(dc::continued, ", extreme: " << acubic.get_extreme());
      Dout(dc::finish, ", starting point: " << x);
//      Debug(libcw_do.off());
    }
  }
//  Debug(libcw_do.on());
  Dout(dc::notice, "Average number of iterations per cubic: " << (total_iterations / number_of_cubics));
}
