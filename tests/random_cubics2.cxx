#include "sys.h"
#include "cairowindow/QuickGraph.h"
#include "math/CubicPolynomial.h"
#include "math/bracket_zero.h"
#include <random>
#include <chrono>
#include <sstream>
#include <array>
#include "debug.h"
#include "cwds/Restart.h"

constexpr int number_of_cubics = 1000000;

int get_roots(math::CubicPolynomial& cubic, std::array<double, 3>& roots_out, int& iterations)
{
  DoutEntering(dc::notice, "get_roots() for " << cubic);

  // Define a coefficients_ for use by math/CubicPolynomial_get_roots.cpp.
  std::array<double, 4>& coefficients_ = reinterpret_cast<std::array<double, 4>&>(cubic[0]);

  auto evaluate = [&](double x){ return coefficients_[0] + (coefficients_[1] + (coefficients_[2] + coefficients_[3] * x) * x) * x; };
  auto derivative = [&](double x){ return coefficients_[1] + (2.0 * coefficients_[2] + 3.0 * coefficients_[3] * x) * x; };
  auto half_second_derivative = [&](double x){ return coefficients_[2] + 3.0 * coefficients_[3] * x; };

  using math::QuadraticPolynomial;
  // Include the body of the function.
# define GETROOTS_ASSIGN_ITERATIONS
# include "math/CubicPolynomial_get_roots.cpp"
# undef GETROOTS_ASSIGN_ITERATIONS
}

void sanity_check(math::CubicPolynomial const& cubic, double const root)
{
  // Course check.
  double fr = cubic(root);

  double root_small = root < 0 ? root * 1.000001 : root * 0.999999;
  double root_large = root < 0 ? root * 0.999999 : root * 1.000001;

  double frm = cubic(root_small);
  double frp = cubic(root_large);

  ASSERT(std::abs(fr) < std::abs(frm) && std::abs(fr) < std::abs(frp));

  try
  {
    double real_root = math::bracket_zero(root_small, root_large, [&](double r){ return cubic(r); });
    Dout(dc::notice, "Real root: " << std::setprecision(18) << real_root << "; root found: " << root);

    int steps = 0;
    double r = root;
    if (r != real_root)
    {
      double direction = real_root > r ? +INFINITY : -INFINITY;
      do
      {
        r = std::nextafter(r, direction);
        ++steps;
      }
      while (r != real_root);
    }

    Dout(dc::notice, std::setprecision(18) << "Cubic: " << cubic << "; root found: " << root << "; actual root: " << real_root << "; offset: " << steps << " steps.");

    double c0 = cubic[0];
    double c1 = cubic[1];
    double c2 = cubic[2];
    double c3 = cubic[3];
    double d = c2 / c3;
    double D = 1.0 - 4.0 * c0 * c2 / utils::square(c1);
    double e = -0.5 * c1 / c2 * (1.0 - std::sqrt(D));
    double epsilon = e / d;
    bool special_case1 = D >= 0.0 && std::abs(epsilon) < 0.1;

    RESTART

    bool special_case2 = false;
    {
      double c2_div_c3 = c2 / c3;
      double c1_div_c3 = c1 / c3;
      double discriminant = utils::square(c2_div_c3) - 3.0 * c1_div_c3;

      if (discriminant >= 0.0)
      {
        double e_minus = (-c2_div_c3 - std::sqrt(discriminant)) / 3.0;
        double e_plus = (-c2_div_c3 + std::sqrt(discriminant)) / 3.0;
        double eval_e_minus = cubic(e_minus);
        double eval_e_plus = cubic(e_plus);
        double extreme = std::abs(eval_e_minus) < std::abs(eval_e_plus) ? e_minus : e_plus;

        double qc2 = 0.5 * cubic.second_derivative(extreme);
        math::QuadraticPolynomial qp(cubic(extreme), cubic.derivative(extreme), qc2);
        std::array<double, 2> roots;
        int n = qp.get_roots(roots);
        if (n > 0)
        {
          if (std::abs(roots[0] * c3 / qc2) < 0.1 && std::abs(roots[1] * c3 / qc2) < 0.1)
            special_case2 = true;
        }
      }
    }

    static int count = 0;
    ++count;
    if (!special_case1 && !special_case2)
    {
      Dout(dc::notice(D >= 0), "epsilon = " << epsilon);
      if (count == 1000 || steps > 9)
      {
        std::ostringstream title;
        title << cubic;
        cairowindow::QuickGraph qg(title.str(), "x", "y", {-10.0, 10.0}, {-10.0, 10.0});
        qg.add_function([&](double x){ return cubic(x); });
        qg.wait_for_keypress();
      }
    }
    else
    {
      double rel_err = (e - real_root) / std::abs(real_root);
      Dout(dc::notice, std::setprecision(18) << "e = " << e << ", epsilon = " << epsilon << ", rel_err = " << rel_err);
      static double abs_max_rel_err = 0.0;
      if (std::abs(rel_err) > abs_max_rel_err)
      {
        abs_max_rel_err = std::abs(rel_err);
        Dout(dc::notice, "abs_max_rel_err = " << abs_max_rel_err);
      }
    }
  }
  catch (std::runtime_error const& error)
  {
    Dout(dc::warning, error.what());
  }
}

void sanity_check(math::CubicPolynomial const& cubic, std::array<double, 3> const& roots, int number_of_roots)
{
  std::array<double, 3> rs = roots;
  if (number_of_roots == 3)
    std::sort(rs.begin(), rs.end(), [](double r1, double r2){ return std::abs(r1) < std::abs(r2); });
//  for (int r = 0; r < number_of_roots; ++r)
  sanity_check(cubic, rs[0]);
}

int main(int argc, char* argv[])
{
  Debug(NAMESPACE_DEBUG::init());

#if 0
  math::CubicPolynomial p(8.59730539601805788, -2.59172811830510774, -8.6459422955246783, 2.09521187806639375);
  int iterations;
  std::array<double, 3> roots;
  int n = get_roots(p, roots, iterations);
  sanity_check(p, roots, n);
  Dout(dc::notice, "Cubic: " << std::setprecision(18) << p << " has " << n << " roots; which we found in " << iterations << " iterations.");
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

  // Generate random cubics.
  std::cout << "Generating random cubic polynomials..." << std::endl;
  std::vector<math::CubicPolynomial> cubics;
  std::uniform_real_distribution<double> dist(-10.0, 10.0);
  for (int i = 0; i < number_of_cubics; ++i)
  {
    cubics.emplace_back(dist(engine), dist(engine), dist(engine), dist(engine));
  }
  std::cout << "Done (" << cubics.size() << " cubics)." << std::endl;

  unsigned long total_number_of_iterations = 0;
  for (int i = 0; i < number_of_cubics; ++i)
  {
    int iterations;
    std::array<double, 3> roots;
    RESTART
    int n = get_roots(cubics[i], roots, iterations);
    sanity_check(cubics[i], roots, n);
    Dout(dc::notice, "Cubic: " << std::setprecision(18) << cubics[i] << " has " << n << " roots; which we found in " << iterations << " iterations.");
    total_number_of_iterations += iterations;
  }

#ifdef CWDEBUG
  Dout(dc::notice, "Average number of iterations: " << (static_cast<double>(total_number_of_iterations) / number_of_cubics));
#else
  std::cout << "Average number of iterations: " << (static_cast<double>(total_number_of_iterations) / number_of_cubics) << std::endl;
#endif
}
