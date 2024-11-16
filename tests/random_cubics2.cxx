#include "sys.h"
#include "math/CubicPolynomial.h"
#include "math/bracket_zero.h"
#include <random>
#include <chrono>
#include <sstream>
#include <array>
#include "debug.h"

constexpr int number_of_cubics = 1000000;

int get_roots(math::CubicPolynomial& cubic, std::array<double, 3>& roots_out, int& iterations)
{
  DoutEntering(dc::notice, "get_roots() for " << cubic);

  // Define a coefficients_ for use by math/CubicPolynomial_get_roots.cpp.
  std::array<double, 4>& coefficients_ = reinterpret_cast<std::array<double, 4>&>(cubic[0]);

  using math::QuadraticPolynomial;
  // Include the body of the function.
# define GETROOTS_ASSIGN_ITERATIONS
# include "math/CubicPolynomial_get_roots.cpp"
# undef GETROOTS_ASSIGN_ITERATIONS
}

void sanity_check(math::CubicPolynomial const& cubic, double root)
{
  // Course check.
  double fr = cubic(root);

  double root_small = root < 0 ? root * 1.000000001 : root * 0.999999999;
  double root_large = root < 0 ? root * 0.999999999 : root * 1.000000001;

  double frm = cubic(root_small);
  double frp = cubic(root_large);

  ASSERT(std::abs(fr) < std::abs(frm) && std::abs(fr) < std::abs(frp));

  try
  {
    double real_root = math::bracket_zero(root_small, root_large, [&](double r){ return cubic(r); });
    Dout(dc::notice, "Real root: " << std::setprecision(18) << real_root << "; root found: " << root);

    int steps = 0;
    if (root != real_root)
    {
      double direction = real_root > root ? +INFINITY : -INFINITY;
      do
      {
        root = std::nextafter(root, direction);
        ++steps;
      }
      while (root != real_root);
    }
    Dout(dc::notice, "Distance: " << steps);
  }
  catch (std::runtime_error const& error)
  {
    Dout(dc::warning, error.what());
  }
}

void sanity_check(math::CubicPolynomial const& cubic, std::array<double, 3> const& roots, int number_of_roots)
{
  for (int r = 0; r < number_of_roots; ++r)
    sanity_check(cubic, roots[0]);
}

int main(int argc, char* argv[])
{
  Debug(NAMESPACE_DEBUG::init());

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
    int n = get_roots(cubics[i], roots, iterations);
    sanity_check(cubics[i], roots, n);
    Dout(dc::notice, "Cubic: " << cubics[i] << " has " << n << " roots; which we found in " << iterations << " iterations.");
    total_number_of_iterations += iterations;
  }

#ifdef CWDEBUG
  Dout(dc::notice, "Average number of iterations: " << (static_cast<double>(total_number_of_iterations) / number_of_cubics));
#else
  std::cout << "Average number of iterations: " << (static_cast<double>(total_number_of_iterations) / number_of_cubics) << std::endl;
#endif
}
