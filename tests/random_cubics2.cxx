#include "sys.h"
#include "math/CubicPolynomial.h"
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
    Dout(dc::notice, "Cubic: " << cubics[i] << " has " << n << " roots; which we found in " << iterations << " iterations.");
    total_number_of_iterations += iterations;
  }

#ifdef CWDEBUG
  Dout(dc::notice, "Average number of iterations: " << (static_cast<double>(total_number_of_iterations) / number_of_cubics));
#else
  std::cout << "Average number of iterations: " << (static_cast<double>(total_number_of_iterations) / number_of_cubics) << std::endl;
#endif
}
