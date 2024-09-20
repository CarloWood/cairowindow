#include "sys.h"
#include "CubicPolynomial.h"
#include "QuadraticPolynomial.h"
#include "debug.h"

// Helper function to create a quadratic with given roots
math::QuadraticPolynomial make_quadratic(double root1, double root2)
{
  double sum = -(root1 + root2);
  double product = root1 * root2;
  return {product, sum, 1.0};  // product + sum*x + x^2
}

// Helper function to create a cubic whose derivative is the given quadratic
math::CubicPolynomial make_cubic(math::QuadraticPolynomial const& quadratic)
{
  return {0.0, quadratic[0], quadratic[1] / 2.0, quadratic[2] / 3.0};
}

int signof(double d)
{
  if (d < 0.0)
    return -1;
  else if (d > 0.0)
    return 1;
  return 0;
}

void test_cubic(math::CubicPolynomial const& poly, double root1, double root2, bool left_most_first)
{
  std::array<double, 2> extrema;
  int num_extrema = poly.get_extrema(extrema, left_most_first);

  std::cout << "Cubic: " << poly[0] << " + " << poly[1] << "x + " << poly[2] << "x^2 + " << poly[3] << "x^3" << std::endl;
  std::cout << "Expected extrema: " << root1 << ", " << root2 << std::endl;
  std::cout << "left_most_first: " << (left_most_first ? "true" : "false") << std::endl;
  std::cout << "Number of extrema: " << num_extrema << std::endl;
  std::cout << "Extrema: ";
  for (int i = 0; i < num_extrema; ++i) {
    std::cout << std::setprecision(6) << extrema[i] << " ";
  }
  std::cout << std::endl;

  // Verify extrema.
  bool correct = true;
  if (num_extrema != (root1 == root2 ? 1 : 2))
  {
    DoutFatal(dc::core, "Expected 2 extrema, got " << num_extrema);
    correct = false;
  }
  else
  {
    double expected_first = left_most_first ? std::min(root1, root2) : (poly[3] < 0 ? std::min(root1, root2) : std::max(root1, root2));
    double expected_second = left_most_first ? std::max(root1, root2) : (poly[3] < 0 ? std::max(root1, root2) : std::min(root1, root2));

    if (std::abs(extrema[0] * extrema[1] - root1 * root2) > 1e-6 || std::abs(extrema[0] + extrema[1] - (root1 + root2)) > 1e-6)
    {
      DoutFatal(dc::core, "Incorrect roots!");
      correct = false;
    }

    if (std::abs(extrema[0] - expected_first) + std::abs(extrema[1] - expected_second) >
        std::abs(extrema[0] - expected_second) + std::abs(extrema[1] - expected_first))
    {
      Dout(dc::notice, "1: " << left_most_first << ", " << signof(poly[3]) << ", " << signof(poly[2]) << ", " << signof(poly[1]));
      DoutFatal(dc::core, "Incorrect order!");
      correct = false;
    }
    else
    {
      Dout(dc::notice, "0: " << left_most_first << ", " << signof(poly[3]) << ", " << signof(poly[2]) << ", " << signof(poly[1]));
    }
  }

  // Verify derivative at extrema
  for (int i = 0; i < num_extrema; ++i)
  {
    double deriv = poly.derivative(extrema[i]);
    if (std::abs(deriv) > 1e-6)
    {
      DoutFatal(dc::core, "Non-zero derivative at extremum " << extrema[i] << ": " << deriv);
      correct = false;
    }
  }

  if (correct)
    std::cout << "Test PASSED" << std::endl;
  else
    std::cout << "Test FAILED" << std::endl;

  std::cout << std::endl;
}

int main()
{
  Debug(NAMESPACE_DEBUG::init());

  for (int root1 = -2; root1 <= 2; ++root1)
    for (int root2 = root1; root2 <= 2; ++root2)
  {
    Dout(dc::notice, "roots: " << root1 << ", " << root2);

    math::QuadraticPolynomial quadratic = make_quadratic(root1, root2);
    Dout(dc::notice, "quadratic = " << quadratic);

    math::CubicPolynomial poly = make_cubic(quadratic);
    Dout(dc::notice, "poly = " << poly << " with derivative: " << poly.derivative());

    // Test with both orderings
    test_cubic(poly, root1, root2, true);
    test_cubic(poly, root1, root2, false);

    // Test with negative cubic term (multiply all coefficients by -1).
    math::CubicPolynomial poly_neg(-poly[0], -poly[1], -poly[2], -poly[3]);
    test_cubic(poly_neg, root1, root2, true);
    test_cubic(poly_neg, root1, root2, false);
  }
}
