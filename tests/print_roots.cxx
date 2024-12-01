#include "sys.h"
#include "utils/square.h"
#include "utils/split.h"
#include "cairowindow/QuickGraph.h"
#include "math/CubicPolynomial.h"
#include "math/bracket_zero.h"
#include "mpreal/mpreal.h"
#include <array>
#include <string>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <concepts>
#include "utils/debug_ostream_operators.h"
#include "debug.h"

using namespace mpfr;

constexpr mp_prec_t precision_in_bits = 256;
constexpr int mpprecision = precision_in_bits * std::log10(2.0);
constexpr double rel_err_threshold = std::pow(10.0, -0.5 * mpprecision);

// Update the concept definition to correctly check the streamability.
template<typename T>
concept StreamExtractable = requires (std::stringstream& ss, T& value) {
  { ss >> value } -> std::same_as<std::istream&>;
};

// Function to parse a string in the desired form and return an array of coefficients.
template<typename T>
std::array<T, 4> parse_string(std::string const& input)
{
  std::regex term_regex(R"(([+-])?\s*((0[.]|[1-9](\d*[.])?)\d*(e[+-]?\d+)?)(\s*(([a-z])(\^([0-9]+))?)?)?)");
  auto terms_begin = std::sregex_iterator(input.begin(), input.end(), term_regex);
  auto terms_end = std::sregex_iterator();

  std::array<T, 4> coefficients{};

  std::string symbol;
  bool first = true;
  for (std::sregex_iterator i = terms_begin; i != terms_end; ++i, first = false)
  {
    std::smatch const& match = *i;
    double sign = match[1].length() == 1 && match[1].str()[0] == '-' ? -1.0 : 1.0;
    if (sign == 1.0 && !first && (match[1].length() == 0 || match[1].str()[0] != '+'))
      THROW_ALERT("Expected leading + or - before \"[TERM]\".", AIArgs("[TERM]", match.str()));

    T value;
    if constexpr (StreamExtractable<T>)
    {
      std::stringstream ss;
      ss << match[2];
      ss >> value;
    }
    else
      value = T{match[2].str()};
    value *= sign;
    if (symbol.empty())
      symbol = match[8].str();
    else if (match[8].length() > 0)
    {
      if (symbol != match[8].str())
        THROW_ALERT("Mismatching symbols [PREV] != [CURRENT].", AIArgs("[PREV]", symbol)("[CURRENT]", match[8]));
    }
    int exponent = symbol.empty() ? 0 : (match[10].length() == 0 ? 1 : match[10].str()[0] - '0');

    if (0 <= exponent && exponent < coefficients.size())
    {
      if (coefficients[exponent] != 0.0)
        THROW_ALERT("Duplicated term with exponent [EXPONENT].", AIArgs("[EXPONENT]", exponent));
      coefficients[exponent] = value;
    }
  }

  return coefficients;
}

mpreal exact_root(double guess, std::array<mpreal, 4> coefficients)
{
  mpreal u = guess;

  auto evaluate = [&](mpreal u) -> mpreal { return coefficients[0] + (coefficients[1] + (coefficients[2] + coefficients[3] * u) * u) * u; };
  auto derivative = [&](mpreal u) -> mpreal { return coefficients[1] + (2 * coefficients[2] + 3 * coefficients[3] * u) * u; };
  auto half_second_derivative = [&](mpreal u) -> mpreal { return coefficients[2] + 3 * coefficients[3] * u; };

  mpreal prev_u = 0;
  mpreal step = 0;
  do
  {
    prev_u = u;
    mpreal Q_u = evaluate(u);
    mpreal dQ_u = derivative(u);
    mpreal half_ddQ_u = half_second_derivative(u);
    step = Q_u * dQ_u / (utils::square(dQ_u) - Q_u * half_ddQ_u);
    u -= step;
  }
  while (abs((prev_u - u) / u) > rel_err_threshold);

  return u;
}

int main(int argc, char* argv[])
{
  Debug(NAMESPACE_DEBUG::init());

  // Set the default floating-point precision.
  mpfr_set_default_prec(precision_in_bits);

  // Print everything with mpprecision digits precision.
  std::streamsize const old_precision = std::cout.precision(mpprecision);
//  std::ios::fmtflags const old_flags = std::cout.setf(std::ios::fixed);

  // Check if there is exactly one argument provided.
  if (argc != 2)
  {
    std::cerr << "Usage: " << argv[0] << " <argument>" << std::endl;
    return 1;
  }

  // Store the argument in a std::string.
  std::string input = argv[1];

  // Print the string.
  std::cout << "The input argument is: \"" << input << "\"." << std::endl;

  try
  {
    math::CubicPolynomial cubic(parse_string<double>(input));
    std::cout << std::setprecision(18) << cubic << '\n';

    double c0 = cubic[0];
    double c1 = cubic[1];
    double c2 = cubic[2];
    double c3 = cubic[3];
    {
      double d = c2 / c3;
      double D = 1.0 - 4.0 * c0 * c2 / utils::square(c1);
      double e = -0.5 * c1 / c2 * (1.0 - std::sqrt(D));
      double epsilon = e / d;
      bool special_case1 = D >= 0.0 && std::abs(epsilon) < 0.1;
      Dout(dc::notice, "special_case1: D = " << D << "; epsilon = " << epsilon << " -->" << (special_case1 ? "" : " NOT") << " special case 1.");
    }
    {
      double c2_div_c3 = c2 / c3;
      double c1_div_c3 = c1 / c3;
      double D = utils::square(c2_div_c3) - 3.0 * c1_div_c3;

      if (D >= 0.0)
      {
        double e_minus = (-c2_div_c3 - std::sqrt(D)) / 3.0;
        double e_plus = (-c2_div_c3 + std::sqrt(D)) / 3.0;
        double eval_e_minus = cubic(e_minus);
        double eval_e_plus = cubic(e_plus);
        double extreme = std::abs(eval_e_minus) < std::abs(eval_e_plus) ? e_minus : e_plus;

        double qc2 = 0.5 * cubic.second_derivative(extreme);
        math::QuadraticPolynomial qp(cubic(extreme), cubic.derivative(extreme), qc2);
        std::array<double, 2> roots;
        int n = qp.get_roots(roots);
        if (n > 0)
        {
          ASSERT(n == 2);
          Dout(dc::notice|continued_cf, "special_case2: qp = " << qp << " where root*c3/qc2 is " << (roots[0] * c3 / qc2) <<
              " and " << (roots[1] * c3 / qc2) << " -->");
          bool special_case2 = false;
          if (std::abs(roots[0] * c3 / qc2) < 0.1 && std::abs(roots[1] * c3 / qc2) < 0.1)
            special_case2 = true;
          Dout(dc::finish, (special_case2 ? "" : " NOT") << " special case 2.");
        }
        else
          Dout(dc::notice, "special_case2: qp = " << qp << " which has no roots! --> NOT special case 2.");
      }
      else
        Dout(dc::notice, "special_case2: D = " << D << " --> NOT special case 2.");
    }

    std::ostringstream title;
    title << cubic;
    cairowindow::QuickGraph qg(title.str(), "x", "y", {-10.0, 10.0}, {-10.0, 10.0});
    qg.add_function([&](double x){ return cubic(x); });

    std::array<mpreal, 4> coefficients = parse_string<mpreal>(input);
    for (int i = 0; i < coefficients.size(); ++i)
      std::cout << "i=" << i << " : " << coefficients[i] << std::endl;

    std::array<double, 3> roots;
    int n = cubic.get_roots(roots);
    std::cout << "This cubic has " << n << " root" << (n != 1 ? "s" : "") << ":\n";
    for (int i = 0; i < n; ++i)
    {
      mpreal root = exact_root(roots[i], coefficients);
      std::cout << "  " << root << " (as double: " << root.toDouble() << ")" << '\n';
    }

    qg.wait_for_keypress();
  }
  catch (AIAlert::Error const& error)
  {
    Dout(dc::warning, error);
  }
}
