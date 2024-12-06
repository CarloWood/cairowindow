#include "mpreal/mpreal.h"
#include "utils/AIAlert.h"
#ifdef CWDEBUG
#include "utils/has_print_on.h"
#include <iomanip>
#endif

#ifdef CWDEBUG
// This class defines a print_on method.
using utils::has_print_on::operator<<;
#endif

class IEEE754
{
 private:
  mpfr::mpreal m_low_value;
  mpfr::mpreal m_high_value;

 public:
  // Default constructor to initialize both low and high to zero.
  IEEE754() : m_low_value(0), m_high_value(0) { }

  // Constructor from double to initialize the value and calculate initial error.
  IEEE754(double value)
  {
    // Calculate round-off error: half the value of the least significant bit.
    mpfr::mpreal mpvalue = value;
    mpfr::mpreal epsilon_plus = std::nextafter(value, +INFINITY);
    mpfr::mpreal epsilon_minus = std::nextafter(value, -INFINITY);

    m_low_value = mpvalue + (epsilon_minus - mpvalue) / 2;
    m_high_value = mpvalue + (epsilon_plus - mpvalue) / 2;
  }

  // Construct an IEEE754 directly from the low and high value.
  IEEE754(mpfr::mpreal low_value, mpfr::mpreal high_value) : m_low_value(std::move(low_value)), m_high_value(std::move(high_value)) { }

  // Copy constructor
  IEEE754(IEEE754 const& other) = default;

  // Assignment operator
  IEEE754& operator=(IEEE754 const& other) = default;

  // Addition operator.
  IEEE754 operator+(IEEE754 const& other) const
  {
    return {m_low_value + other.m_low_value, m_high_value + other.m_high_value};
  }

  // Subtraction operator.
  IEEE754 operator-(IEEE754 const& other) const
  {
    return {m_low_value - other.m_high_value, m_high_value - other.m_low_value};
  }

  // Multiplication operator.
  IEEE754 operator*(IEEE754 const& other) const
  {
    // Calculate product including potential errors, considering all combinations of bounds to handle negative values.
    mpfr::mpreal candidates[4] = {
      m_low_value * other.m_low_value,
      m_low_value * other.m_high_value,
      m_high_value * other.m_low_value,
      m_high_value * other.m_high_value
    };
    return {*std::min_element(std::begin(candidates), std::end(candidates)), *std::max_element(std::begin(candidates), std::end(candidates))};
  }

  // Division operator.
  IEEE754 operator/(IEEE754 const& other) const
  {
    // Ensure no division by zero.
    if (other.m_low_value <= 0 && other.m_high_value >= 0)
      THROW_ALERT("Division by zero error in IEEE754 class: denumerator is in the range [[MIN], [MAX]]",
          AIArgs("[MIN]", other.m_low_value)("[MAX]", other.m_high_value));

    // Calculate division including potential errors, considering all combinations of bounds to handle negative values.
    mpfr::mpreal candidates[4] = {
      m_low_value / other.m_low_value,
      m_low_value / other.m_high_value,
      m_high_value / other.m_low_value,
      m_high_value / other.m_high_value
    };
    return {*std::min_element(std::begin(candidates), std::end(candidates)), *std::max_element(std::begin(candidates), std::end(candidates))};
  }

  // Accessor methods to get the current range.
  mpfr::mpreal const& get_low() const {
    return m_low_value;
  }

  mpfr::mpreal const& get_high() const {
    return m_high_value;
  }

  mpfr::mpreal average() const
  {
    return (m_low_value + m_high_value) / 2;
  }

  // Function to convert back to double (approximate due to precision limits).
  double to_double() const
  {
    return average().toDouble();
  }

#ifdef CWDEBUG
  // Output operator for debugging.
  void print_on(std::ostream& os) const
  {
    mpfr::mpreal avg = average();
    mpfr::mpreal error = (m_high_value - m_low_value) / 2;

    // Calculate the precision required to display two significant digits of the error.
    int error_exponent = std::floor(mpfr::log10(error).toDouble());
    int precision = std::max(0, 1 - error_exponent);

    os << std::setprecision(precision) << avg << " ± " << error << " (" << to_double() << ')';
  }
#endif
};
