#include <iostream>
#include <array>
#include <set>

struct Float
{
  bool negative_;
  int value_;

  Float() : negative_(false), value_(0) { }
  Float(bool negative, int value) : negative_(negative), value_(value) { }

  bool is_negative_zero() const
  {
    return negative_ && value_ == 0;
  }

  bool is_zero() const
  {
    return value_ == 0;
  }

  // For sorting (kinda like the mathematical compare, but -0 < 0).
  bool is_less_than(Float const& val) const
  {
    if (negative_ != val.negative_)
      return negative_;
    return value_ < val.value_;
  }

  // Negation.
  Float operator-() const
  {
    return {!negative_, value_};
  }

  // Mathematical compare (0 == -0).
  friend bool operator==(Float const& val0, Float const& val1)
  {
    return val0.value_ == val1.value_ && (val0.value_ == 0 || val0.negative_ == val1.negative_);
  }

  std::string as_binary(char s, char v) const
  {
    std::string result(1, s);
    result += '=';
    result += negative_ ? '1' : '0';
    result += ' ';
    result += v;
    result += '=';
    result += value_ > 0 ? '1' : '0';
    return result;
  }
};

// Printing.
std::ostream& operator<<(std::ostream& os, Float const& val)
{
  os << (val.negative_ ? '-' : '+') << val.value_;
  return os;
}

// Addition.
Float operator+(Float val0, Float val1)
{
  Float sum{val0};

  if (val0.negative_ == val1.negative_)
    sum.value_ += val1.value_;
  else
  {
    sum.value_ -= val1.value_;
    if (sum.value_ < 0)
    {
      sum.value_ = -sum.value_;
      sum.negative_ = !sum.negative_;
    }
  }

  // Never return -0 unless both values are negative zero.
  if (sum.value_ == 0 && sum.negative_)
    sum.negative_ = val0.is_negative_zero() && val1.is_negative_zero();

  return sum;
}

// Product.
Float operator*(Float val0, Float val1)
{
  return {val0.negative_ != val1.negative_, val0.value_ * val1.value_};
}

std::array<Float, 8> const possible_roots = {{
  { true, 3 },          // -3
  { true, 2 },          // -2
  { true, 1 },          // -1
  { true, 0 },          // -0
  { false, 0 },         // +0
  { false, 1 },         // +1
  { false, 2 },         // +2
  { false, 3 }          // +3
}};

struct Quadratic
{
  // P(x) = coefficients_[0] + coefficients_[1] * x + coefficients_[2] * x^2.
  std::array<Float, 3> coefficients_;

  Float root0_;
  Float root1_;

  // P(x) = a x^2 + b x + c.
  Quadratic(Float a, Float b, Float c, Float root0, Float root1) : coefficients_{c, b, a}, root0_(root0), root1_(root1) { }

  std::string as_binary() const
  {
    std::string result;
    for (int i = 2; i >= 0; --i)
    {
      result += coefficients_[i].as_binary('A' + 2 * (2 - i), 'A' + 2 * (2 - i) + 1);
      if (i != 0)
        result += ' ';
    }
    return result;
  }
};

std::ostream& operator<<(std::ostream& os, Quadratic const& quadratic)
{
  os << quadratic.coefficients_[2] << " * x^2 + " << quadratic.coefficients_[1] << " * x + " << quadratic.coefficients_[0];
  return os;
}

struct Sort
{
  bool operator()(Quadratic const& poly0, Quadratic const& poly1) const
  {
    if (poly0.coefficients_[2].is_less_than(poly1.coefficients_[2]))
      return true;
    else if (poly1.coefficients_[2].is_less_than(poly0.coefficients_[2]))
      return false;
    else if (poly0.coefficients_[1].is_less_than(poly1.coefficients_[1]))
      return true;
    else if (poly1.coefficients_[1].is_less_than(poly0.coefficients_[1]))
      return false;
    else if (poly0.coefficients_[0].is_less_than(poly1.coefficients_[0]))
      return true;
    else
      return false;
  }
};

int output(Quadratic const& quadratic)
{
  if (quadratic.root0_.value_ < quadratic.root1_.value_)
    return 0;
  else if (quadratic.root1_.value_ < quadratic.root0_.value_)
    return 1;

  return quadratic.coefficients_[0].negative_ != quadratic.coefficients_[1].negative_;
}

int formula(Quadratic const& quadratic)
{
  return quadratic.coefficients_[2].negative_ == quadratic.coefficients_[1].negative_;
}

int main()
{
  std::set<Quadratic, Sort> all_quadratics;

  std::cout << "Sum and product of all possible roots:\n";
  // Brute force all the possible combinations of root0 and root1.
  for (int i0 = 0; i0 < possible_roots.size() - 1; ++i0)
    // Make sure that root0 < root1, because if root0 == root1 then we don't care where we put the first calculated root.
    for (int i1 = i0 + 1; i1 < possible_roots.size(); ++i1)
    {
      Float root0 = possible_roots[i0];
      Float root1 = possible_roots[i1];
      Float negative_sum = -root0 + -root1;      // Never -0 unless root0 == root1 == +0, which we don't care about.
      Float product = root0 * root1;

      // Make sure also product is never -0.
      if (product.is_negative_zero())
        product.negative_ = false;

      std::cout << -root0 << " + " << -root1 << " = " << negative_sum <<
        "; " << root0 << " * " << root1 << " = " << product << std::endl;

      if (root0.is_negative_zero())
        root0.negative_ = false;

      if (root1.is_negative_zero())
        root1.negative_ = false;

      if (root0 == root1)
        continue;

      // Store all possible combinations of coefficients.
      all_quadratics.emplace(Float{false, 1}, negative_sum, product, root0, root1);
      all_quadratics.emplace(Float{true, 1}, -negative_sum, -product, root0, root1);
      if (negative_sum.is_zero())
      {
        all_quadratics.emplace(Float{false, 1}, -negative_sum, product, root0, root1);
        all_quadratics.emplace(Float{true, 1}, negative_sum, -product, root0, root1);
      }
      if (product.is_zero())
      {
        all_quadratics.emplace(Float{false, 1}, negative_sum, -product, root0, root1);
        all_quadratics.emplace(Float{true, 1}, -negative_sum, product, root0, root1);
      }
      if (negative_sum.is_zero() && product.is_zero())
      {
        all_quadratics.emplace(Float{false, 1}, -negative_sum, -product, root0, root1);
        all_quadratics.emplace(Float{true, 1}, negative_sum, product, root0, root1);
      }
    }

  std::cout << "All possible coefficients, with their respective roots:\n";
  for (auto&& quadratic : all_quadratics)
  {
    std::cout << quadratic << "; root0 = " << quadratic.root0_ << ", root1 = " << quadratic.root1_ <<
      "; first root: root" << output(quadratic) << '\n';
  }

  std::cout << "Same, but as binary table:\n";
  for (auto&& quadratic : all_quadratics)
  {
    std::cout << output(quadratic) << " <-- " << quadratic.as_binary() << " : " << formula(quadratic) << '\n';
  }
}
