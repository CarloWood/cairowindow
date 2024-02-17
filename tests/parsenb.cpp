#include <fstream>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>
#include <cassert>
#include <stdexcept>
#include <charconv>
#include <cctype>
#include <array>
#include <algorithm>

enum Token
{
  None,
  Keyword,      // RowBox, FractionBox, SqrtBox, SuperscriptBox
  String,       // From opening double quote up till and including closing double quote.
  Other         // Single character that is not one of the above.
};

struct Decoder;
struct GFactor;

enum Operator
{
  Add,
  Subtract,
  Multiply,
  Divide
};

struct Expression
{
  virtual ~Expression() = default;
  virtual void parse(Decoder* decoder) { assert(false); }
  virtual void add(Expression* element) { assert(false); }
  virtual void print_on(std::ostream& os) const = 0;
  virtual void divide_by_sqrt_d0() { assert(false); }
  virtual void multiply_with_sqrt_d0() { assert(false); }
  virtual void multiply_with(int n) { assert(false); }
  virtual void divide_by(int n) { assert(false); }
  virtual void remove_u() { assert(false); }
  virtual bool needs_parens(Operator op) const { return false; }
  virtual void simplify(Expression** self_ptr) = 0;
  virtual void replace_g(Expression** self_ptr, GFactor* gn) = 0;
  virtual bool equals(Expression* expression) const = 0;

  void parse_keyword(Decoder* decoder);
  void parse_string(Decoder* decoder);
  void parse_argument_list(Decoder* decoder, int min = 1, int max = 1000);

  friend std::ostream& operator<<(std::ostream& os, Expression const& expression)
  {
    expression.print_on(os);
    return os;
  }
};

struct StringExpr : Expression
{
  std::string value_;

  StringExpr(std::string_view const& str) : value_(str) { }
  void simplify(Expression** self_ptr) override { assert(false); }
  bool equals(Expression* expression) const override { assert(false); }
  void replace_g(Expression** self_ptr, GFactor* gn) override { assert(false); }

  void print_on(std::ostream& os) const override
  {
    os << '"' << value_ << '"';
  }
};

struct Symbol : Expression
{
  std::string name_;
  int exponent_;

  Symbol(std::string_view name_str, int exponent = 1) : name_(name_str), exponent_(exponent) { }
  void simplify(Expression** self_ptr) override { }
  bool equals(Expression* expression) const override
  {
    Symbol* symbol = dynamic_cast<Symbol*>(expression);
    return symbol && name_ == symbol->name_ && exponent_ == symbol->exponent_;
  }
  void replace_g(Expression** self_ptr, GFactor* gn) override
  {
    // g is never a Symbol.
  }

  void print_on(std::ostream& os) const override
  {
    static std::array<char const*, 10> subscript = {
      "₀", "₁", "₂", "₃", "₄", "₅", "₆", "₇", "₈", "₉"
    };
    static std::array<char const*, 10> superscript = {
      "⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹"
    };

    for (char c : name_)
    {
      if (std::isdigit(c))
        os << subscript[c - '0'];
      else
        os << c;
    }
    int exp = exponent_;
    assert(exp > 0);
    if (exp > 1)
    {
      int n = (exp > 9) ? 10 : 1;
      while (exp)
      {
        os << superscript[exp / n];
        exp %= n;
        n /= 10;
      }
    }
  }
};

struct Integer : Expression
{
  int value_;

  Integer(std::string_view name_str)
  {
    std::from_chars(name_str.begin(), name_str.end(), value_);
  }
  void simplify(Expression** self_ptr) override { }
  bool equals(Expression* expression) const override
  {
    Integer* integer = dynamic_cast<Integer*>(expression);
    return integer && value_ == integer->value_;
  }
  void replace_g(Expression** self_ptr, GFactor* gn) override
  {
    // g is never an Integer.
  }

  void print_on(std::ostream& os) const override
  {
    os << value_;
  }
};

struct Square : Expression
{
  Symbol symbol_;

  Square(std::string_view const& name_str) : symbol_(name_str) { }
  void parse(Decoder* decoder) override;
  bool equals(Expression* expression) const override { assert(false); }

  void simplify(Expression** self_ptr) override
  {
    *self_ptr = new Symbol("d0", 2);
  }
  void replace_g(Expression** self_ptr, GFactor* gn) override { assert(false); }

  void print_on(std::ostream& os) const override
  {
    os << symbol_ << "²";
  }
};

struct Powerd0 : Expression
{
  int exponent_times_two_;
  void parse(Decoder* decoder) override;

  void simplify(Expression** self_ptr) override
  {
    assert((exponent_times_two_ & 1) == 0);
    *self_ptr = new Symbol("d0", exponent_times_two_ / 2);
  }
  bool equals(Expression* expression) const override { assert(false); }
  void replace_g(Expression** self_ptr, GFactor* gn) override { assert(false); }

  void print_on(std::ostream& os) const override
  {
    if ((exponent_times_two_ & 1) == 1)
      os << "d₀^(" << exponent_times_two_ << "/2)";
    else
      os << "d₀^" << (exponent_times_two_ / 2);
  }
};

struct Poweru : Expression
{
  int exponent_;

  void parse(Decoder* decoder) override;
  void simplify(Expression** self_ptr) override { assert(false); }
  bool equals(Expression* expression) const override { assert(false); }
  void replace_g(Expression** self_ptr, GFactor* gn) override { assert(false); }

  void print_on(std::ostream& os) const override
  {
    os << "u^" << exponent_;
  }
};

struct GFactor : Expression
{
  int n_;
  Expression* gn_;

  bool is_equal(Expression* expression) const
  {
    if (expression == gn_)
      return true;
    return gn_->equals(expression);
  }

  GFactor(int n, Expression* gn) : n_(n), gn_(gn) { }
  void simplify(Expression** self_ptr) override { assert(false); }
  bool equals(Expression* expression) const override
  {
    if (this == expression)
      return true;
    return gn_->equals(expression);
  }
  void replace_g(Expression** self_ptr, GFactor* gn) override
  {
    // Nothing to do.
  }

  void print_on(std::ostream& os) const override
  {
    os << "g" << n_; // << " \e[34m[" << *gn_ << "]\e[0m";
  }
};

struct Fraction : Expression
{
  Expression* enumerator_{nullptr};
  Expression* denominator_{nullptr};

  void parse(Decoder* decoder) override;
  void add(Expression* element) override
  {
    if (!enumerator_)
      enumerator_ = element;
    else
    {
      assert(!denominator_);
      denominator_ = element;
    }
  }

  void divide_by_sqrt_d0() override
  {
    denominator_->multiply_with_sqrt_d0();
  }

  void multiply_with(int n) override
  {
    denominator_->divide_by(n);
  }

  void remove_u() override
  {
    enumerator_->remove_u();
  }

  void simplify(Expression** self_ptr) override
  {
    enumerator_->simplify(&enumerator_);
    denominator_->simplify(&denominator_);
  }

  void replace_g(Expression** self_ptr, GFactor* gn) override
  {
    //std::cout << "Fraction::replace_g(" << *this << ", " << *gn << ")" << std::endl;

    enumerator_->replace_g(&enumerator_, gn);
    if (gn->is_equal(this))
      *self_ptr = gn;
  }

  bool equals(Expression* expression) const override
  {
    Fraction* fraction = dynamic_cast<Fraction*>(expression);
    return fraction && enumerator_->equals(fraction->enumerator_) && denominator_->equals(fraction->denominator_);
  }

  bool needs_parens(Operator op) const override { return op == Divide; }

  void print_on(std::ostream& os) const override
  {
    Integer const* enumerator = dynamic_cast<Integer const*>(enumerator_);
    Integer const* denominator = dynamic_cast<Integer const*>(denominator_);
    if (enumerator)
      enumerator->print_on(os);
    else
    {
      bool need_parens = enumerator_->needs_parens(Multiply);
      if (need_parens)
        os << '(';
      enumerator_->print_on(os);
      if (need_parens)
        os << ')';
    }
    if (enumerator && denominator)
      os << "∕";
    else
      os << " ∕ ";
    if (denominator)
      denominator->print_on(os);
    else
    {
      bool need_parens = denominator_->needs_parens(Divide);
      if (need_parens)
        os << '(';
      denominator_->print_on(os);
      if (need_parens)
        os << ')';
    }
  }
};

struct Sqrtd0 : Expression
{
  void parse(Decoder* decoder) override;
  void simplify(Expression** self_ptr) override { assert(false); }
  bool equals(Expression* expression) const override { assert(false); }
  void replace_g(Expression** self_ptr, GFactor* gn) override { assert(false); }

  void print_on(std::ostream& os) const override
  {
    os << "√d₀";
  }
};

struct Product : Expression
{
  std::vector<Expression*> factors_;

  void divide_by_sqrt_d0() override
  {
    for (auto it = factors_.begin(); it != factors_.end(); ++it)
    {
      Sqrtd0* sqrt = dynamic_cast<Sqrtd0*>(*it);
      if (sqrt)
      {
        factors_.erase(it);
        return;
      }
    }
    // Couldn't find it.
    assert(false);
  }

  void multiply_with_sqrt_d0() override
  {
    for (auto it = factors_.begin(); it != factors_.end(); ++it)
    {
      Sqrtd0* sqrt = dynamic_cast<Sqrtd0*>(*it);
      if (sqrt)
      {
        *it = new Symbol("d0");;
        return;
      }
      Powerd0* powerd0 = dynamic_cast<Powerd0*>(*it);
      if (powerd0)
      {
        powerd0->exponent_times_two_++;
        return;
      }
    }
    // Couldn't find it.
    assert(false);
  }

  void multiply_with(int n) override
  {
    for (auto it = factors_.begin(); it != factors_.end(); ++it)
    {
      Fraction* fraction = dynamic_cast<Fraction*>(*it);
      if (fraction)
      {
        Integer* enumerator = dynamic_cast<Integer*>(fraction->enumerator_);
        Integer* denominator = dynamic_cast<Integer*>(fraction->denominator_);
        assert(enumerator && denominator && enumerator->value_ == 1 && denominator->value_ == n);
        factors_.erase(it);
        return;
      }
    }
    // Couldn't find it.
    assert(false);
  }

  void divide_by(int n) override
  {
    for (auto it = factors_.begin(); it != factors_.end(); ++it)
    {
      Integer* i = dynamic_cast<Integer*>(*it);
      if (i)
      {
        assert(i->value_ % n == 0);
        i->value_ /= n;
        return;
      }
    }
    // Couldn't find it.
    assert(false);
  }

  void remove_u() override
  {
    for (auto it = factors_.begin(); it != factors_.end(); ++it)
    {
      Poweru* power_u = dynamic_cast<Poweru*>(*it);
      Symbol* symbol = dynamic_cast<Symbol*>(*it);
      if (power_u || (symbol && symbol->name_ == "u"))
      {
        factors_.erase(it);
        return;
      }
    }
    // Couldn't find it.
    assert(false);
  }

  void simplify(Expression** self_ptr) override
  {
    assert(!factors_.empty());
    for (int n = 0; n < factors_.size(); ++n)
      factors_[n]->simplify(&factors_[n]);
    if (factors_.size() == 1)
      *self_ptr = factors_[0];
  }

  void replace_g(Expression** self_ptr, GFactor* gn) override
  {
    //std::cout << "Product::replace_g(" << *this << ", " << *gn << ")" << std::endl;

    for (int n = 0; n < factors_.size(); ++n)
      factors_[n]->replace_g(&factors_[n], gn);
    if (gn->is_equal(this))
      *self_ptr = gn;
  }

  bool equals(Expression* expression) const override
  {
    Product* product = dynamic_cast<Product*>(expression);

    if (!product || factors_.size() != product->factors_.size())
      return false;

    for (int i = 0; i < factors_.size(); ++i)
      if (!factors_[i]->equals(product->factors_[i]))
        return false;

    return true;
  }

  bool needs_parens(Operator op) const override { return op == Divide; }

  void print_on(std::ostream& os) const override
  {
    char const* separator = "";
    for (Expression* expression : factors_)
    {
      os << separator;
      bool need_parens = expression->needs_parens(Multiply);
      if (need_parens)
        os << '(';
      expression->print_on(os);
      if (need_parens)
        os << ')';
      separator = " ⋆ ";
    }
  }
};

struct Negation : Expression
{
  Expression* expression_;

  Negation(Expression* expression) : expression_(expression) { }

  void simplify(Expression** self_ptr) override
  {
    expression_->simplify(&expression_);
  }

  bool equals(Expression* expression) const override
  {
    Negation* negation = dynamic_cast<Negation*>(expression);
    if (!negation)
      return false;
    return expression_->equals(negation->expression_);
  }

  void replace_g(Expression** self_ptr, GFactor* gn) override
  {
    //std::cout << "Negation::replace_g(" << *this << ", " << *gn << ")" << std::endl;

    expression_->replace_g(&expression_, gn);
    if (gn->is_equal(this))
      *self_ptr = gn;
  }

  void print_on(std::ostream& os) const override
  {
    os << '-';
    bool need_parens = expression_->needs_parens(Subtract);
    if (need_parens)
      os << '(';
    expression_->print_on(os);
    if (need_parens)
      os << ')';
  }
};

struct Sum : Expression
{
  std::vector<Expression*> terms_;

  void simplify(Expression** self_ptr) override
  {
    assert(!terms_.empty());
    for (int n = 0; n < terms_.size(); ++n)
      terms_[n]->simplify(&terms_[n]);
    if (terms_.size() == 1)
      *self_ptr = terms_[0];
  }

  void replace_g(Expression** self_ptr, GFactor* gn) override
  {
    //std::cout << "Sum::replace_g(" << *this << ", " << *gn << ")" << std::endl;

    for (int n = 0; n < terms_.size(); ++n)
      terms_[n]->replace_g(&terms_[n], gn);
    if (gn->is_equal(this))
    {
      *self_ptr = gn;
    }
  }

  bool needs_parens(Operator op) const override { return op != Add; }

  bool equals(Expression* expression) const override
  {
    Sum* sum = dynamic_cast<Sum*>(expression);

    if (!sum || terms_.size() != sum->terms_.size())
      return false;

    for (int i = 0; i < terms_.size(); ++i)
      if (!terms_[i]->equals(sum->terms_[i]))
        return false;

    return true;
  }

  void sort_g()
  {
    std::cout << "Entering: sort_g() for \"" << *this << "\"." << std::endl;
    std::sort(terms_.begin(), terms_.end(), [](Expression* const a, Expression* const b) -> bool {
      Fraction* fraction_a = dynamic_cast<Fraction*>(a);
      if (!fraction_a)
      {
        Negation* negation_a = dynamic_cast<Negation*>(a);
        fraction_a = dynamic_cast<Fraction*>(negation_a->expression_);
      }
      Fraction* fraction_b = dynamic_cast<Fraction*>(b);
      if (!fraction_b)
      {
        Negation* negation_b = dynamic_cast<Negation*>(b);
        fraction_b = dynamic_cast<Fraction*>(negation_b->expression_);
      }
      assert(fraction_a && fraction_b);
      Product* product_a = dynamic_cast<Product*>(fraction_a->enumerator_);
      Product* product_b = dynamic_cast<Product*>(fraction_b->enumerator_);
      if (!product_a)
        std::cout << "Not a product: \"" << *fraction_a->enumerator_ << "\"!" << std::endl;
      if (!product_b)
        std::cout << "Not a product: \"" << *fraction_b->enumerator_ << "\"!" << std::endl;
      return a < b;
    });
  }

  void print_on(std::ostream& os) const override
  {
    char const* separator = "";
    for (Expression* expression : terms_)
    {
      os << separator;
      bool need_parens = expression->needs_parens(Add);
      if (need_parens)
        os << '(';
      expression->print_on(os);
      if (need_parens)
        os << ')';
      separator = " + ";
    }
  }
};

struct RowBox : Expression
{
  std::vector<Expression*> elements_;
  void parse(Decoder* decoder) override;
  void add(Expression* element) override
  {
    elements_.push_back(element);
  }
  void simplify(Expression** self_ptr) override { assert(false); }
  bool equals(Expression* expression) const override { assert(false); }
  void replace_g(Expression** self_ptr, GFactor* gn) override { assert(false); }

  void print_on(std::ostream& os) const override
  {
    char const* separator = "";
    os << '{';
    for (Expression* expression : elements_)
    {
      os << separator;
      expression->print_on(os);
      separator = ", ";
    }
    os << '}';
  }
};

struct Decoder
{
  std::ifstream file_;
  std::string line_;
  std::string_view current_line_;
  std::string_view::const_iterator head_;
  std::string_view current_token_;
  RowBox rowbox_;
  std::vector<Expression*> terms_;

  Decoder(std::string file_name) : current_line_(line_), head_(current_line_.begin())
  {
    file_.open("/home/carlo/Downloads/fadeb220-cb24-4029-af1f-52b6f6ac5085.nb");
  }

  ~Decoder()
  {
    file_.close();
  }

  void putback_other()
  {
    assert(head_ != current_line_.begin());
    --head_;
  }

  bool read_next_line()
  {
    if (!std::getline(file_, line_))
      return false;
    current_line_ = line_;
    head_ = current_line_.begin();
    return true;
  }

  bool ensure_have_data()
  {
    if (head_ != current_line_.end())
      return true;
    if (read_next_line())
      return true;
    throw std::runtime_error("Unexpected end of file");
  }

  Token next_token()
  {
    // Skip any whitespace.
    while (head_ != current_line_.end() && std::isspace(*head_))
      ++head_;

    if (head_ == current_line_.end())
    {
      if (read_next_line())
        return next_token();
      return None;
    }

    auto start = head_;

    // Check for String token.
    if (*head_ == '"')
    {
      start = ++head_; // Skip the opening quote.
      while (ensure_have_data() && *head_ != '"')
        ++head_;
      current_token_ = std::string_view(start, std::distance(start, head_));
      ++head_; // Skip the closing quote.
      return String;
    }

    // Check for Keyword token.
    if (std::isalpha(*head_))
    {
      while (ensure_have_data() && (std::isalnum(*head_) || *head_ == '_'))
        ++head_;
      current_token_ = std::string_view(start, std::distance(start, head_));
      // Check if the token is one of the specific keywords.
      if (current_token_ == "RowBox" || current_token_ == "FractionBox" || current_token_ == "SqrtBox" || current_token_ == "SuperscriptBox")
        return Keyword;
      // If it's not one of the specific keywords.
      assert(false);
    }

    // Handle Other token (single character).
    ++head_;
    current_token_ = std::string_view(start, std::distance(start, head_));
    return Other;
  }

  void parse()
  {
    // Bootstrap.
    auto type = next_token();
    assert(type == Keyword);
    assert(current_token_ == "RowBox");
    type = next_token();
    assert(type == Other);
    assert(current_token_ == "[");
    rowbox_.parse(this);
    type = next_token();
    assert(type == Other);
    assert(current_token_ == "]");

    for (int n = 2; n < rowbox_.elements_.size(); n += 2)
      terms_.push_back(rowbox_.elements_[n]);
  }
};

void Expression::parse_argument_list(Decoder* decoder, int min, int max)
{
  int arguments = 0;
  do
  {
    auto type = decoder->next_token();
    switch (type)
    {
      case Keyword:
        parse_keyword(decoder);
        break;
      case String:
        parse_string(decoder);
        break;
      default:
        // Expected Keyword or String.
        assert(false);
    }
    ++arguments;
    type = decoder->next_token();
    assert(type == Other);
  }
  while (decoder->current_token_ == ",");
  decoder->putback_other();
  assert(min <= arguments && arguments <= max);
}

void RowBox::parse(Decoder* decoder)
{
  auto type = decoder->next_token();
  assert(type == Other);
  assert(decoder->current_token_ == "{");
  parse_argument_list(decoder);
  type = decoder->next_token();
  assert(type == Other);
  assert(decoder->current_token_ == "}");
}

void Sqrtd0::parse(Decoder* decoder)
{
  auto type = decoder->next_token();
  assert(type == String);
  assert(decoder->current_token_ == "d0");
}

void Fraction::parse(Decoder* decoder)
{
  parse_argument_list(decoder, 2, 2);
}

void Expression::parse_keyword(Decoder* decoder)
{
  std::string_view keyword = decoder->current_token_;
  auto type = decoder->next_token();
  assert(type == Other);
  assert(decoder->current_token_ == "[");
  if (keyword == "RowBox")
  {
    RowBox* row_box = new RowBox;
    row_box->parse(decoder);
    assert(row_box->elements_.size() >= 2);
    StringExpr* arg0 = dynamic_cast<StringExpr*>(row_box->elements_[0]);
    int n1 = 1;
    for (;;) // So we can use break.
    {
      if (arg0)
      {
        if (arg0->value_ == "(")
        {
          ++n1;
        }
        else if (arg0->value_ == "-")
        {
          assert(row_box->elements_.size() == 2);
          Negation* negation_expr = new Negation(row_box->elements_[1]);
          add(negation_expr);
          break;
        }
        else
        {
          std::cout << "arg0 is \"" << arg0->value_ << "\"." << std::endl;
          assert(false);
        }
      }
      StringExpr* arg1 = dynamic_cast<StringExpr*>(row_box->elements_[n1]);
      if (arg1)
      {
        if (arg1->value_ == "-" || arg1->value_ == "+")
        {
          Sum* sum_expr = new Sum;
          assert((row_box->elements_.size() & 1) == 1);
          bool saw_minus = false;
          for (int n = n1 - 1; n < row_box->elements_.size(); ++n)
          {
            if (((n - n1) & 1) == 0)
            {
              StringExpr* arg = dynamic_cast<StringExpr*>(row_box->elements_[n]);
              assert(arg && arg->value_ == "+" || arg->value_ == "-" || (arg->value_ == ")" && n == row_box->elements_.size() - 1));
              saw_minus = arg->value_ == "-";
              continue;
            }
            if (saw_minus)
              sum_expr->terms_.push_back(new Negation(row_box->elements_[n]));
            else
              sum_expr->terms_.push_back(row_box->elements_[n]);
          }
          add(sum_expr);
          break;
        }
        else if (arg1->value_ == " ")
        {
          Product* product_expr = new Product;
          assert((row_box->elements_.size() & 1) == 1);
          for (int n = n1 - 1; n < row_box->elements_.size(); ++n)
          {
            if (((n - n1) & 1) == 0)
            {
              StringExpr* arg = dynamic_cast<StringExpr*>(row_box->elements_[n]);
              assert(arg && arg->value_ == " " || (arg->value_ == ")" && n == row_box->elements_.size() - 1));
              continue;
            }
            product_expr->factors_.push_back(row_box->elements_[n]);
          }
          add(product_expr);
          break;
        }
        else if (n1 == 2 && arg1->value_ == ")")
        {
          // ")" must be the last argument.
          assert(row_box->elements_.size() == 3);
          add(row_box->elements_[1]);
          break;
        }
        else
        {
          std::cout << *row_box << std::endl;
          assert(false);
        }
      }
      assert(false);
    }
  }
  else if (keyword == "FractionBox")
  {
    Expression* fraction = new Fraction;
    add(fraction);
    fraction->parse(decoder);
  }
  else if (keyword == "SqrtBox")
  {
    Expression* sqrt = new Sqrtd0;
    add(sqrt);
    sqrt->parse(decoder);
  }
  else if (keyword == "SuperscriptBox")
  {
    auto type = decoder->next_token();
    assert(type == String);
    if (decoder->current_token_ == "u")
    {
      Expression* power_u = new Poweru;
      add(power_u);
      power_u->parse(decoder);
    }
    else
    {
      char n = decoder->current_token_[1];
      assert(decoder->current_token_[0] == 'd' && std::isdigit(n));
      Expression* power;
      if (n == '0')
      {
        power = new Square(decoder->current_token_);
        try
        {
          power->parse(decoder);
        }
        catch (std::domain_error const&)
        {
          assert(n == '0');
          power = new Powerd0;
          power->parse(decoder);
        }
      }
      else
      {
        std::string symbol(decoder->current_token_);
        auto type = decoder->next_token();
        assert(type == Other);
        assert(decoder->current_token_ == ",");
        type = decoder->next_token();
        assert(type == String);
        int exponent;
        std::from_chars(decoder->current_token_.begin(), decoder->current_token_.end(), exponent);
        power = new Symbol(symbol, exponent);
      }
      add(power);
    }
  }
  type = decoder->next_token();
  assert(type == Other);
  assert(decoder->current_token_ == "]");
}

void Expression::parse_string(Decoder* decoder)
{
  if (std::isalpha(decoder->current_token_[0]))
  {
    assert(decoder->current_token_[0] == 'd' || decoder->current_token_[0] == 'u');
    Expression* symbol = new Symbol(decoder->current_token_);
    add(symbol);
    return;
  }
  if (std::isdigit(decoder->current_token_[0]))
  {
    Expression* integer = new Integer(decoder->current_token_);
    add(integer);
    return;
  }
  Expression* string = new StringExpr(decoder->current_token_);
  add(string);
}

void Square::parse(Decoder* decoder)
{
  auto type = decoder->next_token();
  assert(type == Other);
  assert(decoder->current_token_ == ",");
  type = decoder->next_token();
  if (type != String)
  {
    assert(type == Keyword && decoder->current_token_ == "RowBox");
    throw std::domain_error("RowBox");
  }
  assert(decoder->current_token_ == "2");
}

void Powerd0::parse(Decoder* decoder)
{
  auto type = decoder->next_token();
  assert(type == Other);
  assert(decoder->current_token_ == "[");
  type = decoder->next_token();
  assert(type == Other);
  assert(decoder->current_token_ == "{");
  type = decoder->next_token();
  assert(type == String);
  std::from_chars(decoder->current_token_.begin(), decoder->current_token_.end(), exponent_times_two_);
  assert((exponent_times_two_ & 1) == 1);
  type = decoder->next_token();
  assert(type == Other);
  assert(decoder->current_token_ == ",");
  type = decoder->next_token();
  assert(type == String);
  assert(decoder->current_token_ == "/");
  type = decoder->next_token();
  assert(type == Other);
  assert(decoder->current_token_ == ",");
  type = decoder->next_token();
  assert(type == String);
  assert(decoder->current_token_ == "2");
  type = decoder->next_token();
  assert(type == Other);
  assert(decoder->current_token_ == "}");
  type = decoder->next_token();
  assert(type == Other);
  assert(decoder->current_token_ == "]");
}

void Poweru::parse(Decoder* decoder)
{
  auto type = decoder->next_token();
  assert(type == Other);
  assert(decoder->current_token_ == ",");
  type = decoder->next_token();
  assert(type == String);
  std::from_chars(decoder->current_token_.begin(), decoder->current_token_.end(), exponent_);
}

int main()
{
  Decoder decoder("/home/carlo/Downloads/fadeb220-cb24-4029-af1f-52b6f6ac5085.nb");
  decoder.parse();

  std::vector<GFactor*> g;
  for (int n = 1; n < decoder.terms_.size(); ++n)
  {
    Expression*& expr = decoder.terms_[n - 1];
    expr->divide_by_sqrt_d0();
    if (n > 1)
      expr->multiply_with(n);
    expr->remove_u();
    expr->simplify(&expr);
    if (n > 1)
    {
      if (n > 2)
      {
        for (GFactor* gn : g)
          expr->replace_g(&expr, gn);
      }
      g.push_back(new GFactor(n, expr));
    }
    else
    {
      std::cout << "g" << n << " = " << *expr << "\n\n";
    }
  }
  for (int i = 0; i < g.size(); ++i)
  {
    int n = i + 2;
    GFactor* gn = g[i];
    Sum* sum = dynamic_cast<Sum*>(gn->gn_);
    assert(sum);
    sum->sort_g();
    std::cout << "g" << n << " = " << *gn->gn_ << "\n\n";
  }
}
