#include "sys.h"
#include "Algorithm.h"
#include "Function.h"
#include "cairowindow/symbolic/symbolic.h"
#ifdef CWDEBUG
#include "cwds/Restart.h"
#endif

using namespace cairowindow;

struct Function : enable_drawing::Function
{
  symbolic::Function const& function_;
  symbolic::Symbol const& symbol_;
  std::vector<symbolic::Function const*> deps_;

  template<typename... Deps>
  Function(symbolic::Symbol const& symbol, symbolic::Function const& function, Deps&&... deps) :
    function_(function), symbol_(symbol), deps_({&deps...}) { }

  double evaluate(double x) const override
  {
    symbol_ = x;
    reset_evaluation();
    return function_.evaluate();
  }

  double derivative(double x) const
  {
    symbol_ = x;
    reset_evaluation();
    return function_.derivative(symbol_).evaluate();
  }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    function_.print_on(os);
  }
#endif

 private:
  friend class Algorithm;
  void reset_evaluation() const
  {
    function_.reset_evaluation();
    for (symbolic::Function const* dep : deps_)
      dep->reset_evaluation();
  }

#ifdef CWDEBUG
  std::string to_string() const override
  {
    return function_.to_string();
  }
#endif
};

class Algorithm : public enable_drawing::Algorithm
{
 public:
  using enable_drawing::Algorithm::Algorithm;

  void enable_drawing(Function const& L, double w_min, double w_max)
  {
    enable_drawing::Algorithm::enable_drawing(L, w_min, w_max);
    L.reset_evaluation();
  }
};

int main()
{
  Debug(NAMESPACE_DEBUG::init());

  using Sample = gradient_descent::Sample;
  using Scale = gradient_descent::Scale;
  using HorizontalDirection = gradient_descent::HorizontalDirection;
  using ExtremeType = gradient_descent::ExtremeType;
  using CubicToNextSampleType = gradient_descent::CubicToNextSampleType;
  using ExtremeChain = gradient_descent::ExtremeChain;

  // Default values.
  constexpr double L_max = 100.0;
  constexpr double learning_rate = 0.125;

#if 0
  //==========================================================================
  Dout(dc::notice, "*** TEST: starting with a derivative of zero ***");
  {
    // Algorithm: [one sample, derivative is zero, hdirection is unknown]
    constexpr double w0 = 13.0;

    Algorithm gda(learning_rate, L_max);
    // Drop in with a zero derivative.
    double w = w0;
    gda(w, 50.0, 0.0);
    ASSERT(gda.algorithm_str() == "one sample, derivative is zero, hdirection is unknown");

    // EXPECTED: learning_rate was added.
    ASSERT(w == w0 + learning_rate);
  }

  //==========================================================================
  Dout(dc::notice, "*** TEST: starting with a derivative of almost zero (less than one million-th w0) ***");
  {
    // Algorithm: [one sample, derivative is zero, hdirection is unknown]
    constexpr double w0 = 13.0;

    Algorithm gda(learning_rate, L_max);
    // Drop in with an almost zero derivative.
    double w = w0;
    gda(w, 50.0, w0 * 0.999e-6);
    ASSERT(gda.algorithm_str() == "one sample, derivative is zero, hdirection is unknown");

    // EXPECTED: learning_rate was added.
    ASSERT(w == w0 + learning_rate);
  }

  //==========================================================================
  Dout(dc::notice, "*** TEST: starting with a derivative of almost zero (less than epsilon), but w0 is also zero ***");
  {
    // Algorithm: [one sample, derivative is zero, hdirection is unknown]
    constexpr double w0 = 0.0;

    Algorithm gda(learning_rate, L_max);
    // Drop in with an almost zero derivative.
    double w = w0;
    gda(w, 50.0, 0.999 * Algorithm::epsilon);
    ASSERT(gda.algorithm_str() == "one sample, derivative is zero, hdirection is unknown");

    // EXPECTED: learning_rate was added.
    ASSERT(w == w0 + learning_rate);
  }

  //==========================================================================
  Dout(dc::notice, "*** TEST: starting with a small but non-zero derivative while w0 is zero ***");
  {
    // Algorithm: [one sample, gradient descent]
    constexpr double w0 = 0.0;
    constexpr double dLdw = 1.5e-8;

    Algorithm gda(learning_rate, L_max);
    // Drop in with an almost zero derivative.
    double w = w0;
    gda(w, 50.0, dLdw);
    ASSERT(gda.algorithm_str() == "one sample, gradient descent");

    // EXPECTED: learning_rate * derivative was subtracted.
    ASSERT(w == w0 - learning_rate * dLdw);
  }

  //==========================================================================
  Dout(dc::notice, "*** TEST: starting with a derivative of zero (while hdirection is known) ***");
  std::array<HorizontalDirection, 2> hdirections = {
    HorizontalDirection::left,
    HorizontalDirection::right
  };
  for (int hdi = 0; hdi < hdirections.size(); ++hdi)
  {
    // Algorithm: [one sample, derivative is zero]
    constexpr double w0 = 13.0;
    HorizontalDirection const hdirection = hdirections[hdi];

    Dout(dc::notice, "* hdirection = " << hdirection);

    Algorithm gda(learning_rate, L_max);
    ASSERT(gda.debug_small_step() == 0.0);
    gda.debug_set_hdirection_next_extreme_type_small_step(hdirection, gda.debug_next_extreme_type(), 0.0);

    // Drop in with a zero derivative.
    double w = w0;
    gda(w, 50.0, 0.0);
    ASSERT(gda.algorithm_str() == "one sample, derivative is zero");

    // EXPECTED: hdirection * learning_rate was added.
    ASSERT(gda.debug_hdirection() == hdirection);
    ASSERT(w == w0 + static_cast<int>(hdirection) * learning_rate);
  }

  //==========================================================================
  Dout(dc::notice, "*** TEST: starting with a derivative of almost zero (less than one million-th w0) (while hdirection is known) ***");
  for (int hdi = 0; hdi < hdirections.size(); ++hdi)
  {
    // Algorithm: [one sample, derivative is zero]
    constexpr double w0 = 13.0;
    HorizontalDirection const hdirection = hdirections[hdi];

    Dout(dc::notice, "* hdirection = " << hdirection);

    Algorithm gda(learning_rate, L_max);
    gda.debug_set_hdirection_next_extreme_type_small_step(hdirection, gda.debug_next_extreme_type(), 0.0);

    // Drop in with an almost zero derivative.
    double w = w0;
    gda(w, 50.0, w0 * 0.999e-6);
    ASSERT(gda.algorithm_str() == "one sample, derivative is zero");

    // EXPECTED: hdirection * learning_rate was added.
    ASSERT(gda.debug_hdirection() == hdirection);
    ASSERT(w == w0 + static_cast<int>(hdirection) * learning_rate);
  }

  //==========================================================================
  Dout(dc::notice, "*** TEST: starting with a derivative of almost zero (less than epsilon), but w0 is also zero (while hdirection is known) ***");
  for (int hdi = 0; hdi < hdirections.size(); ++hdi)
  {
    // Algorithm: [one sample, derivative is zero]
    constexpr double w0 = 0.0;
    HorizontalDirection const hdirection = hdirections[hdi];

    Dout(dc::notice, "* hdirection = " << hdirection);

    Algorithm gda(learning_rate, L_max);
    gda.debug_set_hdirection_next_extreme_type_small_step(hdirection, gda.debug_next_extreme_type(), 0.0);

    // Drop in with an almost zero derivative.
    double w = w0;
    gda(w, 50.0, 0.999 * Algorithm::epsilon);
    ASSERT(gda.algorithm_str() == "one sample, derivative is zero");

    // EXPECTED: hdirection * learning_rate was added.
    ASSERT(gda.debug_hdirection() == hdirection);
    ASSERT(w == w0 + static_cast<int>(hdirection) * learning_rate);
  }

  //==========================================================================
  Dout(dc::notice, "*** TEST: starting with a small but non-zero derivative while w0 is zero (while hdirection is known) ***");
  for (int sign_of_dLdw = -1; sign_of_dLdw <= 1; sign_of_dLdw +=2)
    for (int hdi = 0; hdi < hdirections.size(); ++hdi)
  {
    // Algorithm: [one sample, same direction]
    constexpr double w0 = 0.0;
    HorizontalDirection const hdirection = hdirections[hdi];
    double const dLdw = sign_of_dLdw * 1.5e-8;

    Dout(dc::notice, "* hdirection = " << hdirection << "; dLdw = " << dLdw);

    Algorithm gda(learning_rate, L_max);
    gda.debug_set_hdirection_next_extreme_type_small_step(hdirection, gda.debug_next_extreme_type(), 0.0);

    // Drop in with an almost zero derivative.
    double w = w0;
    gda(w, 50.0, dLdw);
    ASSERT(gda.algorithm_str() == "one sample, same direction");

    // EXPECTED: hdirection_ * abs(learning_rate_ * dLdW) was added.
    ASSERT(gda.debug_hdirection() == hdirection);
    ASSERT(w == w0 + static_cast<int>(hdirection) * std::abs(learning_rate * dLdw));
  }

  //==========================================================================
  Dout(dc::notice, "*** TEST: starting with zero derivative but known small_step (and thus while hdirection is known) ***");
  for (int hdi = 0; hdi < hdirections.size(); ++hdi)
  {
    // Algorithm: [small step]
    constexpr double w0 = 13.0;
    constexpr double small_step = 0.12345;
    HorizontalDirection const hdirection = hdirections[hdi];

    Dout(dc::notice, "* hdirection = " << hdirection);

    Algorithm gda(learning_rate, L_max);
    gda.debug_set_hdirection_next_extreme_type_small_step(hdirection, ExtremeType::maximum, small_step);

    // Drop in with a zero derivative.
    double w = w0;
    gda(w, 50.0, 0.0);
    ASSERT(gda.algorithm_str() == "small step");

    // EXPECTED: small_step was subtracted (rather randomly, because we are looking for a maximum).
    ASSERT(utils::almost_equal(static_cast<double>(w), w0 - small_step, 10e-15));
  }

  //==========================================================================
  Dout(dc::notice, "*** TEST: first cubic ***");
  {
    symbolic::Symbol const& x = symbolic::Symbol::realize("w");
    symbolic::Symbol const& a = symbolic::Symbol::realize("a");
    symbolic::Symbol const& b = symbolic::Symbol::realize("b");
    symbolic::Symbol const& c = symbolic::Symbol::realize("c");
    symbolic::Symbol const& d = symbolic::Symbol::realize("d");
    symbolic::Function const& sL = symbolic::Function::realize("L", a + b * x + c * (x^2) + d * (x^3) );
    Function L(x, sL);

    // Generate a cubic that has its extremes at -3 and 14.
    constexpr int si = -3;
    constexpr int ti = 14;
    int const ai = 42;

    // By default assume si is the maximum.
    // Lets have two points left of the maximum, three in between, and two on the right of the minimum.
    static_assert(ti > si, "ti must be larger than si");
    static_assert((si + ti) / 2 > si + 1, "ti - si is too small");
    static_assert((si + ti) / 2 < ti - 2, "ti - si is too small");
    std::array<int, 9> xs = { si - 2, si - 1, si, si + 1, (si + ti) / 2, ti - 2, ti, ti + 2, ti + 3 };

    constexpr int i_max = 2;    //    ^^
    constexpr int i_min = 6;    //                                               ^^

    // The inflection point is at
    constexpr double inflection_point = 0.5 * (si + ti);

    for (int no_extremes = 0; no_extremes <= 1; ++no_extremes)
    {
      bool has_extremes = !no_extremes;
      for (int maximum_last = 0; maximum_last <= 1; ++ maximum_last)
      {
        bool inverted = maximum_last;                           // True if ti is the maximum.
        int di = inverted ? -2 : 2;
        int bi = (has_extremes ? 3 : -3) * di * si * ti;
        int two_ci = -3 * di * (si + ti);

        a = ai;
        b = bi;
        c = 0.5 * two_ci;
        d = di;

        for (int i0 = 0; i0 < xs.size() - 1; ++i0)
          for (int i1 = i0 + 1; i1 < xs.size(); ++i1)
          {
            CubicToNextSampleType expected;
            if (i0 < i_max || !has_extremes)
            {
              if (i1 < i_max || !has_extremes)
                expected = inverted ? CubicToNextSampleType::down : CubicToNextSampleType::up;
              else if (i1 == i_max)
                expected = inverted ? CubicToNextSampleType::left_min : CubicToNextSampleType::left_max;
              else if (i1 < i_min)
                expected = inverted ? CubicToNextSampleType::min : CubicToNextSampleType::max;
              else if (i1 == i_min)
                expected = inverted ? CubicToNextSampleType::min_left_max : CubicToNextSampleType::max_left_min;
              else
                expected = inverted ? CubicToNextSampleType::min_max : CubicToNextSampleType::max_min;
            }
            else if (i0 == i_max)
            {
              if (i1 < i_min)
                expected = inverted ? CubicToNextSampleType::right_min : CubicToNextSampleType::right_max;
              else if (i1 == i_min)
                expected = inverted ? CubicToNextSampleType::right_min_left_max : CubicToNextSampleType::right_max_left_min;
              else
                expected = inverted ? CubicToNextSampleType::right_min_max : CubicToNextSampleType::right_max_min;
            }
            else if (i0 < i_min)
            {
              if (i1 < i_min)
                expected = inverted ? CubicToNextSampleType::up : CubicToNextSampleType::down;
              else if (i1 == i_min)
                expected = inverted ? CubicToNextSampleType::left_max : CubicToNextSampleType::left_min;
              else
                expected = inverted ? CubicToNextSampleType::max : CubicToNextSampleType::min;
            }
            else if (i0 == i_min)
              expected = inverted ? CubicToNextSampleType::right_max : CubicToNextSampleType::right_min;
            else
              expected = inverted ? CubicToNextSampleType::down : CubicToNextSampleType::up;

            double w0 = xs[i0];
            double dLdw0 = L.derivative(w0);
            double learning_rate;
            if (dLdw0 == 0)                         // Derivative of first sample is zero. learning rate will be added.
              learning_rate = xs[i1] - xs[i0];
            else
              learning_rate = (xs[i0] - xs[i1]) / dLdw0;

            Algorithm gda(learning_rate, L_max);
            double w = w0;

            Dout(dc::notice, "===========================================");
            //if (!no_extremes && i0 == 2 && i1 == 6)
            //  gda.enable_drawing(L, -6.0, 25.0);

            gda(w, L(w), L.derivative(w));
            ASSERT(gda.algorithm_str() ==
                (dLdw0 == 0.0 ? "one sample, derivative is zero, hdirection is unknown" : "one sample, gradient descent"));
            // Fix floating-point round-off error related to learning_rate step.
            ASSERT(utils::almost_equal(w, (double)xs[i1], 1e-9));
            w = xs[i1];
            gda.debug_chain().find_larger(w);   // Must be called when w was changed.
            ASSERT(gda.debug_hdirection() == HorizontalDirection::undecided);
            ASSERT(gda.debug_next_extreme_type() == ExtremeType::unknown);

            // The center between the two samples is
            double center = 0.5 * (w0 + w);
            // The expected "next extreme" should be set to minimum if the cubic has no extremes,
            // otherwise the extreme that is closest to the center of the two samples.
            ExtremeType expected_next_extreme_type;
            if (!has_extremes)
              expected_next_extreme_type = ExtremeType::minimum;
            else if (inverted)
              expected_next_extreme_type = center > inflection_point ? ExtremeType::maximum : ExtremeType::minimum;
            else
              expected_next_extreme_type = center < inflection_point ? ExtremeType::maximum : ExtremeType::minimum;

            Dout(dc::notice, "-------------------------------------------");
            RESTART
            gda(w, L(w), L.derivative(w));

            bool target_extreme_is_first_extreme = 2 * (xs[i0] + xs[i1]) - maximum_last < 22;
            int target_extreme = target_extreme_is_first_extreme ? si : ti;
            int i = target_extreme_is_first_extreme ? i_max : i_min;
            if (no_extremes || (i0 != i && i1 != i))
              ASSERT(gda.debug_hdirection() == HorizontalDirection::undecided);
            else if (!no_extremes && i0 == 2 && i1 == 6)
            {
              // In this case the algorithm is asking for an extra sample to cut the current local extreme cubic in half.
              ASSERT(gda.debug_hdirection() == HorizontalDirection::undecided);
              ASSERT(w == inflection_point);
              expected = maximum_last ? CubicToNextSampleType::left_max : CubicToNextSampleType::right_max;
              // Call gda again to see what happens next.
              gda(w, L(w), L.derivative(w));
            }
            else
              ASSERT(gda.debug_hdirection() == (i0 == i ? HorizontalDirection::left : HorizontalDirection::right));
            // Because w > w0, the cubic through both is stored in the SampleNode of w0 (last() returns the one of w).
            auto cubic_node = gda.debug_cubic_used();
            ASSERT(cubic_node->type() == expected);
            ASSERT(gda.debug_next_extreme_type() == expected_next_extreme_type);
          }
      }
    }
  }
#endif

#if 0
  //==========================================================================
  Dout(dc::notice, "*** TEST: first local extreme ***");
  {
    constexpr double learning_rate = 0.001;

    symbolic::Symbol const& x = symbolic::Symbol::realize("w");
    symbolic::Symbol const& a = symbolic::Symbol::realize("a");
    symbolic::Symbol const& b = symbolic::Symbol::realize("b");
    symbolic::Symbol const& c = symbolic::Symbol::realize("c");
    symbolic::Symbol const& d = symbolic::Symbol::realize("d");
    symbolic::Function const& sL = symbolic::Function::realize("L", a + b * x + c * (x^2) + d * (x^3) );
    Function L(x, sL);

    // Generate a cubic that has its extremes at -3 and 14.
    constexpr int si = -3;
    constexpr int ti = 14;
    int const ai = 42;

    // By default assume si is the maximum.
    // Lets have two points left of the maximum, three in between, and two on the right of the minimum.
    static_assert(ti > si, "ti must be larger than si");
    static_assert((si + ti) / 2 > si + 1, "ti - si is too small");
    static_assert((si + ti) / 2 < ti - 2, "ti - si is too small");
    std::array<int, 9> xs = { si - 2, si - 1, si, si + 1, (si + ti) / 2, ti - 2, ti, ti + 2, ti + 3 };

    constexpr int i_max = 2;    //    ^^
    constexpr int i_min = 6;    //                                               ^^

    // The inflection point is at
    constexpr double inflection_point = 0.5 * (si + ti);

    for (int no_extremes = 0; no_extremes <= 1; ++no_extremes)
    {
      bool has_extremes = !no_extremes;
      for (int maximum_last = 0; maximum_last <= 1; ++ maximum_last)
      {
        bool inverted = maximum_last;                           // True if ti is the maximum.
        int di = inverted ? -2 : 2;
        int bi = (has_extremes ? 3 : -3) * di * si * ti;
        int two_ci = -3 * di * (si + ti);

        a = ai;
        b = bi;
        c = 0.5 * two_ci;
        d = di;

        for (int i0 = 0; i0 < xs.size() - 1; ++i0)
          for (int i1 = i0 + 1; i1 < xs.size(); ++i1)
          {
            double w0 = xs[i0];

            Algorithm gda(learning_rate, L_max);

            Dout(dc::notice, "===========================================");
            //gda.enable_drawing(L, -15.0, 20.0);

            double w = w0;
            gda(w, L(w), L.derivative(w));

            w = xs[i1];
            gda.debug_chain().find_larger(w);
            Dout(dc::notice, "-------------------------------------------");

            int count = 0;
            [[maybe_unused]] bool not_finished = gda(w, L(w), L.derivative(w));
            ASSERT(not_finished);

            if (no_extremes)
            {
              Dout(dc::notice, maximum_last << ", " << i0 << ", " << i1 << ", w = " << w);
              if (++count == 2)
                break;
            }
            else
            {
              // The extreme picked depends on which extreme is closer to center of the two samples.
              // Normally that would be: 0.5 * (xs[i0] + xs[i1]) < 5.5 (the inflection point), but multiply that by four and add maximum_last into the mix to deal with it being equal.
              bool target_extreme_is_first_extreme = 2 * (xs[i0] + xs[i1]) - maximum_last < 22;
              int target_extreme = target_extreme_is_first_extreme ? si : ti;
              int i = target_extreme_is_first_extreme ? i_max : i_min;

              if (i0 == i || i1 == i)   // Did we reuse one of the first two samples for the extreme?
              {
                // In that case last() was set to the reused sample, which is the local extreme that was found.
                double last_w = gda.debug_chain().last()->w();
                ASSERT(utils::almost_equal(last_w, double(target_extreme), 1e-9));
              }
              else
                ASSERT(utils::almost_equal(w, double(target_extreme), 1e-9));   // w is the local extreme.
            }
          }
      }
    }
  }
#endif

#if 1
  //==========================================================================
  Dout(dc::notice, "*** TEST: parabola connected to dampened sin ***");
  {
    constexpr double w0 = -48.5;//-10.0;
    constexpr double learning_rate = 0.001;
    constexpr double L_max = 2649;

    Algorithm gda(learning_rate, L_max);
    double w = w0;

    //gda.debug_set_hdirection_next_extreme_type_small_step(HorizontalDirection::right, ExtremeType::maximum, 1.0);

    symbolic::Symbol const& x = symbolic::Symbol::realize("w");
    symbolic::Constant const& tp = symbolic::Constant::realize(55);
    symbolic::Function const& sigmoid = symbolic::Function::realize("sigmoid", exp(3 * (x + tp)) / (1 + exp(3 * (x + tp))));
    symbolic::Constant const& tp2 = symbolic::Constant::realize(80);
    symbolic::Function const& sigmoid2 = symbolic::Function::realize("sigmoid", exp(3 * (x + tp2)) / (1 + exp(3 * (x + tp2))));
    symbolic::Constant const& a = symbolic::Constant::realize(146, 10);
    symbolic::Constant const& b = symbolic::Constant::realize(315, 100);
    symbolic::Constant const& c = symbolic::Constant::realize(451, 1000);
    symbolic::Constant const& d = symbolic::Constant::realize(5, 10);
    symbolic::Constant const& e = symbolic::Constant::realize(-6, 1);
    symbolic::Constant const& amplitude = symbolic::Constant::realize(12027, 1000000);
    symbolic::Constant const& level = symbolic::Constant::realize(187838, 100);
    symbolic::Constant const& phase = symbolic::Constant::realize(191892, 100000);
    //symbolic::Function const& sL = symbolic::Function::realize("L",
    //   (1.0 - sigmoid2) * (2800 + e * ((x + 85)^2)) + sigmoid2 * ((1.0 - sigmoid) * (a + b * x + c * (x^2)) + (sigmoid * (amplitude * exp((tp - x) / 10) * sin(d * x + phase) + level))));
    symbolic::Function const& sL = symbolic::Function::realize("L",
       ((1.0 - sigmoid) * (a + b * x + c * (x^2)) + (sigmoid * (amplitude * exp((tp - x) / 10) * sin(d * x + phase) + level))));
    //Function L(x, sL, sigmoid, sigmoid2);
    Function L(x, sL, sigmoid);

    double zoom = 10.0;
    //gda.enable_drawing(L, -51.3629 - zoom, -51.3629 + zoom);
    //gda.enable_drawing(L, -55.7732647023600875968 - zoom, -55.7732647023600875968 + zoom);
    gda.enable_drawing(L, -90.0, 0.0);

#if 0
    gda(w, L(w), L.derivative(w));
    gda(w, L(w), L.derivative(w));
    gda(w, L(w), L.derivative(w));
    w = -78.0;
    gda.debug_chain().find_larger(w);
#endif

    while (gda(w, L(w), L.derivative(w)))
    {
      Dout(dc::notice, "-------------------------------------------");
    }

    ASSERT(gda.success());
    Sample const& result = gda.minimum();
    Dout(dc::notice, "Global minimum: " << result);
  }
#endif

  Dout(dc::notice, "Success!");
}
