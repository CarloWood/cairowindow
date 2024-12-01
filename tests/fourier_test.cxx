#include "sys.h"
#include "Interval.h"
#include "Function.h"
#include "Algorithm.h"
#include <cmath>
#include <random>
#include <vector>
#include <chrono>
#include <utility>
#ifdef CWDEBUG
#include "cwds/Restart.h"
#include "utils/has_print_on.h"
#include "utils/print_using.h"
#include "debug_ostream_operators.h"
#include "cwds/Restart.h"
#endif

#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

using Range = cairowindow::Range;
using Point = cairowindow::Point;
using ExtremeType = gradient_descent::ExtremeType;

class TestFunctionGenerator : public enable_drawing::Function
{
 public:
  using IntervalList = intervallist::IntervalList;
  using Interval = intervallist::Interval;

 private:
  // Filled by constructor.
  unsigned int seed_;
  std::mt19937 generator_;
  double max_frequency_;
  std::vector<double> amplitudes_;
  std::vector<double> frequencies_;
  std::vector<double> phases_;

  // Changed by an optional call to normalize_amplitude.
  double vertical_shift_{0.0};
  double amplitude_scale_{1.0};

  // Filled/used by call to find_extrema.
  utils::UniqueIDContext<int> id_context_;
  ExtremeType first_extreme_;
  ExtremeType last_extreme_;
  std::vector<Point> extrema_;
  double lowest_minimum_;
  double lowest_maximum_;
  double highest_minimum_;
  double highest_maximum_;

  // Filled by call to advanced_normalize.
  double low_low_band_;
  double high_low_band_;
  double low_high_band_;
  double high_high_band_;
  bool advanced_normalization_{false};

 public:
  TestFunctionGenerator(int num_components, Range frequency_range, Range amplitude_range, unsigned int seed = 0) :
    max_frequency_(frequency_range.max()), seed_(seed)
  {
    DoutEntering(dc::notice, "TestFunctionGenerator(" << num_components << ", " << frequency_range << ", " <<
        amplitude_range << ", " << seed << ")");

    if (seed == 0)
      seed_ = std::chrono::system_clock::now().time_since_epoch().count();

    Dout(dc::notice, "seed: " << seed_);
    generator_.seed(seed_);

    std::uniform_real_distribution<> freq_dist(frequency_range.min(), frequency_range.max());
    std::uniform_real_distribution<> amp_dist(amplitude_range.min(), amplitude_range.max());
    std::uniform_real_distribution<> phase_dist(0, 2 * M_PI);

    frequencies_.push_back(max_frequency_);
    amplitudes_.push_back(amp_dist(generator_));
    phases_.push_back(phase_dist(generator_));
    for (int i = 1; i < num_components; ++i)
    {
      frequencies_.push_back(freq_dist(generator_));
      amplitudes_.push_back(amp_dist(generator_));
      phases_.push_back(phase_dist(generator_));
    }
  }

  unsigned int get_seed() const { return seed_; }
  utils::UniqueIDContext<int>& id_context() { return id_context_; }
  double max_frequency() const { return max_frequency_; }

  int sign_of_derivative(double x) const
  {
    return derivative(x) < 0.0 ? -1 : 1;
  }

#ifdef CWDEBUG
  std::string to_string() const override
  {
    return "test function (seed: " + std::to_string(seed_) + ")";
  }

#endif

  void normalize_amplitude(Range desired_amplitude_range, Range x_range, int num_samples = 1000)
  {
    Range current{std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest()};

    for (int i = 0; i < num_samples; ++i)
    {
      double x = x_range.min() + i * x_range.size() / (num_samples - 1);
      double y = evaluate(x);
      current = Range{std::min(current.min(), y), std::max(current.max(), y)};
    }

    amplitude_scale_ = desired_amplitude_range.size() / current.size();
    vertical_shift_ = desired_amplitude_range.center() - amplitude_scale_ * current.center();
  }

  void calculate_alpha_beta(std::vector<Point>::const_iterator right, double& alpha, double& beta) const
  {
    // The extreme on the left of x.
    auto left = std::prev(right);

    // Let f_min be the minimum value of f: f(0).
    // Let f_max be the maximum value of f: f(1).
    double f_min = left->y();
    double f_max = right->y();
    if (f_min > f_max)
      std::swap(f_min, f_max);

    // Linearly transform input_range --> output_range. If the input_range is empty then return output_range.min().
    auto transform_range = [](double input, Range input_range, Range output_range){
      if (input_range.size() == 0.0)
        return output_range.min();
      return output_range.min() + (input - input_range.min()) / input_range.size() * output_range.size();
    };

    // Let g_min be the required output in the minimum at t = 0.
    // Let g_max be the required output in the maximum at t = 1.
    double g_min = transform_range(f_min, {lowest_minimum_, highest_minimum_}, {low_low_band_, high_low_band_});
    double g_max = transform_range(f_max, {highest_maximum_, lowest_maximum_}, {high_high_band_, low_high_band_});
    //
    // g(t) = alpha * f(t) + beta  -->
    //
    //  g(0) = g_min = alpha * f(0) + beta = alpha * f_min + beta
    //  g(1) = g_max = alpha * f(1) + beta = alpha * f_max + beta
    //        --------------------------------------------------- -
    // g_min - g_max = alpha * (f_min - f_max) -->
    //
    alpha = (g_max - g_min) / (f_max - f_min);
    beta = g_min - alpha * f_min;
  }

  double evaluate(double x) const override
  {
    double sum_of_sins = 0.0;
    for (size_t i = 0; i < amplitudes_.size(); ++i)
      sum_of_sins += amplitudes_[i] * std::sin(frequencies_[i] * x + phases_[i]);

    double y = sum_of_sins * amplitude_scale_ + vertical_shift_;

    if (advanced_normalization_)
    {
      // Find the first extreme on the right of x.
      auto right = std::upper_bound(extrema_.begin(), extrema_.end(), x, [](double value, Point const& point) { return value < point.x(); });
      double alpha{}, beta{};

      bool before_first_extreme = right == extrema_.begin();
      // x is before the first extreme. Use a transformation equal to that of the first interval.
      if (before_first_extreme)
        ++right;
      // If x is past the last extreme then leave alpha and beta as they were in the last interval.
      if (right != extrema_.end())
        calculate_alpha_beta(right, alpha, beta);

      // Apply advanced normalization.
      y = alpha * y + beta;
    }

    return y;
  }

  double derivative(double x) const
  {
    double sum_of_cos = 0.0;
    for (size_t i = 0; i < amplitudes_.size(); ++i)
      sum_of_cos += amplitudes_[i] * frequencies_[i] * std::cos(frequencies_[i] * x + phases_[i]);

    double dxdy = amplitude_scale_ * sum_of_cos;

    if (advanced_normalization_)
    {
      // Find the first extreme on the right of x.
      auto right = std::upper_bound(extrema_.begin(), extrema_.end(), x, [](double value, Point const& point) { return value < point.x(); });
      double alpha{}, beta{};

      bool before_first_extreme = right == extrema_.begin();
      // x is before the first extreme. Use a transformation equal to that of the first interval.
      if (before_first_extreme)
        ++right;
      // If x is past the last extreme then leave alpha and beta as they were in the last interval.
      if (right != extrema_.end())
        calculate_alpha_beta(right, alpha, beta);

      // Apply advanced normalization.
      dxdy *= alpha;
    }

    return dxdy;
  }

  // Finds all local extrema in x_range and puts them in extrema_.
  void find_extrema(Range x_range)
  {
    advanced_normalization_ = false;

    IntervalList intervals;
    intervals.push_back(new Interval(
          {x_range.min(), sign_of_derivative(x_range.min())},
          {x_range.max(), sign_of_derivative(x_range.max())},
          intervals.id_context()));
    Interval const* interval = intervals.front();
    first_extreme_ = interval->begin_sign() == -1 ? ExtremeType::minimum : ExtremeType::maximum;
    last_extreme_ = interval->end_sign() == -1 ? ExtremeType::minimum : ExtremeType::maximum;

    // Emperically found constant.
    double const number_of_local_extrema_estimate = 0.255 * x_range.size() * max_frequency();
    // In order to make at least this many cuts we need a minimum depth of,
    int const max_depth = std::ceil(std::log2(number_of_local_extrema_estimate));

    int local_extrema;
    bool had_divide;
    int depth = 0;
    do
    {
      had_divide = false;
      Interval* next;
      for (Interval* interval = intervals.front(); !intervals.is_root(interval); interval = next)
      {
        next = interval->next();
        if (depth < max_depth || interval->must_be_divided(intervals))
        {
          double mid_x = 0.5 * (interval->x_range_begin() + interval->x_range_end());
          interval->split({mid_x, sign_of_derivative(mid_x)}, intervals);
          had_divide = true;
        }
      }
      ++depth;
    }
    while (had_divide);
    // Now intervals contains sign-changing intervals around every local extreme.

    struct Bound {
      double x;
      double y;
      double dxdy;
    };

    double const tolerance = 1e-5 * x_range.size();
    math::CubicPolynomial cubic;

    lowest_minimum_ = std::numeric_limits<double>::max();
    lowest_maximum_ = std::numeric_limits<double>::max();
    highest_minimum_ = std::numeric_limits<double>::lowest();
    highest_maximum_ = std::numeric_limits<double>::lowest();

    // Loop over all intervals with a sign change and find the corresponding local extreme.
    for (Interval const* interval = intervals.front(); !intervals.is_root(interval); interval = interval->next())
    {
      // Skip intervals where the derivative does not change sign.
      if (!interval->contains_sign_change())
        continue;

      Dout(dc::notice, "Determining the extreme between " << *interval << ".");

      using namespace gradient_descent;

      ExtremeType const extreme_type = interval->begin_sign() == -1 ? ExtremeType::minimum : ExtremeType::maximum;

      std::array<Bound, 2> points = {{
        {interval->x_range_begin(), evaluate(interval->x_range_begin()), derivative(interval->x_range_begin())},
        {interval->x_range_end(), evaluate(interval->x_range_end()), derivative(interval->x_range_end())}
      }};

      for (;;)
      {
//        for (int i = 0; i < points.size(); ++i)
//          Dout(dc::notice, "f(" << points[i].x << ") = " << points[i].y << ", f'(" << points[i].x << ") = " << points[i].dxdy);

        // Fit a cubic through both points and find the appropriate extreme.
        cubic.initialize(points[0].x, points[0].y, points[0].dxdy, points[1].x, points[1].y, points[1].dxdy);
//        Dout(dc::notice, "f(x) = " << cubic);
        AnalyzedCubic acubic;
        acubic.initialize(cubic, extreme_type);
        double x_new = acubic.get_extreme();
        double dxdy_new = derivative(x_new);

        // Determine which point is furthest from x_new.
        int furthest = (std::abs(x_new - points[0].x) > std::abs(x_new - points[1].x)) ? 0 : 1;
        int closest = 1 - furthest;

        if (std::abs(x_new - points[closest].x) < tolerance)
        {
          double y_new = evaluate(x_new);
          extrema_.emplace_back(x_new, y_new);
          Dout(dc::notice, "Extreme: x = " << x_new << ", y = " << y_new << ", f'(x) = " << dxdy_new);
          if (extreme_type == ExtremeType::minimum)
          {
            lowest_minimum_ = std::min(lowest_minimum_, y_new);
            highest_minimum_ = std::max(highest_minimum_, y_new);
          }
          else
          {
            lowest_maximum_ = std::min(lowest_maximum_, y_new);
            highest_maximum_ = std::max(highest_maximum_, y_new);
          }
          break;
        }

        // If derivative sign is the same as the closest one, move the new point 1% of the
        // current shift distance (the distance between furthest and new) away from the
        // closest point; this in order to avoid the two points becoming too close.
        if (std::signbit(dxdy_new) == std::signbit(points[closest].dxdy))
        {
          x_new += (x_new > points[closest].x ? 0.01 : -0.01) * std::abs(x_new - points[furthest].x);
          dxdy_new = derivative(x_new);
        }

        // Update the furthest point.
        points[furthest] = {x_new, evaluate(x_new), dxdy_new};
      }
    }

    Dout(dc::notice, "This function has " << extrema_.size() << " local extrema in the interval " << x_range << '.');
  }

  void advanced_normalize(double min_fraction, double max_fraction, Range x_range)
  {
    // Find extrema of the original function.
    find_extrema(x_range);

    // Calculate the bands
    low_low_band_ = lowest_minimum_;
    high_low_band_ = lowest_minimum_ + min_fraction * (highest_maximum_ - lowest_minimum_);
    low_high_band_ = highest_maximum_ - max_fraction * (highest_maximum_ - lowest_minimum_);
    high_high_band_ = highest_maximum_;

    Dout(dc::notice, "Advanced normalization applied with min_fraction = " << min_fraction << " and max_fraction = " << max_fraction);

    // Use advanced normalization.
    advanced_normalization_ = true;
  }
};

using namespace cairowindow;
using Algorithm = enable_drawing::Algorithm;

int main(int argc, char* argv[])
{
  Debug(NAMESPACE_DEBUG::init());

  constexpr double w0 = 5.5;
  constexpr double learning_rate = 0.0001;
  constexpr double L_max = 600;

  unsigned int seed = 0;
  if (argc > 1) {
    seed = std::stoul(argv[1]);
  }

  int const number_of_frequencies = 20;
  Range const frequency_range{0.5, 10.0};
  Range const amplitude_range{25.0, 100.0};
  Range const x_range{0.0, 8.0};

  do
  {
    Algorithm gda(learning_rate, L_max);
    double w = w0;

    TestFunctionGenerator L(number_of_frequencies, frequency_range, amplitude_range, seed);
    L.advanced_normalize(0.1, 0.1, x_range);

    gda.enable_drawing(L, x_range.min(), x_range.max());
//    gda.enable_drawing(L, 7.4622, 7.4623, -50.0, 350.0);        // For seed 2784233688
    while (gda(w, L(w), L.derivative(w)))
    {
      Dout(dc::notice, "-------------------------------------------");
    }

    ASSERT(gda.success());

    gradient_descent::Sample const& result = gda.minimum();
    Dout(dc::notice, "Global minimum: " << result);
  }
  while (seed == 0);
}
