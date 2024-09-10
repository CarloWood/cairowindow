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

class TestFunctionGenerator : public enable_drawing::Function
{
 public:
  using IntervalList = intervallist::IntervalList;
  using Interval = intervallist::Interval;

 private:
  double max_frequency_;
  std::vector<double> amplitudes_;
  std::vector<double> frequencies_;
  std::vector<double> phases_;
  double vertical_shift_;
  double amplitude_scale_;
  unsigned int seed_;
  std::mt19937 generator_;
  std::vector<std::pair<double, double>> minima_;
  std::vector<std::pair<double, double>> maxima_;
  utils::UniqueIDContext<int> id_context_;
  IntervalList intervals_;

 public:
  TestFunctionGenerator(int num_components, Range frequency_range, Range amplitude_range, unsigned int seed = 0) :
    max_frequency_(frequency_range.max()), vertical_shift_(0.0), amplitude_scale_(1.0), seed_(seed)
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
  IntervalList& intervals() { return intervals_; }
  IntervalList const& intervals() const { return intervals_; }
  utils::UniqueIDContext<int>& id_context() { return id_context_; }
  double max_frequency() const { return max_frequency_; }

  int sign_of_derivative(double x) const
  {
    return derivative(x) < 0.0 ? -1 : 1;
  }

  double evaluate(double x) const override
  {
    double result = 0.0;
    for (size_t i = 0; i < amplitudes_.size(); ++i)
      result += amplitudes_[i] * std::sin(frequencies_[i] * x + phases_[i]);
    return result * amplitude_scale_ + vertical_shift_;
  }

  std::string to_string() const override
  {
    return "test function (seed: " + std::to_string(seed_) + ")";
  }

  double derivative(double x) const
  {
    double result = 0;
    for (size_t i = 0; i < amplitudes_.size(); ++i)
      result += amplitudes_[i] * frequencies_[i] * std::cos(frequencies_[i] * x + phases_[i]);
    return amplitude_scale_ * result;
  }

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

  void find_extrema(Range x_range)
  {
    minima_.clear();
    maxima_.clear();

    intervals_.push_back(new Interval(
          {x_range.min(), sign_of_derivative(x_range.min())},
          {x_range.max(), sign_of_derivative(x_range.max())},
          intervals_.id_context()));

    // Emperically found.
    double const number_of_local_extrema_estimate = 0.255 * x_range.size() * max_frequency();
    // In order to make at least this many cuts we need a minimum depth of,
    int const max_depth = std::ceil(std::log2(number_of_local_extrema_estimate));

    int local_extremes;
    bool had_divide;
    int depth = 0;
    do
    {
      had_divide = false;
      Interval* next;
      for (Interval* interval = intervals_.front(); !intervals_.is_root(interval); interval = next)
      {
        next = interval->next();
        if (depth < max_depth || interval->must_be_divided(intervals_))
        {
          double mid_x = 0.5 * (interval->x_range_begin() + interval->x_range_end());
          interval->split({mid_x, sign_of_derivative(mid_x)}, intervals_);
          had_divide = true;
        }
      }
      ++depth;
    }
    while (had_divide);
    // Now intervals contains sign-changing intervals around every local extreme.

    struct Point {
      double x;
      double y;
      double dxdy;
    };

    double const tolerance = 1e-5 * x_range.size();
    math::CubicPolynomial cubic;

    // Loop over all intervals with a sign change and find the corresponding local extreme.
    for (Interval const* interval = intervals_.front(); !intervals_.is_root(interval); interval = interval->next())
    {
      // Skip intervals where the derivative does not change sign.
      if (!interval->contains_sign_change())
        continue;

      Dout(dc::notice, "Determining the extreme between " << *interval << ".");

      using namespace gradient_descent;

      ExtremeType const extreme_type = interval->begin_sign() == -1 ? ExtremeType::minimum : ExtremeType::maximum;

      std::array<Point, 2> points = {{
        {interval->x_range_begin(), evaluate(interval->x_range_begin()), derivative(interval->x_range_begin())},
        {interval->x_range_end(), evaluate(interval->x_range_end()), derivative(interval->x_range_end())}
      }};

      for (;;)
      {
        for (int i = 0; i < points.size(); ++i)
          Dout(dc::notice, "f(" << points[i].x << ") = " << points[i].y << ", f'(" << points[i].x << ") = " << points[i].dxdy);

        // Fit a cubic through both points and find the appropriate extreme.
        cubic.initialize(points[0].x, points[0].y, points[0].dxdy, points[1].x, points[1].y, points[1].dxdy);
        Dout(dc::notice, "f(x) = " << cubic);
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
          if (extreme_type == ExtremeType::minimum)
            minima_.emplace_back(x_new, y_new);
          else
            maxima_.emplace_back(x_new, y_new);
          Dout(dc::notice, "Extreme: x = " << x_new << ", y = " << y_new << ", f'(x) = " << dxdy_new);
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

        // Calculate the new y and dx/dy.
        double y_new = evaluate(x_new);
        Dout(dc::notice, "x_new = " << x_new << ", y_new = " << y_new << ", f'(x_new) = " << dxdy_new);

        // Update the furthest point.
        points[furthest] = {x_new, y_new, dxdy_new};
      }
    }

    Dout(dc::notice, "This function has " << minima_.size() << " minima and " << maxima_.size() << " maxima.");
  }

  void advanced_normalize(double min_percentage, double max_percentage, Range x_range)
  {
    // First, apply regular normalization to [-1, 1].
    //normalize_amplitude({-1, 1}, x_range, 2000);

    // Find extrema.
    find_extrema(x_range);

    // Calculate desired ranges for minima and maxima.
    double min_low = -1.0;
    double min_high = -1.0 + 2.0 * min_percentage;
    double max_low = 1.0 - 2.0 * max_percentage;
    double max_high = 1.0;

    // Create transformation function.
    auto transform = [&](double x) -> double {
      auto it_min = std::lower_bound(minima_.begin(), minima_.end(), std::pair<double, double>(x, 0));
      auto it_max = std::lower_bound(maxima_.begin(), maxima_.end(), std::pair<double, double>(x, 0));

      double dist_to_min_left = (it_min == minima_.begin()) ? std::numeric_limits<double>::max() : x - std::prev(it_min)->first;
      double dist_to_min_right = (it_min == minima_.end()) ? std::numeric_limits<double>::max() : it_min->first - x;
      double dist_to_max_left = (it_max == maxima_.begin()) ? std::numeric_limits<double>::max() : x - std::prev(it_max)->first;
      double dist_to_max_right = (it_max == maxima_.end()) ? std::numeric_limits<double>::max() : it_max->first - x;

      double min_dist_to_min = std::min(dist_to_min_left, dist_to_min_right);
      double min_dist_to_max = std::min(dist_to_max_left, dist_to_max_right);

      bool closer_to_min = min_dist_to_min <= min_dist_to_max;

      double left_x, left_y, right_x, right_y;
      double left_target, right_target;

      if (closer_to_min)
      {
        // Handle the case when x is closer to a minimum.
        if (dist_to_min_right <= dist_to_min_left)
        {
          right_x = it_min->first;
          right_y = it_min->second;
          right_target = min_low + (min_high - min_low) * (right_y - min_low) / 2.0;

          if (it_min != minima_.begin())
          {
            left_x = std::prev(it_min)->first;
            left_y = std::prev(it_min)->second;
            left_target = min_low + (min_high - min_low) * (left_y - min_low) / 2.0;
          }
          else
          {
            left_x = x_range.min();
            left_y = evaluate(left_x);
            left_target = left_y;
          }
        }
        else
        {
          left_x = std::prev(it_min)->first;
          left_y = std::prev(it_min)->second;
          left_target = min_low + (min_high - min_low) * (left_y - min_low) / 2.0;

          if (it_min != minima_.end())
          {
            right_x = it_min->first;
            right_y = it_min->second;
            right_target = min_low + (min_high - min_low) * (right_y - min_low) / 2.0;
          }
          else
          {
            right_x = x_range.max();
            right_y = evaluate(right_x);
            right_target = right_y;
          }
        }
      }
      else
      {
        // Handle the case when x is closer to a maximum.
        if (dist_to_max_right <= dist_to_max_left)
        {
          right_x = it_max->first;
          right_y = it_max->second;
          right_target = max_high - (max_high - max_low) * (max_high - right_y) / 2.0;

          if (it_max != maxima_.begin())
          {
            left_x = std::prev(it_max)->first;
            left_y = std::prev(it_max)->second;
            left_target = max_high - (max_high - max_low) * (max_high - left_y) / 2.0;
          }
          else
          {
            left_x = x_range.min();
            left_y = evaluate(left_x);
            left_target = left_y;
          }
        }
        else
        {
          left_x = std::prev(it_max)->first;
          left_y = std::prev(it_max)->second;
          left_target = max_high - (max_high - max_low) * (max_high - left_y) / 2.0;

          if (it_max != maxima_.end())
          {
            right_x = it_max->first;
            right_y = it_max->second;
            right_target = max_high - (max_high - max_low) * (max_high - right_y) / 2.0;
          }
          else
          {
            right_x = x_range.max();
            right_y = evaluate(right_x);
            right_target = right_y;
          }
        }
      }

      // Linear interpolation between surrounding extrema.
      double t = (x - left_x) / (right_x - left_x);
      double y = evaluate(x);
      double y_normalized = left_target + t * (right_target - left_target);

      return y_normalized;
    };

    // Apply transformation.
    //auto original_function = [this](double x) { return evaluate(x); };
    //*this = TestFunctionGenerator([transform, original_function](double x) { return transform(original_function(x)); });

    Dout(dc::notice, "Advanced normalization applied with min_percentage = " << min_percentage <<
        " and max_percentage = " << max_percentage);
  }
};

using namespace cairowindow;
using Algorithm = enable_drawing::Algorithm;

int main(int argc, char* argv[])
{
  Debug(NAMESPACE_DEBUG::init());

  constexpr double w0 = M_PI;
  constexpr double learning_rate = 0.0001;
  constexpr double L_max = 100;

  Algorithm gda(learning_rate, L_max);
  double w = w0;

  unsigned int seed = 0;
  if (argc > 1) {
    seed = std::stoul(argv[1]);
  }

  int const number_of_frequencies = 20;
  Range const frequency_range{1.0, 20.0};
  Range const amplitude_range{50.0, 100.0};
  Range const x_range{0.0, 2.0 * M_PI};

  TestFunctionGenerator L(number_of_frequencies, frequency_range, amplitude_range, seed);
  L.advanced_normalize(0.1, 0.1, x_range);

  gda.enable_drawing(L, x_range.min(), x_range.max());
//  gda.enable_drawing(L, 2.0, 3.75);
  while (gda(w, L(w), L.derivative(w)))
  {
    Dout(dc::notice, "-------------------------------------------");
  }

  ASSERT(gda.success());
  gradient_descent::Sample const& result = gda.minimum();
  Dout(dc::notice, "Global minimum: " << result);
}
