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

    generator_.seed(seed_);

    std::uniform_real_distribution<> freq_dist(frequency_range.min(), frequency_range.max());
    std::uniform_real_distribution<> amp_dist(amplitude_range.min(), amplitude_range.max());
    std::uniform_real_distribution<> phase_dist(0, 2 * M_PI);

    for (int i = 0; i < num_components; ++i)
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

  int sign_of_derivative(double x) const
  {
    return derivative(x) < 0.0 ? -1 : 1;
  }

  double operator()(double x) const override
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
      double y = operator()(x);
      current = Range{std::min(current.min(), y), std::max(current.max(), y)};
    }

    amplitude_scale_ = desired_amplitude_range.size() / current.size();
    vertical_shift_ = desired_amplitude_range.center() - amplitude_scale_ * current.center();
  }

  int find_extrema(Range x_range)
  {
    minima_.clear();
    maxima_.clear();

    intervals_.push_back(new Interval(
          {x_range.min(), sign_of_derivative(x_range.min())},
          {x_range.max(), sign_of_derivative(x_range.max())},
          intervals_.id_context()));

    int local_extremes;
    bool had_divide;
    do
    {
      had_divide = false;
      Interval* next;
      local_extremes = 0;
      for (Interval* interval = intervals_.front(); !intervals_.is_root(interval); interval = next)
      {
        next = interval->next();
        if (interval->must_be_divided())
        {
          double mid_x = 0.5 * (interval->x_range_begin() + interval->x_range_end());
          interval->split({mid_x, sign_of_derivative(mid_x)}, intervals_);
          had_divide = true;
        }
        if (interval->begin_sign() != interval->end_sign())
          ++local_extremes;
      }
    }
    while (had_divide);

    return local_extremes;
//    Dout(dc::notice, "This function has " << local_extremes << " local extremes.");

//    Dout(dc::notice, "This function has " << minima_.size() << " minima and " << maxima_.size() << " maxima.");
  }

  void advanced_normalize(double min_percentage, double max_percentage, Range x_range)
  {
    // First, apply regular normalization to [-1, 1].
    normalize_amplitude({-1, 1}, x_range, 2000);

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
            left_y = operator()(left_x);
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
            right_y = operator()(right_x);
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
            left_y = operator()(left_x);
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
            right_y = operator()(right_x);
            right_target = right_y;
          }
        }
      }

      // Linear interpolation between surrounding extrema.
      double t = (x - left_x) / (right_x - left_x);
      double y = operator()(x);
      double y_normalized = left_target + t * (right_target - left_target);

      return y_normalized;
    };

    // Apply transformation.
    //auto original_function = [this](double x) { return operator()(x); };
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

  constexpr double w0 = 0.0;
  constexpr double learning_rate = 0.01;
  constexpr double L_max = 100;

  Algorithm gda(learning_rate, L_max);
  double w = w0;

  unsigned int seed = 0;
  if (argc > 1) {
    seed = std::stoul(argv[1]);
  }

  int extremes = 0;
  for (int i = 0; i < 100; ++i)
  {
    Range frequency_range{1.0, 20.0};
    TestFunctionGenerator L(20, frequency_range, {50.0, 100.0}, seed);

    //L.advanced_normalize(0.1, 0.1, {-1.0, 1.0});
    Range x_range{-1.0, 1.0};

//    L.normalize_amplitude({-1, 1}, x_range, 2000);
    extremes += L.find_extrema(x_range);

    if (i == 99)
    {
      Dout(dc::warning, "Prediction: " << 10);
      Dout(dc::warning, "average number of extremes: " << (0.01 * extremes));

#if 0
      // Draw vertical lines.
      gda.enable_drawing(L, -1.0, 1.0);
      std::vector<cairowindow::plot::Line> lines;
      {
        using namespace cairowindow;
        draw::LineStyle line_style({.line_width = 1.0});
        auto& enable_drawing = gda.enable_drawing();
        bool first = true;
        for (Interval* interval = L.intervals().front(); !L.intervals().is_root(interval); interval = interval->next())
        {
          for (int lr = (first ? 0 : 1); lr < 2; ++lr)
          {
            double x = lr == 0 ? interval->x_range_begin() : interval->x_range_end();
            int sign = lr == 0 ? interval->begin_sign() : interval->end_sign();
            Color color = sign == -1 ? color::red : color::blue;

            lines.push_back(enable_drawing.plot().create_line(enable_drawing.second_layer(),
                line_style({.line_color = color}), Point{x, 0.0}, Direction::up));
          }
          first = false;
        }
      }
#endif

#if 0
      while (gda(w, L(w), L.derivative(w)))
      {
        Dout(dc::notice, "-------------------------------------------");
      }
#else
      Debug(dc::notice.off());
//        gda.enable_drawing().wait();
      Debug(dc::notice.on());
#endif
    }
  }

#if 0
  ASSERT(gda.success());
  gradient_descent::Sample const& result = gda.minimum();
  Dout(dc::notice, "Global minimum: " << result);
#endif
}
