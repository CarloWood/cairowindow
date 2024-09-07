class TestFunctionGenerator : public enable_drawing::Function
{
 private:
  std::vector<double> amplitudes_;
  std::vector<double> frequencies_;
  std::vector<double> phases_;
  double vertical_shift_;
  double amplitude_scale_;
  unsigned int seed_;
  std::mt19937 generator_;

 public:
  TestFunctionGenerator(int num_components, Range frequency_range, Range amplitude_range, unsigned int seed = 0) :
    vertical_shift_(0.0), amplitude_scale_(1.0), seed_(seed)
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
    vertical_shift_ = desired_amplitude_range.center() - (current.center() * amplitude_scale_);

    Dout(dc::notice, "Normalization: current range = " << current <<
         ", amplitude_scale = " << amplitude_scale_ <<
         ", vertical_shift = " << vertical_shift_);

    // Verify normalization
    Range normalized{std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest()};
    for (int i = 0; i < num_samples; ++i)
    {
      double x = x_range.min() + i * x_range.size() / (num_samples - 1);
      double y = operator()(x);
      normalized = Range{std::min(normalized.min(), y), std::max(normalized.max(), y)};
    }
    Dout(dc::notice, "After normalization: actual range = " << normalized);
  }
};
