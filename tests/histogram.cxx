#include "sys.h"
#include "HorizontalDirection.h"
#include "VerticalDirection.h"
#include <iostream>
#include <random>
#include <array>
#include <list>
#include <limits>
#include <iomanip>
#include <cassert>
#include "debug.h"

using HorizontalDirection = gradient_descent::HorizontalDirection;
using VerticalDirection = gradient_descent::VerticalDirection;

constexpr int number_of_test_runs = 2; //000;
constexpr int min_number_of_extremes = 1;
constexpr int max_number_of_extremes = 11;
constexpr int max_step_size = 100;
constexpr int lowest_possible_height = ((max_number_of_extremes + 1) / 2) * -max_step_size;
constexpr int highest_possible_height = (2 + max_number_of_extremes / 2) * max_step_size;

std::uniform_int_distribution<int> extremes_distribution((min_number_of_extremes - 1) / 2, (max_number_of_extremes - 1) / 2);
std::uniform_int_distribution<int> height_step_distribution(1, max_step_size);

unsigned long used_height_marker = 0;
std::array<unsigned long, highest_possible_height - lowest_possible_height + 1> used_heights;

class Histogram
{
 private:
  int const number_of_extremes_;                        // The number of valid entries in extremes_.
  std::array<int, max_number_of_extremes> extremes_;    // The height of each extreme (index in [0..number_of_extremes_>).
  int starting_position_;                               // The index of the initial extreme to start with.
  std::array<int, max_number_of_extremes> left_jump_;   // The next index to jump to if going to the left.
  std::array<int, max_number_of_extremes> right_jump_;  // The next index to jump to if going to the right.
  int min_height_;                                      // The lowest value in extremes_.
  int max_height_;                                      // The highest value in extremes_.

 public:
  Histogram(std::mt19937& generator);

  static bool is_minimum(int w);

  // Return the index of the actual global minimum.
  int real_minimum() const;

  // Accessors.
  int number_of_extremes() const { return number_of_extremes_; }
  int operator[](int e) const { return extremes_[e]; }
  int starting_position() const { return starting_position_; }
  int left_jump(int e) const { return left_jump_[e]; }
  int right_jump(int e) const { return right_jump_[e]; }
  int min_height() const { return min_height_; }
  int max_height() const { return max_height_; }

  void print_on(std::ostream& os) const;
};

std::ostream& operator<<(std::ostream& os, Histogram const& histogram)
{
  histogram.print_on(os);
  return os;
}

Histogram::Histogram(std::mt19937& generator) :
  number_of_extremes_(extremes_distribution(generator) * 2 + 1)
{
  // Increment the marker (this is assumed to be single threaded).
  ++used_height_marker;
  // Initial position.
  std::uniform_int_distribution<int> initial_position_distribution(0, number_of_extremes_ - 1);
  starting_position_ = initial_position_distribution(generator);
  // The first extreme is always a minimum.
  int s = 1;
  extremes_[0] = height_step_distribution(generator);                           // [0, 100]
  used_heights[extremes_[0] - lowest_possible_height] = used_height_marker;
  min_height_ = extremes_[0];
  max_height_ = extremes_[0];
  // Fill all other extremes, if any.
  for (int e = 1; e < number_of_extremes_; ++e, s = -s)
  {
    int height;
    do
    {
      height = extremes_[e - 1] + s * height_step_distribution(generator);    // [0, 200], [-100, 200], [-100, 300], [-200, 300], ...
    }
    while (used_heights[height - lowest_possible_height] == used_height_marker);
    extremes_[e] = height;
    used_heights[height - lowest_possible_height] = used_height_marker;
    if (extremes_[e] < min_height_)
      min_height_ = extremes_[e];
    if (extremes_[e] > max_height_)
      max_height_ = extremes_[e];
  }
  // Fill jump tables.
  for (int e = 0; e < number_of_extremes_; ++e)
  {
    if (e == 0)
      left_jump_[0] = -1;
    else
    {
      std::uniform_int_distribution<int> left_jump_distribution(0, e - 1);
      left_jump_[e] = left_jump_distribution(generator);
    }
    if (e == number_of_extremes_ - 1)
      right_jump_[e] = -1;
    else
    {
      std::uniform_int_distribution<int> right_jump_distribution(e + 1, number_of_extremes_ - 1);
      right_jump_[e] = right_jump_distribution(generator);
    }
  }
}

//static
bool Histogram::is_minimum(int w)
{
  // The first extreme is a minimum; then they alternate.
  return (w & 1) == 0;
}

int Histogram::real_minimum() const
{
  int min = extremes_[0];
  int gm = 0;
  for (int e = 0; e < number_of_extremes_; ++e)
  {
    if (extremes_[e] < min)
    {
      min = extremes_[e];
      gm = e;
    }
  }
  return gm;
}

void Histogram::print_on(std::ostream& os) const
{
  os << "Number of extremes: " << number_of_extremes_ << '\n';
  os << "Starting position: " << starting_position_ << '\n';
  // Find lowest and highest values.
  os << "Max height: " << max_height_ << "; min height: " << min_height_ << '\n';
  int total_height = max_height_ - min_height_;
  for (int e = 0; e < number_of_extremes_; ++e)
  {
    os << std::setw(2) << std::setfill(' ') << left_jump_[e] << " <-- ";
    os << std::setw(2) << std::setfill(' ') << e << " --> ";
    os << std::setw(2) << std::setfill(' ') << right_jump_[e] << " : ";
    for (int h = min_height_ - 1; h < extremes_[e]; ++h)
      os << '~';
    if ((e & 1) == 1)
      os << "+\n";
    else
      os << "-\n";
  }
}

constexpr int left = -1;
constexpr int right = 1;

class Weight
{
 protected:
  int w_;

 public:
  Weight(int w) : w_(w) { }

  Weight& operator=(int new_w)
  {
    w_ = new_w;
    return *this;
  }

  operator int() const
  {
    return w_;
  }
};

class LocalExtreme
{
 protected:
  int w_;                               // Extreme index, or -1 if not sampled yet.
  int explored_{0};                     // Bit 1: exploration to the left of this extreme has started.
                                        // Bit 2: same, on the right.

 public:
  LocalExtreme() : w_(-1) {}
  LocalExtreme(int w) : w_(w) { }

  int w() const { return w_; }

  void explored(HorizontalDirection hdirection)
  {
    int explore_flag = hdirection == HorizontalDirection::left ? 1 : 2;
    explored_ |= explore_flag;
  }

  bool done() const
  {
    return explored_ == 3;
  }

  bool is_explored(HorizontalDirection hdirection) const
  {
    int explore_flag = hdirection == HorizontalDirection::left ? 1 : 2;
    return (explored_ & explore_flag) != 0;
  }
};

class AcceleratedGradientDescent
{
 private:
  Histogram const& histogram_;
  HorizontalDirection hdirection_;

  using extremes_type = std::array<LocalExtreme, max_number_of_extremes>;
  extremes_type extremes_;
  extremes_type::iterator best_minimum_;

 public:
  AcceleratedGradientDescent(Histogram const& histogram) :
    histogram_(histogram),
    hdirection_(HorizontalDirection::undecided),
    extremes_{{}},
    best_minimum_(extremes_.end()) { }

  bool operator()(Weight& w, int height);
  bool handle_local_extreme(Weight& w);
  bool handle_abort_hdirection(Weight& w);
  HorizontalDirection lowest_side(Weight const& w);
  int jump(Weight const& w);
  bool do_step(Weight& w);

  bool success() const
  {
    return best_minimum_ != extremes_.end();
  }

  int minimum() const
  {
    // Only call this function when success() returns true.
    ASSERT(success());
    return best_minimum_->w();
  }
};

// lowest_side
//
// Returns HorizontalDirection::undecided, HorizontalDirection::left or HorizontalDirection::right.
//
// In the real code, a fourth degree polynomial is fitted whenever an extreme was found.
// This fourth degree polynomial has three extremes; the one in the middle matches the extreme
// that was just found and the two on the left and right approximate the height and position
// of the extremes of opposite type (minimum <--> maximum). The direction chosen is then the
// direction in which that extreme has the lowest value (for both maximum and minimum).
//
// In this test code we simply use the actual "extremes" left and right of the current one:
// the values of the histogram left and right of the current index w.
//
HorizontalDirection AcceleratedGradientDescent::lowest_side(Weight const& w)
{
  HorizontalDirection result = HorizontalDirection::undecided;
  if (histogram_.number_of_extremes() > 1)
  {
    if (w == histogram_.number_of_extremes() - 1 || (0 < w && histogram_[extremes_[w - 1].w()] < histogram_[extremes_[w + 1].w()]))
      result = HorizontalDirection::left;
    else if (w == 0 || (w < histogram_.number_of_extremes() - 1 && histogram_[extremes_[w - 1].w()] > histogram_[extremes_[w + 1].w()]))
      result = HorizontalDirection::right;
  }
  return result;
}

// jump
//
// Returns the distance that we would jump if going direction hdirection_.
// If hdirection_ is left then a negative value is returned, if hdirection_ is right a positive value,
// so that simply adding the returned value to w gives the new w.
//
// Returns 0 if there is no jump possible in that direction (because the immediate neighbor was already visited).
//
// In the real code a polynomial is fitted and extra polated to predict the next extreme
// to jump to. Since this can completely miss the target, those jumps are generated as
// random values in this test. However, we should never jump past an already sampled point.
// To mimic this here we never jump past an already visited "extreme" (histogram index).
// But, in order not to lose the randomness, the returned value isn't clamped to the nearest
// visited extreme, instead we half the jumping distance until it falls within the allowed
// range.
//
int AcceleratedGradientDescent::jump(Weight const& w)
{
  // Find nearest neighbor that was already visited.
  int const direction = static_cast<int>(hdirection_);
  int const neighbor = w + direction;
  int const target_w = (hdirection_ == HorizontalDirection::left) ? histogram_.left_jump(w) : histogram_.right_jump(w);
  if (target_w == -1)
    return 0;   // Can't go in hdirection_.
  int step = target_w - w;
  int visited = -1;
  for (int nw = neighbor; nw != target_w + direction; nw += direction)
    if (extremes_[nw].w() != -1)
    {
      visited = nw;
      break;
    }
  if (visited != -1)
  {
    if (visited == neighbor)
      return 0; // Abort this direction.
    while ((w + step - visited) * direction >= 0)
      step /= 2;
  }
  return step;
}

// do_step
//
// Returns true upon success, false otherwise.
//
// Execute a 'step': go from the current extreme to the next.
// The algorithm here is as follows:
// - pick a direction: if the current extreme is already explored
//   left and right, then return false. Otherwise if the current
//   extreme is already explored to the left, go right. Otherwise,
//   if the current extreme is already explored to the right, go left.
//   Otherwise, go to the "lowest side".
//
// - call jump. If that returns -1, 0 or 1, mark the extreme as
//   having been explored in this direction and try the other
//   side (if that wasn't already explored) or abort (return false).
//
// The reasoning here is as follows: if jump return -1 or 1 we jump
// to the immediate neighbor, therefore we can never go in that
// direction again. And if it return 0 we can already never go into
// that direction, so don't try again in the future.
//
bool AcceleratedGradientDescent::do_step(Weight& w)
{
  int step;
  do
  {
    LocalExtreme& current = extremes_[w];

    // Was this minimum already explored in both directions?
    if (current.done())
      return false;

    if (current.is_explored(HorizontalDirection::left))
      hdirection_ = HorizontalDirection::right;
    else if (current.is_explored(HorizontalDirection::right))
      hdirection_ = HorizontalDirection::left;
    else
      hdirection_ = lowest_side(w);

    // Happens if there is only a single extreme.
    if (hdirection_ == HorizontalDirection::undecided)
      return false;

    // Return the step that should be made for a jump in hdirection_.
    step = jump(w);

    if (std::abs(step) <= 1)
    {
      // Mark that the current extreme is being explored in the hdirection_.
      current.explored(hdirection_);
    }
  }
  while (step == 0);

  w = w + step;
  return true;
}

// handle_local_extreme
//
// Returns true upon success and false upon failure.
//
// In this test, every sample is an extreme (each column of the
// histogram represents an extreme). Nevertheless we store them
// in the order that we're visiting them.
//
// This function is broken. The idea here we already are exploring
// a certain horizontal direction relative to the best minimum that
// we found so far, and keep going in that direction until we find
// a minimum that is NOT better.
//
// If you assume that we find all extremes exactly in order,
// then this algorithm works. But it is possible that the best
// minimum is jumped over (or is it?).
bool AcceleratedGradientDescent::handle_local_extreme(Weight& w)
{
  // Store it as an extreme.
  extremes_[w] = LocalExtreme{w};
  extremes_type::iterator new_extreme = &extremes_[w];

  // Keep track of the best minimum so far; or abort if this minimum isn't better then one found before.
  if (Histogram::is_minimum(w))
  {
    if (best_minimum_ == extremes_.end() || histogram_[best_minimum_->w()] > histogram_[new_extreme->w()])
      best_minimum_ = new_extreme;
    if (new_extreme != best_minimum_)
      return false;
  }

  // Go to the next extreme.
  return do_step(w);
}

// handle_abort_hdirection
//
// Called when we can't go further in the current hdirection.
// In this case we want to jump to the best minimum so far
// and check that we explored both sides of it already.
bool AcceleratedGradientDescent::handle_abort_hdirection(Weight& w)
{
  // Has a minimum been found at all?
  if (best_minimum_ == extremes_.end())
    return false;

  // Restore the current sample and scale to the values belonging to this minimum.
  w = best_minimum_->w();

  // Go to the next extreme.
  return do_step(w);
}

bool AcceleratedGradientDescent::operator()(Weight& w, int height)
{
  // Every sample is a local extreme here.
  if (!handle_local_extreme(w))
    return handle_abort_hdirection(w);
  return true;        // w was successfully updated by handle_local_extreme.
}

int minimum(Histogram const& histogram)
{
  // Create algorithm object.
  AcceleratedGradientDescent agd(histogram);

  // * start at starting_position_.
  Weight w{histogram.starting_position()};

  // Loop over iterations of w.
  std::cout << w << std::flush;
  for (;;)
  {
    // Replace w with a new point until the global minimum has been reached.
    if (!agd(w, histogram[w]))
      break;
    std::cout << " --> " << w << std::flush;
  }
  std::cout << std::endl;

  ASSERT(agd.success());

  return agd.minimum();
}

int main()
{
  Debug(NAMESPACE_DEBUG::init());

  std::seed_seq seed{1, 2, 3, 4, 10};
  std::mt19937 generator(seed);

  // Generate 1000 tests.
  for (int i = 0; i < number_of_test_runs; ++i)
  {
    Histogram histogram(generator);
    Dout(dc::notice, histogram);
    Dout(dc::notice, "Result: " << minimum(histogram));
  }
}
