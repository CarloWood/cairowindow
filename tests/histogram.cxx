#include "sys.h"
#include "gradient_descent/HorizontalDirection.h"
#include "gradient_descent/ExtremeType.h"
#include <iostream>
#include <random>
#include <array>
#include <list>
#include <limits>
#include <iomanip>
#include <cassert>
#include "debug.h"

constexpr int number_of_test_runs = 1000000;
constexpr int min_number_of_extremes = 1;
constexpr int max_number_of_extremes = 11;
constexpr int max_step_size = 100;
constexpr int lowest_possible_height = ((max_number_of_extremes + 1) / 2) * -max_step_size;
constexpr int highest_possible_height = (2 + max_number_of_extremes / 2) * max_step_size;

std::uniform_int_distribution<int> extremes_distribution((min_number_of_extremes - 1) / 2, (max_number_of_extremes - 1) / 2);
std::uniform_int_distribution<int> height_step_distribution(1, max_step_size);

unsigned long used_height_marker = 0;
std::array<unsigned long, highest_possible_height - lowest_possible_height + 1> used_heights;

using namespace gradient_descent;

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
    os << std::setw(4) << std::setfill(' ') << std::right << extremes_[e] << " ~";
    // Only print one '~' per 10 height diff (so this output fits in discord).
    for (int h = 0; h < (extremes_[e] - min_height_) / 10; ++h)
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

class Algorithm
{
 private:
  Histogram const& histogram_;
  HorizontalDirection hdirection_;

  using extremes_type = std::array<LocalExtreme, max_number_of_extremes>;
  extremes_type extremes_;
  extremes_type::iterator best_minimum_;

  int last_step_ = 0;
  int last_w_;
  bool saw_minimum_ = false;

 public:
  Algorithm(Histogram const& histogram) :
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

  bool is_visited(int w) const
  {
    return extremes_[w].w() != -1;
  }

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

  void reset()
  {
    hdirection_ = HorizontalDirection::undecided;
    last_step_ = 0;
    saw_minimum_ = false;
  }

  void mark_explored(int w, HorizontalDirection hdirection);

  bool is_sane() const;
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
HorizontalDirection Algorithm::lowest_side(Weight const& w)
{
  HorizontalDirection result = HorizontalDirection::undecided;
  if (histogram_.number_of_extremes() > 1)
  {
    if (w == histogram_.number_of_extremes() - 1 || (0 < w && histogram_[w - 1] < histogram_[w + 1]))
      result = HorizontalDirection::left;
    else if (w == 0 || (w < histogram_.number_of_extremes() - 1 && histogram_[w - 1] > histogram_[w + 1]))
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
// In the real code a polynomial is fitted and extrapolated to predict the next extreme
// to jump to. Since this can completely miss the target, those jumps are generated as
// random values in this test. However, we should never jump past an already sampled point.
// To mimic this here we never jump past an already visited "extreme" (histogram index).
// But, in order not to lose the randomness, the returned value isn't clamped to the nearest
// visited extreme, instead we half the jumping distance until it falls within the allowed
// range.
//
int Algorithm::jump(Weight const& w)
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
    if (is_visited(nw))
    {
      visited = nw;
      break;
    }
  if (visited != -1)
  {
    if (visited == neighbor)
    {
      // The immediate neighbor in that direction was already visited.
      // If the neighbor is a maximum that wasn't yet explored yet in this direction, then jump to it!
      if (Histogram::is_minimum(w) && !extremes_[neighbor].is_explored(hdirection_))
        return direction;

      return 0; // Abort this direction.
    }
    while ((w + step - visited) * direction >= 0)
      step /= 2;
  }
  return step;
}

void Algorithm::mark_explored(int w, HorizontalDirection hdirection)
{
  ASSERT(hdirection != HorizontalDirection::undecided);
  std::cout << ":E(" << ((hdirection == HorizontalDirection::left) ? '-' : '+') << ')' << std::flush;
  // Mark that the extreme at w is being explored in hdirection.
  extremes_[w].explored(hdirection);

  // Make sure we don't mark extremes that were not visited yet!
  ASSERT(is_visited(w));
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
// The reasoning here is as follows: if jump returns -1 or 1 we jump
// to the immediate neighbor, therefore we can never go in that
// direction again. And if it returns 0 we can already never go into
// that direction, so don't try again in the future.
//
bool Algorithm::do_step(Weight& w)
{
  int step;
  do
  {
    LocalExtreme& current = extremes_[w];

    // Was this minimum already explored in both directions?
    if (current.done())
      return false;

    if (hdirection_ == HorizontalDirection::undecided)
    {
      if (current.is_explored(HorizontalDirection::left))
        hdirection_ = HorizontalDirection::right;
      else if (current.is_explored(HorizontalDirection::right))
        hdirection_ = HorizontalDirection::left;
      else
        hdirection_ = lowest_side(w);
    }
    else if (current.is_explored(hdirection_))
      return false;

    // Happens if there is only a single extreme.
    if (hdirection_ == HorizontalDirection::undecided)
      return false;

    // Return the step that should be made for a jump in hdirection_.
    step = jump(w);

    if (step == 0)
    {
      // It wasn't possible to go in that direction! This can happen when we reach the edge for example.
      mark_explored(w, hdirection_);
    }
  }
  while (step == 0);

  last_w_ = w;
  last_step_ = step;

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
bool Algorithm::handle_local_extreme(Weight& w)
{
  // Store it as an extreme.
  extremes_[w] = LocalExtreme{w};
  extremes_type::iterator new_extreme = &extremes_[w];

  // If this is an extreme that we found by exploring hdirection_ from a previous exterme,
  // then mark that last extreme as being explored in the hdirection_.
  if (std::abs(last_step_) == 1)
  {
    mark_explored(last_w_, hdirection_);

    // If the current extreme is a maximum, then it can be marked as having been explored to the opposite direction (in that case the previous extreme is a minimum).
    // Or if the current extreme is a minimum and we came from another minimum (with steps of size 1), then mark it as having been explored into the opposite direction.
    if (saw_minimum_)
    {
      // Mark that the current extreme is being explored in the opposite hdirection_.
      mark_explored(w, opposite(hdirection_));
    }
  }

  // Keep track of the best minimum so far; or abort if this minimum isn't better then one found before.
  if (Histogram::is_minimum(w))
  {
    saw_minimum_ = true;
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
bool Algorithm::handle_abort_hdirection(Weight& w)
{
  std::cout << "(a)";
  reset();

  // Has a minimum been found at all?
  if (best_minimum_ == extremes_.end())
    return false;

  // Restore the current sample and scale to the values belonging to this minimum.
  w = best_minimum_->w();
  std::cout << " ==> " << w << std::flush;

  // Go to the next extreme.
  return do_step(w);
}

bool Algorithm::operator()(Weight& w, int height)
{
  // If the step size is larger than one, than we should be able to detect
  // that we skipped an extreme (with a third degree polynomial fit).
  // In that case we want to half the step size *unless* the returned
  // value is lower than the current best minimum.
  if (std::abs(last_step_) > 1)
  {
    // However, if the height here is less than that of the current best minimum then
    // we accept the jump regardless.
    if (best_minimum_ != extremes_.end() && height < histogram_[best_minimum_->w()])
    {
      reset();
    }
    else
    {
      // Reject and half the jump step.
      std::cout << "(r)";
      last_step_ /= 2;
      w = last_w_ + last_step_;
      return true;
    }
  }

  // Every sample is a local extreme here.
  if (!handle_local_extreme(w))
    return handle_abort_hdirection(w);

  return true;        // w was successfully updated by handle_local_extreme.
}

bool Algorithm::is_sane() const
{
  // Get the index of the minimum that was found.
  int result = minimum();
  Dout(dc::finish, result);

  // The result must be a minimum.
  if (!Histogram::is_minimum(result))
  {
    Dout(dc::warning, "The returned result is not a minimum!");
    return false;
  }

  // If there exists a minimum on the left...
  if (result - 2 >= 0)
  {
    // It must be visited.
    if (!is_visited(result - 2))
    {
      Dout(dc::warning, "The minimum on the left wasn't visited!");
      return false;
    }
    // And it must be higher than the found minimum.
    if (!(histogram_[result - 2] > histogram_[result]))
    {
      Dout(dc::warning, "The minimum on the left is lower!");
      return false;
    }
  }
  // If there exists a minimum on the right...
  if (result + 2 < histogram_.number_of_extremes())
  {
    // It must be visited.
    if (!is_visited(result + 2))
    {
      Dout(dc::warning, "The minimum on the right wasn't visited!");
      return false;
    }
    // And it must be higher than the found minimum.
    if (!(histogram_[result + 2] > histogram_[result]))
    {
      Dout(dc::warning, "The minimum on the right is lower!");
      return false;
    }
  }

  return true;
}

int minimum(Histogram const& histogram)
{
  // Create algorithm object.
  Algorithm gda(histogram);

  // * start at starting_position_.
  Weight w{histogram.starting_position()};

  // Loop over iterations of w.
  std::cout << w << std::flush;
  for (;;)
  {
    // Replace w with a new point until the global minimum has been reached.
    if (!gda(w, histogram[w]))
      break;
    std::cout << " --> " << w << std::flush;
  }
  std::cout << std::endl;

  ASSERT(gda.success());
  // Check that the result is sane.
  ASSERT(gda.is_sane());

  return gda.minimum();
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
    Dout(dc::notice|continued_cf, "Result of i = " << i << ": ");       // Finished in is_sane() called from minimum().
    int min = minimum(histogram);
  }
}
