#include "sys.h"
#include "Function.h"
#include "Algorithm.h"
#include "utils/AIRefCount.h"
#include "utils/UniqueID.h"
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

using Range = cairowindow::Range;
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

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
  std::vector<std::pair<double, double>> minima_;
  std::vector<std::pair<double, double>> maxima_;
  utils::UniqueIDContext<int> id_context_;

  int sign_of_derivative(double x) const
  {
    return derivative(x) < 0.0 ? -1 : 1;
  }

  utils::UniqueIDContext<int>& id_context()
  {
    return id_context_;
  }

  // <--------------------------------------------------------------->
  // <--->                                                                    max. stop interval size (1/16 of total interval size)
  // [                                                               ]
  // +                                                               +
  //
  // <--------------------------------------------------------------->
  // <--->                                                                    max. stop interval size
  // [                               |                               ]
  // +                               +                               +
  //
  // .               .               <------------------------------->
  //                                         1               2                max. stop interval sizes
  // [               |               |               |               ]
  // +               -               +               +               +
  //
  // <------->       <------->       <--------------->       .       .
  //     1               1               1       2                            max. stop interval sizes
  // [       |       |       |       |       |       |       |       ]
  // +       +       -       -       +       +       +       -       +
  //
  // <------->   <----------->   <----------->   .   .   <--->   <--->
  //   1   2       1   2   3       1   2   3               1       1          max. stop interval sizes
  // [   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   ]
  // +   +   +   -   -   -   -   +   +   +   +   -   +   -   -   +   +
  //
  // <-------> <-------------> <-------------> <-> <-> <-----> <----->
  //  1 2 3 4   1 2 3 4 5 6 7   1 2 3 4 5 6 7   1   1   1 2 3   1 2 3
  // [ | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | ]
  // + + + + + - - - - - - - - + + + + + + + + - - + + - - - - + + + +
  //
  //                                          ||||||||||||||||||||||||
  // +++++++++----------------++++++++++++++++----+++++--------+++++++

  class IntervalSpan : public AIRefCount
  {
   private:
    Range range_;
    int number_of_intervals_;
    utils::UniqueID<int> id_;

   public:
    IntervalSpan(Range x_range, utils::UniqueIDContext<int>& id_context) :
      range_(x_range), number_of_intervals_(0), id_(id_context.get_id()) { }

    void increment() { ++number_of_intervals_; }
    void decrement() { --number_of_intervals_; }

    void change_begin(double range_begin, int number_of_intervals)
    {
      range_ = Range{range_begin, range_.max()};
      number_of_intervals_ = number_of_intervals;
    }

    void change_begin(double range_begin)
    {
      range_ = Range{range_begin, range_.max()};
    }

    void change_end(double range_end, int number_of_intervals)
    {
      range_ = Range{range_.min(), range_end};
      number_of_intervals_ = number_of_intervals;
    }

    void change_end(double range_end)
    {
      range_ = Range{range_.min(), range_end};
    }

    Range const& range() const { return range_; }
    int number_of_intervals() const { return number_of_intervals_; }
    char id() const { return 'A' + id_; }

    double interval_size() const
    {
      return range_.size() / number_of_intervals_;
    }

#ifdef CWDEBUG
    void print_on(std::ostream& os) const
    {
      os << "{range_:" << range_ << ", number_of_intervals_:" << number_of_intervals_ << "}";
    }
#endif
  };

  class Interval
  {
   public:
    using span_type = std::pair<boost::intrusive_ptr<IntervalSpan>, boost::intrusive_ptr<IntervalSpan>>;

   private:
    Range x_range_;
    int begin_sign_;
    int end_sign_;
    span_type span_;
    Interval* prev_{nullptr};
    Interval* next_{nullptr};

   public:
    bool contains_sign_change() const { return begin_sign_ != end_sign_; }

    // Create a single linked list.
    Interval(Range x_range, TestFunctionGenerator* function_generator) : x_range_(x_range),
      begin_sign_(function_generator->sign_of_derivative(x_range.min())), end_sign_(function_generator->sign_of_derivative(x_range.max()))
    {
      if (!contains_sign_change())
      {
        span_.first = new IntervalSpan(x_range, function_generator->id_context());
        span_.first->increment();
      }
    }

    Interval(IntervalSpan* interval_span, Range x_range, int begin_sign, int end_sign) : x_range_(x_range),
      begin_sign_(begin_sign), end_sign_(end_sign), span_(nullptr, nullptr)
    {
      if (begin_sign == end_sign)
      {
        span_.first = interval_span;
        if (interval_span)
          interval_span->increment();
      }
      else
        span_.second = interval_span;
    }

    double x_range_begin() const { return x_range_.min(); }
    double x_range_end() const { return x_range_.max(); }

    void insert_before(IntervalSpan* interval_span, Range x_range, int begin_sign, int end_sign)
    {
      Interval* new_interval = new Interval(interval_span, x_range, begin_sign, end_sign);
      //                     this
      //                     v
      // .----.              .----.
      // |    |<-------------+prev|
      // |next+------------->|    |
      // `----'              `----'
      //                     |
      //                    -
      //                  -
      //                - ... this is the same border: same x-coordinate and derivative (sign).
      //              -
      //            -
      //           |         this
      //           v         v
      // .----.    .----.    .----.
      // |    |<---+prev|<---+prev|
      // |next+--->|next+--->|    |
      // `----'    `----'    `----'
      //           ^
      //       new_interval

      // This must be the same border.
      ASSERT(begin_sign_ == new_interval->begin_sign_);
      ASSERT(x_range_.min() == new_interval->x_range_.min());

      if (prev_)
        prev_->next_ = new_interval;

      new_interval->next_ = this;
      new_interval->prev_ = prev_;

      prev_ = new_interval;

      x_range_ = Range{new_interval->x_range_.max(), x_range_.max()};
      begin_sign_ = new_interval->end_sign_;
    }

    void insert_after(IntervalSpan* interval_span, Range x_range, int begin_sign, int end_sign)
    {
      Interval* new_interval = new Interval(interval_span, x_range, begin_sign, end_sign);
      // this
      // v
      // .----.              .----.
      // |    |<-------------+prev|
      // |next+------------->|    |
      // `----'              `----'
      //      |
      //       -
      //         -
      //           - ... this is the same border: same x-coordinate and derivative (sign).
      //             -
      //               -
      // this           |
      // v              v
      // .----.    .----.    .----.
      // |    |<---+prev|<---+prev|
      // |next+--->|next+--->|    |
      // `----'    `----'    `----'
      //           ^
      //       new_interval

      // This must be the same border.
      ASSERT(end_sign_ == new_interval->end_sign_);
      ASSERT(x_range_.max() == new_interval->x_range_.max());

      if (next_)
        next_->prev_ = new_interval;

      new_interval->prev_ = this;
      new_interval->next_ = next_;

      next_ = new_interval;

      x_range_ = Range{x_range_.min(), new_interval->x_range_.min()};
      end_sign_ = new_interval->begin_sign_;
    }

    bool must_be_divided() const
    {
      if (!contains_sign_change())
        return span_.first->number_of_intervals() < 4;

      if (!span_.first || !span_.second || span_.first->number_of_intervals() < 4 || span_.second->number_of_intervals() < 2)
        return true;

      return x_range_.size() > std::min(span_.first->interval_size(), span_.second->interval_size());
    }

    bool is_first_interval() const
    {
      // Only call this for intervals that do not change sign.
      ASSERT(!contains_sign_change());
      return span_.first->range().min() == x_range_.min();
    }

    bool is_last_interval() const
    {
      // Only call this for intervals that do not change sign.
      ASSERT(!contains_sign_change());
      return span_.first->range().max() == x_range_.max();
    }

    void split(TestFunctionGenerator* function_generator)
    {
      double mid = x_range_.center();
      int mid_sign = function_generator->sign_of_derivative(mid);

      if (!contains_sign_change())
      {
        IntervalSpan* const P = span_.first.get();
        if (is_first_interval())
        {
          // Case 1 or 3.
          if (is_last_interval())
          {
            // Case 1.
            //                 <-P------------->
            //                         1
            //                 |      P        |
            //                 +               +
            //
            if (mid_sign == begin_sign_)
            {
              //  A) finding a +:
              //                 <-P------------->
              //                     1       2
              //                 |  P    |  P    |
              //                 +       +       +
              insert_after(P, {mid, x_range_.max()}, mid_sign, end_sign_);
            }
            else
            {
              //  B) finding a -:
              //
              //                 |  NN   |  NN   |
              //                 +       -       +
              if (prev_)
                prev_->span_.second.reset();
              if (next_)
                next_->span_.first.reset();
              insert_after(nullptr, {mid, x_range_.max()}, mid_sign, end_sign_);
              // This was the last interval using this IntervalSpan.
              ASSERT(span_.first->unique().is_true());
              span_.first.reset();
            }
          }
          else
          {
            // Case 3.
            //
            //                 <-P-----------------------------?
            //                         1               2
            //                 |      P        |      P        |
            //                 +               +               +
            //
            if (mid_sign == begin_sign_)
            {
              //  A) finding a +:
              //                 <-P-----------------------------?
              //                     1       2           3
              //                 |  P    |  P    |      P        |
              //                 +       +       +               +
              insert_after(P, {mid, x_range_.max()}, mid_sign, end_sign_);
            }
            else
            {
              //  B) finding a -:
              //                                 <-P-------------?
              //                                         1
              //                 |  NN   |  NP   |      P        |
              //                 +       -       +               +
              P->change_begin(x_range_.max(), P->number_of_intervals() - 1);
              insert_after(P, {mid, x_range_.max()}, mid_sign, end_sign_);
              span_.first.reset();
            }
          }
        }
        else
        {
          // Case 2 or 4.
          if (is_last_interval())
          {
            // Case 2.
            //
            // ?-P----------------------------->
            //         n              n+1
            // |      P        |      P        |
            // +               +               +
            //
            if (mid_sign == begin_sign_)
            {
              //  A) finding a +:
              // ?-P----------------------------->
              //         n          n+1     n+2
              // |      P        |   P   |  P    |
              // +               +       +       +
              insert_after(P, {mid, x_range_.max()}, mid_sign, end_sign_);
            }
            else
            {
              //  B) finding a -:
              // ?-P------------->
              //         n
              // |      P        |   PN  |   NN  |
              // +               +       -       +
              P->change_end(x_range_.min(), P->number_of_intervals() - 1);
              insert_after(nullptr, {mid, x_range_.max()}, mid_sign, end_sign_);
            }
          }
          else
          {
            // Case 4.
            //
            // ?-P---------------------------------------------?
            //         n              n+1             n+2
            // |      P        |      P        |      P        |
            // +               +               +               +
            //
            if (mid_sign == begin_sign_)
            {
              //  A) finding a +:
              // ?-P---------------------------------------------?
              //         n          n+1     n+2         n+3
              // |      P        |  P    |  P    |      P        |
              // +               +       +       +               +
              insert_after(P, {mid, x_range_.max()}, mid_sign, end_sign_);
            }
            else
            {
              //  B) finding a -:
              // ?-P------------->               <-Q-------------?
              //         n                               1
              // |      P        |  PN   |  NQ   |      Q!       |  ...Q!
              // +               +       -       +               +

              // Count the number of intervals to the left of 'this'.
              int left_intervals = 1;   // The first interval.
              for (Interval* left = prev_; left && !left->is_first_interval(); left = left->prev_)
                ++left_intervals;

              // Create new IntervalSpan Q.
              auto Q = new IntervalSpan({x_range_.max(), span_.first->range().max()}, function_generator->id_context());

              // Update IntervalSpan P.
              int original_number_of_intervals = P->number_of_intervals();
              P->change_end(x_range_.min(), left_intervals);

              // Create the new NQ interval (sign change interval).
              insert_after(Q, {mid, x_range_.max()}, mid_sign, end_sign_);

              // Update subsequent intervals to use Q instead of P.
              for (Interval* next = next_->next_; next && !next->contains_sign_change(); next = next->next_)
                next->span_.first = Q;

              // Update Q's number of intervals.
              Q->change_begin(x_range_.max(), original_number_of_intervals - left_intervals - 1);

              // Update the current interval.
              x_range_ = Range{x_range_.min(), mid};
              end_sign_ = mid_sign;
            }
          }
        }
      }
      else
      {
        // Get the span on the Left and/or Right, if any.
        IntervalSpan* const L = span_.first.get();
        IntervalSpan* const R = span_.second.get();
        if (!L)
        {
          // Case 1 or 3.
          if (!R)
          {
            // Case 1.
            //
            // |      ?N       |      NN       |      N?       |
            // -               +               -               +
            //
            if (mid_sign == begin_sign_)
            {
              //  A) finding a +:
              //                 <-Q----->
              //                     1
              // |      ?Q       |   Q   |  QN   |      N?       |
              // -               +       +       -               +
              auto Q = new IntervalSpan({x_range_.min(), mid}, function_generator->id_context());
              if (prev_)
              {
                // ?N --> ?Q
                prev_->span_.second = Q;
              }
              insert_before(Q, {x_range_.min(), mid}, begin_sign_, mid_sign);
              // NN --> QN
              span_.first = Q;
            }
            else
            {
              //  B) finding a -:
              //                         <-Q----->
              //                             1
              // |      ?N       |  NQ   |   Q   |      Q?       |
              // -               +       -       -               +
              auto Q = new IntervalSpan({mid, x_range_.max()}, function_generator->id_context());
              if (next_)
              {
                // N? --> Q?
                next_->span_.first = Q;
              }
              insert_after(Q, {mid, x_range_.max()}, mid_sign, end_sign_);
              // QN --> NQ
              span_.second = Q;
              span_.first.reset();
            }
          }
          else
          {
            // Case 3.
            //
            //                                 <-R-------------?
            //                                         1
            // |      ?N       |      NR       |      R        |
            // -               +               -               -
            //
            if (mid_sign == begin_sign_)
            {
              //  A) finding a +:
              //                 <-Q----->       <-R-------------?
              //                     1                   1
              // |      ?Q       |   Q   |  QR   |      R        |
              // -               +       +       -               -
              auto Q = new IntervalSpan({x_range_.min(), mid}, function_generator->id_context());
              if (prev_)
              {
                // ?N --> ?Q
                prev_->span_.second = Q;
              }
              insert_before(Q, {x_range_.min(), mid}, begin_sign_, mid_sign);
              // NR --> QR
              span_.first = Q;
            }
            else
            {
              //  B) finding a -:
              //                         <-R---------------------?
              //                             1           2
              // |      ?N       |   NR  |   R   |      R        |
              // -               +       -       -               -
              insert_after(R, {mid, x_range_.max()}, mid_sign, end_sign_);
              R->change_begin(mid);
            }
          }
        }
        else
        {
          // Case 2 or 4.
          if (!R)
          {
            // Case 2.
            //
            // ?-L------------->
            //         n
            // |      L        |      LN       |      N?       |
            // +               +               -               +
            //
            if (mid_sign == begin_sign_)
            {
              //  A) finding a +:
              // ?-L--------------------->
              //         n          n+1
              // |      L        |   L   |  LN   |      N?       |
              // +               +       +       -               +
              insert_before(L, {x_range_.min(), mid}, begin_sign_, mid_sign);
              L->change_end(mid);
            }
            else
            {
              //  B) finding a -:
              // ?-L------------->       <-Q----->
              //         n                   1
              // |      L        |  LQ   |   Q   |      Q?       |
              // +               +       -       -               +
              auto Q = new IntervalSpan({mid, x_range_.max()}, function_generator->id_context());
              if (next_)
              {
                // N? --> Q?
                next_->span_.first = Q;
              }
              insert_after(Q, {mid, x_range_.max()}, mid_sign, end_sign_);
              // LN --> LQ
              span_.second = Q;
            }
          }
          else
          {
            // Case 4.
            //
            // ?-L------------->               <-R-------------?
            //         n                              1
            // |      L        |      LR       |       R       |
            // +               +               -               -
            //
            if (mid_sign == begin_sign_)
            {
              //  A) finding a +:
              // ?-L--------------------->       <-R-------------?
              //         n          n+1                 1
              // |      L        |   L   |  LR   |       R       |
              // +               +       +       -               -
              insert_before(L, {x_range_.min(), mid}, begin_sign_, mid_sign);
              L->change_end(mid);
            }
            else
            {
              //  B) finding a -:
              // ?-L------------->       <-----------------------?
              //         n                   1          2
              // |      L        |  LR   |   R   |       R       |
              // +               +       -       -               -
              insert_after(R, {mid, x_range_.max()}, mid_sign, end_sign_);
              R->change_begin(mid);
            }
          }
        }
      }
    }

    Range const& x_range() const { return x_range_; }
    int begin_sign() const { return begin_sign_; }
    int end_sign() const { return end_sign_; }
    span_type const& span() const { return span_; }
    Interval* prev() { return prev_; }
    Interval const* prev() const { return prev_; }
    Interval* next() { return next_; }
    Interval const* next() const { return next_; }

    void print_range_begin_on(std::ostream& os) const
    {
      os << x_range_begin() << '(' << (begin_sign() == 1 ? '/' : '\\') << ')';
    }

    void print_range_end_on(std::ostream& os) const
    {
      os << x_range_end() << '(' << (end_sign() == 1 ? '/' : '\\') << ')';
    }

#ifdef CWDEBUG
    mutable bool sanity_check_successful_;
    void fail() const
    {
      sanity_check_successful_ = false;
    }
#endif

    void print_list(std::ostream& os) const
    {
      sanity_check_successful_ = true;
      bool first = true;
      Interval const* prev_interval = nullptr;
      bool printing_end_was_suppressed = false;
      int saw_interval_without_sign_change = 0;
      int number_of_intervals = 0;
      RESTART
      for (Interval const* interval = this; interval != nullptr; interval = interval->next())
      {
        int const interval_contains_sign_change = interval->contains_sign_change();

#ifdef CWDEBUG
        // Sanity check.
        if (interval_contains_sign_change)
        {
          if (interval->span_.first)
          {
            if (!(prev_interval && !prev_interval->contains_sign_change() && prev_interval->span_.first == interval->span_.first))
              fail();
          }
          else
          {
            if (!(!prev_interval || prev_interval->contains_sign_change()))
              fail();
          }
          Interval const* next_interval = interval->next();
          if (interval->span_.second)
          {
            if (!(next_interval && !next_interval->contains_sign_change() && next_interval->span_.first == interval->span_.second))
              fail();
          }
          else
          {
            if (!(!next_interval || next_interval->contains_sign_change()))
              fail();
          }
          if (number_of_intervals > 0)
          {
            if (!interval->span_.first)
              fail();
            else if (!(interval->span_.first->number_of_intervals() == number_of_intervals &&
                  interval->span_.first->range().max() == interval->x_range_begin()))
              fail();
          }
          number_of_intervals = 0;
        }
        else
        {
          if (++number_of_intervals == 1)
          {
            if (!(interval->span_.first->range().min() == interval->x_range_begin()))
              fail();
          }
        }
        if (prev_interval)
        {
          if (!(prev_interval->end_sign_ == interval->begin_sign_ && prev_interval->x_range_end() == interval->x_range_begin()))
            fail();
        }
#endif

        if (!interval_contains_sign_change)
          ++saw_interval_without_sign_change;
        else
          saw_interval_without_sign_change = 0;
        bool interval_is_first_interval = saw_interval_without_sign_change == 1;

        bool print_begin = first || printing_end_was_suppressed;
        if (print_begin && interval_contains_sign_change)
        {
          interval->print_range_begin_on(os);
          print_begin = false;
        }

        bool show_span = interval_contains_sign_change || interval_is_first_interval;
        if (show_span)
        {
          Interval::span_type const& span = interval->span();
          if (!interval_contains_sign_change)
            os << '{' << span.first->id() << ':' << span.first->range() << " <" << span.first->number_of_intervals() << ">}: ";
          else
          {
            os << " <-{";
            if (span.first)
              os << '"' << span.first->id() << '"';
            else
              os << '0';
            os << ", ";
            if (span.second)
              os << '"' << span.second->id() << '"';
            else
              os << '0';
            os << "}-> ";
          }
        }

        if (print_begin)
        {
          interval->print_range_begin_on(os);
          show_span = false;    // No longer the last thing we printed.
        }

        ASSERT(!prev_interval || prev_interval->x_range_end() == interval->x_range_begin());

        if (!interval_contains_sign_change)
        {
          if (!show_span)
            os << ", ";
          interval->print_range_end_on(os);
        }
        printing_end_was_suppressed = interval_contains_sign_change;

        prev_interval = interval;
        first = false;
      }
#ifdef CWDEBUG
      if (number_of_intervals > 0)
      {
        if (prev_interval->span_.first->number_of_intervals() != number_of_intervals ||
            prev_interval->span_.first->range().max() != prev_interval->x_range_end())
          fail();
      }
#endif
      if (printing_end_was_suppressed)
        prev_interval->print_range_end_on(os);
      os << "}";
    }

#ifdef CWDEBUG
    bool sanity_check_successful() const
    {
      return sanity_check_successful_;
    }

    void print_on(std::ostream& os) const
    {
      os << "{x_range_:" << x_range_ <<
          ", begin_sign_:" << begin_sign_ <<
          ", end_sign_:" << end_sign_ <<
          ", span_:first:";
      if (span_.first)
        os << *span_.first;
      else
        os << "null";
      os << ", second:";
      if (span_.second)
        os << *span_.second;
      else
        os << "null";
      os << '}';
    }
#endif
  };

  void find_extrema(Range x_range)
  {
    minima_.clear();
    maxima_.clear();

    Interval intervals(x_range, this);
    Dout(dc::notice, utils::print_using(intervals, &Interval::print_list));

    for (int i = 0; i < 6; ++i)
    {
      Interval* next = nullptr;
      for (Interval* interval = &intervals; interval != nullptr; interval = next)
      {
        next = interval->next();
        if (interval->must_be_divided())
        {
          Dout(dc::notice, "Need to divide " << *interval);
          interval->split(this);
          Dout(dc::notice, utils::print_using(intervals, &Interval::print_list));
          ASSERT(intervals.sanity_check_successful());
        }
        else
          Dout(dc::notice, "There was no need to divide " << *interval);
      }
    }

    Dout(dc::notice, "intervals = {" << utils::print_using(intervals, &Interval::print_list) << '}');
    ASSERT(intervals.sanity_check_successful());

    Dout(dc::notice, "This function has " << minima_.size() << " minima and " << maxima_.size() << " maxima.");
  }

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
    vertical_shift_ = desired_amplitude_range.center() - amplitude_scale_ * current.center();

    Dout(dc::notice, "Normalization: current range = " << current <<
         ", amplitude_scale = " << amplitude_scale_ <<
         ", vertical_shift = " << vertical_shift_);
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

  TestFunctionGenerator L(1000, {1.0, 20.0}, {75.0, 100.0}, seed);
  Dout(dc::notice, "Using seed: " << L.get_seed());

  L.advanced_normalize(0.1, 0.1, {-1.0, 1.0});

  gda.enable_drawing(L, -1.0, 1.0);
  while (gda(w, L(w), L.derivative(w)))
  {
    Dout(dc::notice, "-------------------------------------------");
  }

  ASSERT(gda.success());
  gradient_descent::Sample const& result = gda.minimum();
  Dout(dc::notice, "Global minimum: " << result);
}
