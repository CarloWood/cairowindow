#pragma once

#include "cairowindow/Range.h"
#include "utils/AIRefCount.h"
#include "utils/UniqueID.h"
#ifdef CWDEBUG
#include "utils/has_print_on.h"
#endif

namespace intervallist {
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

using Range = cairowindow::Range;

struct XSignPair
{
  double x_;
  int sign_;          // The sign of the derivative at x_.

  friend bool operator==(XSignPair const& xsp1, XSignPair const& xsp2)
  {
    return xsp1.x_ == xsp2.x_ && xsp1.sign_ == xsp2.sign_;
  }

  void print_on(std::ostream& os) const
  {
    os << "{x_:" << x_ << ", sign_:" << sign_ << '}';
  }
};

class TestFunctionGenerator;
class IntervalList;
class Interval;

class IntervalNode
{
 protected:
  IntervalNode* prev_{nullptr};
  IntervalNode* next_{nullptr};

 public:
  IntervalNode() = default;
  virtual ~IntervalNode() = default;

  IntervalNode* prev() { return prev_; }
  IntervalNode const* prev() const { return prev_; }
  IntervalNode* next() { return next_; }
  IntervalNode const* next() const { return next_; }

  friend class IntervalList;
  friend class Interval;
};

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

class Interval : public IntervalNode
{
 public:
  using span_type = std::pair<boost::intrusive_ptr<IntervalSpan>, boost::intrusive_ptr<IntervalSpan>>;

 private:
  XSignPair begin_;
  XSignPair end_;
  span_type span_;

 public:
  bool contains_sign_change() const { return begin_.sign_ != end_.sign_; }

  // Create a single linked list.
  Interval(XSignPair begin, XSignPair end, utils::UniqueIDContext<int>& id_context);

  Interval(IntervalSpan* interval_span, XSignPair begin, XSignPair end) :
    begin_{begin}, end_{end}, span_(nullptr, nullptr)
  {
    if (begin.sign_ == end.sign_)
    {
      span_.first = interval_span;
      if (interval_span)
        interval_span->increment();
    }
    else
      span_.second = interval_span;
  }

  double x_range_begin() const { return begin_.x_; }
  double x_range_end() const { return end_.x_; }

  void insert_before(IntervalSpan* interval_span, XSignPair begin, XSignPair end)
  {
    Interval* new_interval = new Interval(interval_span, begin, end);
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
    ASSERT(begin_ == new_interval->begin_);

    new_interval->next_ = this;
    new_interval->prev_ = prev_;

    prev_->next_ = new_interval;
    prev_ = new_interval;

    begin_ = new_interval->end_;
  }

  void insert_after(IntervalSpan* interval_span, XSignPair begin, XSignPair end)
  {
    Interval* new_interval = new Interval(interval_span, begin, end);
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
    ASSERT(end_ == new_interval->end_);

    new_interval->prev_ = this;
    new_interval->next_ = next_;

    next_->prev_ = new_interval;
    next_ = new_interval;

    end_ = new_interval->begin_;
  }

  bool must_be_divided(IntervalList const& intervals) const;

  bool is_first_interval() const
  {
    // Only call this for intervals that do not change sign.
    ASSERT(!contains_sign_change());
    return span_.first->range().min() == begin_.x_;
  }

  bool is_last_interval() const
  {
    // Only call this for intervals that do not change sign.
    ASSERT(!contains_sign_change());
    return span_.first->range().max() == end_.x_;
  }

  void split(XSignPair mid, IntervalList& intervals);

  int begin_sign() const { return begin_.sign_; }
  int end_sign() const { return end_.sign_; }
  span_type const& span() const { return span_; }
  Interval* prev() { return static_cast<Interval*>(prev_); }
  Interval const* prev() const { return static_cast<Interval const*>(prev_); }
  Interval* next() { return static_cast<Interval*>(next_); }
  Interval const* next() const { return static_cast<Interval const*>(next_); }

  void print_range_begin_on(std::ostream& os) const
  {
    os << x_range_begin() << '(' << (begin_sign() == 1 ? '/' : '\\') << ')';
  }

  void print_range_end_on(std::ostream& os) const
  {
    os << x_range_end() << '(' << (end_sign() == 1 ? '/' : '\\') << ')';
  }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << "{begin_:" << begin_ <<
        ", end_:" << end_ <<
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

class IntervalList
{
 private:
  IntervalNode root_;
  utils::UniqueIDContext<int> id_context_;

 public:
  IntervalList()
  {
    root_.next_ = &root_;
    root_.prev_ = &root_;
  }

  utils::UniqueIDContext<int>& id_context() { return id_context_; }
  bool empty() const { return root_.next_ == &root_; }

  bool is_root(IntervalNode const* node) const { return node == &root_; }
  bool is_front(Interval const* interval) const { return interval == root_.next_; }
  bool is_last(Interval const* interval) const { return interval == root_.prev_; }

  void push_back(Interval* interval)
  {
    interval->prev_ = root_.prev_;
    interval->next_ = &root_;
    root_.prev_->next_ = interval;
    root_.prev_ = interval;
  }

  Interval* front()
  {
    return empty() ? nullptr : static_cast<Interval*>(root_.next_);
  }

  Interval const* front() const
  {
    return empty() ? nullptr : static_cast<Interval const*>(root_.next_);
  }

  void print_list(std::ostream& os) const;

#ifdef CWDEBUG
  mutable bool sanity_check_successful_;
  void fail() const
  {
    sanity_check_successful_ = false;
  }

  bool sanity_check_successful() const
  {
    return sanity_check_successful_;
  }
#endif
};

} // namespace intervallist
