#include "sys.h"
#include "Interval.h"

namespace intervallist {

Interval::Interval(XSignPair begin, XSignPair end, utils::UniqueIDContext<int>& id_context) : begin_(begin), end_(end)
{
  if (!contains_sign_change())
  {
    span_.first = new IntervalSpan({begin.x_, end.x_}, id_context);
    span_.first->increment();
  }
}

void IntervalList::print_list(std::ostream& os) const
{
  sanity_check_successful_ = true;
  bool first = true;
  Interval const* prev_interval = nullptr;
  bool printing_end_was_suppressed = false;
  int saw_interval_without_sign_change = 0;
  int number_of_intervals = 0;
  for (Interval const* interval = front(); !is_root(interval); interval = interval->next())
  {
    int const interval_contains_sign_change = interval->contains_sign_change();

#ifdef CWDEBUG
    // Sanity check.
    if (interval_contains_sign_change)
    {
      if (interval->span().first)
      {
        if (!(prev_interval && !prev_interval->contains_sign_change() && prev_interval->span().first == interval->span().first))
          fail();
      }
      else
      {
        if (!(!prev_interval || prev_interval->contains_sign_change()))
          fail();
      }
      Interval const* next_interval = interval->next();
      if (interval->span().second)
      {
        if (!(next_interval && !next_interval->contains_sign_change() && next_interval->span().first == interval->span().second))
          fail();
      }
      else
      {
        if (!(!next_interval || next_interval->contains_sign_change()))
          fail();
      }
      if (number_of_intervals > 0)
      {
        if (!interval->span().first)
          fail();
        else if (!(interval->span().first->number_of_intervals() == number_of_intervals &&
              interval->span().first->range().max() == interval->x_range_begin()))
          fail();
      }
      number_of_intervals = 0;
    }
    else
    {
      if (++number_of_intervals == 1)
      {
        if (!(interval->span().first->range().min() == interval->x_range_begin()))
          fail();
      }
    }
    if (prev_interval)
    {
      if (!(prev_interval->end_sign() == interval->begin_sign() && prev_interval->x_range_end() == interval->x_range_begin()))
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
    if (prev_interval->span().first->number_of_intervals() != number_of_intervals ||
        prev_interval->span().first->range().max() != prev_interval->x_range_end())
      fail();
  }
#endif
  if (printing_end_was_suppressed)
    prev_interval->print_range_end_on(os);
  os << "}";
}

void Interval::split(XSignPair mid, IntervalList& intervals)
{
  if (!contains_sign_change())
  {
    IntervalSpan* const P = span_.first.get();
    if (is_first_interval())
    {
      // Case 1 or 3.
      if (is_last_interval())
      {
        // Case 1.
        //     (optional)   <-P------------->  (optional)
        //                          1
        //  |      ?P       |      P        |      P?       |
        //  -               +               +               -
        //                        this
        //
        if (mid.sign_ == begin_.sign_)
        {
          //  A) finding a +:
          //   (unchanged)   <-P------------>  (unchanged)
          //                    1       2
          // |     ?P       |   P   |   P   |      P?       |
          // -              +       +       +               -
          //                   this
          insert_after(P, mid, end_);
        }
        else
        {
          //  B) finding a -:
          //    (optional)                     (optional)
          //
          // |     ?N       |  NN   |  NN   |      N?       |
          // -              +       -       +               -
          //                   this
          if (!intervals.is_front(this))
          {
            // ?P --> ?N
            prev()->span_.second.reset();
          }
          if (!intervals.is_last(this))
          {
            // P? --> N?
            next()->span_.first.reset();
          }
          insert_after(nullptr, mid, end_);     // NN
          // This was the last interval using this IntervalSpan.
          ASSERT(span_.first->unique().is_true());
          span_.first.reset();
        }
      }
      else
      {
        // Case 3.
        //
        //      (optional)  <-P-----------------------------?
        //                          1               2
        //  |      ?P       |      P        |      P        |
        //  -               +               +               +
        //                        this
        //
        if (mid.sign_ == begin_.sign_)
        {
          //  A) finding a +:
          //    (unchanged)  <-P----------------------------?
          //                     1       2          3
          // |     ?P        |   P   |   P  |       P       |
          // -               +       +      +               +
          //                   this
          insert_after(P, mid, end_);   // P (2)
        }
        else
        {
          //  B) finding a -:
          //                                 <-P-------------?
          //                                         1
          // |     ?N        |  NN   |  NP   |      P        |
          // -               +       -       +               +
          //                   this
          if (!intervals.is_front(this))
          {
            // ?P --> ?N
            prev()->span_.second.reset();
          }
          P->change_begin(end_.x_, P->number_of_intervals() - 1);
          insert_after(P, mid, end_);   // NP
          // P --> NN
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
        //  ?-P----------------------------->  (optional)
        //          n              n+1
        //  |      P        |      P        |       P?      |
        //  +               +               +               -
        //                        this
        //
        if (mid.sign_ == begin_.sign_)
        {
          //  A) finding a +:
          // ?-P----------------------------->  (unchanged)
          //         n          n+1     n+2
          // |      P        |   P   |  P    |      P?       |
          // +               +       +       +               -
          //                   this
          insert_after(P, mid, end_);
        }
        else
        {
          //  B) finding a -:
          // ?-P------------->
          //         n
          // |      P        |   PN  |   NN  |      N?       |
          // +               +       -       +               -
          //                   this
          if (!intervals.is_last(this))
          {
            // P? --> N?
            next()->span_.first.reset();
          }
          P->change_end(begin_.x_, P->number_of_intervals() - 1);
          insert_after(nullptr, mid, end_);     // NN
        }
      }
      else
      {
        // Case 4.
        //
        //  ?-P---------------------------------------------?
        //          n              n+1             n+2
        //  |      P        |      P        |      P        |
        //  +               +               +               +
        //                        this
        //
        if (mid.sign_ == begin_.sign_)
        {
          //  A) finding a +:
          // ?-P---------------------------------------------?
          //         n          n+1     n+2         n+3
          // |      P        |  P    |  P    |      P        |
          // +               +       +       +               +
          //                   this
          insert_after(P, mid, end_);
        }
        else
        {
          //  B) finding a -:
          // ?-P------------->               <-Q-------------?
          //         n                               1
          // |      P        |  PN   |  NQ   |      Q!       |  ...Q!
          // +               +       -       +               +
          //                   this

          // Count the number of intervals to the left of 'this'.
          int left_intervals = 1;   // The first interval.
          for (Interval* left = prev(); left && !left->is_first_interval(); left = left->prev())
            ++left_intervals;
          int right_intervals = P->number_of_intervals() - left_intervals - 1;

          auto Q = new IntervalSpan({end_.x_, span_.first->range().max()}, intervals.id_context());

          // Set the correct end/being and update the number of intervals in both.
          P->change_end(begin_.x_, left_intervals);
          Q->change_begin(end_.x_, right_intervals);

          // Create the new NQ interval (sign change interval).
          insert_after(Q, mid, end_);

          // Update all subsequent intervals to use Q instead of P, including the first sign-changing interval that follows.
          for (Interval* next = this->next()->next(); !intervals.is_root(next); next = next->next())
          {
            next->span_.first = Q;
            if (next->contains_sign_change())
              break;
          }
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
        if (mid.sign_ == begin_.sign_)
        {
          //  A) finding a +:
          //                 <-Q----->
          //                     1
          // |      ?Q       |   Q   |  QN   |      N?       |
          // -               +       +       -               +
          auto Q = new IntervalSpan({begin_.x_, mid.x_}, intervals.id_context());
          if (!intervals.is_front(this))
          {
            // ?N --> ?Q
            prev()->span_.second = Q;
          }
          insert_before(Q, begin_, mid);
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
          auto Q = new IntervalSpan({mid.x_, end_.x_}, intervals.id_context());
          if (!intervals.is_last(this))
          {
            // N? --> Q?
            next()->span_.first = Q;
          }
          insert_after(Q, mid, end_);
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
        if (mid.sign_ == begin_.sign_)
        {
          //  A) finding a +:
          //                 <-Q----->       <-R-------------?
          //                     1                   1
          // |      ?Q       |   Q   |  QR   |      R        |
          // -               +       +       -               -
          auto Q = new IntervalSpan({begin_.x_, mid.x_}, intervals.id_context());
          if (!intervals.is_front(this))
          {
            // ?N --> ?Q
            prev()->span_.second = Q;
          }
          insert_before(Q, begin_, mid);
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
          insert_after(R, mid, end_);
          R->change_begin(mid.x_);
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
        if (mid.sign_ == begin_.sign_)
        {
          //  A) finding a +:
          // ?-L--------------------->
          //         n          n+1
          // |      L        |   L   |  LN   |      N?       |
          // +               +       +       -               +
          insert_before(L, begin_, mid);
          L->change_end(mid.x_);
        }
        else
        {
          //  B) finding a -:
          // ?-L------------->       <-Q----->
          //         n                   1
          // |      L        |  LQ   |   Q   |      Q?       |
          // +               +       -       -               +
          auto Q = new IntervalSpan({mid.x_, end_.x_}, intervals.id_context());
          if (!intervals.is_last(this))
          {
            // N? --> Q?
            next()->span_.first = Q;
          }
          insert_after(Q, mid, end_);
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
        if (mid.sign_ == begin_.sign_)
        {
          //  A) finding a +:
          // ?-L--------------------->       <-R-------------?
          //         n          n+1                 1
          // |      L        |   L   |  LR   |       R       |
          // +               +       +       -               -
          insert_before(L, begin_, mid);
          L->change_end(mid.x_);
        }
        else
        {
          //  B) finding a -:
          // ?-L------------->       <-----------------------?
          //         n                   1          2
          // |      L        |  LR   |   R   |       R       |
          // +               +       -       -               -
          insert_after(R, mid, end_);
          R->change_begin(mid.x_);
        }
      }
    }
  }
}

} // namespace intervallist
