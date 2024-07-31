#include "sys.h"
#include "SampleNode.h"
#include "AlgorithmEventType.h"

namespace gradient_descent {

void SampleNode::initialize_cubic(SampleNode const& next
    COMMA_CWDEBUG_ONLY(events::Server<AlgorithmEventType>& event_server, bool this_is_last)) const
{
  DoutEntering(dc::notice, "SampleNode::initialize_cubic(" << next << "..., " << std::boolalpha << this_is_last << ")");

  // next must always be on the right.
  ASSERT(w() < next.w());

  cubic_.initialize(w(), Lw(), dLdw(), next.w(), next.Lw(), next.dLdw());

#ifdef CWDEBUG
  // If the current object (this) was the last sample that was added, then the cubic is on the right
  // of the last sample. We use a different color for that. Otherwise assume next was the last.
  HorizontalDirection side = this_is_last ? HorizontalDirection::right : HorizontalDirection::left;
  event_server.trigger(AlgorithmEventType{cubic_polynomial_event, cubic_, side});
#endif

  // Get the sign of the derivatives.
  int const sign_dLdw_0 = dLdw() < 0.0 ? -1 : dLdw() > 0.0 ? 1 : 0;
  int const sign_dLdw_1 = next.dLdw() < 0.0 ? -1 : next.dLdw() > 0.0 ? 1 : 0;
  bool const neither_derivative_is_zero = (sign_dLdw_0 & sign_dLdw_1) != 0;

  using enum CubicToNextSampleType;

  // The easiest first.
  if (sign_dLdw_0 != sign_dLdw_1)
  {
    // This must be /\, \/ -- or, very unlikely, one of the derivatives is zero.
    if (AI_LIKELY(neither_derivative_is_zero))
    {
      type_ = (sign_dLdw_0 == 1)
          ? max                         // /\.
          : min;                        // \/
    }
    else
    {
      // One of the derivatives is zero!
      if (sign_dLdw_0 == 0)             // The right part of an extreme.
      {
        double d2Ldw2 = cubic_.second_derivative(w());
        if (sign_dLdw_1 == 1)
        {
          type_ = (d2Ldw2 > 0.0)
              ? right_min               // _/
              : right_max_min;          // ‾\/
        }
        else
        {
          type_ = (d2Ldw2 > 0.0)
              ? right_min_max           // _/\.
              : right_max;              // ‾\.
        }
      }
      else                              // The left part of an extreme.
      {
        double d2Ldw2 = cubic_.second_derivative(next.w());
        if (sign_dLdw_0 == 1)
        {
          type_ = (d2Ldw2 > 0.0)
              ? max_left_min          // /\_
              : left_max;             // /‾
        }
        else
        {
          type_ = (d2Ldw2 > 0.0)
              ? left_min              // \_
              : min_left_max;         // \/‾
        }
      }
    }
  }
  else
  {
    // The signs are equal, so this must be /, \, /\/, \/\ -- or, very unlikely, both derivatives are zero.
    if (AI_LIKELY(neither_derivative_is_zero))
    {
      // The signs are equal (and non-zero). Now we have to calculate the extremes of the cubic.
      // Put the minimum (or inflection point) in index 0.
      std::array<double, 2> extremes;
      int number_of_extremes = cubic_.get_extremes(extremes);

      if (number_of_extremes < 2 || extremes[0] < w() || next.w() < extremes[0])
      {
        // If only one extreme would fall in between the samples, then the derivatives can't have the same sign.
        ASSERT(number_of_extremes < 2 || extremes[1] < w() || next.w() < extremes[1]);
        // There are no extremes in between the samples.
        type_ = (sign_dLdw_0 == 1)
            ? up                        // /
            : down;                     // \.
      }
      else
      {
        // Both extremes are between the two samples.
        ASSERT(w() < extremes[1] && extremes[1] < next.w());
        type_ = (sign_dLdw_0 == 1)
            ? max_min                   // /\/
            : min_max;                  // \/\.
      }
    }
    else
    {
      // Both derivatives are zero!
      if (Lw() > next.Lw())             // ‾\_
        type_ = right_max_left_min;
      else if (Lw() < next.Lw())        // _/‾
        type_ = right_min_left_max;
      else
        type_ = flat;                   // __
    }
  }

  Dout(dc::notice, "Result: " << *this);
}

// find_extreme
//
// Return an extreme, if any, of the cubic approximation.
//
// Let `*this` be the left-most sample and `next` the right-most sample.
//
// Input:
//  next:
//      The right-most sample that this cubic is based on.
//
//  extreme_type:
//      maximum: looking for a maximum.
//      minimum: looking for a minimum.
//      unknown: looking for any extreme; the prefered extreme is the one nearest to the center.
//
// Output:
//  extreme_type:
//      maximum: the extreme is a maximum.
//      minimum: the extreme is a minimum.
//      unknown: there is no extreme (of the requested `extreme_type`).
double SampleNode::find_extreme(Sample const& next, ExtremeType& extreme_type) const
{
  DoutEntering(dc::notice, "SampleNode::find_extreme(" << next << ", " << extreme_type << ") [" << *this << "]");

  // Next must be the right-most sample.
  ASSERT(next.w() > w());

#ifdef CWDEBUG
  // The cubic must be based on this sample and next.
  {
    double const Lw_epsilon = std::abs(next.Lw() - Lw()) * 1e-6;
    double const dLdw_epsilon =
      std::max(std::abs((next.Lw() - Lw()) / (next.w() - w())), std::max(std::abs(next.dLdw()), std::abs(dLdw()))) * 1e-6;
    ASSERT(std::abs(cubic_(w()) - Lw()) < Lw_epsilon && std::abs(cubic_.derivative(w()) - dLdw()) < dLdw_epsilon);
    ASSERT(std::abs(cubic_(next.w()) - next.Lw()) < Lw_epsilon && std::abs(cubic_.derivative(next.w()) - next.dLdw()) < dLdw_epsilon);
  }
#endif

  double a = cubic_[0];
  double b = cubic_[1];
  double c = cubic_[2];
  double d = cubic_[3];

  double D = utils::square(c) - 3.0 * b * d;

  // If D is equal to zero we have a point where the derivative is zero, but that isn't a maximum or minimum.
  if (D <= 0.0)
  {
    extreme_type = ExtremeType::unknown;
    Dout(dc::notice, "Returning ExtremeType::unknown because D <= 0.");
#if CW_DEBUG
    return uninitialized_magic;
#else
    return {};  // Not used.
#endif
  }

  double sqrt_D = std::sqrt(D);
  double minimum = (-c + sqrt_D) / (3.0 * d);
  double maximum = (-c - sqrt_D) / (3.0 * d);

  double wl = w();
  double wr = next.w();

  if (std::isnan(minimum) || std::isnan(maximum) || utils::almost_equal(sqrt_D, std::abs(c), 1e-9))
  {
    Dout(dc::notice, "cubic_ = " << cubic_ << "; using parabolic approximation.");
    // The cubic is, almost, a parabola. Use a parabolic approximation.
    double vertex = -0.5 * cubic_[1] / cubic_[2];
    ExtremeType const extreme_type_found = cubic_[2] < 0.0 ? ExtremeType::maximum : ExtremeType::minimum;
    if (extreme_type != ExtremeType::unknown && extreme_type_found != extreme_type)
    {
      Dout(dc::notice, "Returning ExtremeType::unknown because the cubic looks like a parabola with an extreme type (" <<
          extreme_type_found << ") different from what is requested (" << extreme_type << ").");
      extreme_type = ExtremeType::unknown;
#if CW_DEBUG
      return uninitialized_magic;
#else
      return {};  // Not used.
#endif
    }
    extreme_type = extreme_type_found;
    return vertex;
  }

  if (extreme_type == ExtremeType::unknown)
  {
    //      unknown: looking for any extreme; the prefered extreme is the one nearest to the center.
    double center = 0.5 * (wl + wr);
    extreme_type = (std::abs(minimum - center) > std::abs(maximum - center)) ? ExtremeType::maximum : ExtremeType::minimum;
    Dout(dc::notice, "Returning " << extreme_type << " because that extreme is closest to the center of the two given samples.");
  }

  double result = (extreme_type == ExtremeType::maximum) ? maximum : minimum;

  // An extreme_type of unknown means that the returned value should be ignored.
  ASSERT(extreme_type != ExtremeType::unknown);
  return result;
}

} // namespace gradient_descent
