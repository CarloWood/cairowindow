#include "sys.h"
#include "Approximation.h"
#if CW_DEBUG
#include "Algorithm.h"          // Algorithm::uninitialized_magic
#endif

namespace gradient_descent {

void Approximation::add(Sample const* current, bool current_is_replacement, ExtremeType next_extreme_type, bool back_tracking)
{
  DoutEntering(dc::notice, "Approximation::add(" << *current << ", " <<
      std::boolalpha << current_is_replacement << ", " << next_extreme_type << ", " << back_tracking << ")");

//  if (utils::almost_equal(current->w(), -32.4471, 0.001))
//    Debug(attach_gdb());

  Sample const* back_tracking_sample;

  already_had_two_relevant_samples_ = number_of_relevant_samples_ == 2;
  if (current_is_replacement)
  {
    // The sample that was replaced in the history should already be in this approximation.
    ASSERT(number_of_relevant_samples_ > 0);
    // This would be unexpected. Does this happen?
    ASSERT(!back_tracking);
    // Should never end up with twice the same sample.
    ASSERT(!already_had_two_relevant_samples_ || relevant_samples_[1 - current_index_] != current);
    relevant_samples_[current_index_] = current;
  }
  else
  {
    if (back_tracking)
    {
      ASSERT(number_of_relevant_samples_ == 2);
      back_tracking_sample = relevant_samples_[current_index_];
      if (relevant_samples_[1 - current_index_] != &back_tracking_pivot_)
      {
        // Preserve the lifetime of the Sample that is not being replaced during back tracking.
        back_tracking_pivot_ = *relevant_samples_[1 - current_index_];
        relevant_samples_[1 - current_index_] = &back_tracking_pivot_;
      }
    }
    else
      current_index_ = 1 - current_index_;
    // Should never end up with twice the same sample.
    ASSERT(number_of_relevant_samples_ == 0 || relevant_samples_[1 - current_index_] != current);
    relevant_samples_[current_index_] = current;
    number_of_relevant_samples_ = already_had_two_relevant_samples_ ? 2 : number_of_relevant_samples_ + 1;
  }

  // Make sure that the samples are ordered, if we have two of them.
  if (number_of_relevant_samples_ == 2 && relevant_samples_[0]->w() > relevant_samples_[1]->w())
  {
    std::swap(relevant_samples_[0], relevant_samples_[1]);
    current_index_ = 1 - current_index_;
  }

  // LocalExtreme should only store approximation of which the extremes were already determined;
  // hence changing it (by adding a new sample) doesn't make sense.
  ASSERT(!is_extreme_);

  // Can only be back tracking if we already had two relevant samples.
  ASSERT(!back_tracking || (already_had_two_relevant_samples_ && number_of_relevant_samples_ == 2));

  if (number_of_relevant_samples_ == 1)
  {
    // If we have just one point, then the approximation is a linear function:
    //
    // A(w) = current->Lw() + L'(w) * (w - current->w())
    cubic_[3] = 0.0;
    cubic_[2] = 0.0;
    cubic_[1] = current->dLdw();
    cubic_[0] = current->Lw() - current->dLdw() * current->w();
  }
  else
  {
    Sample const* prev = relevant_samples_[1 - current_index_];

    // The theory of this approach is described here:
    // https://math.stackexchange.com/a/4924305/489074
    double w0 = prev->w();
    double Lw0 = prev->Lw();
    double dLdw0 = prev->dLdw();
    double w1 = current->w();
    double Lw1 = current->Lw();
    double dLdw1 = current->dLdw();

    cubic_.initialize(w0, Lw0, dLdw0, w1, Lw1, dLdw1);

    // If a cubic fit is not calculatable then treat it as a straight line through just the current sample.
    if (std::abs(cubic_[3]) < Scale::epsilon && std::abs(cubic_[2]) < Scale::epsilon)
    {
      cubic_[3] = 0.0;
      cubic_[2] = 0.0;
      number_of_relevant_samples_ = 1;        // Ignore the oldest sample.
      scale_.reset();
    }

    Dout(dc::notice, "cubic = " << cubic_);

    if (back_tracking)
    {
      Region region_result;
      ExtremeType extreme_type = next_extreme_type;
      find_extreme(region_result, extreme_type);
      if (extreme_type == ExtremeType::unknown || region_result != Region::inbetween)
      {
        // It shouldn't be on the other side of the two samples than previously predicted.
        ASSERT(extreme_type == ExtremeType::unknown || region_result == (current_index_ == 0 ? Region::left : Region::right));
        // Should never end up with twice the same sample.
        ASSERT(number_of_relevant_samples_ != 2 || relevant_samples_[1 - current_index_] != back_tracking_sample);
        relevant_samples_[current_index_] = back_tracking_sample;
        add(current, current_is_replacement, next_extreme_type, false);
      }
    }
  }
}

ScaleUpdate Approximation::update_scale(bool current_is_replacement, ExtremeType next_extreme_type)
{
  DoutEntering(dc::notice|continued_cf, "Approximation::update_scale(" << std::boolalpha << current_is_replacement << ", " <<
      next_extreme_type << ") --> ");

  ScaleUpdate result = ScaleUpdate::first_sample;
  // Don't update the scale when current was just replaced. In that case we return first_sample,
  // but the caller passed current_is_replacement, so it should know to ignore that.
  if (!current_is_replacement && number_of_relevant_samples_ == 2)
    result = scale_.update(next_extreme_type, relevant_samples_, current_index_, cubic_, already_had_two_relevant_samples_, false);

  Dout(dc::finish, result);
  return result;
}

ScaleUpdate Approximation::update_local_extreme_scale(Sample const& current)
{
  DoutEntering(dc::notice|continued_cf, "Approximation::update_local_extreme_scale(" << current << ") --> ");

  // This function should only be called on Approximation's that are part of a LocalExtreme.
  ASSERT(is_extreme_ && number_of_relevant_samples_ == 2);

  std::array<Sample const*, 2> samples = {{ &current, relevant_samples_[current_index_] }};
  CriticalPointType critical_point_type = scale_.type();
  // Since this Approximation is part of a LocalExtreme, it was used to find an extreme and thus must refer to either a minimum or maximum.
  ASSERT(critical_point_type == CriticalPointType::minimum || critical_point_type == CriticalPointType::maximum);
  ExtremeType current_extreme_type = critical_point_type == CriticalPointType::minimum ? ExtremeType::minimum : ExtremeType::maximum;
  ScaleUpdate result = scale_.update(current_extreme_type, samples, 0, cubic_, true, true);

  Dout(dc::finish, result);
  return result;
}

// find_extreme
//
// Return an extreme, if any, of the cubic approximation.
//
// Let `sl` be the left-most sample and `rl` the right-most sample,
// aka (relevant_samples_[0]) and relevant_samples_[1] respectively,
// because that is how relevant_samples_ is sorted).
//
// Input:
//  extreme_type:
//      maximum: looking for a maximum.
//      minimum: looking for a minimum.
//      unknown: looking for any extreme; the prefered extreme is the one nearest to the center.
//
// Output:
//  region_out:
//         left: the requested extreme is left of `sl`.
//        right: the requested extreme is right of `sr`.
//    inbetween: the requested extreme is in between `sl` and `sr`.
//  extreme_type:
//           up: the extreme is a maximum.
//         down: the extreme is a minimum.
//      unknown: there is no extreme (of the requested `extreme_type`),
//               region_out is set to point in the direction where we go downhill.
Weight Approximation::find_extreme(Region& region_out, ExtremeType& extreme_type) const
{
  DoutEntering(dc::notice, "Approximation::find_extreme(region_out, " << extreme_type << ")");
  // Otherwise the cubic approximation isn't even valid.
  ASSERT(number_of_relevant_samples_ == 2);

  double a = cubic_[0];
  double b = cubic_[1];
  double c = cubic_[2];
  double d = cubic_[3];

  double D = utils::square(c) - 3.0 * b * d;

  // If D is equal to zero we have a point where the derivative is zero, but that isn't a maximum or minimum.
  if (D <= 0.0)
  {
    region_out = (d > 0.0 || (d == 0.0 && b > 0.0)) != (extreme_type == ExtremeType::maximum) ? Region::left : Region::right;
    Dout(dc::notice, "Returning " << region_out << " because that is " << (extreme_type == ExtremeType::maximum ? "uphill" : "downhill") << ".");
    extreme_type = ExtremeType::unknown;
    Dout(dc::notice, "Returning ExtremeType::unknown because D <= 0.");
#if CW_DEBUG
    return Algorithm::uninitialized_magic;
#else
    return {};  // Not used.
#endif
  }

  double sqrt_D = std::sqrt(D);
  double minimum = (-c + sqrt_D) / (3.0 * d);
  double maximum = (-c - sqrt_D) / (3.0 * d);

  double wl = relevant_samples_[0]->w();
  double wr = relevant_samples_[1]->w();

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
      return Algorithm::uninitialized_magic;
#else
      return {};  // Not used.
#endif
    }
    region_out =
      (vertex < wl) ? Region::left : (vertex > wr) ? Region::right : Region::inbetween;
    extreme_type = extreme_type_found;
    Dout(dc::notice, "Returning " << region_out << " because the cubic looks like a parabola and that's where the vertex is.");
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

  //         left: the requested extreme is left of `left`.
  //        right: the requested extreme is right of `right`.
  //    inbetween: the requested extreme is in between `left` and `right`.
  if (result < wl)
    region_out = Region::left;
  else if (result > wr)
    region_out = Region::right;
  else
    region_out = Region::inbetween;
  Dout(dc::notice, "Returning " << region_out << " because this is where the required extreme (" << extreme_type << ") is.");

  // An extreme_type of unknown means that the returned value should be ignored.
  ASSERT(extreme_type != ExtremeType::unknown);
  return result;
}

void Approximation::set_current_index(Region region)
{
  // If region is left/right, then the left/right-most sample must be current.
  current_index_ = ((relevant_samples_[0]->w() < relevant_samples_[1]->w()) == (region == Region::left)) ? 0 : 1;
}

} // namespace gradient_descent
