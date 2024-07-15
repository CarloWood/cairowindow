#include "sys.h"
#include "Approximation.h"
#if CW_DEBUG
#include "Algorithm.h"          // Algorithm::uninitialized_magic
#endif

namespace gradient_descent {

void Approximation::add(Sample const* current, bool current_is_replacement, ExtremeType next_extreme_type)
{
  DoutEntering(dc::notice, "Approximation::add(" << *current << ", " <<
      std::boolalpha << current_is_replacement << ", " << next_extreme_type << ")");

  already_had_two_relevant_samples_ = number_of_relevant_samples_ == 2;
  if (current_is_replacement)
  {
    // The sample that was replaced in the history should already be in this approximation.
    ASSERT(number_of_relevant_samples_ > 0);
    relevant_samples_[current_index_] = current;
  }
  else
  {
    current_index_ = 1 - current_index_;
    relevant_samples_[current_index_] = current;
    number_of_relevant_samples_ = already_had_two_relevant_samples_ ? 2 : number_of_relevant_samples_ + 1;
  }

  // Make sure that the samples are ordered, if we have two of them.
  if (number_of_relevant_samples_ == 2 && relevant_samples_[0]->w() > relevant_samples_[1]->w())
  {
    std::swap(relevant_samples_[0], relevant_samples_[1]);
    current_index_ = 1 - current_index_;
  }

  // LocalExtreme should only store parabolas of which the vertex was already determined;
  // hence changing it (by adding a new sample) doesn't make sense.
  ASSERT(!is_extreme_);

  if (number_of_relevant_samples_ == 1)
  {
    // If we have just one point, then the approximation is a linear function:
    //
    // A(w) = current->Lw() + L'(w) * (w - current->w())
    parabola_[2] = 0.0;
    parabola_[1] = current->dLdw();
    parabola_[0] = current->Lw() - current->dLdw() * current->w();
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

    double delta_w = w0 - w1;
    double delta_w_squared = utils::square(delta_w);
    double delta_Lw = Lw0 - Lw1;
    double delta_dLdw = dLdw0 - dLdw1;
    double sum_w = w0 + w1;
    double sum_dLdw = dLdw0 + dLdw1;
    double product_Lw = dLdw0 * dLdw1;

    parabola_[2] = product_Lw < 0.0 ?
      (delta_Lw * sum_dLdw - 2.0 * delta_w * product_Lw) / (delta_w_squared * delta_dLdw) :
      delta_Lw * delta_dLdw / (delta_w_squared * sum_dLdw);

    // If a parabola fit is not calculatable then treat it as a straight line through just the current sample.
    if (std::isnan(parabola_[2]) || std::isinf(parabola_[2]) ||
        std::abs(parabola_[2]) < Scale::epsilon || std::abs(delta_dLdw) < Scale::epsilon || std::abs(sum_dLdw) < Scale::epsilon)
    {
      parabola_[2] = 0.0;
      number_of_relevant_samples_ = 1;        // Ignore the oldest sample.
      scale_.reset();
      return;
    }

    parabola_[1] = delta_Lw / delta_w - sum_w * parabola_[2];
    parabola_[0] = (w0 * Lw1 - w1 * Lw0) / delta_w + w0 * w1 * parabola_[2];

    Dout(dc::notice, "Parabola fit through (" << w0 << ", " << Lw0 << ") with derivative " << dLdw0);
    Dout(dc::notice, "                 and (" << w1 << ", " << Lw1 << ") with derivative " << dLdw1);
    Dout(dc::notice, "    P(x) = " << parabola_ << " --> P(" << w0 << ") = " << parabola_(w0) << " and P(" << w1 << ") = " << parabola_(w1));
    Dout(dc::notice, "    P'(" << w0 << ") = " << parabola_.derivative(w0) << " (requested: " << dLdw0 << ")");
    Dout(dc::notice, "    P'(" << w1 << ") = " << parabola_.derivative(w1) << " (requested: " << dLdw1 << ")");

    cubic_.initialize(w0, Lw0, dLdw0, w1, Lw1, dLdw1);

    Dout(dc::notice, "cubic = " << cubic_);
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
//  restriction:
//         left: looking for an extreme left of `sr`.
//        right: looking for an extreme right of `sl`.
//         none: looking for an extreme, anywhere.
//  extreme_type:
//      maximum: looking for a maximum.
//      minimum: looking for a minimum.
//      unknown: looking for any extreme; the preference will depend on the value of restriction.
//               If restriction is left, the extreme nearest to, but left of, sr is returned.
//               If restriction is right, the extreme nearest to, but right of, sl is returned.
//               if restriction is none, the prefered extreme is the one nearest to the center.
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
//               in this case, if restriction is none, region_out is set to
//               point in the direction where we go downhill. Otherwise
//               it is set to the direction equivalent to restriction.
Weight Approximation::find_extreme(Region& region_out, ExtremeType& extreme_type, Restriction restriction) const
{
  DoutEntering(dc::notice, "Approximation::find_extreme(region_out, " << extreme_type << ", " << restriction << ")");
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
    Dout(dc::notice, "Returning ExtremeType::unknown because D <= 0.");
    if (restriction == Restriction::none)
    {
      region_out = (d > 0.0 || (d == 0.0 && b > 0.0)) ? Region::left : Region::right;
      Dout(dc::notice, "Returning " << region_out << " because that is downhill.");
    }
    extreme_type = ExtremeType::unknown;
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
    // The cubic is, almost, a parabola. Use our parabolic approximation.
    double vertex = parabola_.vertex_x();
    ExtremeType const extreme_type_found = parabola_has_maximum() ? ExtremeType::maximum : ExtremeType::minimum;
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
    if (region_out != restriction)
    {
      Dout(dc::notice, "Returning ExtremeType::unknown because the cubic looks like a parabola with an extreme type (" <<
          extreme_type_found << " on the other side as what is requested (" << restriction << ").");
      extreme_type = ExtremeType::unknown;
#if CW_DEBUG
      return Algorithm::uninitialized_magic;
#else
      return {};  // Not used.
#endif
    }
    extreme_type = extreme_type_found;
    Dout(dc::notice, "Returning " << region_out << " because the cubic looks like a parabola and that's where the vertex is.");
    return vertex;
  }

  if (extreme_type == ExtremeType::unknown)
  {
    //      unknown: looking for any extreme; the preference will depend on the value of restriction.
    //               If restriction is left, the extreme nearest to, but left of, sr is returned.
    //               If restriction is right, the extreme nearest to, but right of, sl is returned.
    //               if restriction is none, the prefered extreme is the one nearest to the center.
    if (restriction == Restriction::none)
    {
      double center = 0.5 * (wl + wr);
      extreme_type = (std::abs(minimum - center) > std::abs(maximum - center)) ? ExtremeType::maximum : ExtremeType::minimum;
      Dout(dc::notice, "Returning " << extreme_type << " because that extreme is closest to the center of the two given samples.");
    }
    else
    {
      double border = (restriction == Restriction::left) ? wr : wl;
      double rel_pos_min = static_cast<int>(restriction) * (minimum - border);
      double rel_pos_max = static_cast<int>(restriction) * (maximum - border);
      if (rel_pos_min < 0.0 && rel_pos_max < 0.0)
      {
        Dout(dc::notice, "Returning ExtremeType::unknown because both extremes are on the wrong side of border (" << border << ").");
        extreme_type = ExtremeType::unknown;
#if CW_DEBUG
        return Algorithm::uninitialized_magic;
#else
        return {};  // Not used.
#endif
      }
      if (rel_pos_min < 0.0)
      {
        extreme_type = ExtremeType::maximum;
        Dout(dc::notice, "Returning " << extreme_type << " because the minimum is on the wrong side of border (" << border << ").");
      }
      else if (rel_pos_max < 0.0)
      {
        extreme_type = ExtremeType::minimum;
        Dout(dc::notice, "Returning " << extreme_type << " because the maximum is on the wrong side of border (" << border << ").");
      }
      else
      {
        extreme_type = (rel_pos_max < rel_pos_min) ? ExtremeType::maximum : ExtremeType::minimum;
        Dout(dc::notice, "Returning " << extreme_type << " because that extreme is closest to border (" << border << ").");
      }
    }
  }

  double result = (extreme_type == ExtremeType::maximum) ? maximum : minimum;

  //         left: looking for an extreme left of `sr`.
  //        right: looking for an extreme right of `sl`.
  if ((restriction == Restriction::left && result >= wr) ||
      (restriction == Restriction::right && result <= wl))
  {
    Dout(dc::notice, "Returning ExtremeType::unknown because the requested extreme is on the wrong side.");
    extreme_type = ExtremeType::unknown;
#if CW_DEBUG
    return Algorithm::uninitialized_magic;
#else
    return {};  // Not used.
#endif
  }

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
