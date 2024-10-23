#pragma once

#include "Scale.h"
#include "ExtremeType.h"
#include "math/QuadraticPolynomial.h"
#ifdef CWDEBUG
#include "utils/has_print_on.h"
#endif

namespace gradient_descent {
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

class SampleNode;

class LocalExtreme
{
 public:
  using list_type = std::list<SampleNode>;              // The actual list is member of ExtremeChain.
  using iterator = list_type::iterator;
  using const_iterator = list_type::const_iterator;

 public:
  static constexpr int explored_left = 1;               // Bit mask used for explored_.
  static constexpr int explored_right = 2;

 private:
  ExtremeType local_extreme_;                           // Set to minimum or maximum.
  double extreme_Lw_;                                   // The Lw coordinate of the local extreme (w is stored in the scale_).
  mutable int explored_{0};                             // Bit mask 1: exploration to the left of this extreme has started.
                                                        // Bit mask 2: same, on the right.
  math::Polynomial fourth_degree_approximation_;        // The fourth degree approximation around this local extreme.
  std::array<const_iterator, 2> edge_nodes_;            // The left- and right- most samples used for that fit.

  // Opposite direction data.
  bool opposite_direction_is_fourth_degree_extreme_;    // Set iff opposite_direction_w_/opposite_direction_Lw_ refer to the local
                                                        // extreme of a fourth degree approximation.
  double opposite_direction_w_;                         // The w value of the local extreme of the fourth degree approximation,
                                                        // in the opposite direction.
                                                        // Only valid if opposite_direction_is_fourth_degree_extreme_ is set.
  double opposite_direction_Lw_;                        // Same but Lw coordinate.

  const_iterator left_neighbor_;                        // The adjacent local extreme on the left of us, if any.
  const_iterator right_neighbor_;                       // The adjacent local extreme on the left of us, if any.

#ifdef CWDEBUG
  std::string label_;

  void debug_print_label(char const* left_or_right, const_iterator neighbor) const;
#endif

 public:
  LocalExtreme(ExtremeType local_extreme, double extreme_Lw, const_iterator chain_end COMMA_CWDEBUG_ONLY(std::string label)) :
    local_extreme_(local_extreme), extreme_Lw_(extreme_Lw), fourth_degree_approximation_(5 COMMA_CWDEBUG_ONLY("P")),
    left_neighbor_(chain_end), right_neighbor_(chain_end)
  COMMA_CWDEBUG_ONLY(label_(std::move(label))) { }

  void set_left_neighbor(const_iterator left_neighbor)
  {
    Debug(debug_print_label("left", left_neighbor));
    left_neighbor_ = left_neighbor;
  }
  void set_right_neighbor(const_iterator right_neighbor)
  {
    Debug(debug_print_label("right", right_neighbor));
    right_neighbor_ = right_neighbor;
  }

  void set_edge_nodes(math::Polynomial const& fourth_degree_approximation,
      std::array<LocalExtreme::const_iterator, 7> const& samples, int i0, int i1, int i2);
  const_iterator get_edge_node(HorizontalDirectionToInt dir) const { return edge_nodes_[dir.as_index()]; }
  math::Polynomial const& get_fourth_degree_approximation() const { return fourth_degree_approximation_; }

  ExtremeType get_extreme_type() const { return local_extreme_; }
  double extreme_Lw() const { return extreme_Lw_; }
  const_iterator left_neighbor() const { return left_neighbor_; }
  const_iterator right_neighbor() const { return right_neighbor_; }
  std::string const& label() const { return label_; }

  void set_opposite_direction(bool opposite_direction_is_fourth_degree_extreme, double opposite_direction_w, double opposite_direction_Lw)
  {
    opposite_direction_is_fourth_degree_extreme_ = opposite_direction_is_fourth_degree_extreme;
    opposite_direction_w_ = opposite_direction_w;
    opposite_direction_Lw_ = opposite_direction_Lw;
  }

  double opposite_direction_is_fourth_degree_extreme() const
  {
    return opposite_direction_is_fourth_degree_extreme_;
  }
  double opposite_direction_w() const
  {
    // opposite_direction_w is not required if opposite_direction_is_fourth_degree_extreme_ isn't set.
    ASSERT(opposite_direction_is_fourth_degree_extreme_);
    return opposite_direction_w_;
  }
  double opposite_direction_Lw() const
  {
    // opposite_direction_w is not required if opposite_direction_is_fourth_degree_extreme_ isn't set.
    ASSERT(opposite_direction_is_fourth_degree_extreme_);
    return opposite_direction_Lw_;
  }

  void mark_explored(HorizontalDirection hdirection) const
  {
    DoutEntering(dc::notice, "LocalExtreme::mark_explored(" << hdirection << ") \"" << label() << "\"");
    ASSERT(hdirection != HorizontalDirection::undecided);
    int explore_flag = hdirection == HorizontalDirection::left ? explored_left : explored_right;
    explored_ |= explore_flag;
  }

  bool done() const
  {
    return explored_ == (explored_left|explored_right);
  }

  bool is_explored(HorizontalDirection hdirection) const
  {
    int explore_flag = hdirection == HorizontalDirection::left ? explored_left : explored_right;
    return (explored_ & explore_flag) != 0;
  }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << '[' << local_extreme_ << " (" << extreme_Lw_ << ")]";
    if (explored_)
    {
      os << " E:";
      if (explored_ == explored_left)
        os << "←";
      else if (explored_ == explored_right)
        os << "→";
      else
        os << "↔";
    }
  }
#endif
};

} // namespace gradient_descent
