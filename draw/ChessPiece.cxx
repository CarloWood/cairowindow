#include "sys.h"
#include "ChessPiece.h"
#ifdef CWDEBUG
#include "cairowindow/debugcairo.h"
#endif

namespace cairowindow::draw {

static constexpr double black_line_width = 0.026;
static constexpr double white_line_width = 1.5 * black_line_width;

static double snap_bottom(double y, double translate, double scale, double line_width)
{
  if (scale < 27)
    return y;
  return (std::round((y + 0.5 * line_width) * scale - translate) + translate) / scale - 0.5 * line_width;
}

static double snap_top(double y, double translate, double scale, double line_width)
{
  if (scale < 27)
    return y;
  return (std::round((y - 0.5 * line_width) * scale - translate) + translate) / scale + 0.5 * line_width;
}

static double snap_line_width(double line_width, double scale)
{
  if (line_width * scale < 1.0)
    return line_width;
  return std::trunc(line_width * scale + 0.3) / scale;
}

// The fill color of the pieces.
void ChessPiece::set_fill_color(cairo_t* cr, chess::EColor color)
{
#ifdef CWDEBUG
  using namespace debugcairo;
#endif
  if (color == chess::white)
    cairo_set_source_rgb(cr, style_.white_piece_fill_color().red(),
                             style_.white_piece_fill_color().green(),
                             style_.white_piece_fill_color().blue());
  else
    cairo_set_source_rgb(cr, style_.black_piece_fill_color().red(),
                             style_.black_piece_fill_color().green(),
                             style_.black_piece_fill_color().blue());
}

// The line color of the pieces.
void ChessPiece::set_line_color(cairo_t* cr, chess::EColor color)
{
#ifdef CWDEBUG
  using namespace debugcairo;
#endif
  if (color == chess::white)
    cairo_set_source_rgb(cr, style_.white_piece_line_color().red(),
                             style_.white_piece_line_color().green(),
                             style_.white_piece_line_color().blue());
  else
    cairo_set_source_rgb(cr, style_.black_piece_line_color().red(),
                             style_.black_piece_line_color().green(),
                             style_.black_piece_line_color().blue());
}

void ChessPiece::draw_pawn(cairo_t* cr, double x, double y, double scale, chess::EColor color)
{
  DoutEntering(dc::cairowindow, "ChessPiece::draw_pawn(" << cr << ", " << x << ", " << y << ", " << scale << ", " << color << ")");

#ifdef CWDEBUG
  using namespace debugcairo;
#endif
  static double const base_outside_diameter_cm = 3.265;
  static double const width_pawn_cm = 5.31;
  static double const base_radius = 0.5 * (base_outside_diameter_cm / width_pawn_cm - black_line_width);
  static double const mid_outside_diameter_cm = 1.98;
  static double const mid_radius = 0.5 * (mid_outside_diameter_cm / width_pawn_cm - black_line_width);
  static double const head_outside_diameter_cm = 1.12;
  static double const head_radius = 0.5 * (head_outside_diameter_cm / width_pawn_cm - black_line_width);
  static double const height_pawn_cm = 5.43;
  static double const bottom_pawn_cm = 0.58;
  static double const foot_height = 0.0387;
  static double const base_y = 0.5 - bottom_pawn_cm / height_pawn_cm - 0.5 * black_line_width;
  static double const base_scale = 0.931;
  static double const mid_y = -0.0545;
  static double const top_offset_cm = 0.62;
  static double const head_y = -0.5 + top_offset_cm / height_pawn_cm + 0.5 * black_line_width + head_radius;

  static double const base_angle = 1.148;
  static double const mid_angle1 = 0.992;
  static double const inner_neck_width_cm = 0.41;
  static double const neck_right = 0.5 * (inner_neck_width_cm / width_pawn_cm + black_line_width);
  static double const head_angle = asin(neck_right / head_radius);
  static double const mid_scale = (mid_y - (head_y + head_radius * cos(head_angle)) -
      0.1 * black_line_width) / sqrt(mid_radius * mid_radius - neck_right * neck_right);
  static double const mid_angle2 = asin(head_radius * sin(head_angle) / mid_radius);

  double const base_y_sn = snap_bottom(base_y, y, scale, black_line_width);

  cairo_save(cr);
  cairo_translate(cr, x + square_size_ / 2, y + square_size_ / 2);
  cairo_scale(cr, scale, scale);
  cairo_set_line_width(cr, black_line_width);

  // Draw the left-side of the base.
  cairo_move_to(cr, -base_radius, base_y_sn);
  cairo_save(cr);
  cairo_translate(cr, 0.0, base_y_sn - foot_height);
  cairo_scale(cr, 1.0, base_scale);
  cairo_arc(cr, 0.0, 0.0, base_radius, -M_PI, -M_PI + base_angle);
  cairo_restore(cr);

  // Draw the left-side of the mid-section.
  cairo_save(cr);
  cairo_translate(cr, 0.0, mid_y);
  cairo_scale(cr, 1.0, mid_scale);
  cairo_arc(cr, 0.0, 0.0, mid_radius, -M_PI - mid_angle1, -0.5 * M_PI - mid_angle2);
  cairo_restore(cr);

  // Draw the head of the pawn.
  cairo_arc(cr, 0.0, head_y, head_radius, -1.5 * M_PI + head_angle, 0.5 * M_PI - head_angle);

  // Draw the right-side of the mid-section.
  cairo_save(cr);
  cairo_translate(cr, 0.0, mid_y);
  cairo_scale(cr, 1.0, mid_scale);
  cairo_arc(cr, 0.0, 0.0, mid_radius,  -0.5 * M_PI + mid_angle2, mid_angle1);
  cairo_restore(cr);

  // Draw the right-side of the base.
  cairo_save(cr);
  cairo_translate(cr, 0.0, base_y_sn - foot_height);
  cairo_scale(cr, 1.0, base_scale);
  cairo_arc(cr, 0.0, 0.0, base_radius, -base_angle, 0.0);
  cairo_restore(cr);
  cairo_line_to(cr, base_radius, base_y_sn);

  // Draw the base line of the pawn, right to left.
  cairo_close_path(cr);

  set_fill_color(cr, color);
  cairo_fill_preserve(cr);
  if (color == chess::white)
    set_line_color(cr, color);
  double x1, y1, x2, y2;
  cairo_stroke_extents(cr, &x1, &y1, &x2, &y2);
  cairo_stroke(cr);

  cairo_restore(cr);
}

void ChessPiece::draw_king(cairo_t* cr, double x, double y, double scale, chess::EColor color)
{
  DoutEntering(dc::cairowindow, "ChessPiece::draw_king(" << cr << ", " << x << ", " << y << ", " << scale << ", " << color << ")");

#ifdef CWDEBUG
  using namespace debugcairo;
#endif
  // Measurements, in cm.
  static double const blob_left_cm = 1.22;
  static double const band_edge_left_cm = 2.55;
  static double const band_left_cm = 2.67;
  static double const inside_left_cm = 3.06;
  static double const center_blob_left_cm = 4.525;
  static double const cross_left_cm = 4.71;
  static double const width_king_cm = 10.67;
  static double const bottom_king_cm = 1.155;
  static double const band_line_top_cm = 2.95;
  static double const band_top_king_cm = 4.04;
  static double const center_y_cm = 5.02;
  static double const blob_top_cm = 7.4; // 7.28
  static double const center_blob_top_cm = 8.18;
  static double const cross_y_king_cm = 9.17;
  static double const cross_top_cm = 9.86;
  static double const height_king_cm = 10.86;
  // Derived values.
  static double const mid_x_king_cm = width_king_cm / 2;
  static double const mid_y_king_cm = height_king_cm / 2;

  // Same, in coordinates.
  static double const blob_left = (blob_left_cm - mid_x_king_cm) / width_king_cm;
  static double const band_edge_left = (band_edge_left_cm - mid_x_king_cm) / width_king_cm;
  static double const band_left = (band_left_cm - mid_x_king_cm) / width_king_cm;
  static double const inside_left = (inside_left_cm - mid_x_king_cm) / width_king_cm;
  static double const center_blob_left = (center_blob_left_cm - mid_x_king_cm) / width_king_cm;
  static double const cross_left = (cross_left_cm - mid_x_king_cm) / width_king_cm;
  static double const bottom_king = (mid_y_king_cm - bottom_king_cm) / height_king_cm;
  static double const band_line_top = (mid_y_king_cm - band_line_top_cm) / height_king_cm;
  static double const band_top_king = (mid_y_king_cm - band_top_king_cm) / height_king_cm;
  static double const center_y = (mid_y_king_cm - center_y_cm) / height_king_cm;
  static double const blob_top = (mid_y_king_cm - blob_top_cm) / height_king_cm;
  static double const center_blob_top = (mid_y_king_cm - center_blob_top_cm) / height_king_cm;
  static double const cross_y_king = (mid_y_king_cm - cross_y_king_cm) / height_king_cm;
  static double const cross_top = (mid_y_king_cm - cross_top_cm) / height_king_cm;

  // Derived values.
  static double const inside_radius_king = -inside_left;
  static double const inside_scale_king = 0.180132; // Same value as used for the queen.
  static double const band_top_radius = -band_edge_left;
  static double const band_top_scale = inside_scale_king;
  static double const band_top_y = band_top_king + band_top_radius * band_top_scale;
  static double const cos_alpha = band_left / band_edge_left;
  static double const alpha = acos(cos_alpha);
  static double const band_bottom_scale = inside_scale_king;
  static double const band_bottom_radius = band_top_radius;
  static double const band_bottom_y = bottom_king - band_bottom_radius * band_bottom_scale;
  static double const dx = band_top_radius * (1.0 - cos_alpha);
  static double const band_line_scale = band_top_scale;
  static double const band_line_radius = band_top_radius - dx;
  static double const band_line_y = band_line_top + band_line_radius * band_line_scale;
  static double const blob_radius = 0.7071067 * (blob_left + band_top_y - band_left - blob_top);
  static double const blob_x = blob_left + blob_radius;
  static double const blob_y = blob_top + blob_radius;
  static double const center_blob_radius = -center_blob_left;
  static double const center_blob_y = center_blob_top + center_blob_radius;
  // Manual adjustment... looks better.
  static double const adjusted_center_blob_radius = center_blob_radius + 0.01;
  static double const beta_king = asin(adjusted_center_blob_radius / (center_y - center_blob_y));
  static double const center2_y = blob_y - blob_x - 1.4142136 * blob_radius;

  cairo_save(cr);
  cairo_translate(cr, x + square_size_ / 2, y + square_size_ / 2);
  cairo_scale(cr, scale, scale);
  cairo_set_line_width(cr, black_line_width);

  // Draw left blob.
  cairo_move_to(cr, band_left, band_top_y);
  cairo_arc(cr, blob_x, blob_y, blob_radius, 0.75 * M_PI, 1.75 * M_PI);
  cairo_line_to(cr, 0.0, center2_y);

  // Draw right blob.
  cairo_arc(cr, -blob_x, blob_y, blob_radius, -0.75 * M_PI, 0.25 * M_PI);
  cairo_line_to(cr, -band_left, band_top_y);

  set_fill_color(cr, color);
  cairo_fill_preserve(cr);

  // Draw vertical line in the middle.
  cairo_move_to(cr, 0.0, band_top_y);
  cairo_line_to(cr, 0.0, center_y);

  if (color == chess::white)
    set_line_color(cr, color);
  cairo_stroke(cr);

  // Draw center blob.
  cairo_move_to(cr, 0.0, center_y);
  cairo_arc(cr, 0.0, center_blob_y, adjusted_center_blob_radius, M_PI - beta_king, beta_king);
  cairo_close_path(cr);

  if (color == chess::white)
    set_fill_color(cr, color);
  cairo_fill_preserve(cr);
  if (color == chess::white)
    set_line_color(cr, color);
  cairo_stroke(cr);

  // Draw cross.
  cairo_move_to(cr, 0.0, center_blob_y - adjusted_center_blob_radius);
  cairo_line_to(cr, 0.0, cross_top);
  cairo_move_to(cr, cross_left, cross_y_king);
  cairo_line_to(cr, -cross_left, cross_y_king);
  cairo_stroke(cr);

  // Draw half ellipse just below the blobs.
  cairo_save(cr);
  cairo_translate(cr, 0.0, band_top_y);
  cairo_scale(cr, 1.0, band_top_scale);
  cairo_arc(cr, 0.0, 0.0, band_top_radius, M_PI - alpha, 2 * M_PI + alpha);
  cairo_restore(cr);

  // Draw right side of the upper band.
  cairo_line_to(cr, -band_left, band_line_y);

  // Draw right side of lower band and bottom.
  cairo_save(cr);
  cairo_translate(cr, 0.0, band_bottom_y);
  cairo_scale(cr, 1.0, band_bottom_scale);
  cairo_arc(cr, 0.0, 0.0, band_bottom_radius, 0.0, M_PI);
  cairo_restore(cr);

  // Draw left side of lower band.
  cairo_line_to(cr, band_left, band_line_y);

  // Draw left side of upper band.
  cairo_close_path(cr);

  cairo_path_t* path = cairo_copy_path(cr);

  if (color == chess::white)
    set_fill_color(cr, color);
  cairo_fill(cr);

  // Draw lower half ellipse of upper band.
  cairo_save(cr);
  cairo_translate(cr, 0.0, band_line_y);
  cairo_scale(cr, 1.0, band_line_scale);
  cairo_arc(cr, 0.0, 0.0, band_line_radius, -M_PI, 0.0);
  cairo_restore(cr);

  cairo_new_sub_path(cr);

  // Draw opening at bottom, align it with the real bottom.
  cairo_save(cr);
  cairo_translate(cr, 0.0, band_bottom_y + band_bottom_radius * band_bottom_scale - inside_radius_king * inside_scale_king);
  cairo_scale(cr, 1.0, inside_scale_king);
  if (color == chess::white)
    cairo_arc(cr, 0.0, 0.0, inside_radius_king, -M_PI, M_PI);
  else
    cairo_arc(cr, 0.0, 0.0, inside_radius_king, -M_PI - alpha, alpha);
  cairo_restore(cr);

  set_line_color(cr, color);
  if (color == chess::white)
    cairo_stroke(cr);
  else
  {
    cairo_set_line_width(cr, white_line_width);
    cairo_stroke(cr);
    cairo_set_line_width(cr, black_line_width);
  }

  cairo_append_path(cr, path);
  if (color == chess::black)
    set_fill_color(cr, color);
  cairo_stroke(cr);

  cairo_path_destroy(path);

  if (color == chess::black)
  {
    // Draw the white lines inside the blobs.

    static double const av_line_width = 0.5 * (black_line_width + white_line_width);
    static double const da = av_line_width / band_top_radius;
    static double const dy = av_line_width * tan(0.5 * beta_king);

    cairo_save(cr);
    cairo_translate(cr, 0.0, band_top_y);
    cairo_scale(cr, 1.0, band_top_scale);
    cairo_arc_negative(cr, 0.0, 0.0, band_top_radius, -0.5 * M_PI - da, M_PI + alpha + da);
    cairo_restore(cr);

    cairo_arc(cr, blob_x, blob_y, blob_radius - av_line_width, 0.75 * M_PI, 1.75 * M_PI);

    static double const center2b_y = center2_y + av_line_width * 1.4142136;
    static double const sin_beta = adjusted_center_blob_radius / (center_y - center_blob_y);
    static double const x_king = sin_beta * (center_y + av_line_width / sin_beta - center2b_y) / sin(0.25 * M_PI - beta_king);
    static double const y_king = center2b_y - x_king;

    cairo_line_to(cr, -x_king, y_king);
    cairo_line_to(cr, -av_line_width, center_y + dy);

    cairo_close_path(cr);

    cairo_new_sub_path(cr);

    cairo_save(cr);
    cairo_translate(cr, 0.0, band_top_y);
    cairo_scale(cr, 1.0, band_top_scale);
    cairo_arc_negative(cr, 0.0, 0.0, band_top_radius, -alpha - da, -0.5 * M_PI + da);
    cairo_restore(cr);

    cairo_line_to(cr, av_line_width, center_y + dy);
    cairo_line_to(cr, x_king, y_king);

    cairo_arc(cr, -blob_x, blob_y, blob_radius - av_line_width, -0.75 * M_PI, 0.25 * M_PI);

    cairo_close_path(cr);

    cairo_new_sub_path(cr);

    cairo_move_to(cr, 0.0, center_y - av_line_width / sin_beta);
    cairo_arc(cr, 0.0, center_blob_y, adjusted_center_blob_radius - av_line_width, M_PI - beta_king, beta_king);

    cairo_close_path(cr);

    set_line_color(cr, color);
    cairo_set_line_width(cr, white_line_width);
    cairo_stroke(cr);
  }

  cairo_restore(cr);
}

void ChessPiece::draw_queen(cairo_t* cr, double x, double y, double scale, chess::EColor color)
{
  DoutEntering(dc::cairowindow, "ChessPiece::draw_queen(" << cr << ", " << x << ", " << y << ", " << scale << ", " << color << ")");

#ifdef CWDEBUG
  using namespace debugcairo;
#endif
  // Measurements, in cm.
  static double const width_queen_cm = 5.34;
  static double const inside_width_cm = 2.97;
  static double const band1_width_cm = 2.59;
  static double const crown_bottom_width_cm = 3.31;
  static double const height_queen_cm = 5.39;
  static double const bottom_queen_cm = 0.5;
  static double const inside_height_cm = 0.54;
  static double const band1_height_cm = 0.47;
  static double const band2_height_cm = 0.43;
  static double const tooth_outside_cm = 1.83;
  static double const tooth_inside_cm = 2.20;
  static double const tooth_inside2_cm = 2.36;
  static double const ball_outside_diameter_cm = 0.6;
  static double const ball_top1_cm = 4.31;
  static double const ball_right1_cm = 0.90;
  static double const ball_top2_cm = 4.80;
  static double const ball_right2_cm = 1.88;
  static double const tooth3_x_cm = 2.25;
  // Derived values.
  static double const mid_x_queen_cm = width_queen_cm / 2;
  static double const mid_y_queen_cm = height_queen_cm / 2;

  // Same, in coordinates.
  static double const inside_width = inside_width_cm / width_queen_cm;
  static double const band1_width = band1_width_cm / width_queen_cm;
  static double const crown_bottom_width = crown_bottom_width_cm / width_queen_cm;
  static double const bottom_queen = (mid_y_queen_cm - bottom_queen_cm) / height_queen_cm;
  static double const inside_height = inside_height_cm / height_queen_cm;
  static double const band1_height = band1_height_cm / height_queen_cm;
  static double const band2_height = band2_height_cm / height_queen_cm;
  static double const tooth_outside = (mid_y_queen_cm - tooth_outside_cm) / height_queen_cm;
  static double const tooth_inside = (mid_y_queen_cm - tooth_inside_cm) / height_queen_cm;
  static double const tooth_inside2 = (mid_y_queen_cm - tooth_inside2_cm) / height_queen_cm;
  static double const ball_outside_diameter = ball_outside_diameter_cm / height_queen_cm;
  static double const ball_top1 = (mid_y_queen_cm - ball_top1_cm) / height_queen_cm;
  static double const ball_right1 = (ball_right1_cm - mid_x_queen_cm) / width_queen_cm;
  static double const ball_top2 = (mid_y_queen_cm - ball_top2_cm) / height_queen_cm;
  static double const ball_right2 = (ball_right2_cm - mid_x_queen_cm) / width_queen_cm;
  static double const tooth3_x = (tooth3_x_cm - mid_x_queen_cm) / width_queen_cm;

  // Derived values.
  static double const inside_radius_queen = inside_width / 2;
  static double const inside_scale_queen = inside_height / inside_width;
  static double const inside_y_queen = bottom_queen - inside_radius_queen * inside_scale_queen;
  static double const band1_radius = band1_width / 2;
  static double const band1_scale = inside_scale_queen;
  static double const band1_y = bottom_queen - inside_height - band1_height + band1_radius * band1_scale;
  static double const crown_bottom_left = -crown_bottom_width / 2;
  static double const band2_radius = band1_radius + (-band1_radius - crown_bottom_left) * band2_height / (band1_y - tooth_outside);
  static double const band2_scale = band1_scale;
  static double const band2_y = bottom_queen - inside_height - band1_height - band2_height + band2_radius * band2_scale;
  static double const ball1_x = ball_right1 - ball_outside_diameter / 2;
  static double const ball2_x = ball_right2 - ball_outside_diameter / 2;
  static double const ball1_y = ball_top1 + ball_outside_diameter / 2;
  static double const ball2_y = ball_top2 + ball_outside_diameter / 2;
  static double const ball_radius_queen = (ball_outside_diameter - black_line_width) / 2;
  // We calculate ball3_y, so it lays on a perfect circle with the other balls.
  // The distance from ballN to a point (0, ball_center_y) is:
  // sqrt(ballN_x^2 + (ballN_y - ball_center_y)^2), and we find
  // ball_center_y by setting this distance equal for ball1 and 2:
  // ball1_x^2 + ball1_y^2 - 2 ball1_y ball_center_y = ball2_x^2 + ball2_y^2 - 2 ball2_y ball_center_y -->
  static double const ball_center_y = 0.5 * (ball2_x * ball2_x + ball2_y * ball2_y - ball1_x * ball1_x - ball1_y * ball1_y) / (ball2_y - ball1_y);
  static double const ball3_y = ball_center_y - sqrt(ball1_x * ball1_x + (ball1_y - ball_center_y) * (ball1_y - ball_center_y));
  // The tooth points are derived (which turns out better than measuring them).
  static double const ball1_angle = atan((0.5 * (crown_bottom_left + ball2_x) - ball1_x) / (tooth_outside - ball1_y));
  static double const tooth1_x = ball1_x + ball_radius_queen * sin(ball1_angle);
  static double const tooth2_x = ball2_x;
  static double const tooth1_top = ball1_y + ball_radius_queen * cos(ball1_angle);
  static double const tooth2_top = ball2_y + ball_radius_queen;
  static double const tooth3_top = ball3_y + ball_radius_queen;

  cairo_save(cr);
  cairo_translate(cr, x + square_size_ / 2, y + square_size_ / 2);
  cairo_scale(cr, scale, scale);
  cairo_set_line_width(cr, black_line_width);

  // Draw the outline of the queen.
  // First fill, then stroke.
  for (int stroke = 0; stroke < 2; ++stroke)
  {
    // Draw the right side.
    cairo_move_to(cr, -tooth1_x, tooth1_top);
    cairo_line_to(cr, -crown_bottom_left, tooth_outside);
    cairo_line_to(cr, band1_radius, band1_y);
    // The call to arc() draws the last line part, to (inside_radius_queen, inside_y_queen).

    // Draw the half ellipse that makes out the bottom.
    cairo_save(cr);
    cairo_translate(cr, 0.0, inside_y_queen);
    cairo_scale(cr, 1.0, inside_scale_queen);
    cairo_arc(cr, 0.0, 0.0, inside_radius_queen, 0.0, M_PI);
    cairo_restore(cr);

    // Draw the left side.
    cairo_line_to(cr, -band1_radius, band1_y);
    cairo_line_to(cr, crown_bottom_left, tooth_outside);
    cairo_line_to(cr, tooth1_x, tooth1_top);

    // The lines of the teeth should not be 'connected' when we are stroking.
    if (stroke)
    {
      cairo_new_sub_path(cr);
      cairo_move_to(cr, tooth1_x, tooth1_top);
    }

    // Draw right-side of left-most tooth.
    cairo_line_to(cr, tooth2_x, tooth_inside);

    // Draw left-side of second tooth.
    cairo_line_to(cr, tooth2_x, tooth2_top);

    if (stroke)
    {
      cairo_new_sub_path(cr);
      cairo_move_to(cr, tooth2_x, tooth2_top);
    }

    // Draw right-side of second tooth.
    cairo_line_to(cr, tooth3_x, tooth_inside2);

    // Draw left-side of middle tooth.
    cairo_line_to(cr, 0.0, tooth3_top);

    if (stroke)
    {
      cairo_new_sub_path(cr);
      cairo_move_to(cr, 0.0, tooth3_top);
    }

    // Draw right-side of middle tooth.
    cairo_line_to(cr, -tooth3_x, tooth_inside2);

    // Draw left-side of fourth tooth.
    cairo_line_to(cr, -tooth2_x, tooth2_top);

    if (stroke)
    {
      cairo_new_sub_path(cr);
      cairo_move_to(cr, -tooth2_x, tooth2_top);
    }

    // Draw right-side of fourth tooth.
    cairo_line_to(cr, -tooth2_x, tooth_inside);

    // Draw left-side of last tooth.
    cairo_line_to(cr, -tooth1_x, tooth1_top);

    if (stroke)
    {
      if (color == chess::white)
        set_line_color(cr, color);
      else
        set_fill_color(cr, color);
      cairo_stroke(cr);
    }
    else
    {
      set_fill_color(cr, color);
      cairo_fill(cr);

      // Draw the upper side of the bottom ellipse.
      cairo_save(cr);
      cairo_translate(cr, 0.0, inside_y_queen);
      cairo_scale(cr, 1.0, inside_scale_queen);
      cairo_arc(cr, 0.0, 0.0, inside_radius_queen, -M_PI, 0.0);
      cairo_restore(cr);

      cairo_new_sub_path(cr);

      // Draw the half ellipse of band1.
      cairo_save(cr);
      cairo_translate(cr, 0.0, band1_y);
      cairo_scale(cr, 1.0, band1_scale);
      cairo_arc(cr, 0.0, 0.0, band1_radius, -M_PI, 0.0);
      cairo_restore(cr);

      set_line_color(cr, color);
      if (color == chess::white)
	cairo_stroke(cr);
      else
      {
	cairo_set_line_width(cr, white_line_width);
	cairo_stroke(cr);
	cairo_set_line_width(cr, black_line_width);
      }
    }
  }

  // Draw the five balls.

  cairo_arc(cr, ball1_x, ball1_y, ball_radius_queen, -M_PI, M_PI);

  if (color == chess::white)
    set_fill_color(cr, color);
  cairo_fill_preserve(cr);
  if (color == chess::white)
    set_line_color(cr, color);
  cairo_stroke(cr);

  cairo_arc(cr, ball2_x, ball2_y, ball_radius_queen, -M_PI, M_PI);

  if (color == chess::white)
    set_fill_color(cr, color);
  cairo_fill_preserve(cr);
  if (color == chess::white)
    set_line_color(cr, color);
  cairo_stroke(cr);

  cairo_arc(cr, 0.0, ball3_y, ball_radius_queen, -M_PI, M_PI);

  if (color == chess::white)
    set_fill_color(cr, color);
  cairo_fill_preserve(cr);
  if (color == chess::white)
    set_line_color(cr, color);
  cairo_stroke(cr);

  cairo_arc(cr, -ball2_x, ball2_y, ball_radius_queen, -M_PI, M_PI);

  if (color == chess::white)
    set_fill_color(cr, color);
  cairo_fill_preserve(cr);
  if (color == chess::white)
    set_line_color(cr, color);
  cairo_stroke(cr);

  cairo_arc(cr, -ball1_x, ball1_y, ball_radius_queen, -M_PI, M_PI);

  if (color == chess::white)
    set_fill_color(cr, color);
  cairo_fill_preserve(cr);
  if (color == chess::white)
    set_line_color(cr, color);
  cairo_stroke(cr);

  if (color == chess::white)
  {
    // Draw the splines at the bottom of the teeth.
    // The top left y-coordinate.
    static double const y0_queen = 0.0952;
    // The y-coordinate of the middle point.
    static double const ym = 0.0331;
    // The top left y-coordinate lays on the left side of the first tooth.
    // Calculate the x-coordinate:
    static double const x0_queen = tooth1_x + (y0_queen - tooth1_top) * (crown_bottom_left - tooth1_x) / (tooth_outside - tooth1_top);
    // The (apparent) tilting angle.
    static double const tilt_angle = atan((ym - y0_queen) / x0_queen);

    // The angle that the control lines make with the y-axis, before
    // mapping them onto a cylinder and before tilting the cylinder.
    static double const beta_queen = 1.202;
    // The length of the control lines.
    static double const len = 0.1728;
    // The y-value of the control points before tilting (relative to y0_queen).
    static double const py = len * cos(beta_queen);
    static double const y0_plus_py_cos_tilt_angle = y0_queen + py * cos(tilt_angle);
    static double const sin_tilt_angle = sin(tilt_angle);
    // The x-offset of the control points (this is an angle).
    static double px_offset = len * sin(beta_queen);

    cairo_move_to(cr, crown_bottom_left, tooth_outside);
    cairo_line_to(cr, x0_queen, y0_queen);

    // We have to draw N splines.
    int const N = 4;
    for (int i = 0; i < N; ++i)
    {
      double const alpha = i * M_PI / N;
      // The spline points before tilting.
      double px2 = x0_queen * cos(alpha + px_offset);
      double pz2 = - x0_queen * sin(alpha + px_offset);
      double px3 = x0_queen * cos(alpha + M_PI / N - px_offset);
      double pz3 = - x0_queen * sin(alpha + M_PI / N - px_offset);
      double px4 = x0_queen * cos(alpha + M_PI / N);
      double pz4 = - x0_queen * sin(alpha + M_PI / N);
      // Calculate the tilted values. This only has influence on the y value
      // (we rotate around the x-axis, and the resulting z-value doesn't interest us).
      double tpy2 = y0_plus_py_cos_tilt_angle - pz2 * sin_tilt_angle;
      double tpy3 = y0_plus_py_cos_tilt_angle - pz3 * sin_tilt_angle;
      double tpy4 = y0_queen - pz4 * sin_tilt_angle;
      cairo_curve_to(cr, px2, tpy2, px3, tpy3, px4, tpy4);
    }

    cairo_line_to(cr, -crown_bottom_left, tooth_outside);
  }

  // Draw the half ellipse of band2.
  cairo_save(cr);
  cairo_translate(cr, 0.0, band2_y);
  cairo_scale(cr, 1.0, band2_scale);
  cairo_arc_negative(cr, 0.0, 0.0, band2_radius, -0.15, -M_PI + 0.15);
  cairo_restore(cr);

  if (color == chess::white)
  {
    cairo_close_path(cr);
    set_fill_color(cr, color);
    cairo_fill_preserve(cr);
  }
  else
    cairo_set_line_width(cr, white_line_width);
  set_line_color(cr, color);
  cairo_stroke(cr);

  cairo_restore(cr);
}

void ChessPiece::draw_rook(cairo_t* cr, double x, double y, double scale, chess::EColor color)
{
  DoutEntering(dc::cairowindow, "ChessPiece::draw_rook(" << cr << ", " << x << ", " << y << ", " << scale << ", " << color << ")");

#ifdef CWDEBUG
  using namespace debugcairo;
#endif
  // Measurements, in cm.
  static double const width_rook_cm = 5.33;
  static double const foot_left_cm = 0.90;
  static double const base_left_cm = 1.26;
  static double const tower_left_cm = 1.64;
  static double const opening_left_cm = 1.795;
  static double const opening_right_cm = 2.315;
  static double const height_rook_cm = 5.30;
  static double const bottom_rook_cm = 0.58;
  static double const foot_top_cm = 0.95;
  static double const base_top_cm = 1.41;
  static double const tower_bottom_cm = 1.76;
  static double const tower_top_cm = 3.43;
  static double const top_bottom_cm = 3.81;
  static double const opening_bottom_cm = 4.25;
  // static double const top_top_cm = 4.61;

  // In coordinates.
  static double const foot_left = -0.5 + foot_left_cm / width_rook_cm + 0.5 * black_line_width;
  static double const base_left = -0.5 + base_left_cm / width_rook_cm + 0.5 * black_line_width;
  static double const tower_left = -0.5 + tower_left_cm / width_rook_cm + 0.5 * black_line_width;
  static double const opening_left = -0.5 + opening_left_cm / width_rook_cm + 0.5 * black_line_width;
  static double const opening_right = -0.5 + opening_right_cm / width_rook_cm + 0.5 * black_line_width;
  static double const bottom_rook = 0.5 - bottom_rook_cm / height_rook_cm - 0.5 * black_line_width;
  static double const foot_top = 0.5 - foot_top_cm / height_rook_cm - 0.5 * black_line_width;
  static double const base_top = 0.5 - base_top_cm / height_rook_cm - 0.5 * black_line_width;
  static double const tower_bottom = 0.5 - tower_bottom_cm / height_rook_cm - 0.5 * black_line_width;
  static double const tower_top = 0.5 - tower_top_cm / height_rook_cm - 0.5 * black_line_width;
  static double const top_bottom = 0.5 - top_bottom_cm / height_rook_cm - 0.5 * black_line_width;
  static double const opening_bottom = 0.5 - opening_bottom_cm / height_rook_cm - 0.5 * black_line_width;
  // static double const top_top = 0.5 - top_top_cm / height_rook_cm - 0.5 * black_line_width;
  // For alignment purposes, it's better to have the rook *precisely* centered.
  static double const top_top = -bottom_rook;

  // Snapped coordinates.
  double const inner_line_width = (color == chess::white) ? black_line_width : snap_line_width(white_line_width, scale);
  double const bottom_sn = snap_bottom(bottom_rook, y, scale, black_line_width);
  double const foot_top_sn = snap_bottom(foot_top, y, scale, inner_line_width);
  double const base_top_sn = snap_bottom(base_top, y, scale, inner_line_width);
  double const tower_bottom_sn = snap_bottom(tower_bottom, y, scale, inner_line_width);
  double const tower_top_sn = snap_top(tower_top, y, scale, inner_line_width);
  double const top_bottom_sn = snap_top(top_bottom, y, scale, inner_line_width);
  double const opening_bottom_sn = snap_top(opening_bottom, y, scale, black_line_width);
  double const top_top_sn = snap_top(top_top, y, scale, black_line_width);

  cairo_save(cr);
  cairo_translate(cr, x + square_size_ / 2, y + square_size_ / 2);
  cairo_scale(cr, scale, scale);
  cairo_set_line_width(cr, black_line_width);

  // Draw left side.
  cairo_move_to(cr, foot_left, bottom_sn);
  cairo_line_to(cr, foot_left, foot_top_sn);
  cairo_line_to(cr, base_left, foot_top_sn);
  cairo_line_to(cr, base_left, base_top_sn);
  cairo_line_to(cr, tower_left, tower_bottom_sn);
  cairo_line_to(cr, tower_left, tower_top_sn);
  cairo_line_to(cr, base_left, top_bottom_sn);
  cairo_line_to(cr, base_left, top_top_sn);

  // Draw top side.
  cairo_line_to(cr, opening_left, top_top_sn);
  cairo_line_to(cr, opening_left, opening_bottom_sn);
  cairo_line_to(cr, opening_right, opening_bottom_sn);
  cairo_line_to(cr, opening_right, top_top_sn);
  cairo_line_to(cr, -opening_right, top_top_sn);
  cairo_line_to(cr, -opening_right, opening_bottom_sn);
  cairo_line_to(cr, -opening_left, opening_bottom_sn);
  cairo_line_to(cr, -opening_left, top_top_sn);
  cairo_line_to(cr, -base_left, top_top_sn);

  // Draw right side.
  cairo_line_to(cr, -base_left, top_bottom_sn);
  cairo_line_to(cr, -tower_left, tower_top_sn);
  cairo_line_to(cr, -tower_left, tower_bottom_sn);
  cairo_line_to(cr, -base_left, base_top_sn);
  cairo_line_to(cr, -base_left, foot_top_sn);
  cairo_line_to(cr, -foot_left, foot_top_sn);
  cairo_line_to(cr, -foot_left, bottom_sn);

  // Draw bottom line.
  cairo_close_path(cr);
  cairo_path_t* path = cairo_copy_path(cr);

  set_fill_color(cr, color);
  cairo_fill(cr);

  // Draw inner horizontal lines.
  cairo_move_to(cr, base_left + 0.5 * black_line_width, foot_top_sn);
  cairo_line_to(cr, -base_left - 0.5 * black_line_width, foot_top_sn);
  cairo_new_sub_path(cr);
  cairo_move_to(cr, base_left, base_top_sn);
  cairo_line_to(cr, -base_left, base_top_sn);
  cairo_new_sub_path(cr);
  cairo_move_to(cr, tower_left + ((color == chess::white) ? 0.0 : (0.5 * black_line_width)), tower_bottom_sn);
  cairo_line_to(cr, -tower_left - ((color == chess::white) ? 0.0 : (0.5 * black_line_width)), tower_bottom_sn);
  cairo_new_sub_path(cr);
  cairo_move_to(cr, tower_left + ((color == chess::white) ? 0.0 : (0.5 * black_line_width)), tower_top_sn);
  cairo_line_to(cr, -tower_left - ((color == chess::white) ? 0.0 : (0.5 * black_line_width)), tower_top_sn);
  cairo_new_sub_path(cr);
  cairo_move_to(cr, base_left + black_line_width * 0.5, top_bottom_sn);
  cairo_line_to(cr, -base_left - black_line_width * 0.5, top_bottom_sn);

  set_line_color(cr, color);
  if (color == chess::white)
    cairo_stroke(cr);
  else
  {
    cairo_set_line_width(cr, inner_line_width);
    cairo_stroke(cr);
    cairo_set_line_width(cr, black_line_width);
  }

  cairo_append_path(cr, path);
  if (color == chess::black)
    set_fill_color(cr, color);
  cairo_stroke(cr);

  cairo_path_destroy(path);

  cairo_restore(cr);
}

void ChessPiece::draw_bishop(cairo_t* cr, double x, double y, double scale, chess::EColor color)
{
  DoutEntering(dc::cairowindow, "ChessPiece::draw_bishop(" << cr << ", " << x << ", " << y << ", " << scale << ", " << color << ")");

#ifdef CWDEBUG
  using namespace debugcairo;
#endif
  // Measurements, in cm.
  static double const width_bishop_cm = 5.34;
  static double const ribbon_width_cm = 0.49;
  static double const ribbon_bottom_left_cm = 0.72;
  static double const ribbon_top_left_cm = 2.28;
  static double const inside_outer_diameter_cm = 2.0;
  static double const circle_diameter_cm = 2.44;
  static double const cross_width_cm = 0.93;
  static double const ball_outer_diameter_cm = 0.81;
  static double const ball_inner_diameter_cm = 0.41;
  static double const circle_start_angle = 0.767;
  static double const ribbon_end_angle = 1.097;
  static double const height_bishop_cm = 5.44;
  static double const ribbon_bottom_y1_cm = 0.52;
  static double const ribbon_bottom_y2_cm = 0.76;
  static double const ribbon_bottom_y3_cm = 0.55;
  static double const ribbon_top_y1_cm = 0.99;
  static double const ribbon_top_y2_cm = 1.25;
  static double const ribbon_inside_y_cm = 0.93;
  static double const inside_bottom_cm = 1.34;
  static double const inside_top_cm = 1.86;
  static double const band_top_bishop_cm = 2.34;
  static double const circle_y_cm = 3.11;
  static double const cross_y_bishop_cm = 3.24;
  static double const point_y_cm = 4.47;
  static double const ball_y_cm = 4.675;
  static double const sp1_x_cm = 2.1;
  static double const sp1_y_cm = 3.95;
  static double const ribbon_bottom_x1_cm = 3.34;
  static double const ribbon_bottom_x2_cm = 4.1;
  static double const ribbon_top_x1_cm = 3.54;
  static double const ribbon_top_x2_cm = 4.24;

  // Translate to coordinates.
  static double const ribbon_width = ribbon_width_cm / height_bishop_cm;
  static double const ribbon_bottom_left = -0.5 + ribbon_bottom_left_cm / width_bishop_cm;
  static double const ribbon_bottom_x1 = -0.5 + ribbon_bottom_x1_cm / width_bishop_cm;
  static double const ribbon_bottom_x2 = -0.5 + ribbon_bottom_x2_cm / width_bishop_cm;
  static double const ribbon_top_x1 = -0.5 + ribbon_top_x1_cm / width_bishop_cm;
  static double const ribbon_top_x2 = -0.5 + ribbon_top_x2_cm / width_bishop_cm;
  static double const ribbon_top_left = -0.5 + ribbon_top_left_cm / width_bishop_cm;
  static double const inside_radius_bishop = 0.5 * (inside_outer_diameter_cm / width_bishop_cm - black_line_width);
  static double const circle_radius = 0.5 * circle_diameter_cm / width_bishop_cm;
  static double const cross_leg = 0.5 * cross_width_cm / width_bishop_cm;
  static double const ball_radius_bishop  = 0.25 * (ball_outer_diameter_cm + ball_inner_diameter_cm) / width_bishop_cm;
  static double const ball_line_width = black_line_width; // 0.5 * (ball_outer_diameter_cm - ball_inner_diameter_cm) / width_bishop_cm
  static double const ribbon_bottom_y1 = 0.5 - ribbon_bottom_y1_cm / height_bishop_cm - 0.5 * black_line_width;
  static double const ribbon_bottom_y2 = 0.5 - ribbon_bottom_y2_cm / height_bishop_cm + 0.5 * black_line_width;
  static double const ribbon_bottom_y3 = 0.5 - ribbon_bottom_y3_cm / height_bishop_cm;
  static double const ribbon_inside_y = 0.5 - ribbon_inside_y_cm / height_bishop_cm;
  static double const ribbon_top_y1 = 0.5 - ribbon_top_y1_cm / height_bishop_cm - 0.5 * black_line_width;
  static double const ribbon_top_y2 = 0.5 - ribbon_top_y2_cm / height_bishop_cm + 0.5 * black_line_width;
  static double const inside_scale_bishop = ((inside_top_cm - inside_bottom_cm) / height_bishop_cm - black_line_width) / (2 * inside_radius_bishop);
  static double const inside_y_bishop = 0.5 - 0.5 * (inside_top_cm + inside_bottom_cm) / height_bishop_cm;
  static double const inside_bottom = 0.5 - inside_bottom_cm / height_bishop_cm - 0.5 * black_line_width;
  static double const band_top_bishop = 0.5 - band_top_bishop_cm / height_bishop_cm + 0.5 * black_line_width;
  static double const circle_y = 0.5 - circle_y_cm / height_bishop_cm;
  static double const cross_y_bishop = 0.5 - cross_y_bishop_cm / height_bishop_cm;
  static double const point_y = 0.5 - point_y_cm / height_bishop_cm;
  static double const ball_y = 0.5 - ball_y_cm / height_bishop_cm;
  static double const inside_angle = acos(-ribbon_top_left / inside_radius_bishop);
  static double const sp1_x = -0.5 + sp1_x_cm / width_bishop_cm;
  static double const sp1_y = 0.5 - sp1_y_cm / height_bishop_cm;

  // Precalculations for the ribbon.
  static double const spline_magic = 0.551784;
  static double const cp2_x = ribbon_bottom_y1 - ribbon_inside_y;
  static double const sp2_x = spline_magic * cp2_x;
  static double const sp2_y = ribbon_inside_y + spline_magic * (ribbon_bottom_y1 - ribbon_inside_y);
  static double const sp3_x = ribbon_bottom_x1 - spline_magic * (ribbon_bottom_x1 - cp2_x);
  static double const sp3_y = ribbon_bottom_y1;
  static double const sp4_x = ribbon_bottom_x1 + spline_magic * (ribbon_bottom_x2 - ribbon_bottom_x1);
  static double const sp4_y = ribbon_bottom_y1;
  static double const sp5_x = ribbon_bottom_x2 - spline_magic * (ribbon_bottom_x2 - ribbon_bottom_x1);
  static double const sp5_y = ribbon_bottom_y2;
  static double const cp6_x = -ribbon_bottom_left - (ribbon_bottom_y3 - ribbon_bottom_y2) * tan(ribbon_end_angle);
  static double const sp6_x = ribbon_bottom_x2 + spline_magic * (cp6_x - ribbon_bottom_x2);
  static double const sp6_y = ribbon_bottom_y2;
  static double const sp7_x = -ribbon_bottom_left - spline_magic * (-ribbon_bottom_left - cp6_x);
  static double const sp7_y = ribbon_bottom_y3 - spline_magic * (ribbon_bottom_y3 - ribbon_bottom_y2);
  static double const ribbon_end_top_x = -ribbon_bottom_left + ribbon_width * cos(ribbon_end_angle);
  static double const ribbon_end_top_y = ribbon_bottom_y3 - ribbon_width * sin(ribbon_end_angle);
  static double const cp8_x = ribbon_end_top_x - (ribbon_end_top_y - ribbon_top_y2) * tan(ribbon_end_angle);
  static double const sp8_x = ribbon_end_top_x - spline_magic * (ribbon_end_top_x - cp8_x);
  static double const sp8_y = ribbon_end_top_y - spline_magic * (ribbon_end_top_y - ribbon_top_y2);
  static double const sp9_x = ribbon_top_x2 + spline_magic * (cp8_x - ribbon_top_x2);
  static double const sp9_y = ribbon_top_y2;
  static double const sp10_x = ribbon_top_x2 - spline_magic * (ribbon_top_x2 - ribbon_top_x1);
  static double const sp10_y = ribbon_top_y2;
  static double const sp11_x = ribbon_top_x1 + spline_magic * (ribbon_top_x2 - ribbon_top_x1);
  static double const sp11_y = ribbon_top_y1;
  static double const ribbon_top_y3 = 0.2695;
  static double const sp12_x = ribbon_top_x1 - spline_magic * (ribbon_top_x1 + ribbon_top_left);
  static double const sp12_y = ribbon_top_y1;
  static double const sp13_x = -ribbon_top_left;
  static double const sp13_y = ribbon_top_y3 + 0.509 * spline_magic * (ribbon_top_y1 - ribbon_top_y3);

  cairo_save(cr);
  cairo_translate(cr, x + square_size_ / 2, y + square_size_ / 2);
  cairo_scale(cr, scale, scale);
  cairo_set_line_width(cr, black_line_width);

  // Draw the ribbon.

  // Part 14.
  cairo_move_to(cr, -ribbon_top_x1, ribbon_top_y1);
  cairo_curve_to(cr, -sp11_x, sp11_y, -sp10_x, sp10_y, -ribbon_top_x2, ribbon_top_y2);

  // Part 13.
  cairo_curve_to(cr, -sp9_x, sp9_y, -sp8_x, sp8_y, -ribbon_end_top_x, ribbon_end_top_y);

  // Part 12.
  cairo_line_to(cr, ribbon_bottom_left, ribbon_bottom_y3);

  // Part 11.
  cairo_curve_to(cr, -sp7_x, sp7_y, -sp6_x, sp6_y, -ribbon_bottom_x2, ribbon_bottom_y2);

  // Part 10.
  cairo_curve_to(cr, -sp5_x, sp5_y, -sp4_x, sp4_y, -ribbon_bottom_x1, ribbon_bottom_y1);

  // Part 9.
  cairo_curve_to(cr, -sp3_x, sp3_y, -sp2_x, sp2_y, 0.0, ribbon_inside_y);

  // Part 1.
  cairo_curve_to(cr, sp2_x, sp2_y, sp3_x, sp3_y, ribbon_bottom_x1, ribbon_bottom_y1);

  // Part 2.
  cairo_curve_to(cr, sp4_x, sp4_y, sp5_x, sp5_y, ribbon_bottom_x2, ribbon_bottom_y2);

  // Part 3.
  cairo_curve_to(cr, sp6_x, sp6_y, sp7_x, sp7_y, -ribbon_bottom_left, ribbon_bottom_y3);

  // Part 4.
  cairo_line_to(cr, ribbon_end_top_x, ribbon_end_top_y);

  // Part 5.
  cairo_curve_to(cr, sp8_x, sp8_y, sp9_x, sp9_y, ribbon_top_x2, ribbon_top_y2);

  // Part 6.
  cairo_curve_to(cr, sp10_x, sp10_y, sp11_x, sp11_y, ribbon_top_x1, ribbon_top_y1);

  if (color == chess::black)
  {
    set_fill_color(cr, color);
    cairo_fill_preserve(cr);
    cairo_stroke(cr);
    cairo_move_to(cr, ribbon_top_x1, ribbon_top_y1);
  }

  // Part 7.
  cairo_curve_to(cr, sp12_x, sp12_y, sp13_x, sp13_y, -ribbon_top_left, ribbon_top_y3);

  // Part 8 and 17.
  cairo_save(cr);
  cairo_translate(cr, 0.0, inside_y_bishop);
  cairo_scale(cr, 1.0, inside_scale_bishop);
  cairo_arc(cr, 0.0, 0.0, inside_radius_bishop, inside_angle, M_PI - inside_angle);
  cairo_restore(cr);

  // Part 16.
  cairo_line_to(cr, ribbon_top_left, ribbon_top_y3);

  // Part 15.
  cairo_curve_to(cr, -sp13_x, sp13_y, -sp12_x, sp12_y, -ribbon_top_x1 + 0.01 * black_line_width, ribbon_top_y1);
  cairo_close_path(cr);

  if (color == chess::white)
    set_fill_color(cr, color);
  else
    set_line_color(cr, color);
  cairo_fill_preserve(cr);
  if (color == chess::white)
    set_line_color(cr, color);
  else
    set_fill_color(cr, color);
  cairo_stroke(cr);

  // Draw vertical line between left and right ribbon.
  cairo_move_to(cr, 0.0, inside_bottom);
  cairo_line_to(cr, 0.0, ribbon_inside_y);
  cairo_stroke(cr);

  // Draw the outline of the bishop.

  cairo_save(cr);
  cairo_translate(cr, 0.0, inside_y_bishop);
  cairo_scale(cr, 1.0, inside_scale_bishop);
  cairo_arc(cr, 0.0, 0.0, inside_radius_bishop, 0.0, -M_PI);
  cairo_restore(cr);

  cairo_arc(cr, 0.0, circle_y, circle_radius, -M_PI - circle_start_angle, -M_PI);

  cairo_curve_to(cr, -circle_radius, circle_y - 0.0848, sp1_x - 0.02657, sp1_y + 0.01722, sp1_x, sp1_y);
  cairo_curve_to(cr, sp1_x + 0.08845, sp1_y - 0.05733, -0.000333, point_y + 0.000265, 0.0, point_y);
  cairo_curve_to(cr, 0.000333, point_y + 0.000265, -sp1_x - 0.08845, sp1_y - 0.05733, -sp1_x, sp1_y);
  cairo_curve_to(cr, -sp1_x + 0.02657, sp1_y + 0.01722, circle_radius, circle_y - 0.0848, circle_radius, circle_y);

  cairo_arc(cr, 0.0, circle_y, circle_radius, 0.0, circle_start_angle);

  cairo_close_path(cr);

  if (color == chess::white)
    set_fill_color(cr, color);
  cairo_fill_preserve(cr);
  if (color == chess::white)
    set_line_color(cr, color);
  cairo_stroke(cr);

  // Draw inside lines.
  if (color == chess::black)
    set_line_color(cr, color);
  cairo_save(cr);
  if (color == chess::black)
  {
    static double const x2_bishop = -circle_radius * cos(circle_start_angle);
    static double const y2_bishop = (circle_y + circle_radius * sin(circle_start_angle));
    cairo_move_to(cr, -inside_radius_bishop, inside_y_bishop);
    cairo_line_to(cr, x2_bishop, y2_bishop);
    cairo_line_to(cr, -x2_bishop, y2_bishop);
    cairo_line_to(cr, inside_radius_bishop, inside_y_bishop);
    cairo_close_path(cr);
    cairo_clip(cr);
  }
  cairo_save(cr);
  cairo_translate(cr, 0.0, inside_y_bishop);
  cairo_scale(cr, 1.0, inside_scale_bishop);
  cairo_arc(cr, 0.0, 0.0, inside_radius_bishop, -M_PI, 0.0);
  cairo_restore(cr);
  if (color == chess::white)
    cairo_stroke(cr);
  else
  {
    cairo_set_line_width(cr, white_line_width);
    cairo_stroke(cr);
    cairo_set_line_width(cr, black_line_width);
  }
  cairo_restore(cr);

  // Ok, some real math needed here.
  // First, we scale the y-axis - so that the problem changes to a circle instead of an ellipse.
  // Also, flip the y-axis, so that the top of the screen becomes the positive y-axis.
  // This means that instead of using the normal 'y' values, we now should use y / -inside_scale_bishop.
  // Next, calculate the line through (-inside_radius_bishop, inside_y_bishop) and the
  // start of the circle with vector parametric formula: V = X0 + U t1,
  // where U is the unit vector (u1, u2), t1 is the parameter and X0 = (x0, 0), the
  // point where the line crosses the x-axis.
  static double const x1 = -inside_radius_bishop;
  static double const y1 = inside_y_bishop / -inside_scale_bishop;
  static double const x2 = -circle_radius * cos(circle_start_angle);
  static double const y2 = (circle_y + circle_radius * sin(circle_start_angle)) / -inside_scale_bishop;
  static double const d = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
  static double const u1 = (x2 - x1) / d;
  static double const u2 = (y2 - y1) / d;
  static double const x0 = x1 + (x2 - x1) * (0 - y1) / (y2 - y1);
  // A line through the center of the circle, perpendicular to this line is:
  // V = Y0 + W t2, where Y0 = (0, y0) and W = (u2, -u1).
  // Those two lines cross therefore where,
  // x: x0 + u1 t1 = u2 t2
  // y: u2 t1 = y0 - u1 t2
  // Since the t2 parameter gives the distance from the center in terms of a unit vector,
  // the distance from the center to the crossing point is |t2|, which we find by cancelling t1:
  // x0 u2 + u1 (y0 - u1 t2) = u2^2 t2 --> t2 = x0 u2 + u1 y0, since u1^2 + u2^2 = 1.
  // In the crossing point, this value is negative (u2 is positive, and the cross point
  // lays left of the center). Hence we find that the distance is -t2 =
  // - x0 u2 - u1 y0 = band_top_bishop / -inside_scale_bishop - y0.
  static double const y0 = (band_top_bishop / -inside_scale_bishop + x0 * u2) / (1 - u1);
  static double const band_radius = band_top_bishop / -inside_scale_bishop - y0;
  static double const angle = atan(u1 / u2);
  cairo_save(cr);
  cairo_scale(cr, 1.0, -inside_scale_bishop);
  if (color == chess::black)
  {
    static double const t2 = x0 * u2 + u1 * y0;
    static double const t1 = (y0 - u1 * t2) / u2;
    static double const x = x0 + u1 * t1;
    cairo_move_to(cr, x, y0);
    cairo_line_to(cr, x + d * u1, y0 + d * u2);
    cairo_line_to(cr, -x - d * u1, y0 + d * u2);
    cairo_line_to(cr, -x, y0);
    cairo_close_path(cr);
    cairo_clip(cr);
  }
  cairo_arc(cr, 0.0, y0, band_radius, angle, M_PI - angle);
  // Reverse the scale before stroking, without restoring the clipping area.
  cairo_scale(cr, 1.0, -1.0 / inside_scale_bishop);
  if (color == chess::white)
    cairo_stroke(cr);
  else
  {
    cairo_set_line_width(cr, white_line_width);
    cairo_stroke(cr);
    cairo_set_line_width(cr, black_line_width);
  }
  cairo_restore(cr);

  // Draw the cross.
  cairo_move_to(cr, -cross_leg, cross_y_bishop);
  cairo_line_to(cr, cross_leg, cross_y_bishop);
  cairo_move_to(cr, 0.0, cross_y_bishop - cross_leg);
  cairo_line_to(cr, 0.0, cross_y_bishop + cross_leg);
  if (color == chess::white)
    cairo_stroke(cr);
  else
  {
    cairo_set_line_width(cr, white_line_width);
    cairo_stroke(cr);
    cairo_set_line_width(cr, black_line_width);
  }

  if (color == chess::black)
  {
    cairo_move_to(cr, -inside_radius_bishop, inside_y_bishop);
    cairo_arc(cr, 0.0, circle_y, circle_radius, -M_PI - circle_start_angle, -M_PI);
    cairo_move_to(cr, inside_radius_bishop, inside_y_bishop);
    cairo_arc_negative(cr, 0.0, circle_y, circle_radius, circle_start_angle, 0.0);
    set_fill_color(cr, color);
    cairo_stroke(cr);
  }

  // Draw the little ball on the top.
  cairo_set_line_width(cr, ball_line_width);
  cairo_arc(cr, 0.0, ball_y, ball_radius_bishop, -M_PI, M_PI);
  if (color == chess::white)
    set_fill_color(cr, color);
  cairo_fill_preserve(cr);
  if (color == chess::white)
    set_line_color(cr, color);
  cairo_stroke(cr);

  cairo_restore(cr);
}

void ChessPiece::draw_knight(cairo_t* cr, double x, double y, double scale, chess::EColor color)
{
  DoutEntering(dc::cairowindow, "ChessPiece::draw_knight(" << cr << ", " << x << ", " << y << ", " << scale << ", " << color << ")");

#ifdef CWDEBUG
  using namespace debugcairo;
#endif
  // Measurements.
  static double const height_knight_cm = 21.9;
  static double const pixels_per_cm = 32.467;
  static double const min_nose_x_px = 8.0;
  static double const right_ear_y_px = 15.0;		// See 'Draw right ear'.
  static double const bottom_right_x_px = 582.82;	// See 'back' curve.
  static double const bottom_right_y_px = 580.82;
  static double const bottom_left_x_px = 190.00;	// See 'front' curve.
  // Derived.
  static double const pixel_scale = 1.0 / (pixels_per_cm * height_knight_cm);
  static double const knight_black_line_width = 0.95 * black_line_width / pixel_scale;
  static double const knight_white_line_width = 1.3 * knight_black_line_width;
  static double const knight_white_glyp_line_width = knight_white_line_width - knight_black_line_width;

  // The outline of the knight in coordinates, without translation.
  static double const max_y = bottom_right_y_px * pixel_scale;
  static double const min_y = right_ear_y_px * pixel_scale;
  static double const max_x = bottom_right_x_px * pixel_scale;
  static double const min_x = min_nose_x_px * pixel_scale;

  // Calculate coordinate offsets, needed to center the knight.
  static double const pixel_translate_x = -(max_x + min_x) / 2;
  static double const pixel_translate_y = -(max_y + min_y) / 2;

  cairo_save(cr);
  cairo_translate(cr, x + square_size_ / 2, y + square_size_ / 2);
  cairo_scale(cr, scale, scale);

  // At this point, the coordinates run from -0.5 till 0.5.
  // However, we use pixels as coordinates (as measured in gimp)
  // Translate the image so that the pixel-coordinate center falls on (0, 0).
  cairo_translate(cr, pixel_translate_x, pixel_translate_y);
  // Scale, so we can use "pixel-coordinates".
  cairo_scale(cr, pixel_scale, pixel_scale);

  // Fill body.
  cairo_move_to(cr, 319.00, 299.00);
  cairo_curve_to(cr, 322.00, 449.00, 165.00, 445.00, 192.00, 570.00);
  cairo_curve_to(cr, 192.00, 570.00, 568.50, 571.00, 568.50, 571.00);
  cairo_curve_to(cr, 577.00, 426.00, 533.00, 99.00, 340.50, 88.50);
  cairo_curve_to(cr, 245.50, 87.50, 206.00, 86.00, 195.00, 102.00);
  set_fill_color(cr, color);
  cairo_fill(cr);

  // Draw shadow.
  cairo_move_to(cr, 315.00, 300.33);
  cairo_curve_to(cr, 301.43, 300.80, 291.75, 314.52, 282.00, 325.00);
  cairo_curve_to(cr, 298.67, 317.33, 316.33, 325.00, 317.33, 344.33);
  cairo_curve_to(cr, 321.33, 337.33, 326.00, 326.00, 315.00, 300.33);
  if (color == chess::white)
    set_line_color(cr, color);
  cairo_fill(cr);

  // Draw back.
  // Change the thickness of the top of the back to a line width:
  static double const back_top_offset = (93.00 - knight_black_line_width) - 82.00;
  cairo_move_to(cr, 582.82, 580.82);
  cairo_curve_to(cr, 589.00, 359.00, 530.00,85.00, 332.00, 82.00 + back_top_offset);
  cairo_curve_to(cr, 320.87, 82.04 + back_top_offset, 314.25, 82.12 + back_top_offset, 302.50, 82.38 + back_top_offset);
  cairo_curve_to(cr, 302.75, 95.38, 296.22, 93.73, 319.50, 94.00);
  cairo_curve_to(cr, 510.50, 93.00, 556.12, 359.00, 556.12, 563.00);
  cairo_close_path(cr);
  cairo_fill(cr);

  // Draw front.
  cairo_move_to(cr, 190.00, 570.00);
  cairo_curve_to(cr, 190.00, 550.75, 190.00, 549.00, 190.00, 540.00);
  cairo_curve_to(cr, 190.00, 493.25, 210.50, 482.50, 285.00, 409.50);
  cairo_curve_to(cr, 298.25, 391.75, 313.00, 357.50, 317.75, 344.75);
  cairo_curve_to(cr, 320.25, 340.00, 320.25, 330.00, 320.00, 280.00);
  cairo_set_line_width(cr, knight_black_line_width);
  cairo_stroke(cr);

  // Draw head.
  cairo_move_to(cr, 144.00, 31.50);
  cairo_curve_to(cr, 122.50, 67.00, 147.50, 57.50, 146.00, 105.00);
  cairo_curve_to(cr, 112.00, 125.50, 123.00, 140.50, 102.50, 170.00);
  cairo_curve_to(cr, 84.00, 199.50, 128.00, 181.50, 33.50, 313.50);
  cairo_curve_to(cr, -23.00, 414.00, 81.50, 468.00, 130.00, 447.50);
  cairo_curve_to(cr, 182.50, 398.00, 142.50, 427.00, 179.50, 390.00);
  cairo_curve_to(cr, 194.50, 376.50, 212.50, 349.50, 237.50, 347.00);
  cairo_curve_to(cr, 268.00, 344.00, 283.50, 323.50, 306.00, 301.00);
  cairo_curve_to(cr, 327.50, 276.50, 330.00, 264.50, 330.00, 228.50);
  if (color == chess::white)
    set_fill_color(cr, color);
  cairo_fill_preserve(cr);
  cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND);
  if (color == chess::white)
    set_line_color(cr, color);
  cairo_stroke(cr);
  cairo_move_to(cr, 201.00, 94.50);
  cairo_curve_to(cr, 184.50, 54.50, 152.00, 43.50, 144.00, 31.50);
  cairo_stroke(cr);

  // Draw between ears.
  cairo_move_to(cr, 170.50, 136.50);
  cairo_curve_to(cr, 170.00, 129.50, 175.50, 125.00, 183.50, 116.00);
  cairo_curve_to(cr, 204.50, 91.00, 216.00, 94.00, 238.00, 91.00);
  cairo_stroke(cr);

  if (color == chess::black)
  {
    // Draw white hair.
    cairo_move_to(cr, 529.00, 570.00);
    cairo_curve_to(cr, 530.50, 352.00, 476.50, 128.50, 334.00, 121.00);
    cairo_curve_to(cr, 310.50, 118.50, 310.00, 117.50, 296.50, 117.50);
    cairo_curve_to(cr, 291.50, 100.00, 252.50, 95.50, 242.20, 119.35);
    cairo_curve_to(cr, 227.55, 120.95, 212.22, 124.23, 198.50, 130.50);
    cairo_curve_to(cr, 178.00, 137.50, 158.50, 147.50, 154.00, 137.00);
    cairo_curve_to(cr, 149.50, 127.00, 145.50, 121.00, 204.00, 100.00);
    cairo_curve_to(cr, 226.50, 90.00, 276.50, 92.00, 319.50, 94.00);
    cairo_curve_to(cr, 510.50, 93.00, 556.00, 354.00, 556.00, 570.00);
    cairo_curve_to(cr, 548.06, 571.00, 537.73, 569.45, 529.00, 570.00);
    set_line_color(cr, color);
    cairo_fill(cr);
  }

  // Draw bottom.
  double dummy = bottom_left_x_px;
  double bottom_right_y_px_sn = bottom_right_y_px;
  if (scale >= 27)
  {
    // Snap bottom to grid.
    cairo_user_to_device(cr, &dummy, &bottom_right_y_px_sn);
    bottom_right_y_px_sn = round(bottom_right_y_px_sn);
    cairo_device_to_user(cr, &dummy, &bottom_right_y_px_sn);
  }
  cairo_rectangle(cr, bottom_left_x_px - 0.5 * knight_black_line_width,
                      bottom_right_y_px_sn - knight_black_line_width,
                      bottom_right_x_px - (bottom_left_x_px - 0.5 * knight_black_line_width),
		      knight_black_line_width);
  if (color == chess::black)
    set_fill_color(cr, color);
  cairo_fill(cr);

  // Draw mouth.
  cairo_move_to(cr, 113.67, 389.33);
  cairo_curve_to(cr, 121.00, 388.00, 129.33, 406.67, 120.67, 414.67);
  cairo_curve_to(cr, 114.33, 419.33, 104.33, 431.00, 112.67, 444.00);
  //cairo_curve_to(cr, 133.53, 466.89, 104.22, 478.73, 93.00, 446.67);
  cairo_line_to(cr, 93.00, 446.67);
  cairo_curve_to(cr, 89.00, 418.67, 94.67, 417.33, 100.00, 412.67);
  cairo_curve_to(cr, 112.67, 402.00, 100.67, 394.67, 113.67, 389.33);
  if (color == chess::black)
    set_line_color(cr, color);
  if (color == chess::white)
    cairo_fill(cr);
  else
  {
    cairo_fill_preserve(cr);
    cairo_set_line_width(cr, knight_white_glyp_line_width);
    cairo_stroke(cr);
    cairo_set_line_width(cr, knight_black_line_width);
  }

  // Redraw a part of the head.
  cairo_move_to(cr, 33.50, 313.50);
  cairo_curve_to(cr, -23.00, 414.00, 81.50, 468.00, 130.00, 447.50);
  if (color == chess::black)
    set_fill_color(cr, color);
  cairo_stroke(cr);

  if (color == chess::black)
  {
    // Draw jaw.
    cairo_move_to(cr, 312.32, 293.46);
    cairo_curve_to(cr, 328.01, 273.63, 330.00, 260.62, 330.00, 228.50);
    cairo_set_line_width(cr, knight_white_line_width);
    set_line_color(cr, color);
    cairo_stroke(cr);
    cairo_set_line_width(cr, knight_black_line_width);
  }

  // Draw right ear.
  for (int stroke = 0; stroke < 2; ++stroke)
  {
    cairo_move_to(cr, 242.00, 114.00);
    cairo_curve_to(cr, 235.00, 76.00, 235.50, 92.50, 267.00, 15.00);	// 15.00 corresponds to min_y
    if (stroke)
      cairo_move_to(cr, 267.00, 15.00);
    cairo_curve_to(cr, 309.50, 85.50, 312.00, 88.00, 295.00, 117.00);
    if (stroke)
    {
      if (color == chess::white)
        set_line_color(cr, color);
      cairo_stroke(cr);
    }
    else
    {
      if (color == chess::white)
	set_fill_color(cr, color);
      else
        set_fill_color(cr, color);
      cairo_fill(cr);
    }
  }

  if (color == chess::black)
    set_line_color(cr, color);

  // Draw nose.
  cairo_move_to(cr, 76.00, 363.00);
  cairo_curve_to(cr, 66.00, 372.33, 78.33, 379.00, 66.00, 384.00);
  cairo_curve_to(cr, 21.00, 399.00, 61.67, 331.00, 79.67, 341.67);
  cairo_curve_to(cr, 81.00, 342.00, 84.67, 353.33, 76.00, 363.00);
  if (color == chess::white)
    cairo_fill(cr);
  else
  {
    cairo_fill_preserve(cr);
    cairo_set_line_width(cr, knight_white_glyp_line_width);
    cairo_stroke(cr);
    cairo_set_line_width(cr, knight_black_line_width);
  }

  // Draw eye.
  cairo_move_to(cr, 173.33, 208.00);
  cairo_curve_to(cr, 180.67, 207.00, 182.00, 197.67, 182.00, 197.67);
  cairo_curve_to(cr, 184.59, 176.98, 182.28, 177.30, 190.67, 173.00);
  cairo_curve_to(cr, 201.00, 169.33, 198.33, 146.00, 173.33, 161.67);
  cairo_curve_to(cr, 146.00, 181.33, 130.67, 192.00, 128.33, 202.67);
  cairo_curve_to(cr, 124.00, 233.33, 131.00, 227.33, 144.67, 207.00);
  cairo_curve_to(cr, 150.67, 201.00, 158.67, 193.67, 162.33, 203.33);
  cairo_curve_to(cr, 164.67, 206.00, 165.63, 209.29, 173.33, 208.00);
  if (color == chess::white)
    cairo_fill(cr);
  else
  {
    cairo_fill_preserve(cr);
    cairo_set_line_width(cr, knight_white_glyp_line_width);
    cairo_stroke(cr);
  }

  cairo_restore(cr);
}

StrokeExtents ChessPiece::do_draw(cairo_t* cr)
{
  DoutEntering(dc::cairowindow, "draw::ChessPiece::do_draw(cr) [" << this << "]");
  using namespace chess;
  switch (piece_)
  {
    case pawn:
      draw_pawn(cr, x1_, y1_, square_size_, color_);
      break;
    case rook:
      draw_rook(cr, x1_, y1_, square_size_, color_);
      break;
    case knight:
      draw_knight(cr, x1_, y1_, square_size_, color_);
      break;
    case bishop:
      draw_bishop(cr, x1_, y1_, square_size_, color_);
      break;
    case queen:
      draw_queen(cr, x1_, y1_, square_size_, color_);
      break;
    case king:
      draw_king(cr, x1_, y1_, square_size_, color_);
      break;
    default:
      ASSERT(false);
      AI_NEVER_REACHED
  }
  return {x1_ + 1, y1_ + 1, x1_ + square_size_ - 1, y1_ + square_size_ - 1};
}

} // namespace cairowindow::draw
