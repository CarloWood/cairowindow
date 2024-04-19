#include "sys.h"
#include "cairowindow/Window.h"
#include "cairowindow/Plot.h"
#include "cairowindow/BezierFitter.h"
#include "utils/AIAlert.h"
#include "utils/debug_ostream_operators.h"
#include "debug.h"

int main()
{
  Debug(NAMESPACE_DEBUG::init());
  Dout(dc::notice, "Entering main()");

  try
  {
    using namespace cairowindow;
    using Window = cairowindow::Window;

    // Create a window.
    Window window("Arc length", 1200, 900);

    // Create a new layer with a gray background.
    auto background_layer = window.create_background_layer<Layer>(color::white COMMA_DEBUG_ONLY("background_layer"));

    // Create another layer.
    auto second_layer = window.create_layer<Layer>({} COMMA_DEBUG_ONLY("second_layer"));

    // Open the window and start drawing.
    std::thread event_loop([&](){
      // Open window, handle event loop. This must be constructed after the draw stuff, so that it is destructed first!
      // Upon destruction it blocks until the event loop thread finished (aka, the window was closed).
      EventLoop event_loop = window.run();
      event_loop.set_cleanly_terminated();
    });

    // Create and draw plot area.
    plot::Plot plot(window.geometry(), { .grid = {.color = color::orange} },
        "Arc length cubic Bezier", {},
        "x", {},
        "y", {});
    plot.set_xrange({0, 1});
    plot.set_yrange({0, 2});
    plot.add_to(background_layer, true);

    draw::PointStyle point_style({.color_index = 0, .filled_shape = 1});
    draw::TextStyle label_style({.position = draw::centered_above, .font_size = 18.0, .offset = 10});
    draw::LineStyle curve_line_style({.line_width = 1.0});
    draw::LineStyle line_style({.line_color = color::red, .line_width = 1.0});
    draw::LineStyle arrow_style({.line_color = color::green, .line_width = 1.0});
    draw::BezierCurveStyle bezier_curve_style({.line_color = Color::next_color(), .line_width = 1.0});
    draw::TextStyle slider_style({.position = draw::centered_below, .font_size = 18.0, .offset = 10});

//    auto slider_offset = plot.create_slider(second_layer, {978, 83, 7, 400}, -0.205, -M_PI, M_PI);
//    auto slider_offset_label = plot.create_text(second_layer, slider_style, Pixel{978, 483}, "offset");

    auto plot_P0 = plot.create_point(second_layer, point_style, {0.0, 0.0});
    auto plot_P1 = plot.create_point(second_layer, point_style, {1.0, 0.0});
    auto plot_C0 = plot.create_point(second_layer, point_style, {0.2, 0.3});
    auto plot_C1 = plot.create_point(second_layer, point_style, {0.6, 0.2});

    auto P0_label = plot.create_text(second_layer, label_style, plot_P0, "P₀");
    auto P1_label = plot.create_text(second_layer, label_style, plot_P1, "P₁");

    window.register_draggable(plot, &plot_C0);
    window.register_draggable(plot, &plot_C1);

    while (true)
    {
      // Suppress immediate updating of the window for each created item, in order to avoid flickering.
      window.set_send_expose_events(false);

      auto P0C0_arrow = plot.create_connector(second_layer, arrow_style, plot_P0, plot_C0);
      auto P1C1_arrow = plot.create_connector(second_layer, arrow_style, plot_P1, plot_C1);

      auto C0_label = plot.create_text(second_layer, label_style, plot_C0, "C₀");
      auto C1_label = plot.create_text(second_layer, label_style, plot_C1, "C₁");

      Dout(dc::notice, "C0 = " << plot_C0.x() << ", " << plot_C0.y());
      Dout(dc::notice, "C1 = " << plot_C1.x() << ", " << plot_C1.y());

      auto plot_bezier = plot.create_bezier_curve(second_layer, line_style, plot_P0, plot_C0, plot_C1, plot_P1);

#if 0
      double len0 = 0.0;
      double t = 0.0;
      Point prev_P;
      int n = 100000;
      double const dt = 1.0 / n;
      for (int i = 0; i <= n; ++i)
      {
        t = i * dt;
        auto P = plot_bezier.P(t);
        if (i > 0)
          len0 += LinePiece{prev_P, P}.length();
        prev_P = P;
      }
      Dout(dc::notice, "len0 = " << std::setprecision(std::numeric_limits<double>::max_digits10) << len0);
#endif

      auto m = plot_bezier.M();
      double m00 = m.coefficient[0].x();
      double m10 = m.coefficient[0].y();
      double m01 = m.coefficient[1].x();
      double m11 = m.coefficient[1].y();
      double m02 = m.coefficient[2].x();
      double m12 = m.coefficient[2].y();
      double m03 = m.coefficient[3].x();
      double m13 = m.coefficient[3].y();
      double c0 = utils::square(m01) + utils::square(m11);
      double c1 = 4 * m01 * m02 + 4 * m11 * m12;
      double c2 = 4 * utils::square(m02) + 6 * m01 * m03 + 4 * utils::square(m12) + 6 * m11 * m13;
      double c3 = 12 * m02 * m03 + 12 * m12 * m13;
      double c4 = 9 * utils::square(m03) + 9 * utils::square(m13);
      auto f = [c0, c1, c2, c3, c4](double t){
        return std::sqrt(c0 + t * (c1 + t * (c2 + t * (c3 + t * c4))));
      };
      BezierFitter fitter;
      fitter.solve([f](double t) -> Point{ return {t, f(t)}; }, {0.0, 1.0}, {0.0, 0.0, 1.0, 2.0}, 0.001);
      auto plot_curve = plot.create_bezier_fitter(second_layer, curve_line_style, std::move(fitter));

      double d0 = c4 / 16    + c3 / 8     + c2 / 4 + c1 / 2 + c0;
      double d1 = c4 / 2     + c3 * 3 / 4 + c2     + c1;
      double d2 = c4 * 3 / 2 + c3 * 3 / 2 + c2;
      double d3 = c4 * 2     + c3;
      double d4 = c4;
      auto f2 = [d0, d1, d2, d3, d4](double u){
        return std::sqrt(d0 + u * (d1 + u * (d2 + u * (d3 + u * d4))));
      };
      BezierFitter fitter2;
      fitter2.solve([f2](double u) -> Point { return {u + 0.5, f2(u)}; }, {-0.5, 0.5}, {0.0, 0.0, 1.0, 2.0}, 0.001);
//      auto plot_curve2 = plot.create_bezier_fitter(second_layer, curve_line_style({.line_color = color::orange}), std::move(fitter2));

      //   √(d₀ + d₁u + d₂u² + d₃u³ + d₄u⁴) = ∑ₙ hₙ uⁿ    (Taylor series)
      // where
      //   h₀ = √d₀
      //   hₙ = √d₀ / n gₙ  (for n > 0),
      // where
      //   g₁ =                                               +  1d₁/(2d₀)
      //   g₂ =                                + 1d₂/(1d₀)    -  1d₁/(2d₀)  g₁
      //   g₃ =               +  3d₃/(2d₀)     - 0d₂/(1d₀) g₁ -  3d₁/(4d₀)  g₂
      //   g₄ = +2d₄/(1d₀)    +   d₃/(2d₀)  g₁ - 1d₂/(2d₀) g₂ -  5d₁/(6d₀)  g₃
      //   g₅ = +1d₄/(1d₀) g₁ -   d₃/(4d₀)  g₂ - 2d₂/(3d₀) g₃ -  7d₁/( 8d₀) g₄
      //   g₆ =  0d₄/(2d₀) g₂ -  3d₃/(6d₀)  g₃ - 3d₂/(4d₀) g₄ -  9d₁/(10d₀) g₅
      //   g₇ = -1d₄/(3d₀) g₃ -  5d₃/(8d₀)  g₄ - 4d₂/(5d₀) g₅ - 11d₁/(12d₀) g₆
      //   g₈ = -2d₄/(4d₀) g₄ -  7d₃/(10d₀) g₅ - 5d₂/(6d₀) g₆ - 13d₁/(14d₀) g₇
      //   g₉ = -3d₄/(5d₀) g₅ -  9d₃/(12d₀) g₆ - 6d₂/(7d₀) g₇ - 15d₁/(16d₀) g₈
      //   g₁₀= -4d₄/(6d₀) g₆ - 11d₃/(14d₀) g₇ - 7d₂/(8d₀) g₈ - 17d₁/(18d₀) g₉

      double d1_d0 = d1 / d0;
      double d2_d0 = d2 / d0;
      double d3_d0 = d3 / d0;
      double d4_d0 = d4 / d0;

      Dout(dc::notice, "d0 = " << d0 << "; d1 = " << d1 << "; d2 = " << d2 << "; d3 = " << d3 << "; d4 = " << d4);

      std::vector<double> g(1, 1.0);
/* Determined with parsenb.cxx:
g1 =                                                           1d₁∕(2d₀)
g2 =                                         d₂∕d₀         +  -1d₁∕(2d₀)  ⋆ g1
g3 =                    3d₃∕(2d₀)                          +  -3d₁∕(4d₀)  ⋆ g2
g4 =  2d₄∕d₀         +  1d₃∕(2d₀)  ⋆ g1  +  -d₂∕(2d₀) ⋆ g2 +  -5d₁∕(6d₀)  ⋆ g3
*/
      // g1 = d₁∕(2d₀)
      g.push_back(0.5 * d1_d0);
      // g2 = d₂∕d₀ - d₁∕(2d₀) ⋆ g1
      g.push_back(d2_d0 - 0.5 * d1_d0 * g[1]);
      // g3 = 3d₃∕(2d₀) - 3d₁∕(4d₀) ⋆ g2
      g.push_back(1.5 * d3_d0 - 0.75 * d1_d0 * g[2]);
      // g4 = 2d₄∕d₀ + d₃∕(2d₀) ⋆ g1 - d₂∕(2d₀) ⋆ g2 - 5d₁∕(6d₀) ⋆ g3
      g.push_back(2.0 * d4_d0 + 0.5 * d3_d0 * g[1] - 0.5 * d2_d0 * g[2] - 5.0/6.0 * d1_d0 * g[3]);
      // g5 = d₄∕d₀ ⋆ g1 - d₃∕(4d₀) ⋆ g2 - 2d₂∕(3d₀) ⋆ g3 - 7d₁∕(8d₀) ⋆ g4
//      g.push_back(d4_d0 * g[1] - 0.25 * d3_d0 * g[2] - 2.0/3.0 * d2_d0 * g[3] - 7.0/8.0 * d1_d0 * g[4]);
      // g6 = -3d₃∕(6d₀) ⋆ g3 - 3d₂∕(4d₀) ⋆ g4 - 9d₁∕(10d₀) ⋆ g5
//      g.push_back(-0.5 * d3_d0 * g[3] - 3.0 / 4.0 * d2_d0 * g[4] - 9.0/10.0 * d1_d0 * g[5]);

/* Determined with parsenb.cxx:
g5 = +1d₄∕(1d₀) ⋆ g1 - 1d₃∕(4d₀)  ⋆ g2  - 2d₂∕(3d₀) ⋆ g3 -  7d₁∕(8d₀)  ⋆ g4
g6 =  0d₄∕(2d₀) ⋆ g2 - 3d₃∕(6d₀)  ⋆ g3  - 3d₂∕(4d₀) ⋆ g4 -  9d₁∕(10d₀) ⋆ g5
g7 = -1d₄∕(3d₀) ⋆ g3 - 5d₃∕(8d₀)  ⋆ g4  - 4d₂∕(5d₀) ⋆ g5 - 11d₁∕(12d₀) ⋆ g6
g8 = -2d₄∕(4d₀) ⋆ g4 - 7d₃∕(10d₀) ⋆ g5  - 5d₂∕(6d₀) ⋆ g6 - 13d₁∕(14d₀) ⋆ g7
g9 = -3d₄∕(5d₀) ⋆ g5 - 9d₃∕(12d₀) ⋆ g6  - 6d₂∕(7d₀) ⋆ g7 - 15d₁∕(16d₀) ⋆ g8
*/
      for (int n = 5; n < 400; ++n)
      {
        g.push_back(                                            // If i=0
            - (    n - 6) * d4_d0 / (    n - 4) * g[n - 4]        //   d₄/( d₀) g₁
            - (2 * n - 9) * d3_d0 / (2 * n - 6) * g[n - 3]        //  -d₃/(4d₀) g₂
            - (    n - 3) * d2_d0 / (    n - 2) * g[n - 2]        // -2d₂/(3d₀) g₃
            - (2 * n - 3) * d1_d0 / (2 * n - 2) * g[n - 1]);      // -7d₁/(8d₀) g₄
      }

      double sqrt_d0 = std::sqrt(d0);

      double len1 = 0.5;
      double pu = 0.125; // 0.5^(n+1)
      for (int n = 2; n < g.size(); n += 2)
      {
        len1 += pu * g[n] / (n * (n + 1));
        pu *= 0.25;
      }
      len1 *= 2.0 * sqrt_d0;

      pu = 0.5; // 0.5^n
      for (int n = 1; n < g.size(); ++n)
      {
        Dout(dc::notice, "g[" << n << "]/" << n << " * 0.5^" << n << " = " << std::setprecision(15) << (pu * g[n] / n));
        pu *= 0.5;
      }

      Dout(dc::notice, "len1 = " << std::setprecision(std::numeric_limits<double>::max_digits10) << len1);

      auto f3 = [&](double u){
        double result = 1.0;
        double pu = u;
        for (int n = 1; n < g.size(); ++n)
        {
          result += g[n] / n * pu;
          pu *= u;
        }
        return sqrt_d0 * result;
      };
      BezierFitter fitter3;
      fitter3.solve([f,f3](double u) -> Point { return {u + 0.5, f3(u)/* - f(u + 0.5)*/}; }, {-0.5, 0.5}, {0.0, 0.0, 1.0, 2.0}, 0.001);
      auto plot_curve3 = plot.create_bezier_fitter(second_layer, curve_line_style({.line_color = color::orange}), std::move(fitter3));

      double len2 = plot_bezier.arc_length(1e-15);
      Dout(dc::notice, "len2 = " << std::setprecision(std::numeric_limits<double>::max_digits10) << len2);

      Vector K0 = plot_bezier.curvature(0.5);
      Point P05 = plot_bezier.P(0.5);
      double radius = 1.0 / K0.length();
      plot::Circle plot_curvature_circle;
      if (K0.length() > 1e-9)
        plot_curvature_circle = plot.create_circle(second_layer, line_style({.line_color = color::gray}),
            P05 + radius * K0.direction(), radius);

      // Flush all expose events related to the drawing done above.
      window.set_send_expose_events(true);

      // Block until the user moved a draggable object, then go to the top of loop for a redraw.
      window.handle_input_events();
    }

    event_loop.join();
  }
  catch (AIAlert::Error const& error)
  {
    Dout(dc::warning, error);
  }

  Dout(dc::notice, "Leaving main()");
}
