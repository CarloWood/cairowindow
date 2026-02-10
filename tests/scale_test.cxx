#include "sys.h"
#include "math/CubicPolynomial.h"
#include "gradient_descent2/AnalyzedCubic.h"
#include "gradient_descent2/Sample.h"
#include "cairowindow/BezierCurve.h"
#include "cairowindow/Window.h"
#include "cairowindow/Layer.h"
#include "cairowindow/Plot.h"
#include "cairowindow/Direction.h"
#include "cairowindow/draw/Line.h"
#include "cairowindow/draw/Point.h"
#include "utils/AIAlert.h"
#include "utils/debug_ostream_operators.h"
#include "utils/ColorPool.h"
#include "utils/square.h"
#include "utils/UniqueID.h"
#include <thread>
#include "debug.h"

int main()
{
  Debug(NAMESPACE_DEBUG::init());
  Dout(dc::notice, "Entering main()");

  try
  {
    using namespace cairowindow;
    using Window = cairowindow::Window;
    using Sample = gradient_descent::Sample;

    // Create a window.
    Window window("Style Test", 1200, 900);

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

    utils::UniqueIDContext<int> label_context;

    Sample const s1{9.29867, 45.6027, 0.00234467 COMMA_CWDEBUG_ONLY(label_context)};
    Sample const s2{9.73483, 45.8303, 0.951923 COMMA_CWDEBUG_ONLY(label_context)};
    Sample const s3{9.87073, 45.9883, 1.26854 COMMA_CWDEBUG_ONLY(label_context)};
    Sample const s4{10.9106, 47.4864, 2.46341 COMMA_CWDEBUG_ONLY(label_context)};

    double w1 = std::min(s1.w(), s2.w());
    double w2 = std::max(s1.w(), s2.w());
    double dw = w2 - w1;
    double const w_min = 9.0; //w1 - 0.2 * dw;
    double const w_max = 12.0; //w2 + 0.2 * dw;

    double Lw1 = std::min(s1.Lw(), s2.Lw());
    double Lw2 = std::max(s1.Lw(), s2.Lw());
    double dLw = Lw2 - Lw1;
    double const Lw_min = 45.5; //Lw1 - 0.2 * dLw;
    double const Lw_max = 48.0; //Lw2 + 0.2 * dLw;

    // Create and draw plot area.
    plot::Plot plot(window.geometry(), draw::PlotAreaStyle({.color = color::orange}),
        "Scale Test", {},
        "x", {},
        "y", {});
    plot.set_xrange({w_min, w_max});
    plot.set_yrange({Lw_min, Lw_max});
    plot.add_to(background_layer, true);

    utils::ColorPool<32> color_pool;
    draw::PointStyle point_style({.color_index = color_pool.get_and_use_color(), .filled_shape = 1});
    draw::TextStyle label_style({.position = draw::centered_left_of, .font_size = 18.0, .offset = 10});
    draw::TextStyle slider_style({.position = draw::centered_below, .font_size = 18.0, .offset = 10});
    draw::LineStyle curve_line_style({.line_width = 1.0});
    draw::LineStyle solid_line_style({.line_color = color::black, .line_width = 1.0});
    draw::LineStyle line_style({.line_color = color::black, .line_width = 1.0, .dashes = {10.0, 5.0}});
    draw::ArcStyle arc_style({.line_color = color::blue, .line_width = 1.0});

    // Create a point P₀.
    auto plot_P0 = plot.create_point(second_layer, point_style, {s1.w(), s1.Lw()});
    // Create a point P₁.
    auto plot_P1 = plot.create_point(second_layer, point_style, {s2.w(), s2.Lw()});
    // Create a point P.
    auto plot_P2 = plot.create_point(second_layer, point_style, {s3.w(), s3.Lw()});
    // Create a point P.
    auto plot_P3 = plot.create_point(second_layer, point_style, {s4.w(), s4.Lw()});

    double const circle_radius = 0.1 * (w_max - w_min);
    // Create a point Q0 on a circle.
    double phi0 = std::atan(s1.dLdw());
    auto plot_Q0 = plot.create_point(second_layer, point_style({.color_index = 2}), plot_P0 + circle_radius * Direction{phi0});
    // Create a point Q1 on a circle.
    double phi1 = std::atan(s2.dLdw());
    auto plot_Q1 = plot.create_point(second_layer, point_style({.color_index = 2}), plot_P1 + circle_radius * Direction{phi1});
    // Create a point Q2 on a circle.
    double phi2 = std::atan(s3.dLdw());
    auto plot_Q2 = plot.create_point(second_layer, point_style({.color_index = 2}), plot_P2 + circle_radius * Direction{phi2});
    // Create a point Q3 on a circle.
    double phi3 = std::atan(s4.dLdw());
    auto plot_Q3 = plot.create_point(second_layer, point_style({.color_index = 2}), plot_P3 + circle_radius * Direction{phi3});

    // Make all points draggable.
    window.register_draggable(plot, &plot_P0, [&plot, &plot_P0, &plot_Q0](Point const& new_position) -> Point
        {
          auto translation = new_position - plot_P0;
          plot_Q0.move_to(plot_Q0 + translation);
          return new_position;
        }
    );
    window.register_draggable(plot, &plot_P1, [&plot, &plot_P1, &plot_Q1](Point const& new_position) -> Point
        {
          auto translation = new_position - plot_P1;
          plot_Q1.move_to(plot_Q1 + translation);
          return new_position;
        }
    );
    window.register_draggable(plot, &plot_P2, [&plot, &plot_P2, &plot_Q2](Point const& new_position) -> Point
        {
          auto translation = new_position - plot_P2;
          plot_Q2.move_to(plot_Q2 + translation);
          return new_position;
        }
    );
    window.register_draggable(plot, &plot_P3, [&plot, &plot_P3, &plot_Q3](Point const& new_position) -> Point
        {
          auto translation = new_position - plot_P3;
          plot_Q3.move_to(plot_Q3 + translation);
          return new_position;
        }
    );

    window.register_draggable(plot, &plot_Q0, [&plot_P0, circle_radius](Point const& new_position) -> Point
        {
          Direction D0{plot_P0, new_position};
          if (D0.x() < 0.0)
            D0.negate();
          if (D0.x() < 1e-6)
            D0 = Direction{Point{1e-6, 1.0}};
          return plot_P0 + circle_radius * D0;
        }
    );
    window.register_draggable(plot, &plot_Q1, [&plot_P1, circle_radius](Point const& new_position) -> Point
        {
          Direction D1{plot_P1, new_position};
          if (D1.x() < 0.0)
            D1.negate();
          if (D1.x() < 1e-6)
            D1 = Direction{Point{1e-6, 1.0}};
          return plot_P1 + circle_radius * D1;
        }
    );
    window.register_draggable(plot, &plot_Q2, [&plot_P2, circle_radius](Point const& new_position) -> Point
        {
          Direction D2{plot_P2, new_position};
          if (D2.x() < 0.0)
            D2.negate();
          if (D2.x() < 1e-6)
            D2 = Direction{Point{1e-6, 1.0}};
          return plot_P2 + circle_radius * D2;
        }
    );
    window.register_draggable(plot, &plot_Q3, [&plot_P3, circle_radius](Point const& new_position) -> Point
        {
          Direction D3{plot_P3, new_position};
          if (D3.x() < 0.0)
            D3.negate();
          if (D3.x() < 1e-6)
            D3 = Direction{Point{1e-6, 1.0}};
          return plot_P3 + circle_radius * D3;
        }
    );

//    auto slider_d = plot.create_slider(second_layer, {978, 83, 7, 400}, 0.0520743, -0.01, 0.01);
//    auto slider_d_label = plot.create_text(second_layer, slider_style, Pixel{978, 483}, "d");

    while (true)
    {
      // Suppress immediate updating of the window for each created item, in order to avoid flickering.
      window.set_send_expose_events(false);

      // Draw a circle around P₀.
      auto plot_circle0 = plot.create_circle(second_layer, solid_line_style({.line_color = color::gainsboro}), plot_P0, circle_radius);

      // Draw a label for P₀.
      auto P0_label = plot.create_text(second_layer, label_style({.position = draw::centered_right_of}), plot_P0, "P₀");

      // Draw a arrow from P₀ to Q₀.
      auto plot_P0_Q0 = plot.create_connector(second_layer, line_style({.line_color = color::orange}), plot_P0, plot_Q0);

      // Get the direction vector of L'(w₀).
      Direction const D0 = (plot_Q0 - plot_P0).direction();

      // Draw a circle around P₁.
      auto plot_circle1 = plot.create_circle(second_layer, solid_line_style({.line_color = color::gainsboro}), plot_P1, circle_radius);

      // Draw a label for P₁.
      auto P1_label = plot.create_text(second_layer, label_style({.position = draw::centered_right_of}), plot_P1, "P₁");

      // Draw a arrow from P₁ to Q₁.
      auto plot_P1_Q1 = plot.create_connector(second_layer, line_style({.line_color = color::orange}), plot_P1, plot_Q1);

      // Get the direction vector of L'(w₁).
      Direction const D1 = (plot_Q1 - plot_P1).direction();

      // Draw a circle around P₂.
      auto plot_circle2 = plot.create_circle(second_layer, solid_line_style({.line_color = color::gainsboro}), plot_P2, circle_radius);

      // Draw a label for P₂.
      auto P2_label = plot.create_text(second_layer, label_style({.position = draw::centered_right_of}), plot_P2, "P₂");

      // Draw a arrow from P₂ to Q₂.
      auto plot_P2_Q2 = plot.create_connector(second_layer, line_style({.line_color = color::orange}), plot_P2, plot_Q2);

      // Get the direction vector of L'(w₂).
      Direction const D2 = (plot_Q2 - plot_P2).direction();

      // Draw a circle around P₃.
      auto plot_circle3 = plot.create_circle(second_layer, solid_line_style({.line_color = color::gainsboro}), plot_P3, circle_radius);

      // Draw a label for P₃.
      auto P3_label = plot.create_text(second_layer, label_style({.position = draw::centered_right_of}), plot_P3, "P₃");

      // Draw a arrow from P₃ to Q₃.
      auto plot_P3_Q3 = plot.create_connector(second_layer, line_style({.line_color = color::orange}), plot_P3, plot_Q3);

      // Get the direction vector of L'(w₃).
      Direction const D3 = (plot_Q3 - plot_P3).direction();

      double w0 = plot_P0.x();
      double Lw0 = plot_P0.y();
      double dLdw0 = D0.y() / D0.x();
      double w1 = plot_P1.x();
      double Lw1 = plot_P1.y();
      double dLdw1 = D1.y() / D1.x();
      double w2 = plot_P2.x();
      double Lw2 = plot_P2.y();
      double dLdw2 = D2.y() / D2.x();
      double w3 = plot_P3.x();
      double Lw3 = plot_P3.y();
      double dLdw3 = D3.y() / D3.x();
      Dout(dc::notice, "w₀ = " << w0 << "; L(w₀) = " << Lw0 << "; L'(w₀) = " << dLdw0 <<
          "; w₁ = " << w1 << "; L(w₁) = " << Lw1 << "; L'(w₁) = " << dLdw1);
      Dout(dc::notice, "w₂ = " << w2 << "; L(w₂) = " << Lw2 << "; L'(w₂) = " << dLdw2 <<
          "; w₃ = " << w3 << "; L(w₃) = " << Lw3 << "; L'(w₃) = " << dLdw3);

      // Let g(w) = a + b w + c w² + d w³ be a cubic that has the same value and first derivative in w₀ and w₁.
      math::CubicPolynomial<double> g;
      g.initialize(w0, Lw0, dLdw0, w1, Lw1, dLdw1);

      // Plot the cubic.
      BezierFitter fitter0([&](double w) -> Point{ return {w, g(w)}; }, plot.viewport());
      auto plot_g = plot.create_bezier_fitter(second_layer, curve_line_style, std::move(fitter0));

      // The cubic through w2 and w3.
      math::CubicPolynomial<double> h;
      h.initialize(w2, Lw2, dLdw2, w3, Lw3, dLdw3);

      // Plot the cubic.
      BezierFitter fitter1([&](double w) -> Point{ return {w, h(w)}; }, plot.viewport());
      auto plot_h = plot.create_bezier_fitter(second_layer, curve_line_style({.line_color = color::blue}), std::move(fitter1));

      using namespace gradient_descent;

      // Get the minimum of g.
      AnalyzedCubic ag;
      ag.initialize(g, ExtremeType::minimum);
      double const e = ag.get_extreme();
      double const ge = g(e);
      // Draw a horizontal line through the minimum at (e, g(e)).
      plot::Line plot_minimum_line{{Point{e, ge}, Direction::left}};
      plot.add_line(second_layer, line_style, plot_minimum_line);

      // Create and draw 0.9 * g(x) + 0.1 * g(e).
      math::CubicPolynomial<double> g09 = g;
      for (int i = 0; i < 4; ++i)
        g09[i] *= 0.9;
      g09[0] += 0.1 * g(e);
      BezierFitter fitter09([&](double w) -> Point{ return {w, g09(w)}; }, plot.viewport());
      auto plot_g09 = plot.create_bezier_fitter(second_layer, curve_line_style({.line_color = color::green}), std::move(fitter09));

      // Create and draw 1.1 * g(x) - 0.1 * g(e).
      math::CubicPolynomial<double> g11 = g;
      for (int i = 0; i < 4; ++i)
        g11[i] *= 1.1;
      g11[0] -= 0.1 * g(e);
      BezierFitter fitter11([&](double w) -> Point{ return {w, g11(w)}; }, plot.viewport());
      auto plot_g11 = plot.create_bezier_fitter(second_layer, curve_line_style({.line_color = color::red}), std::move(fitter11));

      // Flush all expose events related to the drawing done above.
      window.set_send_expose_events(true);

      // Block until a redraw is necessary (for example because the user moved a draggable object,
      // or wants to print the current drawing) then go to the top of loop for a redraw.
      if (!window.handle_input_events())
        break;          // Program must be terminated.
    }

    event_loop.join();
  }
  catch (AIAlert::Error const& error)
  {
    Dout(dc::warning, error);
  }

  Dout(dc::notice, "Leaving main()");
}
