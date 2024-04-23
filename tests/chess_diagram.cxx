#include "sys.h"
#include "cairowindow/Window.h"
#include "cairowindow/Layer.h"
#include "cairowindow/chess/Diagram.h"
#include "utils/AIAlert.h"
#include "utils/debug_ostream_operators.h"
#include <thread>
#include <iostream>
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
    Window window("A chess diagram", 1200, 900);

    // Create a new layer with a white background.
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

    // Create and draw a chess diagram.
    // ChessDiagramStyleParamsDefault
    chess::Diagram diagram(window.geometry(), {{.coordinate_margin = 50.0, .top_margin = 80.0}}, "My First Chess Diagram", {});

    diagram.create_svg_surface("chess_diagram.svg" COMMA_CWDEBUG_ONLY("diagram"));
    diagram.set_need_print();
    diagram.add_to(background_layer);
    diagram.reset_need_print();

    while (true)
    {
      if (diagram.need_print())
        Dout(dc::notice, "========= STARTING PRINT FRAME ===========");

      // Suppress immediate updating of the window for each created item, in order to avoid flickering.
      window.set_send_expose_events(false);

      // ... stuff here ...

      // Flush all expose events related to the drawing done above.
      window.set_send_expose_events(true);

      if (diagram.need_print())
        diagram.reset_need_print();

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
