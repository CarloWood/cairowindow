#include "sys.h"
#include "cairowindow/Window.h"
#include "cairowindow/Layer.h"
#include "cairowindow/chess/Diagram.h"
#include "cairowindow/draw/ChessPiece.h"
#include "utils/AIAlert.h"
#include "utils/debug_ostream_operators.h"
#include "utils/u8string_to_filename.h"
#include <boost/filesystem/path.hpp>
#include <thread>
#include <iostream>
#include "debug.h"

int main(int argc, char* argv[])
{
  Debug(NAMESPACE_DEBUG::init());
  Dout(dc::notice, "Entering main()");

  // Process command line arguments.
  boost::filesystem::path output_file_name;
  std::string title;
  std::string FEN;
  bool overwrite = false;
  std::vector<std::string> args(argv + 1, argv + argc);
  for (size_t i = 0; i < args.size(); ++i)
  {
    if (args[i] == "-o" && i + 1 < args.size())
      output_file_name = args[++i];
    else if (args[i] == "--title" && i + 1 < args.size())
      title = args[++i];
    else if (args[i] == "--replace")
      overwrite = true;
    else
      FEN = args[i];
  }

  if (FEN.empty())
  {
    std::cerr << "Usage: " << argv[0] << " [-o SVGoutfile] [--replace] [--title \"Diagram header\"] <FEN-code>" << std::endl;
    return 1;
  }
  else if (output_file_name.empty())
  {
    std::string str = FEN + ".svg";
    std::u8string u8_string(str.begin(), str.end());
    output_file_name = utils::u8string_to_filename(u8_string);
  }

  using namespace cairowindow;
  using Window = cairowindow::Window;

  // Create a window.
  Window window("A chess diagram", 200, 200);

  // Create a new layer with a white background.
  auto board_layer = window.create_background_layer<Layer>(color::white COMMA_DEBUG_ONLY("board_layer"));

  // Create another layer.
  auto pieces_layer = window.create_layer<Layer>({} COMMA_DEBUG_ONLY("pieces_layer"));

  // Open the window and start drawing.
  std::thread event_loop([&](){
    // Open window, handle event loop. This must be constructed after the draw stuff, so that it is destructed first!
    // Upon destruction it blocks until the event loop thread finished (aka, the window was closed).
    EventLoop event_loop = window.run();
    event_loop.set_cleanly_terminated();
  });

  try
  {
    // Create and draw a chess diagram.
    chess::Diagram diagram(window.geometry(),
        {{.coordinate_margin = 20.0, .top_margin = !title.empty() ? 40.0 : 10.0, .margin = 10.0}}, title, {});

    diagram.create_svg_surface(output_file_name.string(), overwrite COMMA_CWDEBUG_ONLY("diagram"));
    diagram.set_need_print();
    diagram.add_to(board_layer);

    while (true)
    {
      if (diagram.need_print())
        Dout(dc::notice, "========= STARTING PRINT FRAME ===========");

      // Suppress immediate updating of the window for each created item, in order to avoid flickering.
      window.set_send_expose_events(false);

      // Load a chess position.
      if (!diagram.load_FEN(pieces_layer, FEN, {}))
        THROW_ALERT("Invalid FEN code \"[FENCODE]\"", AIArgs("[FENCODE]", FEN));

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
