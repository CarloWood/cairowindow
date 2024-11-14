#include "sys.h"
#include "AlgorithmEvent.h"

#ifdef CWDEBUG
//static
cairowindow::draw::TextStyle const AlgorithmEvent::s_label_style{{ .position = cairowindow::draw::centered_left_of, .font_size = 18.0,
  .offset = 10}};
//static
cairowindow::draw::ConnectorStyle const AlgorithmEvent::s_difference_expected_style{{.line_color = cairowindow::color::blue, .line_width = 1.0}};
//static
cairowindow::draw::ConnectorStyle const AlgorithmEvent::s_indicator_style{{.line_width = 1}};
#endif
