include(AICxxProject)

# Test for calculating intersection point of two lines (Line::intersection_with).
add_executable(intersecting_lines intersecting_lines.cxx)
target_link_libraries(intersecting_lines AICxx::cairowindow ${AICXX_OBJECTS_LIST})

add_executable(cairo_window cairo_window.cxx)
target_link_libraries(cairo_window AICxx::cairowindow ${AICXX_OBJECTS_LIST})

add_executable(quadratic_bezier quadratic_bezier.cxx)
target_link_libraries(quadratic_bezier AICxx::cairowindow ${AICXX_OBJECTS_LIST})

add_executable(quadratic_bezier2 quadratic_bezier2.cxx)
target_link_libraries(quadratic_bezier2 AICxx::cairowindow ${AICXX_OBJECTS_LIST})
