include(AICxxProject)

# Test for calculating intersection point of two lines (Line::intersection_with).
add_executable(intersecting_lines intersecting_lines.cxx)
target_link_libraries(intersecting_lines AICxx::cairowindow AICxx::symbolic ${AICXX_OBJECTS_LIST})

add_executable(cairo_window cairo_window.cxx)
target_link_libraries(cairo_window AICxx::cairowindow AICxx::symbolic ${AICXX_OBJECTS_LIST})

add_executable(quadratic_bezier quadratic_bezier.cxx)
target_link_libraries(quadratic_bezier AICxx::cairowindow AICxx::symbolic ${AICXX_OBJECTS_LIST})

add_executable(cubic_bezier cubic_bezier.cxx)
target_link_libraries(cubic_bezier AICxx::cairowindow AICxx::symbolic ${AICXX_OBJECTS_LIST})

add_executable(circle_bezier circle_bezier.cxx)
target_link_libraries(circle_bezier AICxx::cairowindow AICxx::symbolic ${AICXX_OBJECTS_LIST})

add_executable(LachanceSchwartz LachanceSchwartz.cxx)
target_link_libraries(LachanceSchwartz AICxx::cairowindow AICxx::symbolic ${AICXX_OBJECTS_LIST} Eigen3::Eigen)

add_executable(style_test style_test.cxx)
target_link_libraries(style_test ${AICXX_OBJECTS_LIST})

add_executable(arc_length_bezier arc_length_bezier.cxx)
target_link_libraries(arc_length_bezier AICxx::cairowindow AICxx::symbolic ${AICXX_OBJECTS_LIST})

add_executable(shortest_arc_length shortest_arc_length.cxx)
target_link_libraries(shortest_arc_length AICxx::cairowindow AICxx::symbolic ${AICXX_OBJECTS_LIST})

add_executable(shortest_arc_length2 shortest_arc_length2.cxx)
target_link_libraries(shortest_arc_length2 AICxx::cairowindow AICxx::symbolic ${AICXX_OBJECTS_LIST})

add_executable(chess_diagram chess_diagram.cxx)
target_link_libraries(chess_diagram AICxx::cairowindow AICxx::symbolic AICxx::chess ${AICXX_OBJECTS_LIST})

add_executable(pango_text pango_text.cxx)
target_link_libraries(pango_text AICxx::cairowindow AICxx::symbolic ${AICXX_OBJECTS_LIST})

add_executable(gradient_descent gradient_descent.cxx Approximation.cxx)
target_link_libraries(gradient_descent AICxx::cairowindow AICxx::symbolic ${AICXX_OBJECTS_LIST} Eigen3::Eigen)

add_executable(parabolic_scale parabolic_scale.cxx)
target_link_libraries(parabolic_scale AICxx::cairowindow AICxx::symbolic ${AICXX_OBJECTS_LIST} Eigen3::Eigen)
