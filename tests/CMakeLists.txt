include(AICxxProject)

# Test for calculating intersection point of two lines (Line::intersection_with).
add_executable(intersecting_lines intersecting_lines.cxx)
target_link_libraries(intersecting_lines AICxx::cairowindow ${AICXX_OBJECTS_LIST})

add_executable(cairo_window cairo_window.cxx)
target_link_libraries(cairo_window AICxx::cairowindow ${AICXX_OBJECTS_LIST})

add_executable(quadratic_bezier quadratic_bezier.cxx)
target_link_libraries(quadratic_bezier AICxx::cairowindow ${AICXX_OBJECTS_LIST})

add_executable(cubic_bezier cubic_bezier.cxx)
target_link_libraries(cubic_bezier AICxx::cairowindow ${AICXX_OBJECTS_LIST})

add_executable(circle_bezier circle_bezier.cxx)
target_link_libraries(circle_bezier AICxx::cairowindow ${AICXX_OBJECTS_LIST})

add_executable(LachanceSchwartz LachanceSchwartz.cxx)
target_link_libraries(LachanceSchwartz AICxx::cairowindow ${AICXX_OBJECTS_LIST} Eigen3::Eigen)

add_executable(style_test style_test.cxx)
target_link_libraries(style_test ${AICXX_OBJECTS_LIST})
