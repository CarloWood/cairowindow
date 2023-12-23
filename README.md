# cairowindow submodule

This repository is a [git submodule](https://git-scm.com/book/en/v2/Git-Tools-Submodules)
providing C++ classes for larger projects, including (all in namespace cairowindow):

* ``Window`` : A new native window, assuming you are on X11. Full support for tripple buffering to allow seamless, non-stop updating.
* ``Layer`` : Abstraction of a layer. Represents a drawing surface, can be smaller than and have an offset relative to the Window.
* ``LayerRegion`` : A rectangular portion of a layer that has a `draw` function associated with it (base class to all drawable objects).
* ``MultiRegion`` : Convenience base class for drawable objects that exist of more than one layer region.
* ``EventLoop`` : Returned by Window::run(). Causes the window to be opened and everything that was added to it to be drawn. Destructing the returned `EventLoop` blocks until the window was closed. You want to destroy the `EventLoop` before destroying any of the above objects, because destroying those will make them disappear (no longer being drawn).
* ``Color`` : A color object.
* ``Rectangle`` : Describes a rectangle with an (optional) offset relative to the main Window.
* ``StrokeExtents`` : Another rectangle object, used to keep track of what area needs to be redrawn for a given `LayerRegion` upon expose.

In namespace `cairowindow::plot` we currently have:

* ``Plot`` : An object representing a plot chart with title, axes, data points etc.
* ``Range`` : Describes the minimum and maximum values of the axes of a plot.

In namespace `cairowindow::draw` we currently have:

* ``Line`` : A line.
* ``Rectangle`` : A rectangle.
* ``Text`` : Some text.
* ``PlotArea`` : The axes and ticks of a plot.

The root project should be using
[cmake](https://cmake.org/overview/)
[cwm4](https://github.com/CarloWood/cwm4) and
[cwds](https://github.com/CarloWood/cwds).

## Checking out a project that uses the cairowindow submodule.

To clone a project example-project that uses cairowindow simply run:

    git clone --recursive <URL-to-project>/example-project.git
    cd example-project
    AUTOGEN_CMAKE_ONLY=1 ./autogen.sh

The ``--recursive`` is optional because ``./autogen.sh`` will fix
it when you forgot it.

In order to use ``cmake`` configure as usual, for example to do a debug build with 16 cores:

    mkdir build_debug
    cmake -S . -B build_debug -DCMAKE_MESSAGE_LOG_LEVEL=DEBUG -DCMAKE_BUILD_TYPE=Debug -DCMAKE_VERBOSE_MAKEFILE=ON -DEnableDebugGlobal:BOOL=OFF
    cmake --build build_debug --config Debug --parallel 16

Or to make a release build:

    mkdir build_release
    cmake -S . -B build_release -DCMAKE_BUILD_TYPE=Release
    cmake --build build_release --config Release --parallel 16

## Adding the cairowindow submodule to a project

To add this submodule to a project, that project should already
be set up to use [cwm4](https://github.com/CarloWood/cwm4).

Simply execute the following in a directory of that project
where you want to have the ``cairowindow`` subdirectory (the
root of the project is recommended as that is the only thing
I've tested so far):

    git submodule add https://github.com/CarloWood/cairowindow.git

This should clone cairowindow into the subdirectory ``cairowindow``, or
if you already cloned it there, it should add it.

### Using cmake

Check out the submodules [cwds](https://github.com/CarloWood/cwds), [cwm4](https://github.com/CarloWood/cwm4)
and [ai-utils](https://github.com/CarloWood/ai-utils) in the root of the project:

    git submodule add https://github.com/CarloWood/cwds.git
    git submodule add https://github.com/CarloWood/cwm4.git
    git submodule add https://github.com/CarloWood/ai-utils.git utils

The easiest way to use libcwd is by using [gitache](https://github.com/CarloWood/gitache).

For that to happen create in the root of the project (that uses utils)
a directory ``cmake/gitache-configs`` and put in it the file ``libcwd_r.cmake``
with the content:

    gitache_config(
      GIT_REPOSITORY
        "https://github.com/CarloWood/libcwd.git"
      GIT_TAG
        "master"
      CMAKE_ARGS
        "-DEnableLibcwdAlloc:BOOL=OFF -DEnableLibcwdLocation:BOOL=ON"
    )

Add the variable ``GITACHE_ROOT`` to your environment,
for example add to your ``~/.bashrc`` the line:

    export GITACHE_ROOT="/opt/gitache"

Add the following lines to the ``CMakeLists.txt`` in the
root of the project (directly under the ``project`` line):

    # Begin of gitache configuration.
    set(GITACHE_PACKAGES libcwd_r)
    include(cwm4/cmake/StableGitache)
    # End of gitache configuration.

    include(cwm4/cmake/AICxxProject)
    include(AICxxSubmodules)

``add_subdirectory`` is not necessary for ``cwds``, ``cwm4``, ``utils`` or ``cairowindow``.

See for example the root [MakeLists.txt](https://github.com/CarloWood/machine-learning/blob/master/CMakeLists.txt) of machine-learning.

Finally, linking is done by adding ``${AICXX_OBJECTS_LIST}`` to
the appropriate ``target_link_libraries``.

For example,

    include(AICxxProject)

    add_executable(register_test register_test.cxx)
    target_link_libraries(register_test PRIVATE ${AICXX_OBJECTS_LIST})

See this [MakeLists.txt](https://github.com/CarloWood/machine-learning/blob/master/src/CMakeLists.txt)
of machine-learning for a complete example.

Finally, run

    ./autogen.sh

to let cwm4 do its magic, and commit all the changes.

Checkout [machine-learning](https://github.com/CarloWood/machine-learning)
for an example of a project that uses this submodule.
