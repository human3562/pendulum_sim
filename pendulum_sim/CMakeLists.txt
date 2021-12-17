cmake_minimum_required (VERSION 3.15)

include(FetchContent) 
FetchContent_Declare(mahi-gui GIT_REPOSITORY https://github.com/mahilab/mahi-gui.git) 
FetchContent_MakeAvailable(mahi-gui)

add_executable (pendulum_sim "pendulum_sim.cpp" "pendulum_sim.h")

target_link_libraries(pendulum_sim mahi::gui)
