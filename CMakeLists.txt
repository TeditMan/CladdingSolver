# CMakeList.txt : CMake project for CladdingSolver, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.10)

# Enable Hot Reload for MSVC compilers if supported.
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

project ("CladdingSolver")

# Add source to this project's executable.
add_executable (CladdingSolver "Source/CladdingSolver.cpp"  "Headers/Grid.h" "Source/Grid.cpp" "Headers/InputData.h" "Source/InputData.cpp" "main.cpp" "Headers/CladdingSolver.h" "Headers/FDSystem.h" "Source/FDSystem.cpp" "Headers/Discretization.h" "Source/Discretization.cpp" "Headers/Triplet.h")

target_include_directories(CladdingSolver PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/Headers
)

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET CladdingSolver PROPERTY CXX_STANDARD 20)
endif()

# TODO: Add tests and install targets if needed.
