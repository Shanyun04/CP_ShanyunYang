# Add test executable for IsingSystemCubic.cpp
enable_testing()
add_executable(Unittest_for_cubic_catch2 Unittest_for_cubic.cpp)
target_include_directories(Unittest_for_cubic_catch2 PRIVATE /home/xuwu/CP2024_code/CP2024-midterm/CP2024-Group1)
#target_include_directories(Unittest_for_cubic_catch2 PRIVATE /home/xuwu/CP2024_code/CP2024-midterm/CP2024-Group1)
target_compile_options(Unittest_for_cubic_catch2 PRIVATE -Wall -Wextra)
target_link_libraries(Unittest_for_cubic_catch2 Catch2::Catch2WithMain)
add_test(Unittest_for_cubic Unittest_for_cubic_catch2)

# Add test executable for IsingSystemCubic.cpp
add_executable(Unittest_for_honeycomb_catch2 Unittest_for_honeycomb.cpp)
#target_include_directories(Unittest_for_cubic_catch2 PRIVATE /home/xuwu/CP2024_code/CP2024-midterm/CP2024-Group1)
target_compile_options(Unittest_for_honeycomb_catch2 PRIVATE -Wall -Wextra)
target_link_libraries(Unittest_for_honeycomb_catch2 Catch2::Catch2WithMain)
add_test(Unittest_for_honeycomb Unittest_for_honeycomb_catch2)

# Add test executable for IsingSystemCubic.cpp
add_executable(Unittest_for_kagome_catch2 Unittest_for_kagome.cpp)
#target_include_directories(Unittest_for_cubic_catch2 PRIVATE /home/xuwu/CP2024_code/CP2024-midterm/CP2024-Group1)
target_compile_options(Unittest_for_kagome_catch2 PRIVATE -Wall -Wextra)
target_link_libraries(Unittest_for_kagome_catch2 Catch2::Catch2WithMain)
add_test(Unittest_for_kagome Unittest_for_kagome_catch2)

# Add test executable for IsingSystem.cpp
# Add test executable for IsingSystemCubic.cpp
add_executable(Unittest_for_Cubic_catch2 Unittest_for_Cubic.cpp)
#target_include_directories(Unittest_for_cubic_catch2 PRIVATE /home/xuwu/CP2024_code/CP2024-midterm/CP2024-Group1)
target_compile_options(Unittest_for_Cubic_catch2 PRIVATE -Wall -Wextra)
target_link_libraries(Unittest_for_Cubic_catch2 Catch2::Catch2WithMain)
add_test(Unittest_for_Cubic Unittest_for_Cubic_catch2)
 #Add test executable for IsingSystem.cpp
#add_executable(IsingSystemtest_catch2 IsingSystemtest.cpp)
#target_include_directories(IsingSystemtest_catch2 PRIVATE /home/xuwu/Example_catch2/)
#target_compile_options(IsingSystemtest_catch2 PRIVATE -Wall -Wextra)
