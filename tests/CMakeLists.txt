#############################################
# Configure test project
#############################################
cmake_minimum_required (VERSION 3.5)
project (atfft_tests)

#############################################
# Define tests
#############################################
add_executable (complex_dft test_functions.c complex_dft.c)
target_link_libraries (complex_dft atfft)
add_test (complex_dft complex_dft)

#############################################
# Add them all to one target
#############################################
add_custom_target (tests DEPENDS complex_dft)
