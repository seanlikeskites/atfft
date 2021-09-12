#############################################
# Configure test project
#############################################
cmake_minimum_required (VERSION 3.5)
project (atfft_tests)

#############################################
# Define tests
#############################################
add_executable (test_dft test_dft.c)
target_link_libraries (test_dft atfft)
add_test (test_dft test_dft)

#############################################
# Add them all to one target
#############################################
add_custom_target (tests DEPENDS test_dft)
