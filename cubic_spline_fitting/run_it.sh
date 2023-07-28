# Linux only... maybe?
# Run this from the normal directory, make sure ./build exists
gcc ./general_cubic_spline_fit.c -o ./build/a -lm && ./build/a && python print_splines.py