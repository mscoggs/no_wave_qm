CC=g++


compile:
	$(CC) -o main main.cpp grid.cpp point.cpp math_functions.cpp input_functions.cpp -std=gnu++11

run:
	$(CC) -o main main.cpp grid.cpp point.cpp math_functions.cpp input_functions.cpp -std=gnu++11
	./main

run_test:
	$(CC) -o math_test math_test.cpp math_functions.cpp
	./math_test

test:
	$(CC) -o test test.cpp math_functions.cpp
	./test
