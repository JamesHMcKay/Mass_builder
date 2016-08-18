calc_amplitudes: calc_amplitudes.o generate_code
	g++ -g -o calc_amplitudes calc_amplitudes.o

calc_amplitudes.o: calc_amplitudes.cpp
	g++ -g -c calc_amplitudes.cpp

generate_code: generate_code.o
	g++ -g -o generate_code generate_code.o

generate_code.o: generate_code.cpp
	g++ -g -c generate_code.cpp