mass_builder: mass_builder.o generate_code generate_figures
	g++ -g -o mass_builder mass_builder.o -I

calc_amplitudes: calc_amplitudes.o generate_code generate_figures
	g++ -g -o calc_amplitudes calc_amplitudes.o

calc_amplitudes.o: calc_amplitudes.cpp
	g++ -g -c calc_amplitudes.cpp

generate_code: generate_code.o
	g++ -g -o generate_code generate_code.o

generate_code.o: generate_code.cpp
	g++ -g -c generate_code.cpp

generate_figures: generate_figures.o
	g++ -g -o generate_figures generate_figures.o

generate_figures.o: generate_figures.cpp
	g++ -g -c generate_figures.cpp