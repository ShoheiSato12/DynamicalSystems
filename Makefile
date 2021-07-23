all:
	g++ -Wall  ./source/biffurcation.cpp  ./source/LinearAlgebra.cpp  ./source/lyapunov.cpp  ./source/main.cpp   ./source/plotting.cpp  ./source/rungekutta4thSquare.cpp  ./source/system.cpp -o main.o -lm
run:
	g++ -Wall  ./source/biffurcation.cpp  ./source/LinearAlgebra.cpp  ./source/lyapunov.cpp  ./source/main.cpp   ./source/plotting.cpp  ./source/rungekutta4thSquare.cpp  ./source/system.cpp -o main.o -lm
	./main.o
clear:
	find . -type f -executable -exec rm '{}' \;