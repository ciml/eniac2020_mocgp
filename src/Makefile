CC = g++
CFLAGS = -g -Wall -O3 -std=c++17
RM = rm ./bin/*

main:
	${CC} ${CFLAGS} main.cpp -o ./bin/cgp -lbdd -DMO
run:
	LD_LIBRARY_PATH=./lib ./bin/cgp 1 circuits/cm42a.ep 1000000 50 SAM NSGA2 100