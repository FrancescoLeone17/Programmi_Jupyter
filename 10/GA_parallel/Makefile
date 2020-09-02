CC = g++
CFLAGS = -Wall -O3 --std=c++11

TSP.exe : TSP.o random.o
	mpicxx random.o TSP.o -o TSP.exe
TSP.o : TSP.cpp
	mpicxx -c TSP.cpp -o TSP.o $(CFLAGS)
random.o : random.cpp random.h
	mpicxx -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o TSP.exe
