CC = g++
CFLAGS = -Wall -O3 --std=c++11

TSP.exe : TSP.o random.o
	$(CC) random.o TSP.o -o TSP.exe
TSP.o : TSP.cpp
	$(CC) -c TSP.cpp -o TSP.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o TSP.exe
