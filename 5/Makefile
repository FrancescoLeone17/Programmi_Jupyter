CC = g++
CFLAGS = -Wall -O3 --std=c++11

1.exe : 1.o random.o
	$(CC) random.o 1.o -o 1.exe
1.o : 1.cpp
	$(CC) -c 1.cpp -o 1.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o 1.exe
