CC = g++
CFLAGS = -Wall -O3 --std=c++11

1_2.exe : 1_2.o random.o
	$(CC) random.o 1_2.o -o 1_2.exe
1_2.o : 1_2.cpp
	$(CC) -c 1_2.cpp -o 1_2.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o 1_2.exe
