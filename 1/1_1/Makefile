CC = g++
CFLAGS = -Wall -O3 --std=c++11

1_1.exe : 1_1.o random.o
	$(CC) random.o 1_1.o -o 1_1.exe
1_1.o : 1_1.cpp
	$(CC) -c 1_1.cpp -o 1_1.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o 1_1.exe
