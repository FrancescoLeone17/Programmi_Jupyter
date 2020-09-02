CC = g++
CFLAGS = -Wall -O3 --std=c++11

2_2.exe : 2_2.o random.o
	$(CC) random.o 2_2.o -o 2_2.exe
2_2.o : 2_2.cpp
	$(CC) -c 2_2.cpp -o 2_2.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o 2_2.exe
