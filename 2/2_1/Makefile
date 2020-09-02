CC = g++
CFLAGS = -Wall -O3 --std=c++11

2_1.exe : 2_1.o random.o
	$(CC) random.o 2_1.o -o 2_1.exe
2_1.o : 2_1.cpp
	$(CC) -c 2_1.cpp -o 2_1.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o 2_1.exe
