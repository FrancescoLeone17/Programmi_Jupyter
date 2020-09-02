CC = g++
CFLAGS = -Wall -O3 --std=c++11

1_3.exe : 1_3.o random.o
	$(CC) random.o 1_3.o -o 1_3.exe
1_3.o : 1_3.cpp
	$(CC) -c 1_3.cpp -o 1_3.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o 3.exe
