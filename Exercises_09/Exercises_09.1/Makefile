CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o City.o
	$(CC) random.o City.o main.o -o main.exe
main.o : main.cpp main.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
City.o : City.cpp City.h
	$(CC) -c City.cpp -o City.o $(CFLAGS)

clean :
	rm *.o main.exe
	rm *.dat
