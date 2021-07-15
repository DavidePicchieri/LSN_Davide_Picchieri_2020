CC = g++
CFLAGS = -Wall -O3 --std=c++11


esercizio3.exe : esercizio3.o random.o functions.o
	$(CC) random.o functions.o esercizio3.o -o esercizio3.exe 
esercizio3.o : esercizio3.cpp
	$(CC) -c esercizio3.cpp -o esercizio3.o $(CFLAGS) 
esercizio2.exe : esercizio2.o random.o
	$(CC) random.o esercizio2.o -o esercizio2.exe 
esercizio2.o : esercizio2.cpp
	$(CC) -c esercizio2.cpp -o esercizio2.o $(CFLAGS) 
esercizio1.exe : esercizio1.o random.o
	$(CC) random.o esercizio1.o -o esercizio1.exe 
esercizio1.o : esercizio1.cpp
	$(CC) -c esercizio1.cpp -o esercizio1.o $(CFLAGS) 
functions.o : functions.cpp functions.h
	$(CC) -c functions.cpp -o functions.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o esercizio1.exe esercizio2.exe esercizio3.exe seed.out
