CC = g++
CFLAGS = -Wall -O3 --std=c++11


esercizio1.exe : esercizio1.o random.o functions.o metropolis.o PDF.o blocking.o
	$(CC) random.o functions.o metropolis.o PDF.o blocking.o esercizio1.o -o esercizio1.exe 
esercizio1.o : esercizio1.cpp
	$(CC) -c esercizio1.cpp -o esercizio1.o $(CFLAGS) 
blocking.o : blocking.cpp blocking.h
	$(CC) -c blocking.cpp -o blocking.o $(CFLAGS)
PDF.o : PDF.cpp PDF.h
	$(CC) -c PDF.cpp -o PDF.o $(CFLAGS)
metropolis.o : metropolis.cpp metropolis.h
	$(CC) -c metropolis.cpp -o metropolis.o $(CFLAGS)
functions.o : functions.cpp functions.h
	$(CC) -c functions.cpp -o functions.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o esercizio1.exe seed.out
