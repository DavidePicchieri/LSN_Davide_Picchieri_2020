CC = mpicxx
CFLAGS = -Wall -O3 --std=c++11


esercizio10.exe : esercizio10.o random.o genetic.o functions.o
	$(CC) random.o genetic.o functions.o esercizio10.o -o esercizio10.exe 
esercizio10.o : esercizio10.cpp
	$(CC) -c esercizio10.cpp -o esercizio10.o $(CFLAGS) 
genetic.o : genetic.cpp genetic.h
	$(CC) -c genetic.cpp -o genetic.o $(CFLAGS)
functions.o : functions.cpp functions.h
	$(CC) -c functions.cpp -o functions.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o esercizio10.exe *.out
