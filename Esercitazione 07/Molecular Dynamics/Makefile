CC = g++
CFLAGS = -Wall -O3 --std=c++11

MolDyn_NVE.exe : MolDyn_NVE.o blocking.o functions.o
	$(CC) functions.o blocking.o MolDyn_NVE.o -o MolDyn_NVE.exe 

MolDyn_NVE.o : MolDyn_NVE.cpp MolDyn_NVE.h
	$(CC) -c MolDyn_NVE.cpp -o MolDyn_NVE.o $(CFLAGS)
blocking.o : blocking.cpp blocking.h
	$(CC) -c blocking.cpp -o blocking.o $(CFLAGS)
functions.o : functions.cpp functions.h
	$(CC) -c functions.cpp -o functions.o $(CFLAGS)

clean :
	rm *.o *.dat MolDyn_NVE.exe 
