# Sample typical Makefile for CSE 597 
#
# Adam Lavely
# Fall 2018

#-#-#-#-#
# Compiler
CC=g++
# Can be any compilers in your path
# Know what modules you have loaded (comment!)
# $ module load gcc/5.3.1

# Flags
CPPFLAGS=-c  -g -std=c++11 -fopenmp
#CFLAGS=-c -Wall -O2
#-#-#-#-#

cg: cg_old.o
	$(CC) -std=c++11 -fopenmp cg_old.o -o cg

cg_old.o: cg_old.cpp
	$(CC) $(CPPFLAGS) cg_old.cpp

cg_omp: cg_omp.o
	$(CC) -std=c++11 -fopenmp cg_omp.o -o cg_omp

cg_omp.o: cg_omp.cpp
	$(CC) $(CPPFLAGS) cg_omp.cpp

cg_omp2: cg_omp2.o
	$(CC) -std=c++11 -fopenmp cg_omp2.o -o cg_omp2

cg_omp2.o: cg_omp2.cpp
	$(CC) $(CPPFLAGS) cg_omp2.cpp

cg_omp3: cg_omp3.o
	$(CC) -std=c++11 -fopenmp cg_omp3.o -o cg_omp3

cg_omp3.o: cg_omp3.cpp
	$(CC) $(CPPFLAGS) cg_omp3.cpp

clean:
	rm *o cg cg_omp

