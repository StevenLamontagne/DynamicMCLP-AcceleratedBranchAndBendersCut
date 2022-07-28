# Exemple de makefile pour un programme cplex en C++
#CFLAGS=-Wall -O -std=c++11
CFLAGS=-Wall -O -std=c++11 -ggdb3
CPPFLAGS=-DNDEBUG -DIL_STD 
CXX=g++
OBJ=main.o Model_BB.o Model_BendersCordeau.o CoverageCallback.o Model_Multicut.o MulticutCallback.o Data.o 
CPLEXLIB=-I/home/ibm/cplex-studio/22.1/CPLEX_Studio/concert/include -I/home/ibm/cplex-studio/22.1/CPLEX_Studio/cplex/include 
CPLEXLIBFLAGS=-lilocplex -lconcert -lcplex -lpthread -lm -ldl

 
main: $(OBJ)
	$(CXX) $(CFLAGS) -o solver $(OBJ) $(CPLEXLIB) $(CPLEXLIBFLAGS) 
 
main.o: main.cpp Model_BB.h Model_BendersCordeau.h Model_Multicut.h Data.h
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c main.cpp $(CPLEXLIB) $(CPLEXLIBFLAGS)

Model_BB.o: Model_BB.cpp Data.h
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c Model_BB.cpp $(CPLEXLIB) $(CPLEXLIBFLAGS)

Model_BendersCordeau.o: Model_BendersCordeau.cpp CoverageCallback.h Data.h
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c Model_BendersCordeau.cpp $(CPLEXLIB) $(CPLEXLIBFLAGS)

CoverageCallback.o: CoverageCallback.cpp Data.h
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c CoverageCallback.cpp $(CPLEXLIB) $(CPLEXLIBFLAGS)

Model_Multicut.o: Model_Multicut.cpp MulticutCallback.h Data.h
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c Model_Multicut.cpp $(CPLEXLIB) $(CPLEXLIBFLAGS)

MulticutCallback.o: MulticutCallback.cpp Data.h
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c MulticutCallback.cpp $(CPLEXLIB) $(CPLEXLIBFLAGS)

Data.o: Data.cpp
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c Data.cpp $(CPLEXLIB) $(CPLEXLIBFLAGS) 


clean:
	rm -f solver main.o Model_BB.o Data.o Model_BendersCordeau.o CoverageCallback.o Model_Multicut.o MulticutCallback.o