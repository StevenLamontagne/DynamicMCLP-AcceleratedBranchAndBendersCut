# Exemple de makefile pour un programme cplex en C++
CFLAGS=-Wall -O
CPPFLAGS=-DNDEBUG -DIL_STD 
CXX=g++
OBJ=main.o Model_BB.o Data.o 
CPLEXLIB=-I/home/ibm/cplex-studio/12.10.0.0/CPLEX_Studio/concert/include -I/home/ibm/cplex-studio/12.10.0.0/CPLEX_Studio/cplex/include 
JSONLIB=-I/local_1/outer/lamste/Librairies/JSON/include
CPLEXLIBFLAGS=-lilocplex -lconcert -lcplex -lpthread -lm -ldl
JSONLIBFLAGS=-lnlohmann
 
main: $(OBJ)
	$(CXX) $(CFLAGS) -o solver $(OBJ) $(CPLEXLIB) $(JSONLIB) $(CPLEXLIBFLAGS) $(JSONLIBFLAG)
 
main.o: main.cpp Model_BB.h Data.h
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c main.cpp $(CPLEXLIB) $(JSONLIB) $(CPLEXLIBFLAGS) $(JSONLIBFLAG)

Model_BB.o: Model_BB.cpp Data.h
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c Model_BB.cpp $(CPLEXLIB) $(JSONLIB) $(CPLEXLIBFLAGS) $(JSONLIBFLAG)

Data.o: Data.cpp
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c Data.cpp $(CPLEXLIB) $(JSONLIB) $(CPLEXLIBFLAGS) $(JSONLIBFLAG)


clean:
	rm -f solver main.o Model_BB.o Data.o