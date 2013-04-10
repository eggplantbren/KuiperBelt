OPTIONS = -O2 -Wall -Wextra -pedantic -ansi -DNDEBUG
LIBS = -ldnest3 -lgsl -lgslcblas -lboost_thread -lboost_system

default:
	g++ $(OPTIONS) -c *.cpp Hyperparameters/*.cpp
	g++ -o main *.o $(LIBS)
	rm -f *.o

