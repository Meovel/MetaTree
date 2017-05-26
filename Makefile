#
#	Makefile for Meta-Tree Labeling
#	- Hao Zhang
#
#	implements:
#		all (default), query, benchmark, output_metatree, clean
#

CC=g++
FLAGS= -std=c++11 -Wall -O3 
SRC = src
MAIN = query
TEST = benchmark
STORE = output_metatree


all: $(MAIN) $(TEST) $(STORE)

# COMPILE
$(MAIN): src/$(MAIN)_main.cpp 
	$(CC) $(FLAGS) -o $(MAIN) $(SRC)/$(MAIN)_main.cpp 

$(TEST): src/$(TEST)_main.cpp 
	$(CC) $(FLAGS) -o $(TEST) $(SRC)/$(TEST)_main.cpp 

$(STORE): src/$(STORE)_main.cpp 
	$(CC) $(FLAGS) -o $(STORE) $(SRC)/$(STORE)_main.cpp 


clean:
	rm -f *.o  $(MAIN)  $(MAIN).exe  $(TEST)  $(TEST).exe  $(STORE)  $(STORE).exe