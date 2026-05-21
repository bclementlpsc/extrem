PYTHIAPATH = ../pythia8317

PYTHIALIB = $(PYTHIAPATH)/lib/libpythia8.a
CXX = g++
CXXFLAGS = -Wall -lm
INCS = -I$(PYTHIAPATH)/include
SRCFILE = $(wildcard *.cpp)
OBJFILE = $(patsubst %.cpp,%_x,$(SRCFILE))

default : all

all : $(OBJFILE)

%_x : %.cpp
	@echo "--> Building " $@
	@$(CXX) $(CXXFLAGS) $(INCS) -o $@ $< $(PYTHIALIB)

clean :
	@echo "--> Removing all executables"
	@rm -rf *_x
	