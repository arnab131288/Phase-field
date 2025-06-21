#Makefile for Phase-Field Simulation

CXX = g++
CXXFLAGS = -std=c++11 -O2 -Wall
LDFLAGS = -lm

#List all source files explicitly

SRCS = \
	src/main.cpp \
	src/write_infile.cpp \
    src/fill_domain.cpp \
    src/respawn.cpp \
	src/memory_allocation.cpp \
    src/boundary.cpp \
	src/free_energy.cpp \
	src/gradient.cpp \
	src/anisotropy.cpp \
    src/phasefield.cpp \
	src/temperature.cpp \
    src/write_output.cpp \
	src/read_infile.cpp 

#Object files

OBJS = $(SRCS:.cpp=.o)

#Executable name

TARGET = src/simulation

.PHONY: all clean

all:$(TARGET)
$(TARGET):$(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS)

#Pattern rule: compile any .cpp to .o

%.o: %.cpp header.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)



