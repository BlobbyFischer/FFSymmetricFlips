CXX = g++

CXXFLAGS = -std=c++17 -O3 -march=native -Wall -Wextra

all: flip extend

flip: main.cpp
	$(CXX) $(CXXFLAGS) -o flip main.cpp

extend: extend.cpp
	$(CXX) $(CXXFLAGS) -o extend extend.cpp

clean:
	rm -f flip extend