CXX = g++

CXXFLAGS = -std=c++17 -O3 -march=native -Wall -Wextra

all: flip reshape verify

flip: main.cpp
	$(CXX) $(CXXFLAGS) -o flip main.cpp

reshape: reshape.cpp
	$(CXX) $(CXXFLAGS) -o reshape reshape.cpp

verify: verify.cpp
	$(CXX) $(CXXFLAGS) -o verify verify.cpp

clean:
	rm -f flip extend