CXX = g++
# -O3 maximizes speed. -march=native allows hardware-specific bitwise ops.
CXXFLAGS = -std=c++17 -O3 -march=native -Wall -Wextra

TARGET = flip
SRCS = main.cpp

all: $(TARGET)

$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRCS)

clean:
	rm -f $(TARGET)