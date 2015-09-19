all: Judge.cpp 	Strategy.cpp
	g++ Judge.cpp Strategy.cpp -m32 -shared -o Strategy.dll