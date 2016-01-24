SHELL   = sh
CXX     = g++
INCLUDE = -I. -I..
CFLAGS = -g $(INCLUDE)
LDFLAGS = -llapack

.SUFFIXES: .cpp .o

.cpp.o:
	$(CXX) $(CFLAGS) -c $<

TESTOBJ = test.o test2.o

OBJ = main.gaussianprocess.o utility.o
HDR = utility.hpp utility.cpp
OBJ2 = scat2.o utility.o mt19937ar.o

SCAT2 : $(OBJ2)
	$(CXX) $(CFLAGS) -o $@ $^ $(LDFLAGS)

test: $(TESTOBJ)
	$(CXX) $(CFLAGS) -o $@ $^ $(LDFLAGS)


scat2.o : scat2.cpp $(HDR)
mt19937ar.o: mt19937ar.c

TAGS: *.hpp *.cpp
	etags $^
tags: *.hpp *.cpp
	ctags $^

# Stripping off the symbols can make the executable smaller.
strip: SCAT
	strip $^

clean:
	rm -rf *.o *.out
ultraclean: clean
	rm -f SCAT
