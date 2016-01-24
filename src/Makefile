SHELL   = sh
CXX     = g++
INCLUDE = -I. -I..
CFLAGS = -g $(INCLUDE)
LDFLAGS = -llapack
OBJ2 = scat2.o utility.o mt19937ar.o

.SUFFIXES: .cpp .o

.cpp.o:
	$(CXX) $(CFLAGS) -c $<

SCAT2 : $(OBJ2)
	$(CXX) $(CFLAGS) -o $@ $^ $(LDFLAGS)

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
