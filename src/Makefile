SHELL   = sh
CXX     = g++
#CXXFLAGS = -g -std=gnu++11 $(INCLUDE)
#CXXFLAGS = -g -std=gnu++11 -D_GLIBCXX_DEBUG $(INCLUDE)
CFLAGS = -O3 -DNDEBUG $(INCLUDE)
CXXFLAGS = -O3 -DNDEBUG -std=gnu++11 $(INCLUDE)

LDFLAGS = -llapack
OBJ2 = scat3.o utility.o mt19937ar.o readboundary.o

.SUFFIXES: .hpp .cpp .o

SCAT3 : $(OBJ2)
	$(CXX) $(CFLAGS) -o $@ $^ $(LDFLAGS)

TAGS: *.hpp *.cpp
	etags $^

tags: *.hpp *.cpp
	ctags $^

# Stripping off the symbols can make the executable smaller.
strip: SCAT3
	strip $^

clean:
	rm -rf *.o *.out

ultraclean: clean
	rm -f SCAT3

mt19937ar.o: mt19937ar.hpp
readboundary.o: readboundary.hpp
utility.o: utility.hpp
scat3.o: scat3.hpp mt19937ar.hpp readboundary.hpp

