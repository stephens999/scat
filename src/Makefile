SHELL   = sh
CXX     = g++
#CFLAGS = -g -std=gnu++11 $(INCLUDE)
#CFLAGS = -g -std=gnu++11 -D_GLIBCXX_DEBUG $(INCLUDE)
#CFLAGS = -O3 -DNDEBUG $(INCLUDE)
CFLAGS = -O3 -DNDEBUG -std=gnu++11 $(INCLUDE)

LDFLAGS = -llapack
OBJ2 = scat3.o utility.o mt19937ar.o readboundary.o

SCAT3 : $(OBJ2) scat3.hpp utility.hpp readboundary.hpp
	$(CXX) $(CFLAGS) -o $@ $^ $(LDFLAGS)

.SUFFIXES: .cpp .o

.cpp.o:
	$(CXX) $(CFLAGS) -c $<

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
