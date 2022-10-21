SHELL   = sh
CXX     = g++
#CXXFLAGS = -g -std=gnu++11 $(INCLUDE)
#CXXFLAGS = -g -std=gnu++11 -D_GLIBCXX_DEBUG $(INCLUDE)
#CXXFLAGS = -O3 -DNDEBUG $(INCLUDE)
CXXFLAGS = -O3 -DNDEBUG -std=gnu++11 $(INCLUDE)

LDFLAGS = -L$(LAPACK_DIR) -llapack
OBJ2 = scat3.o utility.o mt19937ar.o readboundary.o

SCAT3 : $(OBJ2)
	$(CXX) $(CFLAGS) -o $@ $^ $(LDFLAGS)

.SUFFIXES: .hpp .cpp .o

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
