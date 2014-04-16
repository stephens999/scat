SHELL   = sh
CXX     = g++
INCLUDE = -I. -I..

#CFLAGS  = -g -O4 -Wall -W -Wno-sign-compare -pedantic -DDEBUG=5   $(INCLUDE)
#CFLAGS  = -g -O4 -Wall -W -Wno-sign-compare -pedantic -DDEBUG=5 -DCHECK_BOUNDARY  $(INCLUDE)
#CFLAGS  = -g -O4 -Wall -W -Wno-sign-compare -pedantic -DNDEBUG  $(INCLUDE)
#CFLAGS = -O2 -DNDEBUG -DTREAT_SNPS_AS_MS $(INCLUDE)
#CFLAGS = -O2 -DBIGDATASETS $(INCLUDE)
#CFLAGS = -g $(INCLUDE)
#CFLAGS = -O2  $(INCLUDE)
CFLAGS = -g $(INCLUDE)

#LDFLAGS = -llapack -lg2c -lm
#LDFLAGS = -static /usr/lib/libm.a /usr/lib/gcc-lib/i386-redhat-linux/3.2.3/libg2c.a /lib/ssa/gcc-lib/i386-redhat-linux-gnu/3.5-tree-ssa/libg2c.a /usr/lib/liblapack.a
#LDFLAGS = /usr/lib/liblapack-2.a /usr/lib/libg2c.a /usr/lib/libm.a /usr/lib/libblas-2.a
LDFLAGS = -llapack
#LDFLAGS = -llapack /usr/lib/gcc-lib/i386-redhat-linux/3.2.3/libg2c.a /usr/lib/libm.a

SUFFIXES: 
.SUFFIXES: .cpp .o

.cpp.o:
	$(CXX) $(CFLAGS) -c $<

TESTOBJ = test.o test2.o

OBJ = main.gaussianprocess.o utility.o
HDR = utility.hpp utility.cpp
OBJ2 = scat2.o utility.o

SCAT : $(OBJ)
	$(CXX) $(CFLAGS) -o $@ $^ $(LDFLAGS)

SCAT2 : $(OBJ2)
	$(CXX) $(CFLAGS) -o $@ $^ $(LDFLAGS)

test: $(TESTOBJ)
	$(CXX) $(CFLAGS) -o $@ $^ $(LDFLAGS)


main.gaussianprocess.o   : main.gaussianprocess.cpp $(HDR)
scat2.o : scat2.cpp $(HDR)


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
