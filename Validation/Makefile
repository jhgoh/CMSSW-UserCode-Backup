.SUFFIXES : .cpp
CXX = g++
CFLAGS=`root-config --cflags --libs` -l boost_program_options

OBJECTS = lsRoot dumpRoot

all : $(OBJECTS)

lsRoot : lsRoot.cpp
	$(CXX) $(CFLAGS) $< -o $@ 

dumpRoot : dumpRoot.cpp
	$(CXX) $(CFLAGS) $< -o $@ 

clean :
	rm -f core $(OBJECTS)
