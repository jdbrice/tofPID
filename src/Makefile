
# project name and executable name
project = pid
# source suffix
source = .cpp 
# object files to make
objects = pid.o pidHistogramMaker.o TOFrPicoDst.o SimultaneousFit.o\
#			reporter.o 

sources = $(wildcard *.cpp)
headers = $(wildcard ../include/*.h)

# ROOT libs and includes
ROOTCFLAGS    	= $(shell root-config --cflags)
ROOTLDFLAGS    	= $(shell root-config --ldflags)
ROOTLIBS      	= $(shell root-config --libs)
ROOTGLIBS     	= $(shell root-config --glibs)

# Custom Project includes and libs
includes		= -I../include -I. -I$(ROOTSYS)/include -I$(ROOTDEV)/include 
libraries 		= 	$(ROOTLIBS) -lGui -lHtml -lGuiHtml -lRooFitCore -lRooFit $(CERNLIBS) $(SYSLIBS)
staticLibs 		= jdbREF.a


cxx 		= g++ -std=c++0x
flags 		= -Wall -g $(ROOTCFLAGS) $(includes) -Wno-write-strings 
compile 	= $(cxx) $(flags) -c 
ldFlags  	= $(ROOTLDFLAGS) -g

all: $(project)

$(project): $(objects)
		$(cxx) $(ldFlags) -o  $(project) $(objects) $(staticLibs) $(ROOTGLIBS) $(libraries)
		@cp $(project) ../bin

%.o: %$(source)
		$(compile)  $<

clean:
		@rm -f $(objects) $(project)
		@rm -f ../bin/$(project)
		@echo "clean done"

full:
	clean
	$(project)

#doc:
#	@echo "Generating documentation..." 	\					
#	cldoc generate $(CXXFLAGS) --			\
#		--report 							\
#		--merge docs 						\
#		--output html 						\
#		$(sources) $(headers)
