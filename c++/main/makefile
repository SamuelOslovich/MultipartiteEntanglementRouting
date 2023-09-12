################################################################
ROOTLDFLAGS := $(shell root-config --libs --ldflags)
ROOTCFLAGS := $(shell root-config --cflags)
ROOTGLIBS := $(shell root-config --glibs --libs)

MYLIBS = -I. -I../include/

IGRAPH_LIB_DIR = /usr/local/lib
IGRAPH_INCLUDE_DIR = /usr/local/include/igraph
INC_IGRAPH = -I$(IGRAPH_INCLUDE_DIR) -L$(IGRAPH_LIB_DIR) -ligraph

################################################################
PLATFORM = $(shell root-config --platform)
$(info PLATFORM is $(PLATFORM)) 

################################################################
CFLAGS = $(ROOTCFLAGS) $(MYLIBS)
LDFLAGS = $(ROOTLDFLAGS) 
GCCFLAGS = -g -Wall -fPIC -std=c++11
CXX = g++
RM = rm -f


################################################################
SOURCES =$(wildcard *.cpp ../src/*.cpp)
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = main
TESTES = testes
CLEANOBJS = $(wildcard ../src/*.o)

ifeq ($(PLATFORM),macosx)
	EXTRA_CMD = install_name_tool -add_rpath $(shell root-config --libdir) $@
endif

################################################################
all: $(EXECUTABLE) 
run: $(EXECUTABLE) exe

exe:
		@echo "\nrun" 
		@./$(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
		@echo "\nMaking executable " $@ " ..."
		@$(CXX) -o $@ $(GCCFLAGS) $(CFLAGS) $(LDFLAGS) $^
		@echo "\ndone"

.cpp.o:
	@echo  "\nCompiling " $@ 
	$(CXX) -c $(GCCFLAGS) $(CFLAGS) $^ -o $@


clean:          
		@echo "\nclean" 
		@rm -fv $(CLEANOBJS)

distclean:          
		@echo "\nclean" 
		@rm -fv $(EXECUTABLE) $(CLEANOBJS)


