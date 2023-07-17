CC = g++

CFLAGS = -Wall
CFLAGS += -W
CFLAGS += -Werror
CFLAGS += -Winline
CFLAGS += -pedantic
CFLAGS += -std=c++17
CFLAGS += -ffast-math
CFLAGS += -march=native
CFLAGS += -DHAVE_INLINE
#CFLAGS += -Wno-unused-parameter

LIBS = -lgsl
LIBS += -lgslcblas
LIBS += -lstdc++

OBJDIR = ./bin/obj/
BINDIR = ./bin/
SRCDIR = ./src/
INCDIR = ./include/

_OBJ = Atom_info n_lowest_energy RNAData steric_clash_check FragmentAssembly CMB_Manager Combinatorial_Addition DimerLib WatsonCrickPair InputHandler Kabsch OutputString RNA_Math StructureBuilders HopcroftKarp Main Globals
DBUGOBJ = $(patsubst %, %_debug.o,$(_OBJ))
PRFOBJ = $(patsubst %, %_profiled.o,$(_OBJ))
PRODOBJ = $(patsubst %, %_production.o,$(_OBJ))
SRC = $(patsubst %,$(SRCDIR)%.cpp,$(_OBJ))

vpath %.cpp $(SRCDIR)

production: $(BINDIR)TAPERSS
debug: 		$(BINDIR)TAPERSS_debug
profile: 	$(BINDIR)TAPERSS_profile 

$(DBUGOBJ): 
	$(CC) $(CFLAGS) -g -o $(OBJDIR)$@ -c $(SRCDIR)$(subst _debug.o,.cpp,$@) -I $(INCDIR)

$(PRODOBJ): 
	$(CC) $(CFLAGS) -O3 -DGSL_RANGE_CHECK_OFF -o $(OBJDIR)$@ -c $(SRCDIR)$(subst _production.o,.cpp,$@) -I $(INCDIR)

$(PRFOBJ): 
	$(CC) $(CFLAGS) -o $(OBJDIR)$@ -c $(SRCDIR)$(subst _profiled.o,.cpp,$@) -I $(INCDIR) -pg

$(BINDIR)TAPERSS_debug: $(DBUGOBJ)
	$(CC) $(CFLAGS) -g -o $@ $(patsubst %, $(OBJDIR)%, $(DBUGOBJ)) $(LIBS) -I $(INCDIR)

$(BINDIR)TAPERSS_profile: $(PRFOBJ)
	$(CC) $(CFLAGS) -o $@ $(patsubst %, $(OBJDIR)%, $(PRFOBJ)) $(LIBS) -I $(INCDIR) -pg

$(BINDIR)TAPERSS: $(PRODOBJ)
	$(CC) $(CFLAGS) -DGSL_RANGE_CHECK_OFF -O3 -o $@ $(patsubst %, $(OBJDIR)%, $(PRODOBJ)) $(LIBS) -I $(INCDIR)

clean:
	rm -f $(OBJDIR)*.o
