CC = /usr/bin/clang-8

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

_OBJ = Atom_info Attaching_Functions steric_clash_check CMB_Manager Combinatorial_Addition DimerLib Hairpin InputHandler Kabsch output_string RNAData RNADataArrayInternalLoop RNA_Math StructureBuilders HopcroftKarp Main Globals
DBUGOBJ = $(patsubst %, %_debug.o,$(_OBJ))
PRFOBJ = $(patsubst %, %_profiled.o,$(_OBJ))
PRODOBJ = $(patsubst %, %_production.o,$(_OBJ))
SRC = $(patsubst %,$(SRCDIR)%.cpp,$(_OBJ))

vpath %.cpp $(SRCDIR)

debug: 		$(BINDIR)Michael_COMB_hbond_debug
profile: 	$(BINDIR)Michael_COMB_hbond_profile 
production: $(BINDIR)Michael_COMB_hbond

$(DBUGOBJ): 
	$(CC) $(CFLAGS) -g -o $(OBJDIR)$@ -c $(SRCDIR)$(subst _debug.o,.cpp,$@) -I $(INCDIR)

$(PRODOBJ): 
	$(CC) $(CFLAGS) -O3 -DGSL_RANGE_CHECK_OFF -o $(OBJDIR)$@ -c $(SRCDIR)$(subst _production.o,.cpp,$@) -I $(INCDIR)

$(PRFOBJ): 
	$(CC) $(CFLAGS) -o $(OBJDIR)$@ -c $(SRCDIR)$(subst _profiled.o,.cpp,$@) -I $(INCDIR) -pg

$(BINDIR)Michael_COMB_hbond_debug: $(DBUGOBJ)
	$(CC) $(CFLAGS) -g -o $@ $(patsubst %, $(OBJDIR)%, $(DBUGOBJ)) $(LIBS) -I $(INCDIR)

$(BINDIR)Michael_COMB_hbond_profile: $(PRFOBJ)
	$(CC) $(CFLAGS) -o $@ $(patsubst %, $(OBJDIR)%, $(PRFOBJ)) $(LIBS) -I $(INCDIR) -pg

$(BINDIR)Michael_COMB_hbond: $(PRODOBJ)
	$(CC) $(CFLAGS) -DGSL_RANGE_CHECK_OFF -O3 -o $@ $(patsubst %, $(OBJDIR)%, $(PRODOBJ)) $(LIBS) -I $(INCDIR)

clean:
	rm $(OBJDIR)*.o
	rm $(BINDIR)Michael_COMB_hbond_profile $(BINDIR)Michael_COMB_hbond_debug $(BINDIR)Michael_COMB_hbond
