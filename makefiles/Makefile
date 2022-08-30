CC = /usr/bin/g++

CFLAGS = -Wall
CFLAGS += -W
CFLAGS += -Werror
CFLAGS += -pedantic
CFLAGS += -std=c++17 
CFLAGS += -ffast-math

LIBS = -lgsl
LIBS += -lgslcblas

OBJDIR = ./bin/obj/
BINDIR = ./bin/
SRCDIR = ./src/
INCDIR = ./include/

_OBJ = Atom_info Attaching_Functions steric_clash_check CMB_Manager Combinatorial_Addition DimerLib Hairpin InputHandler Kabsch output_string RNAData RNADataArrayInternalLoop RNA_Math StructureBuilders Main Globals
OBJ = $(patsubst %, %.o,$(_OBJ))
PRFOBJ = $(patsubst %, %_profiled.o,$(_OBJ))
SRC = $(patsubst %,$(SRCDIR)%.cpp,$(_OBJ))

vpath %.cpp $(SRCDIR)

debug: 		$(BINDIR)Michael_COMB_hbond_debug
profile: 	$(BINDIR)Michael_COMB_hbond_profile 
production: $(BINDIR)Michael_COMB_hbond

$(OBJ): 
	$(CC) $(CFLAGS) -g -o $(OBJDIR)$@ -c $(SRCDIR)$(subst .o,.cpp,$@) $(LIBS) -I $(INCDIR)

$(PRFOBJ): 
	$(CC) $(CFLAGS) -o $(OBJDIR)$@ -c $(SRCDIR)$(subst _profiled.o,.cpp,$@) $(LIBS) -I $(INCDIR) -pg

$(BINDIR)Michael_COMB_hbond_debug: $(OBJ)
	$(CC) $(CFLAGS) -g -o $@ $(patsubst %, $(OBJDIR)%, $(OBJ)) $(LIBS) -I $(INCDIR)

$(BINDIR)Michael_COMB_hbond_profile: $(PRFOBJ)
	$(CC) $(CFLAGS) -g -o $@ $(patsubst %, $(OBJDIR)%, $(PRFOBJ)) $(LIBS) -I $(INCDIR) -pg

$(BINDIR)Michael_COMB_hbond: $(OBJ)
	$(CC) $(CFLAGS) -O3 -o $@ $(patsubst %, $(OBJDIR)%, $(OBJ)) $(LIBS) -I $(INCDIR)

clean:
	rm $(OBJDIR)*.o
	rm $(BINDIR)Michael_COMB_profile $(BINDIR)Michael_COMB_debug $(BINDIR)Michael_COMB
