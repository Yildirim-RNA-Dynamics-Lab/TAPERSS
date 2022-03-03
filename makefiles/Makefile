CC = /usr/bin/g++

CFLAGS = -Wall
CFLAGS += -W
CFLAGS += -Werror
CFLAGS += -pedantic
CFLAGS += -std=c++11 

LIBS = -lgsl
LIBS += -lgslcblas

OBJDIR = ./bin/obj/
BINDIR = ./bin/
SRCDIR = ./src/
INCDIR = ./include/

_OBJ = Atom_info Attaching_Functions CMB_Manager Combinatorial_Addition DimerLib Hairpin InputHandler Kabsch output_string RNA_data RNA_Math StructureBuilders Main Globals
OBJ = $(patsubst %, %.o,$(_OBJ))
PRFOBJ = $(patsubst %, %_profiled.o,$(_OBJ))
SRC = $(patsubst %,$(SRCDIR)%.cpp,$(_OBJ))

vpath %.cpp $(SRCDIR)

debug: 		$(BINDIR)Ivan_RNACMB_debug
profile: 	$(BINDIR)Ivan_RNACMB_profile
production: $(BINDIR)Ivan_RNACMB

$(OBJ): 
	$(CC) $(CFLAGS) -g -o $(OBJDIR)$@ -c $(SRCDIR)$(subst .o,.cpp,$@) $(LIBS) -I $(INCDIR)

$(PRFOBJ): 
	$(CC) $(CFLAGS) -o $(OBJDIR)$@ -c $(SRCDIR)$(subst _profiled.o,.cpp,$@) $(LIBS) -I $(INCDIR) -pg

$(BINDIR)Ivan_RNACMB_debug: $(OBJ)
	$(CC) $(CFLAGS) -g -o $@ $(patsubst %, $(OBJDIR)%, $(OBJ)) $(LIBS) -I $(INCDIR)

$(BINDIR)Ivan_RNACMB_profile: $(PRFOBJ)
	$(CC) $(CFLAGS) -g -o $@ $(patsubst %, $(OBJDIR)%, $(PRFOBJ)) $(LIBS) -I $(INCDIR) -pg

$(BINDIR)Ivan_RNACMB: $(OBJ)
	$(CC) $(CFLAGS) -O3 -o $@ $(patsubst %, $(OBJDIR)%, $(OBJ)) $(LIBS) -I $(INCDIR)

clean:
	rm $(OBJDIR)*.o
	rm $(BINDIR)Ivan_RNACMB_profile $(BINDIR)Ivan_RNACMB_debug $(BINDIR)Ivan_RNACMB
