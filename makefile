INCLUDE_DIRS =
LIBRARY_DIRS =
STATIC_FLAG = 

OPTIMIZATION_FLAG = -O3 -ffast-math
ARCH = -DNDEBUG

ifeq ($(dbg), 1)
	OPTIMIZATION_FLAG = -O0 -g2 -pg -ggdb3
	ARCH = -DDEBUG
	STATIC_FLAG = -pg
endif

ifeq ($(g), 1)
	OPTIMIZATION_FLAG += -g2 -pg -ggdb3
	ARCH = -DDEBUG
	STATIC_FLAG = -pg
endif

# ----- COMPILER -----#
CC = gcc -Wall -Wshadow

# -- CC FLAGS   ------#
CFLAGS = $(OPTIMIZATION_FLAG) $(ARCH) $(DEF) $(INCLUDE_DIRS)

# --- LINKER FLAGS  -------#
LFLAGS = $(LIBRARY_DIRS) $(STATIC_FLAG)

# ---------  OBJECTS -----------#
OBJ = parse_input.o output.o system.o MC.o utils.o avb.o vmmc.o cells.o

# ---------  EXECUTABLE -----------#
EXE = PatchyParticles

# -------- COMPILATION STEPS ------#

default: PatchyParticles generator

all: PatchyParticles generator

generator: generator.o $(OBJ)
	$(CC) generator.o $(OBJ) $(LIBRARY_DIRS) -lm -o generator

PatchyParticles: main.o $(OBJ)
	$(CC) main.o $(OBJ) $(LIBRARY_DIRS) -lm -o $(EXE)

# -------- SUFFIX RULES ------#
%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<
	
## --- Add dependencies by hand to make the Makefile as generic as possible --- #
avb.o: avb.h MC.h defs.h output.h parse_input.h system.h utils.h
cells.o: cells.h defs.h output.h
generator.o: defs.h output.h utils.h
main.o: defs.h MC.h output.h parse_input.h system.h utils.h
MC.o: avb.h defs.h MC.h output.h utils.h vmmc.h 
output.o: defs.h output.h MC.h parse_input.h utils.h
utils.o: defs.h
vmmc.o: defs.h MC.h defs.h output.h parse_input.h system.h utils.h

clean:
	rm -f $(EXE) $(OBJ) main.o generator generator.o
