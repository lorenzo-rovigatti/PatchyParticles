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

# --- PREPROCESSOR FLAGS --- #
CPPFLAGS = $(INCLUDE_DIRS)

# -- CC FLAGS   ------#
CFLAGS = $(OPTIMIZATION_FLAG) $(ARCH) $(DEF)

# --- LINKER FLAGS  -------#
LFLAGS = $(LIBRARY_DIRS) $(STATIC_FLAG)

# ---------  OBJECTS -----------#
OBJ = parse_input.o output.o system.o MC.o utils.o avb.o vmmc.o

# ---------  EXECUTABLE -----------#
EXE = PatchyParticles

# -------- COMPILATION STEPS ------#

default: PatchyParticles generator

all: PatchyParticles generator

generator: generator.o $(OJB)
	$(CC) generator.o $(OBJ) $(LIBRARY_DIRS) -lm -o generator

PatchyParticles: main.o $(OBJ)
	$(CC) main.o $(OBJ) $(LIBRARY_DIRS) -lm -o $(EXE)

# -------- SUFFIX RULES ------#
%.o: %.c defs.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@


clean:
	rm -f $(EXE) $(OBJ) main.o
