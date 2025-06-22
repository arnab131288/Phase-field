#Makefile for Phase-Field Simulation

# Makefile

# Compiler and flags
CC       := gcc
CFLAGS   := -Wall -Wextra -O2 -g
LDFLAGS  := -lm
INCLUDES := -I.

# Target executable
TARGET   := src/simulation

# Source files
SRCS     := \
	src/main.c \
	src/write_infile.c \
    src/fill_domain.c \
    src/respawn.c \
	src/memory_allocation.c \
    src/boundary.c \
	src/free_energy.c \
	src/gradient.c \
	src/anisotropy.c \
    src/phasefield.c \
	src/temperature.c \
    src/write_output.c \
	src/read_infile.c 
# Automatically compute object files from SRCS
OBJS     := $(SRCS:.c=.o)

.PHONY: all clean

all: $(TARGET)

# Link step
$(TARGET): $(OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

# Compile step
%.o: %.c 
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<

# Clean up artifacts
clean:
	@rm -f $(OBJS) $(TARGET)

# Optional: rebuild everything
rebuild: clean all
