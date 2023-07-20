CC=gcc
CFLAGS=-std=c99 -O2

OBJECTS = mortar

all: $(OBJECTS)

mortar: main.c
	$(CC) $(CFLAGS) main.c paf.c bed.c -o mortar -lz -lm

.PHONY: clean
clean:
	-rm $(OBJECTS)
