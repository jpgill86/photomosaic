# MAKEFILE

CC = gcc
CFLAGS = -I.
LIBS = -lm

TARGET = photomosaic
BUILD = build

SOURCES = main.c antipole.c
_OBJ = main.o antipole.o
OBJ = $(patsubst %,$(BUILD)/%,$(_OBJ))

$(TARGET): $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

$(BUILD)/%.o: %.c
	mkdir -p $(BUILD)
	$(CC) -c -o $@ $< $(CFLAGS)

.PHONY: clean
clean:
	rm -rf $(TARGET) $(BUILD)/ *.o *~


$(BUILD)/main.o: antipole.h
$(BUILD)/antipole.o: antipole.h

# End

