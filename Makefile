TARGET = kClique
SRC = kClique.c
CFLAGS = -Wall -Wextra -O2
LIBS = -lm

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC) $(LIBS)

clean:
	rm -f $(TARGET)
