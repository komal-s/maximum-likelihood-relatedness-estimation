CC	= g++
CFLAGS = -Wall
LIBARG	= -g -std=c++11 
TARGET	= lcMLkin
SRC	= $(addsuffix .cpp, $(TARGET))
CFLAGS=-c -Wall
all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(SRC) $(LIBARG) -o $@

clean:
	rm -f $(TARGET)
