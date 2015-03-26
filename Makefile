CC	= g++
CFLAGS = -c -Wall
LIBARG	= -g -std=c++11 
EIGEN = .
INC = -I $(EIGEN)
TARGET	= lcMLkin
SRC	= $(addsuffix .cpp, $(TARGET))

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(SRC) $(LIBARG) $(INC) -o $@

clean:
	rm -f $(TARGET)
