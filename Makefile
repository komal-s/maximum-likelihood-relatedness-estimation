CC	= g++
CFLAGS = -c -Wall -Werror -fmax-errors=3
LIBARG	= -g -std=c++11 -O3 -fopenmp
EIGEN = .
INC = -I $(EIGEN)
TARGET	= relatedness utils
SRC	= $(addsuffix .cpp, $(TARGET))

$(TARGET): $(SRC)
	$(CC) $(SRC) $(LIBARG) $(INC) -o $@

clean:
	rm -f $(TARGET)
