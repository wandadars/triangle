SRC1 = tri2vtk.c
TARGET1 = tri2vtk
HOME = /home
BIN = bin
CC = gcc
FLAGS = -lm

$(TARGET1) : $(SRC1)
	$(CC) $(FLAGS) $(SRC1) -o $(TARGET1)

install :
	-mv $(TARGET1) $(HOME)/$(BIN)

all : $(TARGET1) 

clean :
	-rm $(TARGET1)  *.o core
 
