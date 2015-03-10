CC		= g++
AR		= ar
FLAG	= -std=c++11 -fext-numeric-literals -Wall -O2 -DDEBUG
INC   = -I/usr/local/include
LIBS	= -L/usr/local/lib64 -lwcs -lm
#LIBSLI= -lsfitsio -lsllib -lreadline -lcurses -lbz2 -lz
OBJ		= opm_generate.o opm_prematch.o opm_final.o

all: libopm.a

%.o:%.cc
	$(CC) $(FLAG) $(INC) -c $< $(LIBS)

libopm.a: $(OBJ) opm.h
	$(AR) rcs $@ $^
