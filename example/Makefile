CC		= g++
FLAG	= -std=c++11 -fext-numeric-literals -Wall -O2
INC   = -I/usr/local/include -I"../include"
LIBS	=  -L"../lib" -L/usr/local/lib64 -lopm -lwcs -lm
LIBSLI= -lsfitsio -lsllib -lreadline -lcurses -lbz2 -lz

all: small medium wcstest

%:%.cc ../include ../lib
	$(CC) $(FLAG) $(INC) -o $@ $< $(LIBS)
