##############################################################################
HOME     = /Users/mdekauwe
CFLAGS   =  -g -Wall #-O3 #-g -Wall
ARCH     =  x86_64
INCLS    = -I./include
LIBS     = -lm
CC       =  gcc
PROGRAM  =  weather_generate_par
SOURCES  =  $(PROGRAM).c
OBJECTS = $(SOURCES:.c=.o)
RM       =  rm -f
##############################################################################

# top level create the program...
all: 		$(PROGRAM)

# Compile the src file...
$(OBJECTS):	$(SOURCES)
		$(CC) ${INCLS} $(CFLAGS) -c $(SOURCES)

# Linking the program...
$(PROGRAM):	$(OBJECTS)
		$(CC) $(OBJECTS) $(LIBS) ${INCLS} $(CFLAGS) -o $(PROGRAM)

# clean up...
clean:
		$(RM) $(OBJECTS) $(PROGRAM)

install:
		cp $(PROGRAM) $(HOME)/bin/$(ARCH)/.
		$(RM) $(OBJECTS)
##############################################################################
