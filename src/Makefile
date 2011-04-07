##############################################################
# photomosaic Makefile                                       #
# ========================================================== #
# Available targets:                                         #
#    photomosaic                                             #
#    photomosaic-debug                                       #
#    all                                                     #
#    profile                                                 #
#    clean                                                   #
# ========================================================== #
#    Copyright (c) 2011, Jeffrey P. Gill                     #
#                                                            #
#    This file is part of photomosaic.                       #
#                                                            #
#    photomosaic is free software: you can redistribute it   #
#    and/or modify it under the terms of the GNU General     #
#    Public License as published by the Free Software        #
#    Foundation, either version 3 of the License, or (at     #
#    your option) any later version.                         #
#                                                            #
#    photomosaic is distributed in the hope that it will be  #
#    useful, but WITHOUT ANY WARRANTY; without even the      #
#    implied warranty of MERCHANTABILITY or FITNESS FOR A    #
#    PARTICULAR PURPOSE.  See the GNU General Public License #
#    for more details.                                       #
#                                                            #
#    You should have received a copy of the GNU General      #
#    Public License along with photomosaic.  If not, see     #
#    <http://www.gnu.org/licenses/>.                         #
##############################################################

# List of executables that can be built
APPS = photomosaic \
		 photomosaic-debug


# When a target is not specified, the default executable is
# built
default: photomosaic


# Specifying the 'all' target will build every executable
all: $(APPS)


# Name of the main build directory
BUILD = build


# Target for running tests and analyzing profiling data to
# assist with optimization
.PHONY: profile
profile: photomosaic-debug
	time -p ./$< > ~/Desktop/out
	gprof $< > gprof.out
	less gprof.out


# Target for removing all files and directories created by
# 'make'
.PHONY: clean
clean:
	rm -rf $(APPS) $(BUILD)/ *.out *~



ifndef OBJDIR
##############################################################
# FIRST PASS:                                                #
# When 'make' is run from the command line, OBJDIR should    #
# not be defined, and the makefile will determine where      #
# object files should be placed before running 'make' again  #
# with OBJDIR defined                                        #
##############################################################

# Set OBJDIR to a unique path for the executable that is
# being compiled, create the build directory if necessary,
# and run 'make' again with OBJDIR defined
$(APPS): OBJDIR=$(BUILD)/$@
$(APPS): FORCE
	-mkdir -p $(OBJDIR)
	make $@ OBJDIR=$(OBJDIR)

FORCE:



else
##############################################################
# SECOND PASS:                                               #
# OBJDIR will be defined when 'make' is called a second time #
# (from the makefile during the first pass) and compilation  #
# can begin                                                  #
##############################################################

# List source code files used
HEADERS = antipole.h
SOURCES = antipole.c \
			 main.c


# Create a list of object files that will be built and
# linked together to construct an execuatable
OBJECTS = $(addprefix $(OBJDIR)/, $(addsuffix .o, $(basename $(notdir $(SOURCES)))))


# Define various compilation flags used by all targets
DEFINES = 
FLAGS   = -pipe -Wall -W
LFLAGS  =
LIBS    = -lm
INCPATH = .


# Define target-specific compilation flags
photomosaic-debug: DEFINES+=DEBUG _GLIBCXX_DEBUG
photomosaic-debug: FLAGS+=-O0 -g -pg
photomosaic-debug: LFLAGS+=-Wl,-O0 -g -pg


# Specify the dependencies and build rules for the
# executables
photomosaic photomosaic-debug: $(OBJECTS)
	gcc $(LFLAGS) $(LIBS) -o $@ $^


# Define the generic rule for building .o object files from
# .c source files
$(OBJDIR)/%.o: %.c
	gcc -c $(FLAGS) $(addprefix -D, $(DEFINES)) $(addprefix -I, $(INCPATH)) -o $@ $<


# Specify dependencies for all object files
$(OBJDIR)/antipole.o: antipole.c \
	antipole.h

$(OBJDIR)/main.o: main.c \
	antipole.h

endif

# End

