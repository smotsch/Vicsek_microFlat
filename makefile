############################
##   Makefile SppFlat2D   ##
############################

##----------------------------------------##
##         main files to compile          ##
##----------------------------------------##
## objects
OBJECTS = toolkit.o             \
  input_output_SppFlat2D.o      \
  boundary_SppFlat2D.o          \
  initial_condition_SppFlat2D.o \
  grid_interaction.o	        \
  interaction_SppFlat2D.o       \
  stat_SppFlat2D.o	        \
  main_SppFlat2D.o
## compiler
CC = ifort
## options for the compiler (-g for debugger, parano: -check all -warn all), 
FLAGS = -assume bscc -msse2 -O2
## name of the executable
EXEC = SppFlat2D


##----------------------------------------##
##            jobs for 'make'             ##
##----------------------------------------##

## Tricks to use make with several folders
##----------------------------------------
# Source tree:
#    Makefile
#    -- src/   *.f90
#    -- build/ *.o *.mod
#    -- bin/   $(EXE)
# 
# Paths
OBJDIR = build
vpath %.f90 src
vpath %.o $(OBJDIR)
# to help make find the object in the folder $(OBJDIR)
OBJECTS2 = $(addprefix $(OBJDIR)/, $(OBJECTS))

## command
all: init_folder $(EXEC) move init_parameters

$(EXEC): $(OBJECTS)
	$(CC) $(FLAGS) -o $(EXEC) $(OBJECTS2)

%.o : %.f90
	$(CC) -c $< -o $(OBJDIR)/$@ -module $(OBJDIR)

.PHONY: init_folder init_parameters move clean
init_folder:
	mkdir -p $(OBJDIR) data images videos
init_parameters:
	@if [ ! -f 'bin/PARAMETER_init.txt' ]; then \
	  cp bin/parameters/PARAMETER_init_bak.txt bin/PARAMETER_init.txt; \
	fi
	@if [ ! -f 'bin/PARAMETER_SppFlat2D.txt' ]; then \
	  cp bin/parameters/PARAMETER_SppFlat2D_bak.txt bin/PARAMETER_SppFlat2D.txt; \
	fi
move:
	mv $(EXEC) bin
clean:
	rm -rf $(OBJDIR)/*.mod
	rm -rf $(OBJDIR)/*.o


##----------------------------------------##
##      remark: changes for gfortan       ##
##----------------------------------------##
# 1. CC = gfortran
# 2. FLAGS2 = -Wall -g -fbackslash -O2
# 3. $(OBJECTS): %.f90
# 	$(CC) $(CFLAGS) -c -o $@ $< -J$(OBJDIR)

