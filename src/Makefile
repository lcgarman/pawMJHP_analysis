# Compiler and flags
CC = gcc
CFLAGS = -Wall -Wextra -g -I/share/apps/dfprograms/fftw3/include -I/home/lcgarman/gsl/include 
LDFLAGS = -Wl,-rpath,/home/lcgarman/gsl/lib -L/share/apps/dfprograms/fftw3/lib -L/home/lcgarman/gsl/lib
LIBS = -lgsl -lgslcblas -lm -lfftw3 -mcmodel=medium 

# Source files for prepare2theta function
SRCS_PREP = prep_mjhp2theta_main.c globals.c structures.c allocate_memory.c cell_parameters.c fftw_functions.c read_binary_v10.c paw_parameters.c hkl_functions.c rwa_functions.c  twotheta_functions.c recgrid_functions.c

# Source files for mjhp2theta function
SRCS_2THETA = mjhp2theta_main.c globals.c structures.c allocate_memory.c cell_parameters.c fftw_functions.c read_binary_v10.c paw_parameters.c rwa_functions.c twotheta_functions.c energy_functions.c extra_functions.c  mjhp_2theta_energy.c hkl_functions.c

# Source files for mjhpHKL function
SRCS_HKL = mjhpHKL_main.c globals.c structures.c allocate_memory.c cell_parameters.c fftw_functions.c paw_parameters.c extra_functions.c mjhp_HKL_energy.c energy_functions.c mjhp_functions.c hkl_functions.c rwa_functions.c read_binary_v10.c

#Source files for mjhpHKL_density function
SRCS_DEN = mjhpHKL_density_main.c globals.c structures.c allocate_memory.c cell_parameters.c fftw_functions.c mjhp_functions.c hkl_functions.c rwa_functions.c xsf_functions.c mjhp_density.c read_binary_v10.c

# Object files
OBJS_PREP = $(SRCS_PREP:.c=.o)
OBJS_2THETA = $(SRCS_2THETA:.c=.o)
OBJS_HKL = $(SRCS_HKL:.c=.o)
OBJS_DEN = $(SRCS_DEN:.c=.o)

# Executable name
TARGET_PREP = prep_paw2theta 
TARGET_2THETA = paw_mjhp2theta 
TARGET_HKL = paw_mjhpHKL 
TARGET_DEN = paw_mjhpHKL_density  

all: $(TARGET_2THETA) $(TARGET_HKL) $(TARGET_PREP)

# Linking rule to create the prep_paw2theta 
$(TARGET_PREP): $(OBJS_PREP)
	$(CC) $(CFLAGS) $(OBJS_PREP) $(LDFLAGS) $(LIBS) -o $(TARGET_PREP)
# Linking rule to create the paw_mjhp2theta 
$(TARGET_2THETA): $(OBJS_2THETA)
	$(CC) $(CFLAGS) $(OBJS_2THETA) $(LDFLAGS) $(LIBS) -o $(TARGET_2THETA)
# Linking rule to create the paw_mjhpHKL
$(TARGET_HKL): $(OBJS_HKL)
	$(CC) $(CFLAGS) $(OBJS_HKL) $(LDFLAGS) $(LIBS) -o $(TARGET_HKL)
# Linking rule to create the paw_mjhpHKL_density 
$(TARGET_DEN): $(OBJS_DEN)
	$(CC) $(CFLAGS) $(OBJS_DEN) $(LDFLAGS) $(LIBS) -o $(TARGET_DEN)

# Compilation rule for object files
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Clean rule
clean:
	rm -f $(OBJS_PREP) $(OBJS_2THETA) $(OBJS_HKL) $(OBJS_DEN)
