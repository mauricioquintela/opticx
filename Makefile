# -----------------------------------------------------------------
# opticx Makefile
# -----------------------------------------------------------------
# Author    : J. J. Esteve-Paredes (JJEP)
# Modified  : D. Hernangómez-Pérez (DH)
# Version   : 1.1
# Date      : 20.06.2025
# -----------------------------------------------------------------

# -----------------------------------------------------------------
# Compiler and flags
# -----------------------------------------------------------------
FC     = gfortran
FFLAGS = -fopt-info-vec -g -fcheck=all -O3
#FFLAGS = -O2-Wall  -g -fcheck=all 

# -----------------------------------------------------------------
# Optional: use MKL  
# Usage: make USE_MKL=1 
# Note: only tested with Ubuntu, install with sudo apt-get install intel-mkl
# -----------------------------------------------------------------
ifeq ($(USE_MKL),1)
LIBS = -lmkl_rt -fopenmp -lpthread -lm -ldl	
# MKLROOT ?= /opt/intel/oneapi/mkl/latest
# LIBS = -L$(MKLROOT)/lib/intel64 \
#       -Wl,--start-group \
#       -lmkl_gf_lp64 -lmkl_core -lmkl_gnu_thread \
#       -Wl,--end-group -fopenmp -lpthread -lm -ldl
else
LIBS   = -lopenblas -fopenmp -lgfortran 
endif

# -----------------------------------------------------------------
# Directories
# -----------------------------------------------------------------
MAINDIR  = main
SRCDIR   = src
BINDIR   = bin
BUILDDIR = build

# -----------------------------------------------------------------
# Files
# -----------------------------------------------------------------
SRC_MODULES = $(wildcard $(SRCDIR)/*.f90)
OBJ_MODULES = $(SRC_MODULES:$(SRCDIR)/%.f90=$(BUILDDIR)/%.o)
SRC_MAIN    = $(MAINDIR)/opticx.f90
OBJ_MAIN    = $(BINDIR)/opticx.o
TARGET      = $(BINDIR)/opticx

# -----------------------------------------------------------------
# Build target
# -----------------------------------------------------------------
all: $(TARGET)
	rm -f $(OBJ_MAIN)

# -----------------------------------------------------------------
# Directory creation rules
# -----------------------------------------------------------------
$(BUILDDIR):
	mkdir -p $(BUILDDIR)

$(BINDIR):
	mkdir -p $(BINDIR)

# -----------------------------------------------------------------
# Linking rules
# -----------------------------------------------------------------
# Use the -J flag to specify the directory for .mod files
$(BUILDDIR)/%.o: $(SRCDIR)/%.f90 | $(BUILDDIR)
	$(FC) -J$(BUILDDIR) -c $< $(FFLAGS) -o $@ $(LIBS)

# -----------------------------------------------------------------------------
# Compilation rules
# -----------------------------------------------------------------------------
# Rule for creating the executable
$(TARGET): $(OBJ_MODULES) $(OBJ_MAIN)
	$(FC) $(FFLAGS) $(OBJ_MODULES) $(OBJ_MAIN) -o $(TARGET) $(LIBS)

$(BINDIR)/opticx.o: $(SRC_MAIN) | $(BINDIR) $(BUILDDIR)
	$(FC) -J$(BUILDDIR) -c $< $(FFLAGS) -o $@

# -----------------------------------------------------------------
# Module dependencies
# -----------------------------------------------------------------
# DH: Add more dependencies as needed for other modules
#     $(BUILDDIR)/some_other_module.o: \
# 	  $(BUILDDIR)/dependency_module.o

$(BUILDDIR)/parser_wannier90_tb.o: \
	$(BUILDDIR)/parser_input_file.o 

$(BUILDDIR)/parser_optics_xatu_dim.o: \
	$(BUILDDIR)/constants_math.o \
	$(BUILDDIR)/parser_wannier90_tb.o \
	$(BUILDDIR)/parser_input_file.o 

$(BUILDDIR)/exciton_envelopes.o: \
    $(BUILDDIR)/constants_math.o \
	$(BUILDDIR)/parser_optics_xatu_dim.o

$(BUILDDIR)/bands.o: \
	$(BUILDDIR)/constants_math.o \
	$(BUILDDIR)/parser_wannier90_tb.o \
	$(BUILDDIR)/parser_optics_xatu_dim.o

$(BUILDDIR)/ome_sp.o: \
	$(BUILDDIR)/constants_math.o \
	$(BUILDDIR)/parser_wannier90_tb.o \
	$(BUILDDIR)/parser_optics_xatu_dim.o  

$(BUILDDIR)/ome_ex.o: \
	$(BUILDDIR)/constants_math.o \
	$(BUILDDIR)/parser_wannier90_tb.o \
	$(BUILDDIR)/parser_optics_xatu_dim.o \
	$(BUILDDIR)/exciton_envelopes.o 

$(BUILDDIR)/ome.o: \
	$(BUILDDIR)/parser_input_file.o \
	$(BUILDDIR)/ome_sp.o \
	$(BUILDDIR)/ome_ex.o

$(BUILDDIR)/optical_response.o: \
	$(BUILDDIR)/parser_input_file.o \
	$(BUILDDIR)/sigma_first_sp.o \
	$(BUILDDIR)/sigma_first_ex.o \
	$(BUILDDIR)/sigma_second_sp.o \
	$(BUILDDIR)/sigma_second_ex.o

$(BUILDDIR)/sigma_first_sp.o: \
	$(BUILDDIR)/constants_math.o \
	$(BUILDDIR)/parser_input_file.o \
	$(BUILDDIR)/parser_optics_xatu_dim.o \
	$(BUILDDIR)/ome.o 

$(BUILDDIR)/sigma_first_ex.o: \
	$(BUILDDIR)/constants_math.o 

$(BUILDDIR)/sigma_second_sp.o: \
	$(BUILDDIR)/constants_math.o \
	$(BUILDDIR)/parser_input_file.o \
	$(BUILDDIR)/parser_optics_xatu_dim.o \
	$(BUILDDIR)/ome.o 

$(BUILDDIR)/sigma_second_ex.o: \
	$(BUILDDIR)/constants_math.o \
	$(BUILDDIR)/ome_ex.o \
	$(BUILDDIR)/sigma_second_sp.o 

# -----------------------------------------------------------------
# Test: shift-kernel equivalence check (fast array path vs. reference
# scalar path). Reuses every module already built for the main target.
# -----------------------------------------------------------------
TESTDIR     = tests
SRC_TEST    = $(TESTDIR)/test_shift_kernel_equivalence.f90
OBJ_TEST    = $(BINDIR)/test_shift_kernel_equivalence.o
TARGET_TEST = $(BINDIR)/test_shift_kernel_equivalence

test: $(TARGET_TEST)
	rm -f $(OBJ_TEST)

$(TARGET_TEST): $(OBJ_MODULES) $(OBJ_TEST)
	$(FC) $(FFLAGS) $(OBJ_MODULES) $(OBJ_TEST) -o $(TARGET_TEST) $(LIBS)

$(OBJ_TEST): $(SRC_TEST) | $(BINDIR) $(BUILDDIR)
	$(FC) -I$(BUILDDIR) -J$(BUILDDIR) -c $< $(FFLAGS) -o $@ $(LIBS)

test_matrix: $(BINDIR)/test_shift_intens_ex_matrix

$(BINDIR)/test_shift_intens_ex_matrix: $(OBJ_MODULES) $(BINDIR)/test_shift_intens_ex_matrix.o
	$(FC) $(FFLAGS) $(OBJ_MODULES) $(BINDIR)/test_shift_intens_ex_matrix.o -o $@ $(LIBS)

$(BINDIR)/test_shift_intens_ex_matrix.o: tests/test_shift_intens_ex_matrix.f90 | $(BUILDDIR) $(BINDIR)
	$(FC) -I$(BUILDDIR) -J$(BUILDDIR) -c $< $(FFLAGS) -o $@ $(LIBS)

test_shift_real_data: $(BINDIR)/test_shift_real_data

$(BINDIR)/test_shift_real_data: $(OBJ_MODULES) $(BINDIR)/test_shift_real_data.o
	$(FC) $(FFLAGS) $(OBJ_MODULES) $(BINDIR)/test_shift_real_data.o -o $@ $(LIBS)

$(BINDIR)/test_shift_real_data.o: tests/test_shift_real_data.f90 | $(BUILDDIR) $(BINDIR)
	$(FC) -I$(BUILDDIR) -J$(BUILDDIR) -c $< $(FFLAGS) -o $@ $(LIBS)

# Runs the real-data shift-current test in bin/test_run_shift (the pipeline writes output files into the cwd)
run_test_shift_real: $(BINDIR)/test_shift_real_data
	mkdir -p $(BINDIR)/test_run_shift
	sed 's|@ROOT@|$(CURDIR)|g' tests/test_shift_real_data.in > $(BINDIR)/test_run_shift/in.txt
	cd $(BINDIR)/test_run_shift && ../test_shift_real_data in.txt
	-python3 tools/plot_test_outputs.py $(BINDIR)/test_run_shift/shift_real_data_spectra.dat

# Regression check of the single-particle shift (sp vs exact IPA and vs the excitonic IPA limit on non-interacting hBN)
check_sp_shift: $(TARGET)
	python3 tools/check_sp_shift.py --opticx $(TARGET) --root $(CURDIR) --workdir $(BINDIR)/check_sp_shift

# Cache_ome_ex: the five modes, asserting both behaviour and what the run says about itself (8.44/8.45)
check_ome_cache: $(TARGET)
	python3 tools/check_ome_cache.py --opticx $(TARGET) --root $(CURDIR) --workdir $(BINDIR)/check_ome_cache

# -----------------------------------------------------------------
# Clean
# -----------------------------------------------------------------
clean:
	rm -f $(BUILDDIR)/*.o $(BUILDDIR)/*.mod $(BINDIR)/opticx.o $(TARGET) $(BINDIR)/test_shift_kernel_equivalence.o $(TARGET_TEST)
