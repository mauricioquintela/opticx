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
FFLAGS = -fopt-info-vec -g -fcheck=all -O3 -ffree-line-length-none
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
# Optional: read Xatu HDF5 exciton archives (xatu -H)
# Usage: make HDF5=1  (after 'make clean' when switching; needs the
# HDF5 Fortran library, e.g. libhdf5-dev / hdf5-fortran)
# -----------------------------------------------------------------
HDF5_INC ?= /usr/include
ifeq ($(HDF5),1)
H5FLAGS = -DOPTICX_HDF5 -I$(HDF5_INC)
LIBS   += -lhdf5_fortran -lhdf5
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

# xatu_h5.f90 holds the HDF5 reader behind #ifdef OPTICX_HDF5
$(BUILDDIR)/xatu_h5.o: FFLAGS += -cpp $(H5FLAGS)

$(BUILDDIR)/parser_optics_xatu_dim.o: \
	$(BUILDDIR)/constants_math.o \
	$(BUILDDIR)/parser_wannier90_tb.o \
	$(BUILDDIR)/parser_input_file.o \
	$(BUILDDIR)/xatu_h5.o

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
	$(BUILDDIR)/exciton_envelopes.o \
	$(BUILDDIR)/ome_sp.o 

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

test_shg_consistency: $(BINDIR)/test_shg_consistency

$(BINDIR)/test_shg_consistency: $(OBJ_MODULES) $(BINDIR)/test_shg_consistency.o
	$(FC) $(FFLAGS) $(OBJ_MODULES) $(BINDIR)/test_shg_consistency.o -o $@ $(LIBS)

$(BINDIR)/test_shg_consistency.o: tests/test_shg_consistency.f90 | $(BUILDDIR) $(BINDIR)
	$(FC) -I$(BUILDDIR) -J$(BUILDDIR) -c $< $(FFLAGS) -o $@ $(LIBS)

test_shg_real_data: $(BINDIR)/test_shg_real_data

$(BINDIR)/test_shg_real_data: $(OBJ_MODULES) $(BINDIR)/test_shg_real_data.o
	$(FC) $(FFLAGS) $(OBJ_MODULES) $(BINDIR)/test_shg_real_data.o -o $@ $(LIBS)

$(BINDIR)/test_shg_real_data.o: tests/test_shg_real_data.f90 | $(BUILDDIR) $(BINDIR)
	$(FC) -I$(BUILDDIR) -J$(BUILDDIR) -c $< $(FFLAGS) -o $@ $(LIBS)

# Runs the real-data SHG test in bin/test_run (the pipeline writes output files into the cwd)
run_test_shg_real: $(BINDIR)/test_shg_real_data
	mkdir -p $(BINDIR)/test_run
	sed 's|@ROOT@|$(CURDIR)|g' tests/test_shg_real_data.in > $(BINDIR)/test_run/in.txt
	cd $(BINDIR)/test_run && ../test_shg_real_data in.txt
	-python3 tools/plot_test_outputs.py $(BINDIR)/test_run/shg_real_data_spectra.dat

test_second_symmetry: $(BINDIR)/test_second_symmetry

$(BINDIR)/test_second_symmetry: $(OBJ_MODULES) $(BINDIR)/test_second_symmetry.o
	$(FC) $(FFLAGS) $(OBJ_MODULES) $(BINDIR)/test_second_symmetry.o -o $@ $(LIBS)

$(BINDIR)/test_second_symmetry.o: tests/test_second_symmetry.f90 | $(BUILDDIR) $(BINDIR)
	$(FC) -I$(BUILDDIR) -J$(BUILDDIR) -c $< $(FFLAGS) -o $@ $(LIBS)

# Symmetry of the general two-frequency branch at r = w_q/w_p other than 1
run_test_second_symmetry: $(BINDIR)/test_second_symmetry
	mkdir -p $(BINDIR)/test_run_sym
	sed 's|@ROOT@|$(CURDIR)|g' tests/test_second_symmetry.in > $(BINDIR)/test_run_sym/in.txt
	cd $(BINDIR)/test_run_sym && ../test_second_symmetry in.txt

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

# Runs the synthetic SHG consistency test in bin/test_run_shg and plots its spectrum
run_test_shg_consistency: $(BINDIR)/test_shg_consistency
	mkdir -p $(BINDIR)/test_run_shg
	cd $(BINDIR)/test_run_shg && ../test_shg_consistency
	-python3 tools/plot_test_outputs.py $(BINDIR)/test_run_shg/shg_consistency_spectra.dat

# Regression check of the single-particle shift (sp vs exact IPA and vs the excitonic IPA limit on non-interacting hBN)
check_sp_shift: $(TARGET)
	python3 tools/check_sp_shift.py --opticx $(TARGET) --root $(CURDIR) --workdir $(BINDIR)/check_sp_shift

check_sp_shg: $(TARGET)
	python3 tools/check_sp_shg.py --opticx $(TARGET) --root $(CURDIR) --workdir $(BINDIR)/check_sp_shg

# Cache_ome_ex: the five modes, asserting both behaviour and what the run says about itself
check_ome_cache: $(TARGET)
	python3 tools/check_ome_cache.py --opticx $(TARGET) --root $(CURDIR) --workdir $(BINDIR)/check_ome_cache

# Gauge covariance: scrambled eigenvector phases must not move X_nm or any second-order output
check_gauge_covariance: $(TARGET)
	python3 tools/check_gauge_covariance.py --opticx $(TARGET) --root $(CURDIR) --workdir $(BINDIR)/check_gauge_covariance

check_shift_covariant: $(TARGET)
	python3 tools/check_shift_covariant.py --opticx $(TARGET) --root $(CURDIR) --workdir $(BINDIR)/check_shift_covariant

check_realtime_sign: $(TARGET)
	python3 tools/check_realtime_sign.py --opticx $(TARGET) --root $(CURDIR) --workdir $(BINDIR)/check_realtime_sign

# Excitonic rectification, the whole causal sigma(0; w, -w): NI limit (Re and injection), reality, selection rule, NumPy
check_ex_rectification: $(TARGET)
	python3 tools/check_ex_rectification.py --opticx $(TARGET) --root $(CURDIR) --workdir $(BINDIR)/check_ex_rectification

# Out-of-plane components and the injection sign vs a real-time propagation (non-interacting buckled hBN)
check_out_of_plane: $(TARGET)
	python3 tools/check_out_of_plane.py --opticx $(TARGET) --root $(CURDIR) --workdir $(BINDIR)/check_out_of_plane

# Band structure along a k-path (Kpath input, default path, Response = bands)
check_bands: $(TARGET)
	python3 tools/check_bands.py --opticx $(TARGET) --root $(CURDIR) --workdir $(BINDIR)/check_bands

# Band window: report and warning for an unusual Bandlist (a list of offsets, not a range)
check_bandlist_guard: $(TARGET)
	python3 tools/check_bandlist_guard.py --opticx $(TARGET) --root $(CURDIR) --workdir $(BINDIR)/check_bandlist_guard

# Xatu HDF5 archive input (needs a 'make HDF5=1' build and h5py): archive = text files, guards.
# XATU=path/to/xatu (built with HDF5=1) adds a real Xatu run written both ways.
check_xatu_h5: $(TARGET)
	python3 tools/check_xatu_h5.py --opticx $(TARGET) --root $(CURDIR) --workdir $(BINDIR)/check_xatu_h5 $(if $(XATU),--xatu $(XATU))

# _tb.dat reader: Hermiticity check and repair of H, S and r
check_tb_hermiticity: $(TARGET)
	python3 tools/check_tb_hermiticity.py --opticx $(TARGET) --root $(CURDIR) --workdir $(BINDIR)/check_tb_hermiticity

# Eq. (A4) basis guard: OME_sp = none must refuse the excitonic path unless a cache covers it
check_a4_basis_guard: $(TARGET)
	python3 tools/check_a4_basis_guard.py --opticx $(TARGET) --root $(CURDIR) --workdir $(BINDIR)/check_a4

# -----------------------------------------------------------------
# Clean
# -----------------------------------------------------------------
clean:
	rm -f $(BUILDDIR)/*.o $(BUILDDIR)/*.mod $(BINDIR)/opticx.o $(TARGET) $(BINDIR)/test_shift_kernel_equivalence.o $(TARGET_TEST)
