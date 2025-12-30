# --- CONFIGURATION ---
FC = gfortran
FLAGS_COMMON = -std=f2008 -Wall -O2 -J. -cpp

# specific flags:
FLAGS_SERIAL = $(FLAGS_COMMON)
FLAGS_OPENMP = $(FLAGS_COMMON) -fopenmp -DUSE_OMP

# common source files
SOURCES = src/geometry.f90 src/particle.f90 src/io_module.f90 src/barnes_hut.f90 src/main.f90

# executable names
EXEC_SERIAL = barnes_hut_serial
EXEC_OMP = barnes_hut_omp

# --- EXECUTION RULES ---

help:
	@echo "Available options:"
	@echo "  make serial            -> Compile serial version"
	@echo "  make omp               -> Compile parallel OpenMP version"
	@echo "  make input             -> Compile the input generator"
	@echo "  make run_serial N=XX   -> Run serial version with input_disk_XX.dat"
	@echo "  make run_omp N=XX      -> Run OpenMP version with input_disk_XX.dat"
	@echo "  make plot N=XX         -> Run the Python animatiom script for N=XX"
	@echo "  make clean             -> Clean all generated files"

# default target: build both versions
all: serial omp

# --- SERIAL VERSION ---
serial: $(SOURCES)
	@echo "Compiling SERIAL version..."
	$(FC) $(FLAGS_SERIAL) -o $(EXEC_SERIAL) $(SOURCES)
	@echo "Done: ./$(EXEC_SERIAL)"

# --- OPENMP VERSION ---
omp: $(SOURCES)
	@echo "Compiling OPENMP version..."
	$(FC) $(FLAGS_OPENMP) -o $(EXEC_OMP) $(SOURCES)
	@echo "Done: ./$(EXEC_OMP)"

# --- INPUT GENERATOR ---
input_gen: generate_input_disk.f90
	$(FC) $(FLAGS_SERIAL) -o generate_input_disk generate_input_disk.f90

input: input_gen
	./generate_input_disk

# --- PARAMETRIZED EXECUTION ---
# we set a default value if N is not specified
N ?= 50

run_serial: serial
	@if [ ! -f input_disk_$(N).dat ]; then \
		echo "Error: 'input_disk_$(N).dat' not found."; \
		exit 1; \
	fi
	./$(EXEC_SERIAL) < input_disk_$(N).dat

run_omp: omp
	@if [ ! -f input_disk_$(N).dat ]; then \
		echo "Error: 'input_disk_$(N).dat' not found."; \
		exit 1; \
	fi
	export OMP_NUM_THREADS=8 && ./$(EXEC_OMP) < input_disk_$(N).dat 
# 4 threads is the optimal choice for my laptop as it has 4 physical cores; 
#using 8 threads causes contention due to "hyperthreading".

# --- ANIMATION ---
plot:
	python3 animation.py $(N)

# --- CLEANUP ---
clean:
	rm -f *.o *.mod $(EXEC_SERIAL) $(EXEC_OMP) generate_input_disk

.PHONY: help all serial omp input run_serial run_omp plot clean



