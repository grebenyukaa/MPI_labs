MPICC ?= mpicc
MPICC-TV ?= mpicc-TV

CFLAGS = -Wall -std=c99
LDFLAGS = -lm

CMN_SOURCE = utils.c
INCLUDE = .
OUTPUT_DIR = /gpfs/home/iu7/stud051/labs

.PHONY: lab1 lab2 lab3 lab4 lab5 clean all

lab5:
	$(MPICC-TV) $(CFLAGS) $(LDFLAGS) $(CMN_SOURCE) -I$(INCLUDE) lab5/*.c -o $(OUTPUT_DIR)/lab5

lab4:
	$(MPICC) $(CFLAGS) $(LDFLAGS) $(CMN_SOURCE) -I$(INCLUDE) lab4/*.c -o $(OUTPUT_DIR)/lab4

lab3:
	$(MPICC) $(CFLAGS) $(LDFLAGS) $(CMN_SOURCE) -I$(INCLUDE) lab3/*.c -o $(OUTPUT_DIR)/lab3

lab2:
	$(MPICC) $(CFLAGS) $(LDFLAGS) $(CMN_SOURCE) -I$(INCLUDE) lab2/*.c -o $(OUTPUT_DIR)/lab2

lab1:
	$(MPICC) $(CFLAGS) $(LDFLAGS) $(CMN_SOURCE) -I$(INCLUDE) lab1/*.c -o $(OUTPUT_DIR)/lab1

all: lab1 lab2 lab3 lab4 lab5

clean:
	@rm -f $(OUTPUT_DIR)/labs/*
