##########################################
# Makefile                               #
# Makefile for the code developed in TP1 #
#                                        #
# T. Dufaud                              #
##########################################
################################
# Variables for this makefile
################################
# 
# -- option to dedicated machine
#
# Rk: You should create a file such as your-machineName.mk
# Follow the example of the machine called "ambre" in the 
# file ambre.mk
#
HOSTNAME?=$(shell hostname)
include $(HOSTNAME).mk

OPTC=${OPTCLOCAL}

TPDIR=.
TPDIRSRC=$(TPDIR)/src

LIBS=${LIBSLOCAL}

INCLATLAS=${INCLUDEBLASLOCAL}
INCL= -I $(TPDIR)/include $(INCLATLAS)

OBJENV= tp_env.o
OBJLIBPOISSON= lib_poisson1D.o lib_poisson1D_writers.o lib_poisson1D_richardson.o
OBJTP2ITER= $(OBJLIBPOISSON) tp_poisson1D_iter.o
OBJTP2DIRECT= $(OBJLIBPOISSON) tp_poisson1D_direct.o

.PHONY: all data_dir img_dir

all: data_dir img_dir bin/tp_testenv bin/tpPoisson1D_iter bin/tpPoisson1D_direct

data_dir:
	mkdir -p data

img_dir:
	mkdir -p img

%.o : $(TPDIRSRC)/%.c
	$(CC) $(OPTC) -c $(INCL) $<

bin/tp_testenv: $(OBJENV)
	$(CC) -o $@ $(OPTC) $(OBJENV) $(LIBS)

bin/tpPoisson1D_iter: data_dir img_dir $(OBJTP2ITER)
	$(CC) -o $@ $(OPTC) $(OBJTP2ITER) $(LIBS)

bin/tpPoisson1D_direct: data_dir img_dir $(OBJTP2DIRECT)
	$(CC) -o $@ $(OPTC) $(OBJTP2DIRECT) $(LIBS)

run_testenv:
	bin/tp_testenv

run_tpPoisson1D_iter:
	bin/tpPoisson1D_iter

run_tpPoisson1D_direct:
	bin/tpPoisson1D_direct

run: run_testenv run_tpPoisson1D_iter run_tpPoisson1D_direct

clean:
	rm -f *.o bin/*
