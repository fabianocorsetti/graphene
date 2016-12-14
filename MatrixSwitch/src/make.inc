ROOTDIR = ../..
include $(dir $(lastword $(MAKEFILE_LIST)))/$(ROOTDIR)/make.inc

LIBNAME = $(MSLIBNAME)
PSPINCLUDEPATH = ./
