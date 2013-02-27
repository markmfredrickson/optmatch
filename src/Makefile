OS := $(shell uname -s)

ifeq ($(OS), Linux)
	include makefile.linux
endif

ifeq ($(OS), Darwin)
	include makefile.osx
endif

#include makefile.windows
