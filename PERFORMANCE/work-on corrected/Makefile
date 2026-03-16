# Student's Makefile for the CS:APP Performance Lab
TEAM = bovik
VERSION = 1
HANDINDIR = 

# Detect operating system
UNAME_S := $(shell uname -s 2>/dev/null || echo "Windows")
UNAME_M := $(shell uname -m 2>/dev/null || echo "x86_64")

# Default compiler
CC = gcc
CFLAGS = -Wall -O2
LIBS = -lm

# Platform-specific settings
ifeq ($(UNAME_S),Darwin)
    # macOS
    CFLAGS += -DMACOS
    ifeq ($(UNAME_M),arm64)
        # Apple Silicon (M1/M2)
        CFLAGS += -DARM64
    endif
else ifeq ($(UNAME_S),Linux)
    # Linux/WSL
    CFLAGS += -DLINUX
else ifeq ($(OS),Windows_NT)
    # Native Windows (MinGW/MSYS)
    CC = gcc
    CFLAGS += -DWINDOWS
    LIBS += -lws2_32
else ifneq ($(findstring MINGW,$(UNAME_S)),)
    # MinGW/MSYS2
    CC = gcc
    CFLAGS += -DWINDOWS
    LIBS += -lws2_32
else ifneq ($(findstring MSYS,$(UNAME_S)),)
    # MSYS
    CC = gcc
    CFLAGS += -DWINDOWS
    LIBS += -lws2_32
endif

OBJS = driver.o kernels.o fcyc.o clock.o

all: driver

driver: $(OBJS) fcyc.h clock.h defs.h config.h
	$(CC) $(CFLAGS) $(OBJS) $(LIBS) -o driver

handin:
	cp kernels.c $(HANDINDIR)/$(TEAM)-$(VERSION)-kernels.c

clean: 
	-rm -f $(OBJS) driver core *~ *.o


