# Makefile for Glimmer3/src/ICM

LOCAL_WORK = $(shell cd ../..; pwd)

ICM_SRCS = icm.cc em_icm.cc build-icm.cc build-fixed.cc score-fixed.cc simple-score.cc em_build-icm.cc
ICM_OBJS = $(ICM_SRCS:.cc=.o)

SOURCES = $(ICM_SRCS)
OBJECTS = $(ICM_OBJS)

PROGS = em_build-icm build-icm build-fixed score-fixed simple-score

LIBRARIES = libGLMicm.a

include  $(LOCAL_WORK)/src/c_make.glm

all:    $(OBJECTS) $(LIBRARIES) $(PROGS)

em_build-icm: em_build-icm.o libGLMicm.a libGLMcommon.a

build-icm:  build-icm.o libGLMicm.a libGLMcommon.a

build-fixed:  build-fixed.o libGLMicm.a libGLMcommon.a

score-fixed:  score-fixed.o libGLMicm.a libGLMcommon.a

simple-score: simple-score.o libGLMicm.a libGLMcommon.a

libGLMicm.a:  $(ICM_OBJS)
