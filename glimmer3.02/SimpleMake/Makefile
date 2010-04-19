# Simple Makefile for Glimmer3

BINDIR = ../bin
OBJDIR = ../obj
LIBDIR = ../lib

VPATH = $(BINDIR):$(OBJDIR):$(LIBDIR)

COMMON_SRCS = delcher.cc fasta.cc gene.cc
COMMON_OBJS = $(COMMON_SRCS:.cc=.o)

GLIMMER_SRCS = anomaly.cc glimmer3.cc long-orfs.cc test.cc
GLIMMER_OBJS = $(GLIMMER_SRCS:.cc=.o)
GLIMMER_PROGS = $(GLIMMER_SRCS:.cc=)

ICM_SRCS = icm.cc build-icm.cc build-fixed.cc score-fixed.cc
ICM_OBJS = $(ICM_SRCS:.cc=.o)
ICM_PROGS = build-icm build-fixed score-fixed

UTIL_SRCS = entropy-profile.cc entropy-score.cc extract.cc multi-extract.cc \
  start-codon-distrib.cc uncovered.cc window-acgt.cc
UTIL_OBJS = $(UTIL_SRCS:.cc=.o)
UTIL_PROGS = $(UTIL_SRCS:.cc=)

SOURCES = $(COMMON_SRCS) $(GLIMMER_SRCS) $(ICM_SRCS) $(UTIL_SRCS)
OBJECTS = $(COMMON_OBJS) $(GLIMMER_OBJS) $(ICM_OBJS) $(UTIL_OBJS)
PROGS = $(COMMON_PROGS) $(GLIMMER_PROGS) $(ICM_PROGS) $(UTIL_PROGS)

LIBRARIES = libGLMcommon.a libGLMicm.a

AR = ar
ARFLAGS = -rsv
CPPC = g++
CFLAGS = -g
LDFLAGS = -lm -L $(LIBDIR)

DEPEND_FILES = *.cc *.hh
CLEANABLE_FILES = $(OBJDIR)/*.o *~

ALL = $(LIBRARIES) $(PROGS)

.SUFFIXES: .cc

.cc.o:
	$(CPPC) $(CFLAGS) -c $*.cc -o $(OBJDIR)/$*.o 

all: $(ALL)

anomaly:  anomaly.o libGLMcommon.a
	$(CPPC) -o $(BINDIR)/$@ $(OBJDIR)/anomaly.o $(LDFLAGS) -lGLMcommon

build-fixed:  build-fixed.o libGLMicm.a libGLMcommon.a
	$(CPPC) -o $(BINDIR)/$@ $(OBJDIR)/build-fixed.o $(LDFLAGS) -lGLMicm -lGLMcommon

build-icm:  build-icm.o libGLMicm.a libGLMcommon.a
	$(CPPC) -o $(BINDIR)/$@ $(OBJDIR)/build-icm.o $(LDFLAGS) -lGLMicm -lGLMcommon

entropy-profile:  entropy-profile.o libGLMcommon.a
	$(CPPC) -o $(BINDIR)/$@ $(OBJDIR)/entropy-profile.o $(LDFLAGS) -lGLMcommon

entropy-score:  entropy-score.o libGLMcommon.a
	$(CPPC) -o $(BINDIR)/$@ $(OBJDIR)/entropy-score.o $(LDFLAGS) -lGLMcommon

extract:  extract.o libGLMcommon.a
	$(CPPC) -o $(BINDIR)/$@ $(OBJDIR)/extract.o $(LDFLAGS) -lGLMcommon

glimmer3:  glimmer3.o libGLMicm.a libGLMcommon.a
	$(CPPC) -o $(BINDIR)/$@ $(OBJDIR)/glimmer3.o $(LDFLAGS) -lGLMicm -lGLMcommon

long-orfs:  long-orfs.o libGLMcommon.a
	$(CPPC) -o $(BINDIR)/$@ $(OBJDIR)/long-orfs.o $(LDFLAGS) -lGLMcommon

multi-extract:  multi-extract.o libGLMcommon.a
	$(CPPC) -o $(BINDIR)/$@ $(OBJDIR)/multi-extract.o $(LDFLAGS) -lGLMcommon

score-fixed:  score-fixed.o libGLMicm.a libGLMcommon.a
	$(CPPC) -o $(BINDIR)/$@ $(OBJDIR)/score-fixed.o $(LDFLAGS) -lGLMicm -lGLMcommon

start-codon-distrib:  start-codon-distrib.o libGLMcommon.a
	$(CPPC) -o $(BINDIR)/$@ $(OBJDIR)/start-codon-distrib.o $(LDFLAGS) -lGLMcommon

test:  test.o libGLMcommon.a
	$(CPPC) -o $(BINDIR)/$@ $(OBJDIR)/test.o $(LDFLAGS) -lGLMcommon

uncovered:  uncovered.o libGLMcommon.a
	$(CPPC) -o $(BINDIR)/$@ $(OBJDIR)/uncovered.o $(LDFLAGS) -lGLMcommon

window-acgt:  window-acgt.o libGLMcommon.a
	$(CPPC) -o $(BINDIR)/$@ $(OBJDIR)/window-acgt.o $(LDFLAGS) -lGLMcommon


libGLMcommon.a:  $(COMMON_OBJS)
	$(AR) $(ARFLAGS) $(LIBDIR)/$@ $(OBJDIR)/delcher.o $(OBJDIR)/fasta.o $(OBJDIR)/gene.o

libGLMicm.a:  icm.o
	$(AR) $(ARFLAGS) $(LIBDIR)/$@ $(OBJDIR)/icm.o

depend:
	makedepend $(DEPEND_FILES)

clean:
	/bin/rm -f $(CLEANABLE_FILES)

# DO NOT DELETE THIS LINE -- make depend depends on it.

