PROGNAME = plot

SOURCEFILES = plot.C \
              src/binHistos.C \
              src/dpdxGraph.C \
              src/elPlot.C \
              src/sT.C \
              src/splitTreeForPlot.C

ROOTCFLAGS   = $(shell root-config --cflags)
ROOTLIBS     = $(shell root-config --libs)
ROOTGLIBS    = $(shell root-config --glibs)
ROOTINC      = $(shell root-config --incdir)

CC      = g++
OBJS    = $(patsubst %.C, %.o, $(SOURCEFILES))
INCPATH = -I$(ROOTINC) -Iinc -I$(ROOFITSYS)/include
CFLAGS  = -Wall 
#-Werror
#prova a mettere -O2
#cmq prova su lxplus
LIBS    = $(ROOTLIBS) -L. -L$(ROOFITSYS)/lib -lRooFit -lRooFitCore

%.o: %.C
	${CC} ${INCPATH} ${CFLAGS} -c -o $@ $<

${PROGNAME}: ${OBJS}
	${CC} -o $@ ${INCPATH} ${OBJS} ${LIBS}

clean:	
	-rm -f ${PROGNAME} ${OBJS}

