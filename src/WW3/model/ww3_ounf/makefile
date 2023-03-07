########################################################
PROGRAM  = ww3_ounf
########################################################

SHELL    = /bin/sh

CC       = 
CFLAGS   = 
CPPFLAGS = 
FC       = 
FFLAGS   = 
LDFLAGS  = 
LIBS     = 

INCS     = w3macros.h
SRCS     = ww3_ounf.F90 wmmdatmd.F90 w3gdatmd.F90 w3wdatmd.F90 w3adatmd.F90 w3idatmd.F90 w3odatmd.F90 w3initmd.F90 w3sdb1md.F90 w3sbt1md.F90 w3parall.F90 w3triamd.F90 w3src6md.F90 w3snl1md.F90 w3sis1md.F90 w3iogrmd.F90 w3iogomd.F90 w3iorsmd.F90 w3iopomd.F90 constants.F90 w3servmd.F90 w3timemd.F90 w3arrymd.F90 w3dispmd.F90 w3gsrumd.F90 w3nmlounfmd.F90
OBJS     = $(addsuffix .o, $(basename $(SRCS)))

SWITCHES = F90 NCO NOGRB DIST MPI PR3 UQ ST6 STAB0 FLX4 LN1 NL1 IC3 IS1 REF0 BT1 DB1 MLIM TR0 BS0 XX0 WNX1 WNT1 RWND CRX1 CRT1 O0 O1 O2 O4 O5 O6 O7

########################################################
# targets
########################################################
default : $(PROGRAM)

$(PROGRAM) : $(OBJS)
	$(FC) $(LDFLAGS) -o $(PROGRAM) $(OBJS) $(LIBS)

clean :
	\rm -f *.o *.mod

########################################################
# compile rules
########################################################
.SUFFIXES : .f .f90 .F .F90 .c

%.o : %.c
	$(CC) -c $(CFLAGS) $<

%.o : %.f90
	$(FC) -c $(FFLAGS) $<

%.o : %.F90
	$(FC) -c $(FFLAGS) $<

%.o : %.f
	$(FC) -c $(FFLAGS) $<

%.o : %.F
	$(FC) -c $(FFLAGS) $<

########################################################
# dependencies
########################################################

ww3_ounf.o : ww3_ounf.F90 w3macros.h constants.o w3adatmd.o w3arrymd.o w3gdatmd.o w3initmd.o w3iogomd.o w3iogrmd.o w3nmlounfmd.o w3odatmd.o w3servmd.o w3timemd.o w3wdatmd.o

wmmdatmd.o : wmmdatmd.F90 w3macros.h constants.o w3gdatmd.o w3odatmd.o w3servmd.o

w3gdatmd.o : w3gdatmd.F90 w3macros.h constants.o w3gsrumd.o w3servmd.o

w3wdatmd.o : w3wdatmd.F90 w3macros.h constants.o w3gdatmd.o w3odatmd.o w3parall.o w3servmd.o

w3adatmd.o : w3adatmd.F90 w3macros.h constants.o w3gdatmd.o w3idatmd.o w3odatmd.o w3servmd.o

w3idatmd.o : w3idatmd.F90 w3macros.h w3gdatmd.o w3servmd.o

w3odatmd.o : w3odatmd.F90 w3macros.h constants.o w3gdatmd.o w3servmd.o

w3initmd.o : w3initmd.F90 w3macros.h constants.o w3adatmd.o w3arrymd.o w3dispmd.o w3gdatmd.o w3idatmd.o w3iogomd.o w3iogrmd.o w3iopomd.o w3iorsmd.o w3odatmd.o w3parall.o w3servmd.o w3timemd.o w3wdatmd.o

w3sdb1md.o : w3sdb1md.F90 w3macros.h constants.o w3gdatmd.o w3odatmd.o

w3sbt1md.o : w3sbt1md.F90 w3macros.h w3gdatmd.o

w3parall.o : w3parall.F90 w3macros.h constants.o w3adatmd.o w3gdatmd.o w3idatmd.o w3odatmd.o

w3triamd.o : w3triamd.F90 w3macros.h constants.o w3adatmd.o w3gdatmd.o w3idatmd.o w3odatmd.o w3servmd.o

w3src6md.o : w3src6md.F90 w3macros.h constants.o w3gdatmd.o w3odatmd.o w3servmd.o w3timemd.o w3wdatmd.o

w3snl1md.o : w3snl1md.F90 w3macros.h constants.o w3adatmd.o w3gdatmd.o w3odatmd.o

w3sis1md.o : w3sis1md.F90 w3macros.h w3gdatmd.o w3odatmd.o w3servmd.o

w3iogrmd.o : w3iogrmd.F90 w3macros.h constants.o w3adatmd.o w3dispmd.o w3gdatmd.o w3odatmd.o w3servmd.o w3snl1md.o w3timemd.o

w3iogomd.o : w3iogomd.F90 w3macros.h constants.o w3adatmd.o w3dispmd.o w3gdatmd.o w3odatmd.o w3parall.o w3servmd.o w3wdatmd.o

w3iorsmd.o : w3iorsmd.F90 w3macros.h constants.o w3adatmd.o w3gdatmd.o w3odatmd.o w3parall.o w3servmd.o w3triamd.o w3wdatmd.o

w3iopomd.o : w3iopomd.F90 w3macros.h constants.o w3adatmd.o w3gdatmd.o w3gsrumd.o w3odatmd.o w3servmd.o w3triamd.o w3wdatmd.o

constants.o : constants.F90 w3macros.h

w3servmd.o : w3servmd.F90 w3macros.h constants.o

w3timemd.o : w3timemd.F90 w3macros.h w3odatmd.o w3servmd.o

w3arrymd.o : w3arrymd.F90 w3macros.h w3servmd.o

w3dispmd.o : w3dispmd.F90 w3macros.h constants.o w3gdatmd.o w3odatmd.o w3servmd.o

w3gsrumd.o : w3gsrumd.F90 w3macros.h w3servmd.o

w3nmlounfmd.o : w3nmlounfmd.F90 w3macros.h constants.o w3odatmd.o w3servmd.o
