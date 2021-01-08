#RSS = -DUSE_RSS
RSS = 

CFLAGS =    $(RSS) -DNDEBUG -O5  
#CFLAGS =   $(RSS) -D__cplusplus -g    

LIBS = 	-lstdc++ -lm 

CC       = g++

current: mcs

MCS_C =	slist.C chain.C eef1.C rss.C crmsd.C spheres.C energy.cpp \
		node.C leaf.C simulation.C pairtree.C pdb.C 

MCS_OBJS	=	$(MCS_C:.C=.o)

mcs:	$(MCS_OBJS) Makefile
	rm -f mcs	
	$(CC)  $(CFLAGS) $(MCS_OBJS) $(LIBS) -o mcs

.C.o:
	$(CC) $(CFLAGS) -c $<

clean:
	rm *.o

depend:
	makedepend $(INCL) $(MCS_C)


# DO NOT DELETE THIS LINE -- make depend depends on it.
