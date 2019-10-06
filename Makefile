CC = g++
CFLAGS = -g -Wall
TARGET = analysis
VERSION = 1.0
TARGETV = $(TARGET)-$(VERSION)
INC=constants.h util_func.h container.h message.h
SRC=analysis.cpp
OBJS = analysis.o
DOCS = $(OBJS:%.o=%.txt)


SAVE =  Makefile $(SRC) $(INC)

%.o: %.cpp $(INC)
	$(CC) $(CFLAGS) -o $@ -c $<

ALL: $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS)

clean:
	rm -rf *.o *~ $(TARGET)


tar:$(TARGETV).tar
$(TARGETV).tar: $(SAVE)
	tar -cvf $(TARGETV).tar $(SAVE)
gzip: $(TARGETV).tar.gz
$(TARGETV).tar.gz: $(TARGETV).tar       
	gzip $(TARGETV).tar

save: $(SAVE) $(TARGETV)
	tar -cvf $(TARGETV).tar $(SAVE)
	cd $(TARGETV); tar -xvf ../$(TARGETV).tar
	tar -cvzf $(TARGETV).tar.gz $(TARGETV)
	rm -rf $(TARGETV) $(TARGETV).tar
$(TARGETV):
	mkdir $(TARGETV)
