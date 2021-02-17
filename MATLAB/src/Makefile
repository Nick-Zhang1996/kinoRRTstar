#########################
# customise these
CFILES := $(wildcard ./*.c)

PROG := benchmarkRRT
CFLAGS := -Wall -Wextra -g
LDFLAGS := 
LDLIBS := -lm
########################

# -MMD generates dependencies while compiling
CFLAGS += -MMD
CC := gcc

OBJFILES := $(CFILES:.c=.o)
DEPFILES := $(CFILES:.c=.d)

$(PROG) : $(OBJFILES)
	$(LINK.o) $(LDFLAGS) -o $@ $^ $(LDLIBS)

clean :
	rm -f $(PROG) $(OBJFILES) $(DEPFILES)

-include $(DEPFILES)
