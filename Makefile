# Dynamic Power Capping
# - depends on PoLiMEr and PAPI
# - so set POLILIB and POLIINC to correct paths

CFLAGS   = -DTHETA -fopenmp
PAPIPATH = 
POLIPATH = 

LIBDIR  = lib
OBJDIR  = bin
INCDIR  = include
TESTDIR = test

CC = cc
AR = ar

CFLAGS += -Wall -fPIC
TESTCFLAGS = -g -Wall -Werror -fopenmp
ifeq ($(NOMPI), yes)
CFLAGS += -D_NOMPI
endif
ifeq ($(NOOMP), yes)
CFLAGS += -D_NOOMP
endif

INCLUDE =\
  -I$(PAPIPATH)/include\
  -I$(POLIPATH)/include\
  -I./$(INCDIR)

DYNAMICLFLAGS =\
  -L$(PAPIPATH)/lib\
  -L$(POLIPATH)/lib\
  -lpapi\
  -lpolimer\

STATICLFLAGS =\
  -L$(PAPIPATH)/libpapi.a\
  -L$(POLIPATH)/libpolimer.a


SRCS = dynamicap.c
OBJS = $(SRCS:%.c=$(OBJDIR)/%.o)

.PHONY: all lib clean


all: lib

#lib: $(LIBDIR)/libdynamicap.a $(LIBDIR)/libdynamicap.so
lib: $(LIBDIR)/libdynamicap.so $(LIBDIR)/libdynamicap.a
	
$(OBJDIR)/%.o: ./%.c 
		$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@ $(STATICLFLAGS)

$(LIBDIR)/libdynamicap.a: $(OBJS)
		$(AR) rcs $@ $(OBJS) 

$(LIBDIR)/libdynamicap.so: $(OBJS)
		$(CC) -shared -o $@ $(OBJS) 

clean:
		rm -f $(OBJDIR)/*.o
		rm -f $(LIBDIR)/*.a
		rm -f $(LIBDIR)/*.so
		rm -f $(TESTDIR)/test
