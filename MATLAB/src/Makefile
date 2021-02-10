#########################
# customise these
CFILES := anyNonFinite.c main.c sqrt.c benchmarkRRT.c norm.c state_traj.c benchmarkRRT_data.c opt_time.c xdlapy3.c benchmarkRRT_emxutil.c rand.c xgehrd.c benchmarkRRT_initialize.c recip.c xnrm2.c benchmarkRRT_rtwutil.c roots.c xzgeev.c benchmarkRRT_terminate.c RRTstar3D.c xzhgeqz.c collision.c rtGetInf.c xzhseqr.c collisionKnowTf.c rtGetNaN.c xzlarfg.c cost_eval.c rt_nonfinite.c xzlartg.c eml_rand_mt19937ar_stateful.c segment_cost.c linspace.c sqr_eucl_dist.c

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
