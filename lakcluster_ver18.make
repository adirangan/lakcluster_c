
# Compiler
CC=gcc
CCFLAGS_dev= -w -O2 -lpthread -fopenmp -lm -lgslcblas -I./dir_h -I./dir_c 
CCFLAGS=$(CCFLAGS_dev) -D_MONOLITH

vpath %.c ./dir_c

project=lakcluster_ver18
project_dev=lakcluster_ver18_dev

CLINK=$(CC) -o $(project_dev)

header_dev = ./dir_h/S_handle.h \
	./dir_h/P_handle.h \
	./dir_h/R_handle.h \
	./dir_h/L_handle.h \
	./dir_h/M_handle.h \
	./dir_h/bcc.h \
	./dir_h/dcc.h \
	./dir_h/dcg.h \
	./dir_h/lakcluster_header.h

sources_dev = boxmuller.c \
	An_ajdk_v.c \
	AnAt_vv.c \
	An_v.c \
	AnZt_mm.c \
	An_ZtSn_ww.c \
	An_ZtSWn_Yt_vv.c \
	AnZt_S_WnYt_vv.c \
	An_ZtSWn_Yt_ww.c \
	AnZt_vv.c \
	AttAn_vv.c \
	AtTAn_vv.c \
	AttYn_vv.c \
	AtTYn_vv.c \
	AttYn____WtsZn_vv.c \
	At_T_YnWt_S_Zn_vv.c \
	AtTYn____WtSZn_vv.c \
	At_T_YnWt_S_Zn_ww.c \
	bcc.c \
	bcc_flattenloop.c \
	bcc_lf.c \
	bcc_lrup.c \
	bcc_scorebox.c \
	bcc_sumscores.c \
	binary_read.c \
	D_AtTn_ZtSn.c \
	dcc.c \
	dcc_lf.c \
	dcc_scorebox.c \
	dcc_sumscores.c \
	dexcluster_driver.c \
	dexcluster_scorebox.c \
	dcg.c \
	dcg_lf.c \
	dcg_scorebox.c \
	dcg_sumscores.c \
	dcgxpander_driver.c \
	fillzero.c \
	GLOBAL_pthread.c \
	lakcluster_driver.c \
	lakcluster_scorebox.c \
	L_init.c \
	M_Ax_to_L2.c \
	mda_io.c \
	M_init.c \
	PNMfile.c \
	popcount.c \
	Quicksort.c \
	raprintf.c \
	ra_stats.c \
	S_init.c \
	TAnZtS.c \
	updateglobals.c \
	wkspace.c \
	xcalc.c \
	xdrop.c \
	At_T_Xn_ww.c \
	An_Xn_ww.c \
	P_init.c \
	pca_driver.c \
	R_init.c \
	lakcluster_ver18.c

objects_dev = $(patsubst %.c,./dir_o/%.o,$(sources_dev)) 

all: $(objects_dev) $(project_dev)

$(objects_dev): | ./dir_o

./dir_o:
	@mkdir -p $@

$(project_dev): $(objects_dev) $(header_dev)
	rm -f $(project_dev)
	$(CLINK) $(objects_dev) $(CCFLAGS_dev)

$(project): $(objects_dev) $(header_dev)
	rm -f $(project)
	gcc -w -O2 -lpthread -fopenmp ./dir_c/lakcluster_ver18.c -o $(project) -L./dir_h -L./dir_c -I./dir_h -I./dir_c -D_MONOLITH -lgslcblas -lm

./dir_o/%.o : %.c ./dir_h/lakcluster_header.h
	@echo $< 
	@$(CC) $(CCFLAGS_dev) -c $< -o $@

.c.o: ./dir_h/lakcluster_header.h
	$(CC) -c $(CCFLAGS_dev) $<

clean: 
	rm -f $(objects_dev)

list_sources: $(sources_dev)
	echo $^

list_objects: $(objects_dev)
	echo $^

