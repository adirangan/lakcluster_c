# to redirect use: (in bash);
# make lakcluster_ver18 > error.log 2>&1 ; grep "error" error.log;
# make lakcluster_ver18 > error.log 2>&1 ; grep "Wunused-variable" error.log;

CC = gcc
CFLAGS= -O2 -lpthread #-pg
#CFLAGS= -lpthread #-pg
#CFLAGS= -Wall -O2 -lpthread #-pg
#CFLAGS= -Wall -lpthread #-pg

VER18FILES= dir_c/An_ajdk_v.c dir_c/AnAt_vv.c dir_c/An_v.c dir_c/An_ZtSn_ww.c dir_c/An_ZtSWn_Yt_vv.c dir_c/AnZt_S_WnYt_vv.c dir_c/An_ZtSWn_Yt_ww.c dir_c/AnZt_vv.c dir_c/AttAn_vv.c dir_c/AtTAn_vv.c dir_c/AttYn_vv.c dir_c/AtTYn_vv.c dir_c/AttYn____WtsZn_vv.c dir_c/At_T_YnWt_S_Zn_vv.c dir_c/AtTYn____WtSZn_vv.c dir_c/At_T_YnWt_S_Zn_ww.c dir_c/bcc.c dir_c/bcc_flattenloop.c dir_h/bcc.h dir_c/bcc_lf.c dir_c/bcc_lrup.c dir_c/bcc_sumscores.c dir_c/binary_read.c dir_c/D_AtTn_ZtSn.c dir_c/dcc.c dir_h/dcc.h dir_c/dcc_lf.c dir_c/dcc_sumscores.c dir_c/dexcluster_driver.c dir_c/fillzero.c dir_c/GLOBAL_pthread.c dir_c/lakcluster_driver.c dir_c/bcc_scorebox.c dir_c/dcc_scorebox.c dir_c/lakcluster_scorebox.c dir_c/dexcluster_scorebox.c dir_c/mda_io.c dir_c/S_init.c dir_c/L_init.c dir_c/M_init.c dir_c/PNMfile.c dir_c/popcount.c dir_c/Quicksort.c dir_c/raprintf.c dir_c/rastats.c dir_c/TAnZtS.c dir_c/updateglobals.c dir_c/wkspace.c dir_c/xcalc.c dir_c/xdrop.c

lakcluster_ver18: lakcluster_ver18.c $(VER18FILES)
	gcc $(CFLAGS) -fopenmp -lm lakcluster_ver18.c -o lakcluster_ver18 -L.
lakcluster_ver18.o : lakcluster_ver18.c 

clean: lakcluster_ver18.c $(VER18FILES)
	rm lakcluster_ver18.o lakcluster_ver18
