#ifndef TRANSFER_H
#define TRANSFER_H

#include <QMutex>

extern int showingVTK;

#define MAX_CELLS 200000
#define MAX_BC 100000
#define MAX_DC 500
#define MAX_BOND 1
#define NINFO 12

extern int VTKbuffer[100];
extern double ecell_list[NINFO*MAX_CELLS];  // pos[3],tensor[9]
extern int necell_list;

extern int summaryData[100];
extern double concData[4000];
//extern bool concUsed[16];
extern int conc_nc;
extern double conc_dx;
extern int NX, NY, NZ;
extern int MAX_CHEMO;
extern int nvariables;
extern int ncomponents;
extern double volProb[100];
extern int vol_nv;
extern double vol_v0;
extern double vol_dv;
extern double oxyProb[100];
extern int oxy_nv;
extern double oxy_dv;
extern int nt_vtk;
extern int istep;
extern bool leftb;
extern double DELTA_T;
extern double dfraction;
extern int i_hypoxia_cutoff;
extern int i_growth_cutoff;
extern bool goflag;

#endif // TRANSFER_H
