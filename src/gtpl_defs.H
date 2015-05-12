//
// file gtpl_defs.H
// David Cosgrove
// AstraZeneca
// 31st July 2007
//
// Some definitions used by Triphic, Plurality and Loob, in their
// own namespace.

#ifndef DAC_GTPL_DEFS
#define DAC_GTPL_DEFS

namespace GtplDefs {

typedef enum { RMS_AND_SIZE , GRID_SHAPE_TANI , GAUSS_SHAPE_TANI ,
               SURFACE_OVLP_VOLUME , PROTEIN_CLASH , MMFF_NRG ,
               ROBINS_SCORE_PARETO , OVERALL_SCORE_PARETO } SCORE_METHOD;
typedef enum { IGNORE , SCREEN , ALIGN } DIRS_USAGE;
typedef enum { RING_NORMAL , H_VECTOR , LP_VECTOR , UNKNOWN } DIRS_TYPE;
typedef enum { NONE , FREE , FLIP } DIR_MOVES;
typedef enum { OUTSIDE = 0 , OUTSIDE_SHELL = 1 , INSIDE_SHELL = 2 ,
               CORE = 4 } VOL_MARKER;

static const float FLT_TOL = 1.0e-7;

};

#endif
