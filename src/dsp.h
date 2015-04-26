/*

dsp.h

typesdef and constant for libdsp.c

(C) Philippe Strauss <philippe.strauss@urbanet.ch>, 1993, 1995, 1996.

this program is licensed under the term of the GNU Library General Public
License 2, which is included in the file LGPL2.

*/

#ifndef _dsp_h
#define _dsp_h

#include <math.h>

#define DIRECT   (int) 1
#define REVERSE  (int) -1 

#define OK         (int) 0
#define ARRAY_SIZE (int) -1
#define MEM_ALLOC  (int) -2

typedef int tInt;
typedef tInt *ptInt;

typedef float tFloat;
typedef tFloat *ptFloatPt;

typedef int tError;

typedef struct
{
    ptFloatPt pfReal, pfImag;
} tsCmplxRect;

typedef struct
{
    ptFloatPt pfMagn, pfPhase;
} tsCmplxPol;

/*
Exported funcs
*/

void   InitR2SwapTable (ptInt piSwapTable, tInt iN);
void   InitR4SwapTable (ptInt piSwapTable, tInt iLog4N, tInt iN);
void   InitTwidTable   (tsCmplxRect sTwidTable, int iDirectOrReverse, tInt iN);
void   InitWinBlackman (ptFloatPt pfWindow, tInt iN);
void   InitWinKaiser   (ptFloatPt pfWindow, tInt iBeta, tInt iN);
tFloat Bessel_Io       (tFloat fX, tInt iAccuracy);
tFloat Fact            (tInt iX);
void   WindowingCR     (tsCmplxRect sFFT, ptFloatPt pfWindow, tInt iN);
void   WindowingRR     (ptFloatPt pfData, ptFloatPt pfWindow, tInt iN);
void   SwapFFTC        (tsCmplxRect sFFT, ptInt piSwapTable, tInt iN);
void   SwapFFTR        (ptFloatPt pfData, ptInt piSwapTable, tInt iN);
void   DFT             (tsCmplxRect sIn, tsCmplxRect sOut,
								tsCmplxRect sTwidTable, tInt iN);
void   R2FFTdif        (tsCmplxRect sFFT, tsCmplxRect sTwidTable,
                        tInt iLog2N, tInt iN);
void   R2FFTdit        (tsCmplxRect sFFT, tsCmplxRect sTwidTable,
                        tInt iLog2N, tInt iN);
void   R4FFTdif        (tsCmplxRect sFFT, tsCmplxRect sTwidTable,
                        tInt iLog4N, tInt iN);
void   R4FFTdit        (tsCmplxRect sFFT, tsCmplxRect sTwidTable,
                        tInt iLog4N, tInt iN);
void   Magnitude       (tsCmplxRect sIn, ptFloatPt pfOut, tInt iN);
void   Magnitude2      (tsCmplxRect sIn, ptFloatPt pfOut, tInt iN);
void   RectToPol       (tsCmplxRect sIn, tsCmplxPol sOut, tInt iN);
tInt   Sign            (tFloat fIn);
void   Scaling         (tsCmplxRect sFFT, tFloat fCste, tInt iN);
void   AddCmplxRect    (tsCmplxRect sA, tsCmplxRect sB, tInt iBOffset,
                        tInt iNA, tInt iNB);
void   MultCmplxRect   (tsCmplxRect sA, tsCmplxRect sB, tInt iBOffset,
                        tInt iNA, tInt iNB);
tError BlockConvolve   (ptFloatPt pfX, ptFloatPt pfH, ptFloatPt pfY, tInt iNX,
                        tInt iNH, tInt iNY);

#endif /* _dsp_h */
