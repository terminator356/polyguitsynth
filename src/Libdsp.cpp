/*

libdsp.c

Low-level library of signal processing functions. Memory is never allocated
inside one of these functions, the user must deal with the OS himself.
Calculus are done in place if possible.
The goal was to provide efficient routines in C to be used on a computer, the
memory requirement is rather large.

(C) Philippe Strauss <philippe.strauss@urbanet.ch>, 1993, 1995, 1996.

this program is licensed under the term of the GNU General Public License 2,
which is included in the file GPL2.

*/

#include "dsp.h"

/*
InitR2SwapTable : Fill in an array piSwapTable with index for reordering a
                  shuffled input or output of a R2 FFT.
Input           : Array pointed by piSwapTable
Output          : Initialized array pointed by piSwapTable
Comment         : Quiet efficient algorithm (for processor without hardwired
                  bit reversed addressing), quiet original, too, there'not any
                  >> or << operator, its pure arithmetic. 100% P.Strauss
                  design.  
*/

void InitR2SwapTable (ptInt piSwapTable, tInt iN)
{
    tInt iCnt1, iL, iM;

    iL = iN / 2;
    iM = 1;
    piSwapTable[0] = 0;

    while (iL >= 1)
    {
        for (iCnt1 = 0; iCnt1 < iM; ++iCnt1)
            piSwapTable[iCnt1 + iM] = piSwapTable[iCnt1] + iL;

        iL /= 2;
        iM *= 2;
    }
}

/*
InitR$SwapTable : Fill in an array piSwapTable with index for reordering a
                  shuffled input or output of a R4 FFT.
Input           : Array pointed by piSwapTable, base 4 log of N in iLog4N
Output          : Initialized array pointed by piSwapTable
Comment         : Quiet efficient algorithm (for processor without hardwired
                  bit reversed addressing), quiet original, too, there'not any
                  >> or << operator, its pure arithmetic. 100% P.Strauss
                  design.  
*/

void InitR4SwapTable (ptInt piSwapTable, tInt iLog4N, tInt iN)
{
    tInt iL, iM, iCnt1, iCnt2;

    iL = iN / 4;
    iM = 1;
    piSwapTable[0] = 0;

    for (iCnt1 = 0; iCnt1 < iLog4N; ++iCnt1)
    {
        for (iCnt2 = 0; iCnt2 < iM; ++iCnt2)
        {
            piSwapTable[    iM + iCnt2] = piSwapTable[iCnt2] +     iL;
            piSwapTable[2 * iM + iCnt2] = piSwapTable[iCnt2] + 2 * iL;
            piSwapTable[3 * iM + iCnt2] = piSwapTable[iCnt2] + 3 * iL;
        }
        iL /= 4;
        iM *= 4;
    }
}

/*
InitTwidTable : Initialize sin and cos table for later use in FFT routines.
Input         : Array pointed by sTwidTable.{piReal, piImag}
Output        : Initialized arrays
Comment       : Memory requirement is large
*/

void InitTwidTable (tsCmplxRect sTwidTable, int iDirectOrReverse, tInt iN)
{
    tInt iP;
    tFloat fTheta, fPhi;

    fTheta = (2 * M_PI) / iN;
    if (iDirectOrReverse == DIRECT)
    { 
        for (iP = 0; iP < iN; ++iP)
        {
            fPhi = fTheta * iP;
            sTwidTable.pfReal[iP] = cos (-fPhi);
            sTwidTable.pfImag[iP] = sin (-fPhi);
        }
    }
    if (iDirectOrReverse == REVERSE)
    {
        for (iP = 0; iP < iN; ++iP)
        {
            fPhi = fTheta * iP;
            sTwidTable.pfReal[iP] = cos (fPhi);
            sTwidTable.pfImag[iP] = sin (fPhi);
        }
    }
}

/*
InitWinBlackman : Initialize an array pfWindow with a Blackman window
*/

void InitWinBlackman (ptFloatPt pfWindow, tInt iN)
{
    tInt iCnt1;

    for (iCnt1 = 0; iCnt1 < iN; ++iCnt1)
    {
        pfWindow[iCnt1] = 0.42 - 0.50 * cos ((2 * M_PI * iCnt1) / (iN - 1))
                               + 0.08 * cos ((4 * M_PI * iCnt1) / (iN - 1));
    }
}

/*
InitWinKaiser : Initialize an array pfWindow with a Kaiser window
iBeta         : Low value: good spectral resolution, small gibbson effects
                attenuation.
Comment       : Due to the way Bessel_Io is calculated, this algorithm sucks!
*/

void InitWinKaiser (ptFloatPt pfWindow, tInt iBeta, tInt iN)
{
    tInt iCnt1;

    for (iCnt1 = 0; iCnt1 < iN; ++iCnt1)
    {
        pfWindow[iCnt1] = Bessel_Io (2 * iBeta * sqrt (iCnt1 /
                                     (tFloat) (iN - 1) - pow (iCnt1 /
                                     (tFloat) (iN - 1), 2)), 10) /
                                     Bessel_Io ((tFloat) (iBeta), 10);
    }
}

/*
Bessel_Io : Calculation of the zero order Bessel function.
Comment   : This one really sucks!
*/

tFloat Bessel_Io (tFloat fX, tInt iAccuracy)
{
    tInt iCnt1;
    tFloat fOut = 0, fX_2 = fX / 2;

    for (iCnt1 = 1; iCnt1 < iAccuracy; ++iCnt1)
        fOut = fOut + pow ((pow (fX_2, iCnt1) / Fact (iCnt1)), 2);

    fOut += 1.0;
    return fOut;
}

/*
Fact : rather explicit
*/

tFloat Fact (tInt iX)
{
    tInt iCnt1;
    tFloat fFact = iX;

    for (iCnt1 = iX - 1; iCnt1 > 0; --iCnt1)
        fFact *= iCnt1;

    return fFact;
}

/*
WindowingCR : Point - point multiplication between an array of complex and an
              array of real
*/

void WindowingCR (tsCmplxRect sFFT, ptFloatPt pfWindow, tInt iN)
{
    tInt iCnt1;
    tFloat WinTemp;

    for (iCnt1 = 0; iCnt1 < iN; ++iCnt1)
    {
        WinTemp = pfWindow[iCnt1]; /* a little faster with GCC 2.6.3 */
        sFFT.pfReal[iCnt1] = sFFT.pfReal[iCnt1] * WinTemp;
        sFFT.pfImag[iCnt1] = sFFT.pfReal[iCnt1] * WinTemp;
        /*
           Using *= for the two lines above slow down the Windowing by a
           4.5 ratio with GCC !! dont ask why.
           (2.6.3 -O2 -m486)
        */
    }
}

/*
WindowingRR : Point - point multiplication between two arrays of real.
*/

void WindowingRR (ptFloatPt pfData, ptFloatPt pfWindow, tInt iN)
{
    tInt iCnt1;

    for (iCnt1 = 0; iCnt1 < iN; ++iCnt1)
        pfData[iCnt1] = pfData[iCnt1] * pfWindow[iCnt1];
        /*
           Using *= for the line above slow down the Windowing by a
           4.5 ratio with GCC !! dont ask why.
           (2.6.3 -O2 -m486)
        */
}

/*
SwapFFTC : Do the (in place) swapping for reordering the in/out coeffs
           in bit or digit4 reversed order. piSwapTable [] must be
           initialized first.
*/

void SwapFFTC (tsCmplxRect sFFT, ptInt piSwapTable, tInt iN)
{
    tInt iSwapAddr, iCnt1;
    tFloat fRealOutSwap, fImagOutSwap;

    for (iCnt1 = 0; iCnt1 < iN; ++iCnt1)
    {
        fRealOutSwap = sFFT.pfReal[iCnt1];
        fImagOutSwap = sFFT.pfImag[iCnt1];
        iSwapAddr    = piSwapTable[iCnt1];
        if (iCnt1 < iSwapAddr)
        {
            sFFT.pfReal[iCnt1]     = sFFT.pfReal[iSwapAddr];
            sFFT.pfImag[iCnt1]     = sFFT.pfImag[iSwapAddr];
            sFFT.pfReal[iSwapAddr] = fRealOutSwap;
            sFFT.pfImag[iSwapAddr] = fImagOutSwap;
        }
    }
}

/*
SwapFFTC : Do the (in place) swapping for reordering the in/out coeffs
           in bit or digit4 reversed order. piSwapTable [] must be
           initialized first.
Comment  : This one is for real only coeffs.
*/

void SwapFFTR (ptFloatPt pfData, ptInt piSwapTable, tInt iN)
{
    tInt iSwapAddr, iCnt1;
    tFloat fOutSwap;

    for (iCnt1 = 0; iCnt1 < iN; ++iCnt1)
    {
        fOutSwap  = pfData     [iCnt1];
        iSwapAddr = piSwapTable[iCnt1];
        if (iCnt1 < iSwapAddr)
        {
            pfData[iCnt1]     = pfData[iSwapAddr];
            pfData[iSwapAddr] = fOutSwap;
        }
    }
}

/*
DFT   : Discrete Fourier Transform. Plain simple coding of

         N-1
         ___
         \
          \
  X(k) =  /   x(n) * e ^ (-2*pi*j*k*n/N)
         /___

         n=0
*/

void DFT (tsCmplxRect sIn, tsCmplxRect sOut, tsCmplxRect sTwidTable, tInt iN)
{
    tInt iIdx, iK, iKN;
    tFloat fRealSum = 0.0, fImagSum = 0.0;

    for (iK = 0; iK < iN; ++iK)
    {
        iKN = 0;
        for (iIdx = 0; iIdx < iN; ++iIdx)
        {
            fRealSum += sIn.pfReal[iIdx] * sTwidTable.pfReal[iKN]
                      - sIn.pfImag[iIdx] * sTwidTable.pfImag[iKN];
            fImagSum += sIn.pfReal[iIdx] * sTwidTable.pfImag[iKN]
                      + sIn.pfImag[iIdx] * sTwidTable.pfReal[iKN];
            iKN += iK;
            if (iKN >= iN)
                iKN = iKN - iN; /* circular addressing */
        }
        sOut.pfReal[iK] = fRealSum;
        sOut.pfImag[iK] = fImagSum;
        fRealSum = fImagSum = 0.0;
    }
}

/*

R2FFTdif : Radix 2, in place, complex in & out, decimation in frequency FFT
           algorithm.
Required : A tsCmplxRect array with data to transform, another tsCmplxRect
           array with data twiddle factors (sine & cosine table),
           base 2 log of N, N
*/

void R2FFTdif (tsCmplxRect sFFT, tsCmplxRect sTwidTable, tInt iLog2N, tInt iN)
{
  tInt
      iCnt1, iCnt2, iCnt3,
      iQ,    iL,    iM,
      iA,    iB;
  tFloat
      fRealTemp, fImagTemp,
      fReal_Wq,  fImag_Wq;

    iL = 1;
    iM = iN / 2;

    for (iCnt1 = 0; iCnt1 < iLog2N; ++iCnt1)
    {
        iQ = 0;
        for (iCnt2 = 0; iCnt2 < iM; ++iCnt2)
        {
            iA = iCnt2;
            fReal_Wq = sTwidTable.pfReal[iQ];
            fImag_Wq = sTwidTable.pfImag[iQ];

            for (iCnt3 = 0; iCnt3 < iL; ++iCnt3)
            {
                iB = iA + iM;

                /* Butterfly: 10 FOP, 4 FMUL, 6 FADD */

                fRealTemp        = sFFT.pfReal[iA] - sFFT.pfReal[iB];
                sFFT.pfReal[iA] +=                   sFFT.pfReal[iB];
                fImagTemp        = sFFT.pfImag[iA] - sFFT.pfImag[iB];
                sFFT.pfImag[iA] +=                   sFFT.pfImag[iB];

                sFFT.pfReal[iB]  = fRealTemp * fReal_Wq - fImagTemp * fImag_Wq;
                sFFT.pfImag[iB]  = fImagTemp * fReal_Wq + fRealTemp * fImag_Wq;

                iA = iA + 2 * iM;
            }
            iQ += iL;
        }
        iL *= 2;
        iM /= 2;
    }
}

/*

R2FFTdit : Radix 2, in place, complex in & out, decimation in time FFT
           algorithm.
Required : A tsCmplxRect array with data to transform, another tsCmplxRect
           array with data twiddle factors (sine & cosine table),
           base 2 log of N, N
*/

void R2FFTdit (tsCmplxRect sFFT, tsCmplxRect sTwidTable, tInt iLog2N, tInt iN)
{
    tInt
        iCnt1, iCnt2,iCnt3,
        iQ,    iL,   iM,
        iA,    iB;
    tFloat
        fRealTemp, fImagTemp,
        fReal_Wq,  fImag_Wq;

    iL = iN / 2;
    iM = 1;

    for (iCnt1 = 0; iCnt1 < iLog2N; ++iCnt1)
    {
        iQ = 0;
        for (iCnt2 = 0; iCnt2 < iM; ++iCnt2)
        {
            iA = iCnt2;
            fReal_Wq = sTwidTable.pfReal[iQ];
            fImag_Wq = sTwidTable.pfImag[iQ];

            for (iCnt3 = 0; iCnt3 < iL; ++iCnt3)
            {
                iB = iA + iM;

                /* Butterfly: 10 FOP, 4 FMUL, 6 FADD */

                fRealTemp = sFFT.pfReal[iB] * fReal_Wq
                          - sFFT.pfImag[iB] * fImag_Wq;
                fImagTemp = sFFT.pfReal[iB] * fImag_Wq
                          + sFFT.pfImag[iB] * fReal_Wq;

                sFFT.pfReal[iB]  = sFFT.pfReal[iA] - fRealTemp;
                sFFT.pfReal[iA] +=                   fRealTemp;
                sFFT.pfImag[iB]  = sFFT.pfImag[iA] - fImagTemp;
                sFFT.pfImag[iA] +=                   fImagTemp;

                iA = iA + 2 * iM;
            }
            iQ += iL;
        }
        iL /= 2;
        iM *= 2;
    }
}

/*

R4FFTdif : Radix 4, in place, complex in & out, decimation in frequency FFT
           algorithm.
Required : A tsCmplxRect array with data to transform, another tsCmplxRect
           array with data twiddle factors (sine & cosine table),
           base 4 log of N, N
*/

void R4FFTdif (tsCmplxRect sFFT, tsCmplxRect sTwidTable, tInt iLog4N, tInt iN)
{
    tInt
        iCnt1, iCnt2, iCnt3,
        iL,    iM,    iQ,
        iA,    iB,    iC,     iD;
    tFloat
        fRealA,   fRealB,    fRealC,    fRealD,
        fImagA,   fImagB,    fImagC,    fImagD,
        fReal_Wq, fReal_W2q, fReal_W3q,
        fImag_Wq, fImag_W2q, fImag_W3q;

    iL = 1;
    iM = iN / 4;

    for (iCnt1 = 0; iCnt1 < iLog4N; ++iCnt1)
    {
        iQ = 0;
        for (iCnt2 = 0; iCnt2 < iM; ++iCnt2)
        {
            iA = iCnt2;
            fReal_Wq  = sTwidTable.pfReal[    iQ];
            fImag_Wq  = sTwidTable.pfImag[    iQ];
            fReal_W2q = sTwidTable.pfReal[2 * iQ];
            fImag_W2q = sTwidTable.pfImag[2 * iQ];
            fReal_W3q = sTwidTable.pfReal[3 * iQ];
            fImag_W3q = sTwidTable.pfImag[3 * iQ];

            for (iCnt3 = 0; iCnt3 < iL; ++iCnt3)
            {
                iB = iA +     iM;
                iC = iA + 2 * iM;
                iD = iA + 3 * iM;

                /* Butterfly: 42 FOP, 12 FMUL, 30 FADD */

                fRealA = sFFT.pfReal[iA] + sFFT.pfReal[iB]
                       + sFFT.pfReal[iC] + sFFT.pfReal[iD];
                fImagA = sFFT.pfImag[iA] + sFFT.pfImag[iB]
                       + sFFT.pfImag[iC] + sFFT.pfImag[iD];
                fRealB = sFFT.pfReal[iA] + sFFT.pfImag[iB]
                       - sFFT.pfReal[iC] - sFFT.pfImag[iD];
                fImagB = sFFT.pfImag[iA] - sFFT.pfReal[iB]
                       - sFFT.pfImag[iC] + sFFT.pfReal[iD];
                fRealC = sFFT.pfReal[iA] - sFFT.pfReal[iB]
                       + sFFT.pfReal[iC] - sFFT.pfReal[iD];
                fImagC = sFFT.pfImag[iA] - sFFT.pfImag[iB]
                       + sFFT.pfImag[iC] - sFFT.pfImag[iD];
                fRealD = sFFT.pfReal[iA] - sFFT.pfImag[iB]
                       - sFFT.pfReal[iC] + sFFT.pfImag[iD];
                fImagD = sFFT.pfImag[iA] + sFFT.pfReal[iB]
                       - sFFT.pfImag[iC] - sFFT.pfReal[iD];

                sFFT.pfReal [iA] = fRealA;
                sFFT.pfImag [iA] = fImagA;
                sFFT.pfReal [iB] = fRealB * fReal_Wq  - fImagB * fImag_Wq;
                sFFT.pfImag [iB] = fRealB * fImag_Wq  + fImagB * fReal_Wq;
                sFFT.pfReal [iC] = fRealC * fReal_W2q - fImagC * fImag_W2q;
                sFFT.pfImag [iC] = fRealC * fImag_W2q + fImagC * fReal_W2q;
                sFFT.pfReal [iD] = fRealD * fReal_W3q - fImagD * fImag_W3q;
                sFFT.pfImag [iD] = fRealD * fImag_W3q + fImagD * fReal_W3q;

                iA = iA + 4 * iM;
            }
            iQ += iL;
        }
        iL *= 4;
        iM /= 4;
    }
}

/*

R4FFTdit : Radix 4, in place, complex in & out, decimation in time FFT
           algorithm.
Required : A tsCmplxRect array with data to transform, another tsCmplxRect
           array with data twiddle factors (sine & cosine table),
           base 4 log of N, N
*/

void R4FFTdit (tsCmplxRect sFFT, tsCmplxRect sTwidTable, tInt iLog4N, tInt iN)
{
    tInt
        iCnt1, iCnt2, iCnt3,
        iL,    iM,    iQ,
        iA,    iB,    iC,     iD;
    tFloat
        fRealA,   fRealB,    fRealC,    fRealD,
        fImagA,   fImagB,    fImagC,    fImagD,
        fReal_Wq, fReal_W2q, fReal_W3q,
        fImag_Wq, fImag_W2q, fImag_W3q;

    iL = iN / 4;
    iM = 1;

    for (iCnt1 = 0; iCnt1 < iLog4N; ++iCnt1)
    {
        iQ = 0;
        for (iCnt2 = 0; iCnt2 < iM; ++iCnt2)
        {
            iA = iCnt2;
            fReal_Wq  = sTwidTable.pfReal[    iQ];
            fImag_Wq  = sTwidTable.pfImag[    iQ];
            fReal_W2q = sTwidTable.pfReal[2 * iQ];
            fImag_W2q = sTwidTable.pfImag[2 * iQ];
            fReal_W3q = sTwidTable.pfReal[3 * iQ];
            fImag_W3q = sTwidTable.pfImag[3 * iQ];

            for (iCnt3 = 0; iCnt3 < iL; ++iCnt3)
            {
                iB = iA +     iM;
                iC = iA + 2 * iM;
                iD = iA + 3 * iM;

                /* Butterfly: 42 FOP, 12 FMUL, 30 FADD */

                fRealA = sFFT.pfReal[iA];
                fImagA = sFFT.pfImag[iA];
                fRealB = sFFT.pfReal[iB] * fReal_Wq
                       - sFFT.pfImag[iB] * fImag_Wq;
                fImagB = sFFT.pfReal[iB] * fImag_Wq
                       + sFFT.pfImag[iB] * fReal_Wq;
                fRealC = sFFT.pfReal[iC] * fReal_W2q
                       - sFFT.pfImag[iC] * fImag_W2q;
                fImagC = sFFT.pfReal[iC] * fImag_W2q
                       + sFFT.pfImag[iC] * fReal_W2q;
                fRealD = sFFT.pfReal[iD] * fReal_W3q
                       - sFFT.pfImag[iD] * fImag_W3q;
                fImagD = sFFT.pfReal[iD] * fImag_W3q
                       + sFFT.pfImag[iD] * fReal_W3q;

                sFFT.pfReal[iA] = fRealA + fRealB + fRealC + fRealD;
                sFFT.pfImag[iA] = fImagA + fImagB + fImagC + fImagD;
                sFFT.pfReal[iB] = fRealA + fImagB - fRealC - fImagD;
                sFFT.pfImag[iB] = fImagA - fRealB - fImagC + fRealD;
                sFFT.pfReal[iC] = fRealA - fRealB + fRealC - fRealD;
                sFFT.pfImag[iC] = fImagA - fImagB + fImagC - fImagD;
                sFFT.pfReal[iD] = fRealA - fImagB - fRealC + fImagD;
                sFFT.pfImag[iD] = fImagA + fRealB - fImagC - fRealD;

                iA = iA + 4 * iM;
            }
            iQ += iL;
        }
        iL /= 4;
        iM *= 4;
    }
}

/*
Magnitude : Point - point modulus of the tsCmplxRect array.
*/

void Magnitude (tsCmplxRect sIn, ptFloatPt pfOut, tInt iN)
{
    tInt iCnt1;

    for (iCnt1 = 0; iCnt1 < iN; ++iCnt1)
    {
        pfOut [iCnt1] = sqrt (pow (sIn.pfReal [iCnt1], 2)
                            + pow (sIn.pfImag [iCnt1], 2));
    }
}

void Magnitude2 (tsCmplxRect sIn, ptFloatPt pfOut, tInt iN)
{
    tInt iCnt1;

    for (iCnt1 = 0; iCnt1 < iN; ++iCnt1)
    {
        pfOut [iCnt1] = pow (sIn.pfReal [iCnt1], 2)
                      + pow (sIn.pfImag [iCnt1], 2);
    }
}

/*
RectToPol : Convert from rectangular coord to polar.
*/

void RectToPol (tsCmplxRect sIn, tsCmplxPol sOut, tInt iN)
{
    tInt iCnt1;

    for (iCnt1 = 0; iCnt1 < iN; ++iCnt1)
    {
        sOut.pfMagn [iCnt1]  = sqrt (pow (sIn.pfReal [iCnt1], 2)
                                   + pow (sIn.pfImag [iCnt1], 2));
        sOut.pfPhase [iCnt1] = atan2 (sIn.pfImag [iCnt1], sIn.pfReal [iCnt1]);
    }
}

/*
Sign
*/

tInt Sign (tFloat fIn)
{
    if (fIn > 0)
        return 1;
    else
        return -1;
}

/*
Scaling : multiply tsCmplxRect with a constant
*/

void Scaling (tsCmplxRect sFFT, tFloat fCste, tInt iN)
{
    tInt iCnt1;

    for (iCnt1 = 0; iCnt1 < iN; ++iCnt1)
    {
        sFFT.pfReal [iCnt1] = sFFT.pfReal [iCnt1] * fCste;
        sFFT.pfImag [iCnt1] = sFFT.pfImag [iCnt1] * fCste; 
    }   
}

/*
AddCmplxRect : A' = A + B;
Ex.     : A=[0,2,3,1,3,5,4],
          B=[2,3,1,9],
          BOffset=3 =>
          A'=[0,2,3,3,6,6,13].
NB. if sA.iN < sB.iN + iBOffset, result is truncated.
*/

void AddCmplxRect (tsCmplxRect sA, tsCmplxRect sB, tInt iBOffset, tInt iNA,
                   tInt iNB)
{
    tInt iCnt1, iB = 0;

    for (iCnt1 = iBOffset; iCnt1 < iNA; ++iCnt1)
    {
        sA.pfReal [iCnt1] += sB.pfReal [iB];
        sA.pfImag [iCnt1] += sB.pfImag [iB];
        ++iB;
    }
}

void MultCmplxRect (tsCmplxRect sA, tsCmplxRect sB, tInt iBOffset, tInt iNA,
                    tInt iNB)
{
    tInt iCnt1, iB = 0;
    tFloat fRealA, fImagA;

    for (iCnt1 = iBOffset; iCnt1 < iNA; ++iCnt1)
    {
        fRealA = sA.pfReal [iCnt1];
        fImagA = sA.pfImag [iCnt1];
        sA.pfReal [iCnt1] = fRealA * sB.pfReal [iB] - fImagA * sB.pfImag [iB];
        sA.pfImag [iCnt1] = fRealA * sB.pfImag [iB] + fImagA * sB.pfReal [iB];
        ++iB;
    }
}

/*
BlockConvolve : compute the convolution y(t) = x(t) * h(t)
                Return defErrArraySize if the length of Y != iNX + iNH - 1
                STILL UNTESTED
*/

tError BlockConvolve (ptFloatPt pfX, ptFloatPt pfH, ptFloatPt pfY, tInt iNX,
                      tInt iNH, tInt iNY)
{
    tInt iOffset, iIndex, iUBound, iLBound;
    tFloat fSum;

    if ((iNX + iNH - 1) != (iNY))
        return ARRAY_SIZE;

    fSum = 0.0;
  
    for (iOffset = 0; iOffset < iNY; ++iOffset)
    {
        if (iOffset < iNH)
            iUBound = iOffset;
        else
            iUBound = iNH;

        if (iOffset < iNX)
            iLBound = 0;
        else
            iLBound = iOffset - iNX;

        for (iIndex = iUBound; iIndex >= iLBound; --iIndex)
            fSum += pfX [iOffset - iIndex] * pfH [iIndex];

        pfY [iOffset] = fSum;
        fSum = 0.0;
    }
    return OK;
}
















