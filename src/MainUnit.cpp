/*
Edited from \GuitSyn\1 :
For 22 fret guitar, the lowest note, E, is 82.40688922822Hz(Windows Calc).
The highest note, D, is 1174.65907167Hz.

To allow for bends I will plan for low D at 73.41619197935Hz,
  and high E at 1318.510227651Hz.
That's a span of 1245.094035672Hz.
The difference between low E and low D is 8.99069724887Hz.
In a most precise implementation (my first) the peak detector needs
  one more frequency below and above those limits.
That was way too picky. If the guitar notes were slightly out of tune,
  the detector didn't play the note.
To allow for imprecision could be difficult.
Detector note bandwidth will have to be introduced, and it will
  depend on frequency, maybe using a table.
But the problem is, with such low note frequencies and
  increments, the low (card's minimum) sample rate and larger buffer
  mean higher latency. -- Actually the sample rate is not the problem...
Specifically, this simple relation was found: delay = 1 / (f1 - f0).
So for the 9 Hz difference at low E to low D, this unfortunately
means a delay of 1/9 second (111 ms).
*/

//---------------------------------------------------------------------------
#include <vcl\vcl.h>
#pragma hdrstop

//#include "WaveUnit.h"
#include "MainUnit.h"
#pragma link "Grids"
//#pragma link "sampreg"
//#pragma link "CSPIN"

//#pragma link "zPanel"
//#pragma link "zPanel2"

#pragma link "RollEdit"
#pragma link "knob"
#pragma resource "*.dfm"

TMainForm *MainForm;

// R2_ORDER = 8 for 256 samples, 9 for 512 samples, 10 for 1024 samples, etc.
int R2_ORDER = 10;
int WAVEIN_BUFSIZE = 1 << R2_ORDER;
int FREQUENCY_BINS = WAVEIN_BUFSIZE >> 1;
int WAVEIN_MAXBUFS = 64;
int WAVEIN_SAMPRATE = 11025;
// July 4, 2002 Test: To make it think it's sampling at a lower frequency than it is.
double WAVEIN_SAMPRATE_FACTOR = 1.0;
double WAVEIN_SAMPRATE_OFFSET = 0.0;

bool UseGraph;
bool GraphDone = true;
DWORD GraphAvgCount = 0;
tFloat *GraphAvgVals = NULL;
bool GraphUsePeakElseAvg = true;
double SenseLevelExp;
tFloat GraphYScale;


TDSPData DSPData;
TDSPControl DSPControl;

HMIDIIN CurrMidiInHandle;
HMIDIOUT CurrMidiOutHandle;
UINT NumMidiInDevs;
UINT NumMidiOutDevs;
MIDIINCAPS MidiInCapArray[MAX_MIDIIN_DEVS];
MIDIOUTCAPS MidiOutCapArray[MAX_MIDIOUT_DEVS];
TStringList *MidiInDevNames = NULL;
TStringList *MidiOutDevNames = NULL;
BYTE KeyToNoteTable[256];
BYTE KeyRepeatTable[256];
int CurrPlayPrgs[MAX_PLAY_CHN];
int CurrPlayOct = 1;
int CurrPlayChn;
int  *NoteToSpecMap = NULL;
BYTE *SpecToNoteMap = NULL;
BYTE *NoteRepeatTable = NULL;
int   TotalNotesPlaying;

CRITICAL_SECTION CriticalSection;
HWAVEIN CurrWaveInHandle;
//WaveInProcAddress;
int NumWaveInDevs;
PWAVEFORMATEX pWaveInFormat;
HGLOBAL  hWaveInFormat;
WAVEINCAPS WaveInCapArray[MAX_WAVEIN_DEVS];

//PWAVEHDR pWaveInHeaders[WAVEIN_MAXBUFS];
PWAVEHDR *pWaveInHeaders = NULL;
//HGLOBAL hWaveInHeaders[WAVEIN_MAXBUFS];
HGLOBAL *hWaveInHeaders = NULL;
//HGLOBAL hWaveInBuffers[WAVEIN_MAXBUFS];
HGLOBAL *hWaveInBuffers = NULL;
//SAMPLE_TYPE *pWaveInBuffers[WAVEIN_MAXBUFS];
SAMPLE_TYPE **pWaveInBuffers = NULL;

int WaveInBufIndex;
int WaveInBufPend;
int WaveInCurrBufs;
int WaveInBufsMissed;
//UINT PrivMsgNumberIn;
//bool WaveInPrivMsgRdy;

tsCmplxRect        dspIn, dspOut, dspTwidTable;
ptFloatPt          pfdspOutMag;
ptInt              pidspR2SwapTable;

TStringList *WaveInDevNames = NULL;
int StopWaveInFlag;
int Recording;

/*const tInt         dspiN         =  512; // N Samples
tsCmplxRect        dspIn, dspOut, dspTwidTable;
ptFloat            pfdspOutMag;
ptInt              pidspR2SwapTable; */


//---------------------------------------------------------------------------
void DestroyMemory()
{
  DestroyWaveMemory();

  if(SpecToNoteMap)
    delete SpecToNoteMap;
  SpecToNoteMap = NULL;
  if(GraphAvgVals)
    delete GraphAvgVals;
  GraphAvgVals = NULL;  

  if(dspIn.pfReal)
    delete dspIn.pfReal;
  dspIn.pfReal = NULL;
  if(dspIn.pfImag)
    delete dspIn.pfImag;
  dspIn.pfImag = NULL;
  if(dspOut.pfReal)
    delete dspOut.pfReal;
  dspOut.pfReal = NULL;
  if(dspOut.pfImag)
    delete dspOut.pfImag;
  dspOut.pfImag = NULL;
  if(pfdspOutMag)
    delete pfdspOutMag;
  pfdspOutMag = NULL;
  if(dspTwidTable.pfReal)
    delete dspTwidTable.pfReal;
  dspTwidTable.pfReal = NULL;
  if(dspTwidTable.pfImag)
    delete dspTwidTable.pfImag;
  dspTwidTable.pfImag = NULL;
  if(pidspR2SwapTable)
    delete pidspR2SwapTable;
  pidspR2SwapTable = NULL;
}
//---------------------------------------------------------------------------
bool InitMemory()
{
  //DestroyMemory();

  dspIn.pfReal        = new tFloat[WAVEIN_BUFSIZE];
  dspIn.pfImag        = new tFloat[WAVEIN_BUFSIZE];
  dspOut.pfReal       = new tFloat[WAVEIN_BUFSIZE];
  dspOut.pfImag       = new tFloat[WAVEIN_BUFSIZE];
  pfdspOutMag         = new tFloat[WAVEIN_BUFSIZE];
  dspTwidTable.pfReal = new tFloat[WAVEIN_BUFSIZE];
  dspTwidTable.pfImag = new tFloat[WAVEIN_BUFSIZE];
  pidspR2SwapTable    = new tInt[WAVEIN_BUFSIZE];

  SpecToNoteMap = new BYTE[WAVEIN_BUFSIZE];//really only half is used,but be safe
  GraphAvgVals = new tFloat[WAVEIN_BUFSIZE];

  AllocateWaveMemory();
  return true;
}
//---------------------------------------------------------------------------
void CALLBACK WaveInFunc(HWAVEIN InHandle, UINT Msg,
    DWORD dwInstance, DWORD dwParam1, DWORD dwParam2)
{
  LARGE_INTEGER PerfFreq;
  LARGE_INTEGER PerfCnt;
  float ElapsedTime;
  float Freq;

//  if(QueryPerformanceCounter(&PerfCnt))
//  {
    // This line will only be executed once...
//    static LARGE_INTEGER LastPerfCnt = PerfCnt;
    // Get the frequency in counts per second
//    QueryPerformanceFrequency(&PerfFreq);
    // Calculate the elapsed time in seconds
    // This part ONLY runs on 64-bit processors.
    // Rewrite to deal with split LONGLONG later...
//    ElapsedTime = float(PerfCnt.QuadPart - LastPerfCnt.QuadPart) / float(PerfFreq.QuadPart);
//    Freq = float(PerfFreq.QuadPart) / float(PerfCnt.QuadPart - LastPerfCnt.QuadPart);
//    LastPerfCnt = PerfCnt;
//    MainForm->Label20->Caption = String(Freq);
//  }
//  else
    // Get the relative time, in seconds
//    ElapsedTime = (float(timeGetTime()) - GlobalTime) * 0.001f;

  static FLOAT fFPS      = 0.0f;
  static FLOAT fLastTime = 0.0f;
  static DWORD dwFrames  = 0L;

  // Keep track of the time lapse and frame count
  FLOAT fTime = timeGetTime() * 0.001f; // Get current time in seconds
  ++dwFrames;

  // Update the frame rate once per second
  if( fTime - fLastTime > 4.00f )
  {
      fFPS      = dwFrames / (fTime - fLastTime);
      fLastTime = fTime;
      dwFrames  = 0L;
  }

  MainForm->Label20->Caption = String(fFPS);


  int i, j;
//  int k, l, m;
  bool Found;
  LPWAVEHDR HdrPtr;
  int BufNum;
  //tFloat SmpTmp;
  //tFloat ScaledSenseLevel;
  //BYTE NoteTmp;
  DWORD VolTmp;
  const int HalfBufSizeM1 = (WAVEIN_BUFSIZE >> 1) -1;

  switch(Msg)
  {
    case
    WIM_DATA :
    {
      /*
      if(WaveInPrivMsgRdy)
      {
        WaveInPrivMsgRdy = false;
        PostMessage(Application->Handle,
                  PrivMsgNumberIn,
                  WaveInBufIndex, 0);
      }
      else
        WaveInBufsMissed++;
      */

     if(StopWaveInFlag == 0)
     {
      EnterCriticalSection(&CriticalSection);
      HdrPtr = (LPWAVEHDR) dwParam1;
      BufNum = (int) HdrPtr->dwUser;
      //if(HdrPtr->dwBytesRecorded
      for(i = 0; i < WAVEIN_BUFSIZE; i++)
      {
        if(sizeof(SAMPLE_TYPE) == 2)
          dspIn.pfReal[i] =
            (((SAMPLE_TYPE)pWaveInBuffers[BufNum][i]) / 65536.0);
        else
          dspIn.pfReal[i] =
            (((BYTE)pWaveInBuffers[BufNum][i] - 128) / 256.0);
        dspIn.pfImag[i] = 0.0;
      }
      waveInAddBuffer(InHandle, (LPWAVEHDR)dwParam1, sizeof(WAVEHDR));


      //for(i = 0; i < WAVEIN_BUFSIZE; i++)

      R2FFTdif(dspIn, dspTwidTable, R2_ORDER, WAVEIN_BUFSIZE); // 2^9 = 512
      SwapFFTC(dspIn, pidspR2SwapTable, WAVEIN_BUFSIZE);
      Magnitude(dspIn, pfdspOutMag, WAVEIN_BUFSIZE);

      for(i = 0; i < FREQUENCY_BINS; i++)
        pfdspOutMag[i] =
          pow(pfdspOutMag[i] * (DSPControl.GainPerFreq * (float)i +1.0), SenseLevelExp);


      // 09/29/01
      // Disabled for now, because Zeigler collection removed - sharware.
      // Must make our own spectral display.
      /*
      if(UseGraph)
      {
        if(!GraphDone)
        {
          for(i = 0; i < FREQUENCY_BINS; i++)
          {
            if(GraphUsePeakElseAvg)
            {
              if(pfdspOutMag[i] > GraphAvgVals[i])
                GraphAvgVals[i] = pfdspOutMag[i];
            }
            else
              GraphAvgVals[i] += pfdspOutMag[i];
          }
          // Increment the counter
          GraphAvgCount++;
        }
        else
        {
          GraphAvgCount++;
  //        for(i = 0; i < WAVEIN_BUFSIZE; i++)
  //          GraphAvgVals[i] = ((tFloat)GraphAvgVals[i] + pfdspOutMag[i]) / (tFloat)GraphAvgCount;
          GraphDone = false;
  //        MainForm->zcSpectrum1->SetAllChannel(GraphAvgVals);

          int Cnt = MainForm->zcSpectrum1->NumberOfChannels;
          if(Cnt > FREQUENCY_BINS)
            Cnt = FREQUENCY_BINS;
          for(i = 0; i < Cnt; i++)
          {
            if(GraphUsePeakElseAvg)
            {
              if(pfdspOutMag[i] > GraphAvgVals[i])
                GraphAvgVals[i] = pfdspOutMag[i];
              //MainForm->zcSpectrum1->SetChannel(i, (int)(GraphYScale * pow(GraphAvgVals[i], SenseLevelExp)));
              MainForm->zcSpectrum1->SetChannel(i, (int)(GraphYScale * GraphAvgVals[i]));
            }
            else
              //MainForm->zcSpectrum1->SetChannel(i, (int)(GraphYScale * pow(((GraphAvgVals[i] + pfdspOutMag[i]) / (tFloat)GraphAvgCount), SenseLevelExp)));
              MainForm->zcSpectrum1->SetChannel(i, (int)(GraphYScale * ((GraphAvgVals[i] + pfdspOutMag[i]) / (tFloat)GraphAvgCount)));
            GraphAvgVals[i] = 0.0f;
          }
          // Reset the counter
          GraphAvgCount = 0;
          GraphDone = true;
        }
      }
      */


      LeaveCriticalSection(&CriticalSection);

      /*
      for(i = 1; i < HalfBufSizeM1; i++)
      {
        NoteTmp = SpecToNoteMap[i] & 0x7F; // Be safe by and-ing
        if(pfdspOutMag[i] >= ScaledSenseLevel)
        {
          if( ((pfdspOutMag[i + 1] - pfdspOutMag[  i  ]) < 0.0) &&
              ((pfdspOutMag[  i  ] - pfdspOutMag[i - 1]) > 0.0) )
          {
            if(NoteRepeatTable[NoteTmp] == 0)
            {
              NoteRepeatTable[NoteTmp] = 1;
              midiOutShortMsg(CurrMidiOutHandle,
                   0x00600090 |
                   (((DSPControl.NoteOffset + NoteTmp) & 0x7F) << 8) |
                   (CurrPlayChn & 0x0F));
            }
          }
          else
          {
            if(NoteRepeatTable[NoteTmp] == 1)
            {
              NoteRepeatTable[NoteTmp] = 0;
              midiOutShortMsg(CurrMidiOutHandle,
                   0x00600080 |
                   (((DSPControl.NoteOffset + NoteTmp) & 0x7F) << 8) |
                   (CurrPlayChn & 0x0F));
            }
          }
        }
        else
        {
          if(NoteRepeatTable[NoteTmp] == 1)
          {
            NoteRepeatTable[NoteTmp] = 0;
            midiOutShortMsg(CurrMidiOutHandle,
                 0x00600080 |
                 (((DSPControl.NoteOffset + NoteTmp) & 0x7F) << 8) |
                 (CurrPlayChn & 0x0F));
          }
        }
      }   */


      // Don't do first and last notes because of +1/-1 logic below...
      for(i = 1; i < MAX_POLY_NOTES -1; i++)
      {
        // Get the 'frequency bin' index (j) which corresponds to this note (i).
        j = NoteToSpecMap[i];
        // If the bin is marked as not used (because it's too high - out of range), just continue.
//        if(j == -1)
//          continue;

        //if((j > 0) && (j < WAVEIN_BUFSIZE -1))    PUT IN INIT ROUTINE!

//        l = (NoteToSpecMap[i] -
//            (NoteToSpecMap[i] - NoteToSpecMap[i - 1]) /2) + 1;
//        m = (NoteToSpecMap[i] +
//            (NoteToSpecMap[i + 1] - NoteToSpecMap[i]) /2) - 1;
        Found = false;
//        for(k = l; k < m; k++) // check this
        if(j != -1)
        {

        //}

//        if(pfdspOutMag[k] > DSPControl.ScaledSenseLevel)
        if(pow(pfdspOutMag[j], SenseLevelExp) >= DSPControl.ScaledSenseLevel)
        {
//          if( ((pfdspOutMag[k + 1] - pfdspOutMag[  k  ]) < 0.0) &&
//              ((pfdspOutMag[  k  ] - pfdspOutMag[k - 1]) > 0.0) )
//          if( ((pfdspOutMag[j + 1] - pfdspOutMag[  j  ]) < 0.0) &&
//              ((pfdspOutMag[  j  ] - pfdspOutMag[j - 1]) > 0.0) )
          {
            Found = true;
            if(NoteRepeatTable[i] == 0)
            {
              NoteRepeatTable[i] = 1;
              // hiword
//              VolTmp = (DWORD)(pfdspOutMag[k] / DSPControl.ScaledSenseLevel);
              VolTmp = (DWORD)(pow(pfdspOutMag[j], 2) - DSPControl.ScaledSenseLevel);
              if(VolTmp > 0x7F)
                VolTmp = 0x7F;
              VolTmp = (VolTmp << 16) & 0x007F0000;

              midiOutShortMsg(CurrMidiOutHandle,
                0x00600090 | VolTmp |
                 (((DSPControl.NoteOffset + i + CurrPlayOct * 12) & 0x7F) << 8) |
                   (CurrPlayChn & 0x0F));
//              break;
//              continue;
            }
          }
          /*
          else
          {

            if(NoteRepeatTable[i] == 1)
            {
              NoteRepeatTable[i] = 0;
              midiOutShortMsg(CurrMidiOutHandle,
                0x00600080 | (((DSPControl.NoteOffset + i) & 0x7F) << 8) | (CurrPlayChn & 0x0F));
            }
          }
          */
        }
        /*
        else
        {
          if(NoteRepeatTable[i] == 1)
          {
            NoteRepeatTable[i] = 0;
            midiOutShortMsg(CurrMidiOutHandle,
              0x00600080 | (((DSPControl.NoteOffset + i) & 0x7F) << 8) | (CurrPlayChn & 0x0F));
          }
        }
        */

        }
        if(!Found)
        {
          if(NoteRepeatTable[i] == 1)
          {
            NoteRepeatTable[i] = 0;
            midiOutShortMsg(CurrMidiOutHandle,
              0x00000080 | (((DSPControl.NoteOffset + i + CurrPlayOct * 12) & 0x7F) << 8) | (CurrPlayChn & 0x0F));
            //  ..vvnn..

          }
        }   // ..from for

      } // ..for
     } // if(StopWaveInFlag...
    }
  }
}

      /*
        (cos((FLOAT_TYPE)i * SinGen1Freq * TwoPiTScale) *
	  					SinGen1Amp + SinGen1Bias) +
	  				 (cos((FLOAT_TYPE)i * SinGen2Freq * TwoPiTScale) *
	  					SinGen2Amp + SinGen2Bias);
	    //dspIn.pfImag[i] = 0.0;
        */


//---------------------------------------------------------------------------
/*
void __fastcall TMainForm::MyAppOnMessage(MSG &Msg, bool &Handled)
{
  int i, ii;

  if(Msg.message == PrivMsgNumberIn)
  {

      Handled = true;
      if(StopWaveInFlag == 0)
      {
        ii = WaveInBufsMissed;
        WaveInBufsMissed = 0;
        WaveInPrivMsgRdy = true; // Notify the callback
        // Limit the number of bufs to send
        if(ii > WAVEIN_MAXBUFS)
             ii = WAVEIN_MAXBUFS;

        Loop1:

        waveInAddBuffer(CurrWaveInHandle,
          pWaveInHeaders[WaveInBufPend], sizeof(WAVEHDR));
        WaveInBufPend++;
        if(WaveInBufPend >= WAVEIN_MAXBUFS)
          WaveInBufPend = WaveInBufPend - WAVEIN_MAXBUFS;
        if(ii != 0)
        {
          ii--;
          goto Loop1;  // To Loop1 >>>>>>
        }
      }

  }
}
*/

//---------------------------------------------------------------------------
__fastcall TMainForm::TMainForm(TComponent* Owner)
	: TForm(Owner)
{
}
//---------------------------------------------------------------------------
void __fastcall TMainForm::FormCreate(TObject *Sender)
{
  dspIn.pfReal        = NULL;
  dspIn.pfImag        = NULL;
  dspOut.pfReal       = NULL;
  dspOut.pfImag       = NULL;
  pfdspOutMag         = NULL;
  dspTwidTable.pfReal = NULL;
  dspTwidTable.pfImag = NULL;
  pidspR2SwapTable    = NULL;

  NoteToSpecMap       = NULL;
  SpecToNoteMap       = NULL;
  NoteRepeatTable     = NULL;
  MidiInDevNames      = NULL;
  MidiOutDevNames     = NULL;
  GraphAvgVals        = NULL;

  int i;

  Recording = 0;
  StopWaveInFlag = 1;
  WaveInBufsMissed = 0;
  WaveInBufPend = 0;
  WaveInBufIndex = 0;

  InitializeCriticalSection(&CriticalSection);

  /*
  PrivMsgNumberIn = RegisterWindowMessage(PRIV_MSG_NAME_IN);
  WaveInPrivMsgRdy = False;
  Application->OnMessage = MyAppOnMessage;
  */

//  dspIn.pfReal        = new tFloat[WAVEIN_BUFSIZE];
//  dspIn.pfImag        = new tFloat[WAVEIN_BUFSIZE];
//  dspOut.pfReal       = new tFloat[WAVEIN_BUFSIZE];
//  dspOut.pfImag       = new tFloat[WAVEIN_BUFSIZE];
//  pfdspOutMag         = new tFloat[WAVEIN_BUFSIZE];
//  dspTwidTable.pfReal = new tFloat[WAVEIN_BUFSIZE];
//  dspTwidTable.pfImag = new tFloat[WAVEIN_BUFSIZE];
//  pidspR2SwapTable    = new tInt[WAVEIN_BUFSIZE];
  InitMemory();

  InitTwidTable(dspTwidTable, DIRECT, WAVEIN_BUFSIZE);
  InitR2SwapTable(pidspR2SwapTable, WAVEIN_BUFSIZE);
  for(i = 0; i < WAVEIN_BUFSIZE; i++)
  {
	 dspIn.pfReal[i] = 0.0;
	 dspIn.pfImag[i] = 0.0;
  };
                   /* (cos((FLOAT_TYPE)i * SinGen1Freq * TwoPiTScale) *
			 			SinGen1Amp + SinGen1Bias) +
			 			 (cos((FLOAT_TYPE)i * SinGen2Freq * TwoPiTScale) *
			 				SinGen2Amp + SinGen2Bias);    */

  int Val;

  NoteToSpecMap = new int[MAX_POLY_NOTES];
//  SpecToNoteMap = new BYTE[WAVEIN_BUFSIZE];//really only half is used,but be safe
  NoteRepeatTable = new BYTE[MAX_POLY_NOTES];
  MidiInDevNames = new TStringList();
  MidiOutDevNames = new TStringList();
//  GraphAvgVals = new tFloat[WAVEIN_BUFSIZE];

  for(i = 0; i < WAVEIN_BUFSIZE; i++);
    GraphAvgVals[i] = 0.0f;

  ClearRepeatTables();

  // Do this first before init
  DSPControl.TuneOffset = TuneOffsetZKnob->Value;
  InitNoteToSpecMap();

  InitKeyToNoteTable();
  //InitSpecToNoteMap(); // Not used
  TotalNotesPlaying = 0; // also cleared in Turn Off All Notes button event

  QueryMidi();
  ComboBox1->Items->Assign(MidiInDevNames);
  ComboBox2->Items->Assign(MidiOutDevNames);

  for(i = 0; i < MAX_PLAY_CHN; i++)
    CurrPlayPrgs[i] = MIN_PLAY_PRG;

//  CurrPlayOct = (MAX_PLAY_OCT - MIN_PLAY_OCT) >> 1;
//  CurrPlayChn = MIN_PLAY_CHN;

  Val = MidiChannelRollEdit->Value -1;
  if((Val >= MIN_PLAY_CHN) || (Val <= MAX_PLAY_CHN))
    CurrPlayChn = Val;

//  MaskEdit1->Text = IntToStr(CurrPlayOct + 1);
//  MaskEdit3->Text = IntToStr(CurrPlayChn + 1);
//  MaskEdit2->Text = IntToStr(CurrPlayPrgs[CurrPlayChn] + 1);

  DSPData.PrevEdgeOffset = 0;
  DSPData.BufsToPrevEdge = 0;
  DSPData.SmpToPrevEdge = 0;
  DSPData.AvgPer = 0;
  DSPData.NoteFreq = 0;
  DSPData.PrevNoteOn = 0;
  DSPData.InSquared = 0;

  DSPControl.SenseLevel = SenseLevelZKnob->Value;
  DSPControl.ScaledSenseLevel =
    (float)(DSPControl.SenseLevel) / (float)(MAG_OUT_SCALE);

  DSPControl.NoteOffset = (int)NoteOffsetRollEdit->Value & 0xFF;

  SenseLevelExp = (double)SenseLevelExpZKnob->Value / 100.0;

  DSPControl.GainPerFreq = (float)FreqGainZKnob->Value / 10000.0f;

  GraphYScale = (tFloat)GraphYScaleZKnob->Value / 10.0f;

  if(UseGraphCheckBox->Checked)
    UseGraph = true;
  else
    UseGraph = false;

//  DSPControl.TuneOffset = TuneOffsetZKnob->Value;

  Val = (int)MidiOctaveRollEdit->Value;
  if((Val >= MIN_PLAY_OCT) || (Val <= MAX_PLAY_OCT))
    CurrPlayOct = Val;




  WaveInDevNames = new TStringList;
  WaveInDevNames->Sorted = False;
//  AllocateWaveMemory();
  SetWaveInFormat();
  InitWaveBufs();
  QueryWave();
  WaveDevComboBox->Items->Assign(WaveInDevNames);
  // Preset these to be safe. They're also set in StartWaveInRec
  WaveDevComboBox->Text = ""; //WaveInDevNames[0];
  /*OpenWaveIn(Handle,*/ /*WaveDevComboBox->ItemIndex*/ /*0);*/

  for(i = 0; i < 16; i++)
    StringGrid1->Cells[i + 1][0] = IntToHex(i, 1);
  for(i = 0; i < 16; i++)
    StringGrid1->Cells[0][i + 1] = IntToHex((i << 4) & 240, 2);
}
//---------------------------------------------------------------------------
void __fastcall TMainForm::FormDestroy(TObject *Sender)
{
  CloseWaveIn();
  CloseMidiOut();

  DeleteCriticalSection(&CriticalSection);
  //DestroyWaveMemory();
  DestroyMemory();
  delete WaveInDevNames;

//  delete GraphAvgVals;

  delete MidiInDevNames;
  delete MidiOutDevNames;

  delete NoteRepeatTable;
//  delete SpecToNoteMap;
  delete NoteToSpecMap;

//  delete pidspR2SwapTable;
//  delete dspTwidTable.pfImag;
//  delete dspTwidTable.pfReal;
//  delete pfdspOutMag;
//  delete dspOut.pfImag;
//  delete dspOut.pfReal;
//  delete dspIn.pfImag;
//  delete dspIn.pfReal;
}
//---------------------------------------------------------------------------
void __fastcall TMainForm::FormKeyDown(TObject *Sender, WORD &Key,
	TShiftState Shift)
{
  if(!StringGrid1->Focused())
    return;

  int i;

  if(KeyRepeatTable[Key] == 0)
  {
    KeyRepeatTable[Key] = 1;
    StringGrid1->Cells[(Key & 15) + 1][(Key >> 4) + 1] = 'X';

    if(KeyToNoteTable[Key] != 0xFF)
    {
      i = (KeyToNoteTable[Key] +
            CurrPlayOct * 12) & 0x7F;
      midiOutShortMsg(CurrMidiOutHandle,
        0x00400090 | (i << 8) | (CurrPlayChn & 0x0F));
    }
  }
  Key = 0;
}
//---------------------------------------------------------------------------
void __fastcall TMainForm::FormKeyUp(TObject *Sender, WORD &Key,
	TShiftState Shift)
{
  if(!StringGrid1->Focused())
    return;

  int i;
  StringGrid1->Cells[(Key & 15) + 1][(Key >> 4) + 1] = "";

  if(KeyToNoteTable[Key] != 0xFF)
  {
    i = (KeyToNoteTable[Key] +
            CurrPlayOct * 12) & 0x7F;
    midiOutShortMsg(CurrMidiOutHandle,
      0x00400080 | (i << 8) | (CurrPlayChn & 0x0F));
  }
  KeyRepeatTable[Key] = 0;
  Key = 0;
}
//---------------------------------------------------------------------------
void __fastcall TMainForm::NoteResetButtonClick(TObject *Sender)
{
  AnsiString s;

  s = Button2->Caption;
  Button2->Caption = "Wait...";
  midiOutReset(CurrMidiOutHandle);
  TotalNotesPlaying = 0;
  Button2->Caption = s;
  Button2->SetFocus();
}
//---------------------------------------------------------------------------
void __fastcall TMainForm::ComboBox2Change(TObject *Sender)
{
  int i;
  if(ComboBox2->ItemIndex != -1)
  {
    CloseMidiOut();
    OpenMidiOut(Handle, ComboBox2->ItemIndex);

//    MaskEdit1->Text = IntToStr(CurrPlayOct + 1);
//    MaskEdit3->Text = IntToStr(CurrPlayChn + 1);
//    MaskEdit2->Text = IntToStr(CurrPlayPrgs[CurrPlayChn] + 1);

    for(i = 0; i < MAX_PLAY_CHN; i++)
    {
      /*CurrPlayPrgs[i] := MIN_PLAY_PRG;*/
      midiOutShortMsg(CurrMidiOutHandle,
        0x000000C0 | (CurrPlayPrgs[i] << 8) | (i & 0x0F));
    }
    midiOutShortMsg(CurrMidiOutHandle,
      0x000000C0 | (CurrPlayPrgs[CurrPlayChn] << 8) | (CurrPlayChn & 0x0F));
  }
  Button2->SetFocus();
}
//---------------------------------------------------------------------------
void __fastcall TMainForm::WaveDevComboBoxChange(TObject *Sender)
{
  if(WaveDevComboBox->ItemIndex != -1)
  {
    /*waveInReset(CurrWaveInHandle);*/
    CloseWaveIn();
    OpenWaveIn(Handle, WaveDevComboBox->ItemIndex);
  }
}
//---------------------------------------------------------------------------
void __fastcall TMainForm::StartWaveButtonClick(TObject *Sender)
{
  if(Recording == 1)
    StopWaveInRec();
  StartWaveInRec();
}
//---------------------------------------------------------------------------
void __fastcall TMainForm::StopWaveButtonClick(TObject *Sender)
{
// Unconditional:
  StopWaveInRec();
}
//---------------------------------------------------------------------------


void __fastcall TMainForm::RadioButton1Click(TObject *Sender)
{
  midiOutShortMsg(CurrMidiOutHandle, 0x00007DB0 | (CurrPlayChn & 0x0F));
}
//---------------------------------------------------------------------------
void __fastcall TMainForm::RadioButton2Click(TObject *Sender)
{
  midiOutShortMsg(CurrMidiOutHandle, 0x00007CB0 | (CurrPlayChn & 0x0F));
}
//---------------------------------------------------------------------------
void QueryMidi()
{
  UINT i;
  NumMidiInDevs = midiInGetNumDevs();
  NumMidiOutDevs = midiOutGetNumDevs();
  for(i = 0; i < NumMidiInDevs; i++)
  {
    midiInGetDevCaps(i, &MidiInCapArray[i], sizeof(MIDIINCAPS));
    MidiInDevNames->Add(MidiInCapArray[i].szPname);
  }
  for(i = 0; i < NumMidiOutDevs; i++)
  {
    midiOutGetDevCaps(i, &MidiOutCapArray[i], sizeof(MIDIINCAPS));
    MidiOutDevNames->Add(MidiOutCapArray[i].szPname);
  }
}
//---------------------------------------------------------------------------
void OpenMidiOut(HWND WindowHandle, UINT ID)
{
  if(midiOutOpen(&CurrMidiOutHandle,
                 ID, (DWORD)(WindowHandle),
                 0, CALLBACK_WINDOW) != MMSYSERR_NOERROR)
     Application->MessageBox("Error opening midi out", "Error", MB_OK);
}
//---------------------------------------------------------------------------
void CloseMidiOut()
{
  midiOutClose(CurrMidiOutHandle);
}
//---------------------------------------------------------------------------
void InitKeyToNoteTable()
{
  int i;
  for(i = 0; i < 256; i++) KeyToNoteTable[(BYTE)(i)] = (BYTE)(0xFF);
KeyToNoteTable[(BYTE)('Z')] = 0;
KeyToNoteTable[(BYTE)('S')] = 1;
KeyToNoteTable[(BYTE)('X')] = 2;
KeyToNoteTable[(BYTE)('D')] = 3;
KeyToNoteTable[(BYTE)('C')] = 4;
KeyToNoteTable[(BYTE)('V')] = 5;
KeyToNoteTable[(BYTE)('G')] = 6;
KeyToNoteTable[(BYTE)('B')] = 7;
KeyToNoteTable[(BYTE)('H')] = 8;
KeyToNoteTable[(BYTE)('N')] = 9;
KeyToNoteTable[(BYTE)('J')] = 10;
KeyToNoteTable[(BYTE)('M')] = 11;
KeyToNoteTable[(BYTE)(0xBC)] = 12;
KeyToNoteTable[(BYTE)('L')] = 13;
KeyToNoteTable[(BYTE)(0xBE)] = 14;
KeyToNoteTable[(BYTE)(0xBA)] = 15;
KeyToNoteTable[(BYTE)(0xBF)] = 16;
KeyToNoteTable[(BYTE)(0x10)] = 17;
KeyToNoteTable[(BYTE)('Q')] = 12;
KeyToNoteTable[(BYTE)('2')] = 13;
KeyToNoteTable[(BYTE)('W')] = 14;
KeyToNoteTable[(BYTE)('3')] = 15;
KeyToNoteTable[(BYTE)('E')] = 16;
KeyToNoteTable[(BYTE)('R')] = 17;
KeyToNoteTable[(BYTE)('5')] = 18;
KeyToNoteTable[(BYTE)('T')] = 19;
KeyToNoteTable[(BYTE)('6')] = 20;
KeyToNoteTable[(BYTE)('Y')] = 21;
KeyToNoteTable[(BYTE)('7')] = 22;
KeyToNoteTable[(BYTE)('U')] = 23;
KeyToNoteTable[(BYTE)('I')] = 24;
KeyToNoteTable[(BYTE)('9')] = 25;
KeyToNoteTable[(BYTE)('O')] = 26;
KeyToNoteTable[(BYTE)('0')] = 27;
KeyToNoteTable[(BYTE)('P')] = 28;
KeyToNoteTable[(BYTE)(0xDB)] = 29;
KeyToNoteTable[(BYTE)(0xBB)] = 30;
KeyToNoteTable[(BYTE)(0xDD)] = 31;
KeyToNoteTable[(BYTE)(0x08)] = 32;
}
//---------------------------------------------------------------------------
void ClearRepeatTables()
{
  int i;

  for(i = 0; i < 256; i++)
    KeyRepeatTable[(BYTE)(i)] = (BYTE)(0x00);
  for(i = 0; i < MAX_POLY_NOTES; i++)
    NoteRepeatTable[i] = 0;
}
//---------------------------------------------------------------------------
/*   Function not used
void InitSpecToNoteMap()
{
  const tFloat EqTempStep = 1.059463094359;
  const tFloat SpecFreqStep = (tFloat)(WAVEIN_SAMPRATE) /
                              (tFloat)(WAVEIN_BUFSIZE);
  int i;

  // Handle 0th element to prevent error in log() function below:
  //SpecToNoteMap[0] = 0;


  for(i = 1; i < (WAVEIN_BUFSIZE >> 1); i++)
  {
    SpecToNoteMap[i] =  log((tFloat)(i) * (tFloat)(DETECTOR_BASE_F) /
                          SpecFreqStep) /
                         log(EqTempStep) ;
  }
}
*/
//---------------------------------------------------------------------------
void InitNoteToSpecMap()
{
  int i;
//  const tFloat EqTempStep = 1.059463094359;
//  const tFloat SpecFreqStep = (tFloat)(WAVEIN_SAMPRATE) / (tFloat)(WAVEIN_BUFSIZE);
  const double EqTempStep = pow(2.0, 1.0/12.0);
  double SpecFreqStep =
    // July 4, 2002 Test: Added factor and offset to make it think it's sampling
    //  at a lower frequency than it is.
    ((double)WAVEIN_SAMPRATE * WAVEIN_SAMPRATE_FACTOR +
     (((WAVEIN_SAMPRATE+WAVEIN_SAMPRATE_OFFSET)>0.0)?WAVEIN_SAMPRATE_OFFSET:0.0)) /
    (double)WAVEIN_BUFSIZE;
  // Safety.
  if(SpecFreqStep < 1.0e-5)
    SpecFreqStep = 1.0e-5;

  for(i = 0; i < MAX_POLY_NOTES; i++)
  {
    NoteToSpecMap[i] =
      int( ( pow(EqTempStep, (double)i) *
             ((double)DETECTOR_BASE_F + (double)DSPControl.TuneOffset / 100.0) ) / SpecFreqStep );
 // Limit the result, to be safe!
 // No... mark as unused (-1)... we don't want a frequency present in this bin to
 //  trigger this and all the rest of the notes!!
    if(NoteToSpecMap[i] >= (FREQUENCY_BINS))
//       NoteToSpecMap[i] = (WAVEIN_BUFSIZE /2) -1;
       NoteToSpecMap[i] = -1;
  }
}
//---------------------------------------------------------------------------
void AllocateWaveMemory()
{
  //DestroyWaveMemory();

  pWaveInHeaders = new PWAVEHDR[WAVEIN_MAXBUFS];
  hWaveInHeaders = new HGLOBAL[WAVEIN_MAXBUFS];
  hWaveInBuffers = new HGLOBAL[WAVEIN_MAXBUFS];
  pWaveInBuffers = new SAMPLE_TYPE*[WAVEIN_MAXBUFS];

  int i;
  hWaveInFormat = GlobalAlloc(GMEM_MOVEABLE, sizeof(WAVEFORMATEX));
  if(!hWaveInFormat)
  {
    Application->MessageBox("Error allocating wave in format", "Error", MB_OK);
    Application->Terminate();
  }

  pWaveInFormat = (WAVEFORMATEX *)(GlobalLock(hWaveInFormat));
  if(!pWaveInFormat)
  {
    Application->MessageBox("Error locking wave in format", "Error", MB_OK);
    Application->Terminate();
  }

  for(i = 0; i < WAVEIN_MAXBUFS; i++)
  {
    hWaveInHeaders[i] =
      GlobalAlloc(GMEM_MOVEABLE | GMEM_SHARE | GMEM_ZEROINIT, sizeof(WAVEHDR));
    if(!hWaveInHeaders[i])
    {
      Application->MessageBox("Error allocating wave in headers", "Error", MB_OK);
      Application->Terminate();
    }

    pWaveInHeaders[i] = (WAVEHDR *)(GlobalLock(hWaveInHeaders[i]));
    if(!pWaveInHeaders[i])
    {
      Application->MessageBox("Error locking wave in headers", "Error", MB_OK);
      Application->Terminate();
    }
  }

  for(i = 0; i < WAVEIN_MAXBUFS; i++)
  {
    hWaveInBuffers[i] =
      GlobalAlloc(GMEM_MOVEABLE | GMEM_SHARE, WAVEIN_BUFSIZE * sizeof(SAMPLE_TYPE));
    if(!hWaveInBuffers[i])
    {
      Application->MessageBox("Error allocating wave in buffers", "Error", MB_OK);
      Application->Terminate();
    }

    pWaveInBuffers[i] = (SAMPLE_TYPE *)(GlobalLock(hWaveInBuffers[i]));
    if(!pWaveInBuffers[i])
    {
      Application->MessageBox("Error locking wave in buffers", "Error", MB_OK);
      Application->Terminate();
    }
  }
}
//---------------------------------------------------------------------------
void DestroyWaveMemory()
{
  int i;

  if(hWaveInBuffers)
  {
    for(i = 0; i < WAVEIN_MAXBUFS; i++)
    {
      if(hWaveInBuffers[i])
      {
        GlobalUnlock(hWaveInBuffers[i]);
        GlobalFree(hWaveInBuffers[i]);
      }
      hWaveInBuffers[i] = NULL;
    }
    delete hWaveInBuffers;
  }
  hWaveInBuffers = NULL;
  if(pWaveInBuffers)
    delete pWaveInBuffers;
  pWaveInBuffers = NULL;

  if(hWaveInHeaders)
  {
    for(i = 0; i < WAVEIN_MAXBUFS; i++)
    {
      if(hWaveInHeaders[i])
      {
        GlobalUnlock(hWaveInHeaders[i]);
        GlobalFree(hWaveInHeaders[i]);
      }
      hWaveInHeaders[i] = NULL;
    }
    delete hWaveInHeaders;
  }
  hWaveInHeaders = NULL;
  if(pWaveInHeaders)
    delete pWaveInHeaders;
  pWaveInHeaders = NULL;

  if(hWaveInFormat)
  {
    GlobalUnlock(hWaveInFormat);
    GlobalFree(hWaveInFormat);
  }
  hWaveInFormat = NULL;  
}
//---------------------------------------------------------------------------
void QueryWave()
{
  int i;
  NumWaveInDevs = waveInGetNumDevs();
  for(i = 0; i < NumWaveInDevs; i++)
  {
    waveInGetDevCaps(i, &WaveInCapArray[i], sizeof(WAVEINCAPS));
    WaveInDevNames->Add(WaveInCapArray[i].szPname);
    //WaveInDevNames.Add(StrPas(WaveInCapArray[i].szPname));
  }
}
//---------------------------------------------------------------------------
void SetWaveInFormat()
{
  int i;
  pWaveInFormat->wFormatTag = WAVE_FORMAT_PCM;
  pWaveInFormat->nChannels = 1;
  pWaveInFormat->nSamplesPerSec = WAVEIN_SAMPRATE;
  if(sizeof(SAMPLE_TYPE) == 2)
    pWaveInFormat->wBitsPerSample = 16;
  else
    pWaveInFormat->wBitsPerSample = 8;
  //pWaveInFormat->nBlockAlign = 1;
  pWaveInFormat->nBlockAlign = pWaveInFormat->nChannels * pWaveInFormat->wBitsPerSample / 8;
  //pWaveInFormat->nAvgBytesPerSec = WAVEIN_SAMPRATE;
  pWaveInFormat->nAvgBytesPerSec = pWaveInFormat->nSamplesPerSec * pWaveInFormat->nBlockAlign;
  pWaveInFormat->cbSize = 0;

  for(i = 0; i < WAVEIN_MAXBUFS; i++)
  {
    pWaveInHeaders[i]->lpData = (char*)pWaveInBuffers[i];
    pWaveInHeaders[i]->dwBufferLength = WAVEIN_BUFSIZE * sizeof(SAMPLE_TYPE);
    pWaveInHeaders[i]->dwBytesRecorded = 0;
    pWaveInHeaders[i]->dwUser = i;
    pWaveInHeaders[i]->dwFlags = 0;
    pWaveInHeaders[i]->dwLoops = 0;
  }
}
//---------------------------------------------------------------------------
void PrepareWaveInHeaders()
{
  int i;
  for(i = 0; i < WAVEIN_MAXBUFS; i++)
    waveInPrepareHeader(CurrWaveInHandle, pWaveInHeaders[i], sizeof(WAVEHDR));
}
//---------------------------------------------------------------------------
void UnPrepareWaveInHeaders()
{
  int i;
  for(i = 0; i < WAVEIN_MAXBUFS; i++)
  waveInUnprepareHeader(CurrWaveInHandle, pWaveInHeaders[i], sizeof(WAVEHDR));
}
//---------------------------------------------------------------------------
void InitWaveBufs()
{
  int i, j;

  for(j = 0; j < WAVEIN_MAXBUFS; j++)
  {
    for(i = 0; i < WAVEIN_BUFSIZE; i++)
      pWaveInBuffers[j][i] = SAMPLE_TYPE(0);
  }
}
//---------------------------------------------------------------------------
void OpenWaveIn(HWND WindowHandle, UINT ID)
{
  MMRESULT r;
  r = waveInOpen(&CurrWaveInHandle,
                 ID, pWaveInFormat,
                 DWORD( /*WindowHandle*/ WaveInFunc /*WaveOutProcAddress*/),
                 0, /*CALLBACK_WINDOW*/ CALLBACK_FUNCTION );
  if(r != MMSYSERR_NOERROR)
    Application->MessageBox("Error opening Wave in", "Error", MB_OK);

  /*if waveOutSetPitch(CurrWaveOutHandle, $00020000) =
      MMSYSERR_NOTSUPPORTED then
    Application.MessageBox('Pitch change not supported', 'Error', MB_OK);
  if waveOutSetPlaybackRate(CurrWaveOutHandle, $00030000) =
      MMSYSERR_NOTSUPPORTED then
    Application.MessageBox('PBRate change not supported', 'Error', MB_OK);
   */
  PrepareWaveInHeaders();
}
//---------------------------------------------------------------------------
void CloseWaveIn()
{
  if(Recording == 1)
    StopWaveInRec();
  //waveInReset(CurrWaveInHandle);
  UnPrepareWaveInHeaders;
  //waveInClose(CurrWaveInHandle);
}
//---------------------------------------------------------------------------
void StartWaveInRec()
{
  int i;
  if(Recording == 0)
  {
    StopWaveInFlag = 1;
    WaveInBufsMissed = 0;  //Set these to 0 to be safe
    WaveInBufPend = 0;
    WaveInBufIndex = 0;
    ClearRepeatTables();
    for(i = 0; i < WAVEIN_MAXBUFS; i++)
    {
      if(waveInAddBuffer(CurrWaveInHandle,
          pWaveInHeaders[i], sizeof(WAVEHDR)) != MMSYSERR_NOERROR)
      {
        Application->MessageBox("Error adding Wave in buffers", "Error", MB_OK);
        Application->Terminate();
      }
    }
    //WaveInPrivMsgRdy = True;
    Recording = 1;
    StopWaveInFlag = 0;
    if(waveInStart(CurrWaveInHandle) != MMSYSERR_NOERROR)
      Application->MessageBox("Error starting Wave in", "Error", MB_OK);
  }
}
//---------------------------------------------------------------------------
void StopWaveInRec()
{
  StopWaveInFlag = 1;
  waveInReset(CurrWaveInHandle);
  ClearRepeatTables();
  //WaveInPrivMsgRdy = False;
  WaveInBufsMissed = 0;  //Set these to 0 to be safe
  WaveInBufPend = 0;
  WaveInBufIndex = 0;
  Recording = 0;
}
//---------------------------------------------------------------------------


/*
  //2 down
  if(CurrPlayPrgs[CurrPlayChn] > MIN_PLAY_PRG)
  {
    CurrPlayPrgs[CurrPlayChn]--;
    MaskEdit2->Text = IntToStr(CurrPlayPrgs[CurrPlayChn] + 1);
  }

  // 2 up
  if(CurrPlayPrgs[CurrPlayChn] < MAX_PLAY_PRG)
  {
    CurrPlayPrgs[CurrPlayChn]++;
    MaskEdit2->Text = IntToStr(CurrPlayPrgs[CurrPlayChn] + 1);
  }

  // 1 down
  if(CurrPlayOct > MIN_PLAY_OCT)
  {
    CurrPlayOct--;
    MaskEdit1->Text = IntToStr(CurrPlayOct + 1);
  }

  // 1 up
  if(CurrPlayOct < MAX_PLAY_OCT)
  {
    CurrPlayOct++;
    MaskEdit1->Text = IntToStr(CurrPlayOct + 1);
  }

  // 3 down
  if(CurrPlayChn > MIN_PLAY_CHN)
  {
    CurrPlayChn--;
    MaskEdit3->Text = IntToStr(CurrPlayChn + 1);
  }

  // 3 up
  if(CurrPlayChn < MAX_PLAY_CHN)
  {
    CurrPlayChn++;
    MaskEdit3->Text = IntToStr(CurrPlayChn + 1);
  }
*/

/*
void __fastcall TMainForm::CSpinButton2DownClick(TObject *Sender)
{
  //2 down
  if(CurrPlayPrgs[CurrPlayChn] > MIN_PLAY_PRG)
  {
    CurrPlayPrgs[CurrPlayChn]--;
    MaskEdit2->Text = IntToStr(CurrPlayPrgs[CurrPlayChn] + 1);
  }
}
//---------------------------------------------------------------------------
void __fastcall TMainForm::CSpinButton2UpClick(TObject *Sender)
{
  // 2 up
  if(CurrPlayPrgs[CurrPlayChn] < MAX_PLAY_PRG)
  {
    CurrPlayPrgs[CurrPlayChn]++;
    MaskEdit2->Text = IntToStr(CurrPlayPrgs[CurrPlayChn] + 1);
  }
}
//---------------------------------------------------------------------------
void __fastcall TMainForm::CSpinButton1DownClick(TObject *Sender)
{
  // 1 down
  if(CurrPlayOct > MIN_PLAY_OCT)
  {
    CurrPlayOct--;
    MaskEdit1->Text = IntToStr(CurrPlayOct + 1);
  }
}
//---------------------------------------------------------------------------
void __fastcall TMainForm::CSpinButton1UpClick(TObject *Sender)
{
  // 1 up
  if(CurrPlayOct < MAX_PLAY_OCT)
  {
    CurrPlayOct++;
    MaskEdit1->Text = IntToStr(CurrPlayOct + 1);
  }
}
//---------------------------------------------------------------------------
void __fastcall TMainForm::CSpinButton3DownClick(TObject *Sender)
{
  // 3 down
  if(CurrPlayChn > MIN_PLAY_CHN)
  {
    CurrPlayChn--;
    MaskEdit3->Text = IntToStr(CurrPlayChn + 1);
  }
}
//---------------------------------------------------------------------------
void __fastcall TMainForm::CSpinButton3UpClick(TObject *Sender)
{
  // 3 up
  if(CurrPlayChn < MAX_PLAY_CHN)
  {
    CurrPlayChn++;
    MaskEdit3->Text = IntToStr(CurrPlayChn + 1);
  }
}     */
//---------------------------------------------------------------------------

void __fastcall TMainForm::SenseLevelZKnobChange(TObject *Sender)
{
  DSPControl.SenseLevel = int(RollEdit1->Max - RollEdit1->Value) & 0xFF;
  DSPControl.ScaledSenseLevel =
    (float)(DSPControl.SenseLevel) / (float)(MAG_OUT_SCALE);
}
//---------------------------------------------------------------------------
void __fastcall TMainForm::TuneOffsetZKnobChange(TObject *Sender)
{
  DSPControl.TuneOffset = (int)RollEdit3->Value;
  // Call this to re-init the table:
  InitNoteToSpecMap();
}
//---------------------------------------------------------------------------
void __fastcall TMainForm::MidiProgramRollEditChange(TObject *Sender)
{
//  if(CurrPlayPrgs[CurrPlayChn] > MIN_PLAY_PRG)
//  {
//    CurrPlayPrgs[CurrPlayChn]--;
//    MaskEdit2->Text = IntToStr(CurrPlayPrgs[CurrPlayChn] + 1);
//  }
  if(((int)MidiProgramRollEdit->Value >= (int)MIN_PLAY_PRG) || ((int)MidiProgramRollEdit->Value <= (int)MAX_PLAY_PRG))  {
    CurrPlayPrgs[CurrPlayChn] = (int)MidiProgramRollEdit->Value;
    midiOutShortMsg(CurrMidiOutHandle,
      0x000000C0 | (CurrPlayPrgs[CurrPlayChn] << 8) | (CurrPlayChn & 0x0F));
  }
}
//---------------------------------------------------------------------------
void __fastcall TMainForm::MidiOctaveRollEditChange(TObject *Sender)
{
  if(((int)MidiOctaveRollEdit->Value >= (int)MIN_PLAY_OCT) || ((int)MidiOctaveRollEdit->Value <= (int)MAX_PLAY_OCT))
    CurrPlayOct = (int)MidiOctaveRollEdit->Value;
}
//---------------------------------------------------------------------------
void __fastcall TMainForm::MidiChannelRollEditChange(TObject *Sender)
{
  int Val = (int)MidiChannelRollEdit->Value -1;
  if((Val >= MIN_PLAY_CHN) || (Val <= MAX_PLAY_CHN))
    CurrPlayChn = Val;
}
//---------------------------------------------------------------------------
void __fastcall TMainForm::NoteOffsetRollEditChange(TObject *Sender)
{
  DSPControl.NoteOffset = (int)NoteOffsetRollEdit->Value & 0xFF;
}
//---------------------------------------------------------------------------
void __fastcall TMainForm::RadioGroup1Click(TObject *Sender)
{
  if(RadioGroup1->ItemIndex == 1)
    GraphUsePeakElseAvg = false;
  else
    GraphUsePeakElseAvg = true;
}
//---------------------------------------------------------------------------

void __fastcall TMainForm::SenseLevelExpZKnobChange(TObject *Sender)
{
  SenseLevelExp = (double)RollEdit2->Value / RollEdit2->Max;
}
//---------------------------------------------------------------------------

void __fastcall TMainForm::GraphYScaleZKnobChange(TObject *Sender)
{
  GraphYScale = (tFloat)GraphYScaleZKnob->Value / 10.0f;
}
//---------------------------------------------------------------------------
void __fastcall TMainForm::UseGraphCheckBoxClick(TObject *Sender)
{
  if(UseGraphCheckBox->Checked)
    UseGraph = true;
  else
    UseGraph = false;
}
//---------------------------------------------------------------------------

void __fastcall TMainForm::FreqGainZKnobChange(TObject *Sender)
{
  DSPControl.GainPerFreq = (float)FreqGainZKnob->Value / 10000.0f;
}
//---------------------------------------------------------------------------

void __fastcall TMainForm::Knob3Change(TObject *Sender)
{
  // July 4, 2002 Test: To make it think it's sampling at a lower frequency than it is.
  WAVEIN_SAMPRATE_FACTOR = (double)RollEdit7->Value / RollEdit7->Max;
  InitNoteToSpecMap();
}
//---------------------------------------------------------------------------

void __fastcall TMainForm::Knob4Change(TObject *Sender)
{
  // July 4, 2002 Test: To make it think it's sampling at a lower frequency than it is.
  WAVEIN_SAMPRATE_OFFSET =  WAVEIN_SAMPRATE * (double)RollEdit8->Value / (double)RollEdit8->Max;
  InitNoteToSpecMap();
}
//---------------------------------------------------------------------------


