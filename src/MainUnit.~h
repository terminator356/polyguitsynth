//---------------------------------------------------------------------------
#ifndef MainUnitH
#define MainUnitH
//---------------------------------------------------------------------------
#include <vcl\Classes.hpp>
#include <vcl\Controls.hpp>
#include <vcl\StdCtrls.hpp>
#include <vcl\Forms.hpp>
#include "Grids.hpp"
//#include <sampreg.h>
//#include <vcl\cspin.h>
#include <vcl\Mask.hpp>

#include <mmsystem.h>
#include "dsp.h"
//#include "CSPIN.h"

// Removed Zeigler collection - sharware! Using my own now.
//#include "zPanel.hpp"
//#include "zPanel2.hpp"

#include "RollEdit.h"
#include <ExtCtrls.hpp>
#include "knob.h"

#define   MAX_WAVEIN_DEVS       20
//#define   R2_ORDER              9      // 2^9 = 512
//#define   WAVEIN_BUFSIZE        512
//#define   FREQUENCY_BINS       (WAVEIN_BUFSIZE >> 1)
//#define   WAVEIN_MAXBUFS        64
//#define   WAVEIN_SAMPRATE       8000
typedef   SHORT SAMPLE_TYPE;   // BYTE or WORD
#define   PRIV_MSG_NAME_IN      "FX_WIM_DATA"

#define   MAX_MIDIIN_DEVS   20
#define   MAX_MIDIOUT_DEVS  20
#define   MIN_PLAY_PRG      0
#define   MAX_PLAY_PRG      127
#define   MIN_PLAY_OCT      0
#define   MAX_PLAY_OCT      9
#define   MIN_PLAY_CHN      0
#define   MAX_PLAY_CHN      15
#define   MAX_POLY_NOTES    128  // Never set above 128
#define   MAX_NOTES         128  // for guitar - 4 octaves, but use 128 anyways

//#define   DETECTOR_BASE_F   82.40688922822
// The lowest recognizable note given the most acceptable delay (.09s)
// Make it second lowest 'E' ...
//  Acceptable!!! - Can play any four note chord
//  from fret position 1 on up etc...
//#define   DETECTOR_BASE_F   (82.40688922822 * 2.0)
#define   DETECTOR_BASE_F   82.40688922822
//#define   DETECTOR_BASE_F_ADD
#define   MAG_OUT_SCALE     50

typedef struct
{
  int InSquared;
  /*int PrevVal;
  int PrevMax;*/
  int PrevEdgeOffset;
  int BufsToPrevEdge;
  int SmpToPrevEdge;
  int AvgPer;
  int NoteFreq;
  int PrevNoteOn;
}TDSPData;

typedef struct
{
  int         SenseLevel;
  tFloat      ScaledSenseLevel;
  int         NoteOffset;
  int         TuneOffset;
  float       GainPerFreq;
}TDSPControl;


void QueryMidi();
void OpenMidiOut(HWND WindowHandle, UINT ID);
void CloseMidiOut();
void InitKeyToNoteTable();
//void InitSpecToNoteMap();   Not used
void InitNoteToSpecMap();
void ClearRepeatTables();

void AllocateWaveMemory();
void DestroyWaveMemory();
void QueryWave();
void SetWaveInFormat();
void PrepareWaveInHeaders();
void UnPrepareWaveInHeaders();
void OpenWaveIn(HWND WindowHandle, UINT ID);
void CloseWaveIn();
void StartWaveInRec();
void StopWaveInRec();
void InitWaveBufs();
void CALLBACK WaveInFunc(HWAVEIN InHandle, UINT Msg,
    DWORD dwInstance, DWORD dwParam1, DWORD dwParam2);


/*
extern TDSPData DSPData;
extern TDSPControl DSPControl;

extern HMIDIIN CurrMidiInHandle;
extern HMIDIOUT CurrMidiOutHandle;
extern UINT NumMidiInDevs;
extern UINT NumMidiOutDevs;
extern MIDIINCAPS *MidiInCapArray;
extern MIDIOUTCAPS *MidiOutCapArray;
extern TStringList *MidiInDevNames;
extern TStringList *MidiOutDevNames;
extern BYTE *KeyToNoteTable;
extern BYTE *KeyRepeatTable;
extern int *CurrPlayPrgs;
extern int CurrPlayOct;
extern int CurrPlayChn;

extern HWAVEIN CurrWaveInHandle;
//WaveInProcAddress;
extern int NumWaveInDevs;
extern WAVEINCAPS *WaveInCapArray;
extern PWAVEFORMATEX pWaveInFormat;
extern HGLOBAL  hWaveInFormat;
extern PWAVEHDR *pWaveInHeaders;
extern HGLOBAL *hWaveInHeaders;
extern HGLOBAL *hWaveInBuffers;
extern char **pWaveInBuffers;
extern int WaveInBufIndex;
extern TStringList *WaveInDevNames;
extern int StopWaveInFlag;
extern int Recording;      */

/*const tInt         dspiN         =  512; // N Samples
tsCmplxRect        dspIn, dspOut, dspTwidTable;
ptFloat            pfdspOutMag;
ptInt              pidspT2SwapTable; */

//---------------------------------------------------------------------------
class TMainForm : public TForm
{
__published:	// IDE-managed Components
	TLabel *Label1;
	TLabel *Label2;
	TLabel *Label3;
	TLabel *Label4;
	TLabel *Label5;
	TLabel *Label6;
	TLabel *Label7;
	TLabel *Label8;
	TLabel *Label9;
	TLabel *Label10;
	TStringGrid *StringGrid1;
	TComboBox *ComboBox1;
	TComboBox *ComboBox2;
	TButton *NoteResetButton;
	TButton *Button2;
	TRadioButton *RadioButton1;
	TRadioButton *RadioButton2;
	TComboBox *WaveDevComboBox;
	TButton *StartWaveButton;
	TButton *StopWaveButton;
        TRollEdit *MidiProgramRollEdit;
        TRollEdit *MidiOctaveRollEdit;
        TRollEdit *MidiChannelRollEdit;
        TRollEdit *NoteOffsetRollEdit;
        TRadioGroup *RadioGroup1;
        TLabel *Label11;
        TLabel *Label12;
        TCheckBox *UseGraphCheckBox;
        TLabel *Label13;
        TKnob *SenseLevelZKnob;
        TKnob *SenseLevelExpZKnob;
        TKnob *TuneOffsetZKnob;
        TKnob *FreqGainZKnob;
        TKnob *GraphYScaleZKnob;
        TComboBox *ComboBox3;
        TLabel *Label15;
        TKnob *Knob1;
        TLabel *Label14;
        TKnob *Knob2;
        TLabel *Label16;
        TLabel *Label17;
        TLabel *Label18;
        TKnob *Knob4;
        TKnob *Knob3;
        TLabel *Label19;
        TLabel *Label20;
        TRollEdit *RollEdit1;
        TRollEdit *RollEdit2;
        TRollEdit *RollEdit3;
        TRollEdit *RollEdit4;
        TRollEdit *RollEdit5;
        TRollEdit *RollEdit6;
        TRollEdit *RollEdit7;
        TRollEdit *RollEdit8;
        //TzcSpectrum *zcSpectrum1;
        //TzKnob *SenseLevelZKnob;
        //TzKnob *TuneOffsetZKnob;
        //TzKnob *FreqGainZKnob;
        //TzKnob *GraphYScaleZKnob;
        //TzKnob *SenseLevelExpZKnob;
    //void __fastcall MyAppOnMessage(MSG &Msg, bool &Handled);
	void __fastcall NoteResetButtonClick(TObject *Sender);
	void __fastcall FormCreate(TObject *Sender);
	void __fastcall FormDestroy(TObject *Sender);
	void __fastcall FormKeyDown(TObject *Sender, WORD &Key, TShiftState Shift);
	void __fastcall FormKeyUp(TObject *Sender, WORD &Key, TShiftState Shift);
	void __fastcall ComboBox2Change(TObject *Sender);
	void __fastcall WaveDevComboBoxChange(TObject *Sender);
	void __fastcall StartWaveButtonClick(TObject *Sender);
	void __fastcall StopWaveButtonClick(TObject *Sender);
	void __fastcall RadioButton1Click(TObject *Sender);
	void __fastcall RadioButton2Click(TObject *Sender);
        void __fastcall MidiProgramRollEditChange(TObject *Sender);
        void __fastcall MidiOctaveRollEditChange(TObject *Sender);
        void __fastcall MidiChannelRollEditChange(TObject *Sender);
        void __fastcall NoteOffsetRollEditChange(TObject *Sender);
        void __fastcall RadioGroup1Click(TObject *Sender);
        void __fastcall UseGraphCheckBoxClick(TObject *Sender);
        void __fastcall SenseLevelZKnobChange(TObject *Sender);
        void __fastcall TuneOffsetZKnobChange(TObject *Sender);
        void __fastcall SenseLevelExpZKnobChange(TObject *Sender);
        void __fastcall GraphYScaleZKnobChange(TObject *Sender);
        void __fastcall FreqGainZKnobChange(TObject *Sender);
        void __fastcall Knob3Change(TObject *Sender);
        void __fastcall Knob4Change(TObject *Sender);
private:	// User declarations
public:		// User declarations
	__fastcall TMainForm(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern TMainForm *MainForm;
//---------------------------------------------------------------------------
#endif
