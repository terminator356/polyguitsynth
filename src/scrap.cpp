Triggered = false;
      TotalTriggers = 0;
      DSPData.AvgPer = 0;
      for(i = 0; i < WAVEIN_BUFSIZE; i++)
      {
        k = (int)(pWaveInBuffers[((LPWAVEHDR)dwParam1)->dwUser][i]) ;
        if(DSPData.InSquared == 0)
        {
          if(k >= DSPControl.SenseLevel)
          {
            DSPData.SmpToPrevEdge =
              DSPData.BufsToPrevEdge * WAVEIN_BUFSIZE +
                (i - DSPData.PrevEdgeOffset);
            DSPData.AvgPer = DSPData.AvgPer + DSPData.SmpToPrevEdge;
            TotalTriggers++;

            Triggered = true;
            DSPData.InSquared = 1;
            DSPData.PrevEdgeOffset = i;
            DSPData.BufsToPrevEdge = 0;
          }
        }
        else
        if(k < DSPControl.NegSenseLevel)
        {
          DSPData.InSquared = 0;
        }
      }
      if(Triggered)
      {
        DSPData.AvgPer = DSPData.AvgPer / TotalTriggers; //div <<<<<<<<<<<<<<<<<<<<
        DSPData.NoteFreq =
          DSPControl.NoteOffset - ((int)((float)DSPData.AvgPer / 12.0) & 0x7F);
        midiOutShortMsg(CurrMidiOutHandle,
          0x00400080 | (DSPData.PrevNoteOn << 8) | (CurrPlayChn & 0x0F));
        midiOutShortMsg(CurrMidiOutHandle,
          0x00400090 | (DSPData.NoteFreq << 8) | (CurrPlayChn & 0x0F));
        DSPData.PrevNoteOn = DSPData.NoteFreq;
      }
      DSPData.BufsToPrevEdge++;

