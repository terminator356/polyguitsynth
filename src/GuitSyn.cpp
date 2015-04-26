//---------------------------------------------------------------------------
#include <vcl.h>
#pragma hdrstop
USERES("GuitSyn.res");
USEFORM("MainUnit.cpp", MainForm);
USEUNIT("Libdsp.cpp");
//---------------------------------------------------------------------------
WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int)
{
        try
        {
                 Application->Initialize();
                 Application->CreateForm(__classid(TMainForm), &MainForm);
                 Application->Run();
        }
        catch (Exception &exception)
        {
                 Application->ShowException(&exception);
        }
        return 0;
}
//---------------------------------------------------------------------------
