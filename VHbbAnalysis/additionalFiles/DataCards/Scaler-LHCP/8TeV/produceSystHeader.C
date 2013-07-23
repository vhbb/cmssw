void produceSystHeader(int Bin = 0)
{
TString FileName = "XXX";
TString WSName;

TFile *f = new TFile(FileName.Data(),"read");
if(FileName.Contains("WmnH")) WSName= "WmnHighPt_8TeV";
if(FileName.Contains("WenH")) WSName= "WenHighPt_8TeV";
if(FileName.Contains("WmnM")) WSName= "WmnMidPt_8TeV";
if(FileName.Contains("WenM")) WSName= "WenMidPt_8TeV";
if(FileName.Contains("WmnL")) WSName= "WmnLowPt_8TeV";
if(FileName.Contains("WenL")) WSName= "WenLowPt_8TeV";
if(FileName.Contains ("Zn") && (Bin == 0)) WSName = "ZnunuLowPt_8TeV";
if(FileName.Contains ("Zn") && (Bin == 1)) WSName = "ZnunuMedPt_8TeV";
if(FileName.Contains ("Zn") && (Bin == 2)) WSName = "ZnunuHighPt_8TeV";
if(FileName.Contains ("ZmmL")) WSName = "ZmmLowPt_8TeV";
if(FileName.Contains ("ZmmH")) WSName = "ZmmHighPt_8TeV";
if(FileName.Contains ("ZeeL")) WSName = "ZeeLowPt_8TeV";
if(FileName.Contains ("ZeeH")) WSName = "ZeeHighPt_8TeV";
if(FileName.Contains ("Wt")) WSName = "Wtn";




RooWorkspace *myWS = (RooWorkspace*) f->Get(WSName.Data());


myWS->Print();
}

