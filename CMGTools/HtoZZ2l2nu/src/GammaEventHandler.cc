#include "CMGTools/HtoZZ2l2nu/interface/GammaEventHandler.h"

using namespace std;

//
GammaEventHandler::GammaEventHandler(const edm::ParameterSet &runProcess)
  : isGoodEvent_(false)
{
  //trigger thresholds to consider
  bool isMC = runProcess.getParameter<bool>("isMC");
  gammaCats_ = runProcess.getParameter<std::vector<int> >("gammaCategories"); 
  if(!isMC) gammaTriggerRenWeights_ = runProcess.getParameter<std::vector<double> >("gammaTriggerRenWeights");
  else      gammaTriggerRenWeights_.resize(gammaCats_.size(),1.0);
 
  //open file and retrieve weights + mass shapes
  std::vector<string> weightsFiles =  runProcess.getParameter<std::vector<std::string> >("weightsFile");  
  for(std::vector<string>::iterator it = weightsFiles.begin();  it != weightsFiles.end(); it++)
    {
      TString gammaPtWeightsFile(*it); 
      gSystem->ExpandPathName(gammaPtWeightsFile);
      fwgt_=TFile::Open(gammaPtWeightsFile);
      if(fwgt_==0 || fwgt_->IsZombie()) return;
      cout << gammaPtWeightsFile << endl;
      TString cats[]   =  {"eq0jets","eq1jets","eq2jets","geq3jets","vbf"};
      TString dilCats[] = {"ee","mumu","ll"};
      TString wgtType( isMC ? "mcwgts" : "datawgts");
      TString wgtName(gammaPtWeightsFile.Contains("nvtx")? "nvtx":"qt");
      for(size_t ic=0; ic<sizeof(cats)/sizeof(TString); ic++)
	{
	  for(size_t id=0; id<sizeof(dilCats)/sizeof(TString); id++)
	    {
	      //event weights
	      TString key = dilCats[id] + "_" + wgtName + "_" + cats[ic] + "_" + wgtType;
	      TH1 *h = (TH1 *) fwgt_->Get(key);
	      if(h!=0)
		{
		  key = dilCats[id] + "_"+wgtName+"_" + cats[ic];
		  TGraph *iwgtgr=new TGraph(h);
		  iwgtgr->SetName(key);
		  wgtsH_[key] = iwgtgr;
		}
	
	      //mass shape 
	      key = dilCats[id];
	      if(zmassH_.find(key)!=zmassH_.end()) continue;
	      key += (isMC ? "_mczmass" : "_datazmass");
	      h = (TH1 *) fwgt_->Get(key);
	      if(h!=0)
		{
		  key = dilCats[id];
		  zmassH_[key]= h;
		  zmassH_[key]->SetDirectory(0);
		}
	    }
	}

      fwgt_->Close();
    }
  
  if(wgtsH_.size()) std::cout << "[GammaEventHandler] gamma spectrum will be reweighted using distributions found in "  << weightsFiles.size() << " files" << std::endl;
}


//
bool GammaEventHandler::isGood(PhysicsEvent_t &phys)
{
  //reset
  isGoodEvent_=false;
  triggerThr_=0;
  photonCategory_="";
  massiveGamma_.clear();
  evWeights_.clear();
  triggerPrescaleWeight_=1;

  //check if it is a gamma event
  if( phys.cat<22) return isGoodEvent_;
  triggerThr_ = (phys.cat-22)/1000;
  if(triggerThr_==0) return isGoodEvent_;

  //check which category this event belongs to (use the trigger)
  int eventTriggerCat(-1);
  for(size_t icat=0; icat<gammaCats_.size()-1; icat++)
    {
      if(triggerThr_<gammaCats_[icat])    break;
      eventTriggerCat=icat;   
    }
  if(eventTriggerCat<0) return isGoodEvent_;
  triggerPrescaleWeight_ = gammaTriggerRenWeights_[eventTriggerCat];

  photonCategory_="photon";  photonCategory_ += triggerThr_; 
  
  //all done here
  isGoodEvent_=true;
  return isGoodEvent_;
}

//
std::map<TString,float> GammaEventHandler::getWeights(PhysicsEvent_t &phys, TString evCategoryLabel)
{
  //loop over categories
  LorentzVector gamma=phys.gammas[0];
  TString dilCats[]={"ee","mumu","ll"};
  for(size_t id=0; id<sizeof(dilCats)/sizeof(TString); id++)
    {
      //generate a mass
      TString key = dilCats[id];
      float mass(0);
      if(zmassH_.find(key)!=zmassH_.end())
	{
	  if(zmassH_[key]->Integral())
	    while(fabs(mass-91)>15) 
	      mass = zmassH_[key]->GetRandom();
	}
      massiveGamma_[dilCats[id]]=LorentzVector(gamma.px(),gamma.py(),gamma.pz(),sqrt(pow(mass,2)+pow(gamma.energy(),2)));
      
      //get event weight 
      float weight(1.0);
      if(wgtsH_.size())
	{
	  key = dilCats[id] + "_qt_" + evCategoryLabel;
	  if(wgtsH_.find(key) != wgtsH_.end()) 
	    {
	      TGraph *iwgtgr = wgtsH_[key];
	      weight *= iwgtgr->Eval(gamma.pt());
	    }

	  key = dilCats[id] + "_nvtx_" + evCategoryLabel;
	  if(wgtsH_.find(key) != wgtsH_.end()) 
	    {
	      TGraph *iwgtgr = wgtsH_[key];
	      weight *= iwgtgr->Eval(phys.nvtx);
	    }

	}
      evWeights_[dilCats[id]]=weight;
    }
  
  return evWeights_;
}

//
GammaEventHandler::~GammaEventHandler()
{
}
