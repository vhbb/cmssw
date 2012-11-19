#ifndef JETSUBSTRUCTURETOOLS_H
#define JETSUBSTRUCTURETOOLS_H


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h" 

#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbEventAuxInfo.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbEvent.h"

#include <fastjet/ClusterSequence.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/GhostedAreaSpec.hh"

#include "VHbbAnalysis/VHbbDataFormats/interface/Nsubjettiness.h"

#include <iostream>

class JetSubstructureTools {
    
    // member functions
 public:
    
    // constructor
    JetSubstructureTools( std::vector<float> c_px, std::vector<float> c_py, std::vector<float> c_pz, std::vector<float> c_e, std::vector<float> c_pdgId );    
    ~JetSubstructureTools(){}
    
    fastjet::PseudoJet getPrunedJet(){ return prunedJet_; }
    fastjet::PseudoJet getTrimmedJet(){ return trimmedJet_; }
    fastjet::PseudoJet getFilteredJet(){ return filteredJet_; }    
    float getPrunedJetArea(){ return prunedJetArea_; }
    float getTrimmedJetArea(){ return trimmedJetArea_; }
    float getFilteredJetArea(){ return filteredJetArea_; }
    
    float getTau1(){ return tau1_; }
    float getTau2(){ return tau2_; }
    float getTau3(){ return tau3_; }
    float getTau4(){ return tau4_; }    
    
    fastjet::PseudoJet getPrunedSubJet(int i){ 
        if (i == 0) return prunedSubJet1_;
        else if (i == 1) return prunedSubJet2_;
        else throw cms::Exception("JetSubstructureTools") << "Too many subjets..." << std::endl;
    }

    // data members
 public:
    
    std::vector<fastjet::PseudoJet> FJconstituents_;
    std::vector<float> c_pdgIds_; 
    
    fastjet::PseudoJet prunedJet_;
    fastjet::PseudoJet trimmedJet_;
    fastjet::PseudoJet filteredJet_;
    
    float prunedJetMass_;
    float trimmedJetMass_;
    float filteredJetMass_;
    float ungroomedJetMass_;
    
    float prunedJetArea_;
    float trimmedJetArea_;
    float filteredJetArea_;
    
    float tau1_;
    float tau2_;
    float tau3_;
    float tau4_;
    
    float rcore01_;
    float rcore02_;
    float rcore03_;
    float rcore04_;
    float rcore05_;
    float rcore06_;
    float rcore07_;
    float rcore08_;
    float rcore09_;
    float rcore10_;
    float rcore11_;
    float rcore12_;

    fastjet::PseudoJet prunedSubJet1_;
    fastjet::PseudoJet prunedSubJet2_;    
    
    float QJetVolatility_;
};

    // -------------------------------------------
    // -------------------------------------------
    // -------------------------------------------
    // -------------------------------------------

    // constructor
JetSubstructureTools::JetSubstructureTools( std::vector<float> c_px, std::vector<float> c_py, std::vector<float> c_pz, std::vector<float> c_e, std::vector<float> c_pdgId )   
{
    
    std::cout << "hey y'all!" << std::endl;
    
        // check that they are all the same size
    if ((c_px.size() == c_py.size())&&(c_py.size() == c_pz.size())&&(c_pz.size() == c_e.size())&&(c_e.size() == c_pdgId.size())){
        for (unsigned int i = 0; i < c_px.size(); i++){
            
            FJconstituents_.push_back( fastjet::PseudoJet( c_px[i], c_py[i], c_pz[i], c_e[i] ) );
            c_pdgIds_.push_back( c_pdgId[i] );
            
        }
    }
    else throw cms::Exception("JetSubstructureTools") << "Constituent size mismatch..." << std::endl;
    
        // -------------------------------------------
        // recluster on the fly....
    double mJetRadius = 1.2;
    fastjet::JetDefinition jetDef(fastjet::cambridge_algorithm, mJetRadius);
    
    int activeAreaRepeats = 1;
    double ghostArea = 0.01;
    double ghostEtaMax = 5.0;

    fastjet::GhostedAreaSpec fjActiveArea(ghostEtaMax,activeAreaRepeats,ghostArea);
    fjActiveArea.set_fj2_placement(true);
    fastjet::AreaDefinition fjAreaDefinition( fastjet::active_area_explicit_ghosts, fjActiveArea );

    fastjet::ClusterSequenceArea thisClustering(FJconstituents_, jetDef, fjAreaDefinition);
    std::vector<fastjet::PseudoJet> out_jets = sorted_by_pt(thisClustering.inclusive_jets(50.0));

    fastjet::ClusterSequence thisClustering_basic(FJconstituents_, jetDef);
    std::vector<fastjet::PseudoJet> out_jets_basic = sorted_by_pt(thisClustering_basic.inclusive_jets(50.0));

        // -------------------------------------------
        // define groomers
    fastjet::Filter trimmer( fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.2), fastjet::SelectorPtFractionMin(0.03)));
    fastjet::Filter filter( fastjet::Filter(fastjet::JetDefinition(fastjet::cambridge_algorithm, 0.3), fastjet::SelectorNHardest(3)));
    fastjet::Pruner pruner(fastjet::cambridge_algorithm, 0.1, 0.5);
    
    std::vector<fastjet::Transformer const *> transformers;
    transformers.push_back(&trimmer);
    transformers.push_back(&filter);
    transformers.push_back(&pruner);

        // define n-subjettiness
    float mNsubjettinessKappa = 1.;
    NsubParameters paraNsub = NsubParameters(mNsubjettinessKappa, mJetRadius);   
    Nsubjettiness routine(nsub_kt_axes, paraNsub);
    
    
    ungroomedJetMass_ = out_jets.at(0).m();
        // -------------------------------------------    
        // compute pruning, trimming, filtering  -------------    
    int transctr = 0;
    for ( std::vector<fastjet::Transformer const *>::const_iterator 
         itransf = transformers.begin(), itransfEnd = transformers.end(); 
         itransf != itransfEnd; ++itransf ) {  
        
        fastjet::PseudoJet transformedJet = out_jets.at(0);
        transformedJet = (**itransf)(transformedJet);
        
        if (transctr == 0){ // trimmed
            trimmedJet_ = transformedJet;
            trimmedJetMass_ = transformedJet.m();            
            trimmedJetArea_ = transformedJet.area();
        }
        else if (transctr == 1){ // filtered
            filteredJet_ = transformedJet;
            filteredJetMass_ = transformedJet.m();            
            filteredJetArea_ = transformedJet.area();
        }
        else if (transctr == 2){ // pruned
            prunedJet_ = transformedJet;
            prunedJetMass_ = transformedJet.m();            
            prunedJetArea_ = transformedJet.area();
            
                //decompose into requested number of subjets:
            if (transformedJet.constituents().size() > 1){
                
                int nsubjetstokeep = 2;
                std::vector<fastjet::PseudoJet> prunedSubjets = transformedJet.associated_cluster_sequence()->exclusive_subjets(transformedJet,nsubjetstokeep);    
                
                prunedSubJet1_ = prunedSubjets.at(0);
                prunedSubJet2_ = prunedSubjets.at(1);
                
            }
        }
        else{ std::cout << "error in number of transformers" << std::endl;}                    
        transctr++;
    }        
    
        // -------------------------------------------    
        // compute n-subjettiness  -------------
    tau1_ = routine.getTau(1, out_jets.at(0).constituents()); 
    tau2_ = routine.getTau(2, out_jets.at(0).constituents());
    tau3_ = routine.getTau(3, out_jets.at(0).constituents());
    tau4_ = routine.getTau(4, out_jets.at(0).constituents());
    
}
    // -------------------------------------------
    // -------------------------------------------
    // -------------------------------------------
    // -------------------------------------------


#endif






