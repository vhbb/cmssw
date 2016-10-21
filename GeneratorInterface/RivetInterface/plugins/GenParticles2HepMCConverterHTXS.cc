#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

#include <iostream>
#include <map>

using namespace std;

class GenParticles2HepMCConverterHTXS : public edm::stream::EDProducer<>
{
public:
  explicit GenParticles2HepMCConverterHTXS(const edm::ParameterSet& pset);
  ~GenParticles2HepMCConverterHTXS() {};

  void produce(edm::Event& event, const edm::EventSetup& eventSetup) override;
  bool checkForAncestor(const reco::Candidate * particle, int pdgId);

private:
  edm::EDGetTokenT<reco::CandidateView> prunedGenParticlesToken_;
  edm::EDGetTokenT<reco::CandidateView> packedGenParticlesToken_;
  edm::EDGetTokenT<GenEventInfoProduct> genEventInfoToken_;
  edm::ESHandle<ParticleDataTable> pTable_;

private:
  inline HepMC::FourVector FourVector(const reco::Candidate::Point& point)
  {
    return HepMC::FourVector(10*point.x(), 10*point.y(), 10*point.z(), 0);
  };

  inline HepMC::FourVector FourVector(const reco::Candidate::LorentzVector& lvec)
  {
    // Avoid negative mass, set minimum m^2 = 0
    return HepMC::FourVector(lvec.px(), lvec.py(), lvec.pz(), std::hypot(lvec.P(), std::max(0., lvec.mass())));
  };


};

GenParticles2HepMCConverterHTXS::GenParticles2HepMCConverterHTXS(const edm::ParameterSet& pset)
{
  prunedGenParticlesToken_ = consumes<reco::CandidateView>(pset.getParameter<edm::InputTag>("prunedGenParticles")); // contain hard scattering info
  packedGenParticlesToken_ = consumes<reco::CandidateView>(pset.getParameter<edm::InputTag>("packedGenParticles")); // contain status 1 particle info
  genEventInfoToken_ = consumes<GenEventInfoProduct>(pset.getParameter<edm::InputTag>("genEventInfo"));

  produces<edm::HepMCProduct>("unsmeared");
}

bool GenParticles2HepMCConverterHTXS::checkForAncestor(const reco::Candidate * particle, int pdgId)
{

    if ((int)particle->pdgId()==pdgId) return true;

    for(size_t i=0;i< particle->numberOfMothers();i++)
    {
        if (checkForAncestor(particle->mother(i),pdgId)) return true;
    }

    return false;

}

void GenParticles2HepMCConverterHTXS::produce(edm::Event& event, const edm::EventSetup& eventSetup)
{

  edm::Handle<reco::CandidateView> prunedGenParticlesHandle;
  event.getByToken(prunedGenParticlesToken_, prunedGenParticlesHandle);

  edm::Handle<reco::CandidateView> packedGenParticlesHandle;
  event.getByToken(packedGenParticlesToken_, packedGenParticlesHandle);

  edm::Handle<GenEventInfoProduct> genEventInfoHandle;
  event.getByToken(genEventInfoToken_, genEventInfoHandle);

  eventSetup.getData(pTable_);

  HepMC::GenEvent* hepmc_event = new HepMC::GenEvent();
  hepmc_event->set_event_number(event.id().event());
  hepmc_event->set_signal_process_id(genEventInfoHandle->signalProcessID());
  hepmc_event->set_event_scale(genEventInfoHandle->qScale());
  hepmc_event->set_alphaQED(genEventInfoHandle->alphaQED());
  hepmc_event->set_alphaQCD(genEventInfoHandle->alphaQCD());

  hepmc_event->weights() = genEventInfoHandle->weights();

  // Set PDF
  const gen::PdfInfo* pdf = genEventInfoHandle->pdf();
  const int pdf_id1 = pdf->id.first, pdf_id2 = pdf->id.second;
  const double pdf_x1 = pdf->x.first, pdf_x2 = pdf->x.second;
  const double pdf_scalePDF = pdf->scalePDF;
  const double pdf_xPDF1 = pdf->xPDF.first, pdf_xPDF2 = pdf->xPDF.second;
  HepMC::PdfInfo hepmc_pdfInfo(pdf_id1, pdf_id2, pdf_x1, pdf_x2, pdf_scalePDF, pdf_xPDF1, pdf_xPDF2);
  hepmc_event->set_pdf_info(hepmc_pdfInfo);

  /////////////////////////////
  // SETUP hard scattering info from PRUNED GEN PARTICLES
  /////////////////////////////

  // Prepare list of HepMC::GenParticles from prunedGenParticles
  std::map<const reco::Candidate*, HepMC::GenParticle*> prunedGenCandToHepMCMap;
  std::vector<HepMC::GenParticle*> hepmc_prunedParticles;
  for ( unsigned int i=0, n=prunedGenParticlesHandle->size(); i<n; ++i )
  {
    const reco::Candidate* p_pruned = &prunedGenParticlesHandle->at(i);
    //cout << "pruned genparticle with pdgid " << p_pruned->pdgId() << " and status " << p_pruned->status() << " and pt,eta " << p_pruned->pt() << "," << p_pruned->eta() << " Ndaughters: " <<p_pruned->numberOfDaughters() << endl;
    HepMC::GenParticle* hepmc_prunedParticle = new HepMC::GenParticle(FourVector(p_pruned->p4()), p_pruned->pdgId(), p_pruned->status());
    hepmc_prunedParticle->suggest_barcode(i+1);

    // Assign particle's generated mass from the standard particle data table
    double particleMass;
    if ( pTable_->particle(p_pruned->pdgId()) ) particleMass = pTable_->particle(p_pruned->pdgId())->mass();
    else particleMass = p_pruned->mass();
    hepmc_prunedParticle->set_generated_mass(particleMass);

    hepmc_prunedParticles.push_back(hepmc_prunedParticle);
    prunedGenCandToHepMCMap[p_pruned] = hepmc_prunedParticle;
  }
  

  // Put incident beam particles : proton -> parton vertex
  const reco::Candidate* parton1_pruned = prunedGenParticlesHandle->at(0).daughter(0);
  const reco::Candidate* parton2_pruned = prunedGenParticlesHandle->at(1).daughter(0);
  HepMC::GenVertex* vertex1_pruned = new HepMC::GenVertex(FourVector(parton1_pruned->vertex()));
  HepMC::GenVertex* vertex2_pruned = new HepMC::GenVertex(FourVector(parton2_pruned->vertex()));
  hepmc_event->add_vertex(vertex1_pruned);
  hepmc_event->add_vertex(vertex2_pruned);
  vertex1_pruned->add_particle_in(hepmc_prunedParticles[0]);
  vertex2_pruned->add_particle_in(hepmc_prunedParticles[1]);
  hepmc_event->set_beam_particles(hepmc_prunedParticles[0], hepmc_prunedParticles[1]);

  // STORE THE FIRST (RANDOM) VERTEX IN THE hepmc_event
  HepMC::GenVertex* vertex3 = new HepMC::GenVertex(FourVector(parton2_pruned->vertex()));
  hepmc_event->add_vertex(vertex3);
  vertex3->add_particle_in(hepmc_prunedParticles[0]);
  vertex3->add_particle_in(hepmc_prunedParticles[1]);

  // Prepare vertex list
  typedef std::map<const reco::Candidate*, HepMC::GenVertex*> ParticleToVertexMap;
  ParticleToVertexMap particleToVertexMap;
  particleToVertexMap[parton1_pruned] = vertex1_pruned;
  particleToVertexMap[parton2_pruned] = vertex2_pruned;
  for ( unsigned int i=2, n=prunedGenParticlesHandle->size(); i<n; ++i )
  {
    const reco::Candidate* p = &prunedGenParticlesHandle->at(i);

    if (p->status()==1) continue;

    //cout<<"adding pruned id: "<<p->pdgId()<<" pt: "<<p->pt()<<" status: "<<p->status()<<endl;

    // Connect mother-daughters for the other cases
    for ( unsigned int j=0, nMothers=p->numberOfMothers(); j<nMothers; ++j )
    {
      // Mother-daughter hierarchy defines vertex
      const reco::Candidate* elder = p->mother(j)->daughter(0);
      HepMC::GenVertex* vertex;
      if ( particleToVertexMap.find(elder) == particleToVertexMap.end() )
      {
        vertex = new HepMC::GenVertex(FourVector(elder->vertex()));
        hepmc_event->add_vertex(vertex);
        particleToVertexMap[elder] = vertex;
      }
      else
      {
        vertex = particleToVertexMap[elder];
      }

      // Vertex is found. Now connect each other
      const reco::Candidate* mother = p->mother(j);
      vertex->add_particle_in(prunedGenCandToHepMCMap[mother]);
      vertex->add_particle_out(hepmc_prunedParticles[i]);
    }
  }

  // Finalize HepMC event record
  // hepmc_event->set_signal_process_vertex(*(vertex1_pruned->vertices_begin()));

  //Loop over all vertices 
  int particle_pruned = 0;
  int vtx = 0;
  bool endloop = false;
  for (HepMC::GenEvent::vertex_iterator ver = hepmc_event->vertices_begin(); ver != hepmc_event->vertices_end(); ver++){    
    if(endloop || hepmc_event->signal_process_vertex()) break;
    vtx++;
    //cout << "vtx= " << vtx ;
    particle_pruned = 0;
    HepMC::GenVertex::particle_iterator par = (*ver)->particles_begin(HepMC::children);    
    for (; par != (*ver)->particles_end(HepMC::children); ++par){
      particle_pruned++;
      //cout << " particle_pruned= " << particle_pruned << " pdgid= " << (*par)->pdg_id() << " ; ";       
      if ((*par)->pdg_id() == 25){
        hepmc_event->set_signal_process_vertex(*ver);
        endloop=true;
      }   
    }
    //cout<<endl;
    //if(endloop) cout << " assigned as signal_process_vertex!" << endl;   
  }

  /////////////////////////////
  // SETUP status 1 particles info from PACKED GEN PARTICLES
  /////////////////////////////

  // HepMC::GenEvent* packed_hepmc_event = new HepMC::GenEvent();
  
  // Prepare list of HepMC::GenParticles from packedGenParticles
  std::map<const reco::Candidate*, HepMC::GenParticle*> packedGenCandToHepMCMap;
  std::vector<HepMC::GenParticle*> hepmc_packedParticles;
  for ( unsigned int i=0, n=packedGenParticlesHandle->size(); i<n; ++i )
  {
    const reco::Candidate* p_packed = &packedGenParticlesHandle->at(i);

    /*
    //check that the particle is not already added to the event  
    bool matched=false;
    for ( unsigned int j=0, nj=prunedGenParticlesHandle->size(); j<nj; ++j )
    {
        const reco::Candidate* p_pruned = &prunedGenParticlesHandle->at(j);
        if (p_pruned->status()!=1) continue;
        if (p_pruned->pdgId()!=p_packed->pdgId()) continue;
        //cout<<"inside packed loop, found pruned id: "<<p_pruned->pdgId()<<" pt: "<<p_pruned->pt()<<" eta: "<<p_pruned->eta()<<endl;
        if (abs(p_pruned->pt()-p_packed->pt())/p_packed->pt()<0.01 && abs(p_pruned->eta()-p_packed->eta())/p_packed->eta()<0.05) matched=true;
    }        
    if (matched) {
        //cout<<"skipping packed id: "<<p_packed->pdgId()<<" pt: "<<p_packed->pt()<< " eta: "<< p_packed->eta()<<endl;
        continue;
    }
    */

    //cout << "packed id: " << p_packed->pdgId() << " pt: "<< p_packed->pt()<< " eta: " << p_packed->eta()<<endl;
    HepMC::GenParticle* hepmc_packedParticle = new HepMC::GenParticle(FourVector(p_packed->p4()), p_packed->pdgId(), p_packed->status());
    hepmc_packedParticle->suggest_barcode(i+1);

    // Assign particle's generated mass from the standard particle data table
    double particleMass;
    if ( pTable_->particle(p_packed->pdgId()) ) particleMass = pTable_->particle(p_packed->pdgId())->mass();
    else particleMass = p_packed->mass();
    hepmc_packedParticle->set_generated_mass(particleMass);

    hepmc_packedParticles.push_back(hepmc_packedParticle);
    packedGenCandToHepMCMap[p_packed] = hepmc_packedParticle;

    const reco::Candidate* motherInPruned = p_packed->mother(0);
    if (checkForAncestor(motherInPruned,25)==true) 
    {
        //cout<<"from Higgs, skipping id:"<<p_packed->pdgId()<<" pt: "<<p_packed->pt()<<endl;
        continue;
    }
    // STORE THE PACKED HEPMC PARTICLES INTO THE FIRST PRUNED VERTEX STORED IN THE hepmc_event
    vertex3->add_particle_out(hepmc_packedParticle);

  }

  /////////////////////////////
  // STORE HEPMC RESULT
  /////////////////////////////

  std::auto_ptr<edm::HepMCProduct> hepmc_product(new edm::HepMCProduct());
  hepmc_product->addHepMCData(hepmc_event);
  event.put(hepmc_product, "unsmeared");

}

DEFINE_FWK_MODULE(GenParticles2HepMCConverterHTXS);
