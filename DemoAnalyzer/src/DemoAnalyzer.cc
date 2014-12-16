// -*- C++ -*-
//
// Package:    DemoAnalyzer
// Class:      DemoAnalyzer
// 
/**\class DemoAnalyzer DemoAnalyzer.cc Demo/DemoAnalyzer/src/DemoAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Dominik Haitz, R. 8-22, 47243
//         Created:  Wed Nov  5 15:49:26 CET 2014
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TProfile.h"
#include "TNtuple.h"


#include <DataFormats/Math/interface/deltaR.h>
#include <DataFormats/JetReco/interface/PFJet.h>
#include <DataFormats/JetReco/interface/GenJet.h>
#include <DataFormats/JetReco/interface/GenJetCollection.h>

#include <DataFormats/VertexReco/interface/Vertex.h>
#include <DataFormats/VertexReco/interface/VertexFwd.h>
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include <DataFormats/JetReco/interface/Jet.h>

//
// class declaration
//

class DemoAnalyzer : public edm::EDAnalyzer {
   public:
      explicit DemoAnalyzer(const edm::ParameterSet&);
      ~DemoAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual int matchParton(reco::GenJet const&, edm::Handle<reco::GenParticleCollection> const&);
      virtual bool matchRecoGen(reco::PFJet const &, reco::GenJet const&);
      // ----------member data ---------------------------

  TProfile *hasmatch_genpt;
  TNtuple *ntuple;

  std::string jetCollection;
  std::string folder;
  bool m_leptonIso;


};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
DemoAnalyzer::DemoAnalyzer(const edm::ParameterSet& iConfig):
jetCollection(iConfig.getParameter<edm::InputTag>("jetCollection").label()),
folder(iConfig.getParameter<edm::InputTag>("folder").label()),
m_leptonIso(iConfig.getUntrackedParameter<bool>("leptonIso"))
{
    edm::Service<TFileService> fs;
    int min = 20;
    int max = 100;
    int nbins = 8;

    ntuple = fs->make<TNtuple>(folder.c_str(), folder.c_str(), 
    "genpt:recopt:flavour:rho:recoarea:nconst:spread:zmass:g:r:npv:iso:ptdist:maxdist");

    hasmatch_genpt = fs->make<TProfile>("hasmatch_genpt", "hasmatch_genpt", nbins, min, max);
}


DemoAnalyzer::~DemoAnalyzer()
{
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------
void
DemoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;


Handle<reco::VertexCollection> vertices;
iEvent.getByLabel("offlinePrimaryVertices", vertices);

Handle<double> rho;
iEvent.getByLabel("kt6PFJets", "rho", rho);

Handle<reco::PFJetCollection> recojets;
iEvent.getByLabel(jetCollection.c_str(), recojets);

Handle<reco::GenJetCollection> genjets;
iEvent.getByLabel("ak5GenJetsNoNu", genjets);


Handle<reco::GenParticleCollection> genparticles;
iEvent.getByLabel("genParticles", genparticles);



bool hasmatch;
int flavour;
float zmass;
int iso = 1;

if (genjets->size() < 1)
    return;

std::vector<reco::GenParticle> leptons;

//std::cout << "\n\n\n" << std::endl;/////////
for (unsigned int i=0; i < genparticles->size(); i++)
{
    //std::cout << genparticles->at(i).pdgId() << " ";/////////
    if (leptons.size() <2 &&  i>5 && std::abs(genparticles->at(i).pdgId()) <20 && std::abs(genparticles->at(i).pdgId())>10)
        leptons.push_back(genparticles->at(i));
}
//std::cout << "\n";////////

//for (unsigned int i=0; i < leptons.size(); i++){
//std::cout << leptons.at(i).pdgId()<< std::endl;////////
//}


for (unsigned int g = 0; g < genjets->size(); g++)
{

    flavour = DemoAnalyzer::matchParton(genjets->at(g), genparticles);


    hasmatch = false;
    for (unsigned int r = 0; r < recojets->size(); r++)
    {
        if (DemoAnalyzer::matchRecoGen(recojets->at(r), genjets->at(g)))
        {
            if (m_leptonIso && leptons.size() > 1)
            {
                if ((deltaR(genjets->at(g), leptons.at(0)) < 0.5) || (deltaR(genjets->at(g), leptons.at(1)) < 0.5))
                    iso=0;
            }

            //match
            hasmatch = true;

            //iterate over gen particles to get Z mass
            zmass = 0;
            for (unsigned int p = 0; p < genparticles->size(); p++)
            {
                if (genparticles->at(p).pdgId() == 23)
                    zmass = genparticles->at(p).mass();
            }

            ntuple->Fill(
                genjets->at(g).pt(),
                recojets->at(r).pt(),
                flavour,
                (*rho),
                recojets->at(r).jetArea(),
                recojets->at(r).nConstituents(),
                recojets->at(r).constituentEtaPhiSpread(),
                zmass,
                g,
                r,
                vertices->size(),
                iso,
                recojets->at(r).constituentPtDistribution(),
                recojets->at(r).maxDistance()
            );
        }
    }
    hasmatch_genpt->Fill(genjets->at(g).pt(), hasmatch);
}



#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
   
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


//this function matches a parton to a genjet
int DemoAnalyzer::matchParton(const reco::GenJet& genjet, const edm::Handle<reco::GenParticleCollection>& genparticles)
{
    std::vector<reco::GenParticle> matching_partons;
    //check for parton flavour   phys def
    for (unsigned int p = 0; p < genparticles->size(); p++)
    {
        if (
            (
                std::abs(genparticles->at(p).pdgId()) > 0
                && (genparticles->at(p).status() == 3)
                && (
                    (std::abs(genparticles->at(p).pdgId()) < 6)
                    || (std::abs(genparticles->at(p).pdgId()) == 21)
                )
            )
            && (deltaR(genjet, genparticles->at(p))<0.3)
        )
            matching_partons.push_back(genparticles->at(p));
    }
    if (matching_partons.size() == 1)
        return std::abs(matching_partons.at(0).pdgId());
    else
        return 0;
}

bool DemoAnalyzer::matchRecoGen(const reco::PFJet& recojet, const reco::GenJet& genjet)
{
    return (
        genjet.pt() > 20
        && (recojet.pt() > 12)
        //&& (std::abs(recojet.eta()) < 1.3)
        && (deltaR(genjet, recojet) < 0.25)
    );
}


// ------------ method called once each job just before starting event loop  ------------
void 
DemoAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DemoAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
DemoAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
DemoAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
DemoAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
DemoAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DemoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DemoAnalyzer);
