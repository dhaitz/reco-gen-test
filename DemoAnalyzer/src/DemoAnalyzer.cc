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

      // ----------member data ---------------------------

  TProfile *recogen_genpt;
  TProfile *npv_genpt;
  TH1D *recogen3040;
  TH1D *recogen5080;
  TH1D *recogen100;
  TH1D *ids2030;
  TH1D *recopt;
  TProfile *recogen_npv;

  TH1D *ngenjets;
  TH1D *nrecojets;
  TH1D *njetsRGratio;

  //reco jet properties
  TProfile *mu_genpt;
  TProfile *e_genpt;
  TProfile *ph_genpt;
  TProfile *chargedHad_genpt;
  TProfile *neutralHad_genpt;

  TProfile *nCharged_genpt;
  TProfile *nConst_genpt;

  //gen
  TProfile *emGen_genpt;
  TProfile *hadGen_genpt;
  TProfile *recogen_emGen;
  TProfile *recogen_hadGen;

  std::string jetCollection;

  std::vector<int> m_idvector;

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
jetCollection(iConfig.getParameter<edm::InputTag>("jetCollection").label())
{
    m_idvector =   {11, 13, 111, 130, 211, 2112, 22, 2212, 310, 3112, 3122, 321, 3222, 3312, 3322};

    int min = 20;
    int max = 100;
    int nbins = 8;

    edm::Service<TFileService> fs;
    recogen_genpt = fs->make<TProfile>("recogen_genpt", "recogen_genpt", nbins, min, max);
    recopt = fs->make<TH1D>("recopt", "recopt", 20, 0, 100);
    npv_genpt = fs->make<TProfile>("npv", "npv", nbins, min, max);

    recogen3040 = fs->make<TH1D>("recogen3040", "recogen3040" , 30 , 0 , 3);
    recogen5080 = fs->make<TH1D>("recogen5080", "recogen5080" , 30 , 0 , 3);
    recogen100 = fs->make<TH1D>("recogen100", "recogen100" , 30 , 0 , 3);

    ids2030 = fs->make<TH1D>("ids2030", "ids2030", 20, -0.5, 19.5);
    recogen_npv = fs->make<TProfile>("recogen_npv", "recogen_npv", 7, -0.5, 34.5);

    ngenjets = fs->make<TH1D>("ngenjets", "ngenjets", 100, -0.5, 99.5);
    nrecojets = fs->make<TH1D>("nrecojets", "nrecojets", 100, -0.5, 99.5);
    njetsRGratio = fs->make<TH1D>("njetsratio", "njetsratio", 100, -0.5, 99.5);

    //gen jet
    emGen_genpt = fs->make<TProfile>("emGen", "emGen" , nbins, min, max);
    hadGen_genpt = fs->make<TProfile>("hadGen", "hadGen", nbins, min, max);

    // recogen as function of gen composition
    recogen_emGen = fs->make<TProfile>("recogen_emGen", "recogen_emGen" , 10 , 0 , 1);
    recogen_hadGen = fs->make<TProfile>("recogen_hadGen", "recogen_hadGen" , 10 , 0 , 1);

    //reco jet properties
    mu_genpt = fs->make<TProfile>("mu", "mu", nbins, min, max);
    chargedHad_genpt = fs->make<TProfile>("chargedHad", "chargedHad", nbins, min, max);
    neutralHad_genpt = fs->make<TProfile>("neutralHad", "neutralHad", nbins, min, max);
    e_genpt = fs->make<TProfile>("e", "e", nbins, min, max);
    ph_genpt = fs->make<TProfile>("ph", "ph", nbins, min, max);

    nConst_genpt = fs->make<TProfile>("nConst", "nConst", nbins, min, max);
    nCharged_genpt = fs->make<TProfile>("nCharged", "nCharged", nbins, min, max);
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

Handle<reco::PFJetCollection> recojets;
iEvent.getByLabel(jetCollection.c_str(), recojets);

Handle<reco::GenJetCollection> genjets;
iEvent.getByLabel("ak5GenJetsNoNu", genjets);


Handle<reco::GenParticleCollection> genparticles;
iEvent.getByLabel("genParticles", genparticles);
std::vector<reco::GenParticle> leptons;

for (unsigned int i=0; i < genparticles->size(); i++)
{
    if (leptons.size() <3 &&  i>5 && std::abs(genparticles->at(i).pdgId()) <20 && std::abs(genparticles->at(i).pdgId())>10)
        leptons.push_back(genparticles->at(i));
}

ngenjets->Fill(genjets->size());
nrecojets->Fill(recojets->size());

if (genjets->size() > 0)
njetsRGratio->Fill(recojets->size()/genjets->size());

for (unsigned int g = 0; g < genjets->size(); g++)
{
    for (unsigned int r = 0; r < recojets->size(); r++)
    {
        if (genjets->at(g).pt() > 20
            && (recojets->at(r).pt() > 12)
            && (std::abs(recojets->at(r).eta()) < 1.3)
            && deltaR(genjets->at(g), recojets->at(r)) < 0.25
            //&& (deltaR(genjets->at(g), leptons.at(0)) > 0.5 && deltaR(genjets->at(g), leptons.at(1)) > 0.5)
        )
        {
            recogen_genpt->Fill(genjets->at(g).pt(), recojets->at(r).pt()/genjets->at(g).pt());
            npv_genpt->Fill(genjets->at(g).pt(), vertices->size());
            recogen_npv->Fill(vertices->size(), recojets->at(r).pt()/genjets->at(g).pt());

            recopt->Fill(recojets->at(r).pt());

            //reco jet properties
            mu_genpt->Fill(genjets->at(g).pt(), recojets->at(r).muonEnergyFraction());
            e_genpt->Fill(genjets->at(g).pt(), recojets->at(r).electronEnergyFraction());
            ph_genpt->Fill(genjets->at(g).pt(), recojets->at(r).photonEnergyFraction());
            chargedHad_genpt->Fill(genjets->at(g).pt(), recojets->at(r).chargedHadronEnergyFraction());
            neutralHad_genpt->Fill(genjets->at(g).pt(), recojets->at(r).neutralHadronEnergyFraction());
            nConst_genpt->Fill(genjets->at(g).pt(), recojets->at(r).nConstituents());
            nCharged_genpt->Fill(genjets->at(g).pt(), recojets->at(r).chargedMultiplicity());

            //gen jet properties
            emGen_genpt->Fill(genjets->at(g).pt(), genjets->at(g).emEnergy()/genjets->at(g).energy());
            hadGen_genpt->Fill(genjets->at(g).pt(), genjets->at(g).hadEnergy()/genjets->at(g).energy());

            recogen_emGen->Fill(genjets->at(g).emEnergy()/genjets->at(g).energy(), recojets->at(r).pt()/genjets->at(g).pt());
            recogen_hadGen->Fill(genjets->at(g).hadEnergy()/genjets->at(g).energy(), recojets->at(r).pt()/genjets->at(g).pt());

            //only fill these histograms if gen jet pt between 20 and 30
            if (genjets->at(g).pt() > 30 && genjets->at(g).pt() < 40)
            {
                recogen3040->Fill(recojets->at(r).pt()/genjets->at(g).pt());

                // ids
                std::vector< const reco::GenParticle * > gens = genjets->at(g).getGenConstituents();
                for (unsigned int i = 0; i < gens.size(); i++)
                {
                    std::vector<int>::iterator it = find(m_idvector.begin(), m_idvector.end(), std::abs(gens.at(i)->pdgId()));
                    unsigned int pos = it - m_idvector.begin();
                    if (pos >= m_idvector.size())
                        std::cout << std::abs(gens.at(i)->pdgId()) << std::endl;
                    else
                        ids2030->Fill(pos);
                }
            }
            else if (genjets->at(g).pt() > 50 && genjets->at(g).pt() < 80)
                recogen5080->Fill(recojets->at(r).pt()/genjets->at(g).pt());
            else if (genjets->at(g).pt() > 100)
                recogen100->Fill(recojets->at(r).pt()/genjets->at(g).pt());
        }
    }
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
