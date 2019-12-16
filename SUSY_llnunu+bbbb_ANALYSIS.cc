// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"

namespace Rivet {

  using namespace Cuts;

  class SUSY_HH_ANALYSIS : public Analysis {
  public:

    /// Constructor
    SUSY_HH_ANALYSIS()
      : Analysis("SUSY_HH_ANALYSIS")
    {    }


    /// @name Analysis methods
    //@{

    // Book histograms
    void init() {
      // Basic final state
      const FinalState fs(-5.0, 5.0, 0.5*GeV);

      nvetoed1 = 0;
      nvetoed2 = 0;
      nvetoed3 = 0;
      nvetoed4 = 0;
      nvetoed5 = 0;
      nvetoed6 = 0;
      nvetoed7 = 0;
      nvetoed = 0;
      
      // Tracks and jets
      addProjection(ChargedFinalState(fs), "Tracks");
      VetoedFinalState vfs;
      vfs.addVetoPairId(PID::MUON);
      vfs.addVetoPairId(PID::ELECTRON);
      addProjection(FastJets(vfs, FastJets::ANTIKT, 0.4), "Jets");

    
      IdentifiedFinalState muelfs(etaIn(-5, 5) & (pT >= 0.5*GeV));
      muelfs.acceptIdPair(PID::MUON);
      muelfs.acceptIdPair(PID::ELECTRON);
      addProjection(muelfs, "ElecMuons");


      MissingMomentum missing(fs);
      addProjection(missing, "MET");

      // for HT
      addProjection(VisibleFinalState(etaIn(-5, 5) & (pT >= 0.5*GeV)),"vfs");


      _hist_n_trk   = bookHisto1D("n-trk", 50, 0.5, 300.5);
      _hist_phi_trk = bookHisto1D("phi-trk", 50, -PI, PI);
      _hist_eta_trk = bookHisto1D("eta-trk", 50, -4, 4);
      _hist_pt_trk  = bookHisto1D("pt-trk", 100, 0.0, 1500);

      _hist_n_jet   = bookHisto1D("n-jet", 21, -0.5, 20.5);
      _hist_phi_jet = bookHisto1D("phi-jet", 50, -PI, PI);
      _hist_eta_jet = bookHisto1D("eta-jet", 50, -4, 4);
      _hist_pt_jet  = bookHisto1D("pt-jet", 100, 0.0, 1500);

      _hist_n_bjet   = bookHisto1D("n-bjet", 21, -0.5, 20.5);
      _hist_phi_bjet = bookHisto1D("phi-bjet", 50, -PI, PI);
      _hist_eta_bjet = bookHisto1D("eta-bjet", 50, -4, 4);
      _hist_pt_bjet  = bookHisto1D("pt-bjet", 100, 0.0, 1500);

      _hist_n_l   = bookHisto1D("n-l", 11, -0.5, 10.5);
      _hist_phi_l = bookHisto1D("phi-l", 50, -PI, PI);
      _hist_eta_l = bookHisto1D("eta-l", 50, -4, 4);
      _hist_pt_l  = bookHisto1D("pt-l", 100, 0.0, 500);

      _hist_met = bookHisto1D("Etmiss", 100, 0.0, 800);
      
      _hist_met_incl = bookHisto1D("Etmiss_incl", 100, 0.0, 1200);

      _hist_mll = bookHisto1D("mll", 50, 0.0, 500);
     
      _hist_ht = bookHisto1D("HT", 100, 0, 1500);

      _hist_deltaPhi = bookHisto1D("deltaPhi", 50, 0, PI);
      _hist_deltaPhi_cross = bookHisto1D("deltaPhi_cross", 50, 0, PI);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();
      
      const Particles lpfs = applyProjection<FinalState>(event, "ElecMuons").particlesByPt(pT>=15*GeV);
      const FastJets& jetpro = applyProjection<FastJets>(event, "Jets");
      const Jets jets = jetpro.jetsByPt(20*GeV);
      
      //cout << "n leptons before isolation: " << lpfs.size() << endl;
      ParticleVector leptons;
      foreach( const Particle& l, lpfs) { 
		int nnotiso=0;
	    foreach (const Jet& j, jets) {
          if (deltaR(l.momentum(),j.momentum())<=0.1) nnotiso++;
	    }
        if (nnotiso==0)  leptons.push_back(l);
	    
	  }
	  
	   
      if (leptons.size() !=2) {
        nvetoed1++;
        vetoEvent;
      }

      if (leptons[0].charge() * leptons[1].charge() != -1){
        nvetoed2++;
        vetoEvent;
      } 

                        
      // Calculate missing Et
      const MissingMomentum& met = applyProjection<MissingMomentum>(event, "MET");
      
      // met inclusive
      _hist_met_incl->fill(met.vectorEt().mod()/GeV, weight);
      
      if (met.vectorEt().mod()<=10*GeV) {
		nvetoed3++;
		vetoEvent;
	  }

      int bjets=0;
      foreach (const Jet& j, jets) {
		  if (j.bTagged());
			bjets++;
	  }
	  
      //cout << bjets <<endl;
      if (bjets!=4) {
		  nvetoed4++;
		  vetoEvent;
	  }

      //mll cut
      double mll = (leptons[0].momentum()+leptons[1].momentum()).mass();
      if (mll<10*GeV) {
		  nvetoed5++;
		  vetoEvent;
	  }

      
      if (deltaPhi(leptons[0].momentum(), leptons[1].momentum())>2) {
        nvetoed6++;
		vetoEvent;
	  }
      

      // HT calculation and histo
      Particles visible = applyProjection<VisibleFinalState>(event, "vfs").particles();
      double HT=0.0;
      foreach(Particle &p, visible) {
        HT += p.perp();
      }

      if (HT>800*GeV) {
	nvetoed7++;
	vetoEvent;
	}

      nvetoed=nvetoed1+nvetoed2+nvetoed3+nvetoed4+nvetoed5+nvetoed6+nvetoed7;
	  
	  //met histo  
      _hist_met->fill(met.vectorEt().mod()/GeV, weight);

	  //deltaPhi histo
      _hist_deltaPhi -> fill(deltaPhi(leptons[0].momentum(), leptons[1].momentum()), weight);
      _hist_deltaPhi_cross -> fill(deltaPhi(leptons[0].momentum(), leptons[1].momentum()), weight);

      
      const FinalState& tracks = applyProjection<FinalState>(event, "Tracks");

      // Fill track histos
      _hist_n_trk->fill(tracks.size(), weight);
      foreach (const Particle& t, tracks.particles()) {
        const FourMomentum& p = t.momentum();
        _hist_phi_trk->fill(mapAngleMPiToPi(p.phi()), weight);
        _hist_eta_trk->fill(p.eta(), weight);
        _hist_pt_trk->fill(p.pT()/GeV, weight);
      }

      // Get jets and fill jet histos
      _hist_n_jet->fill(jets.size(), weight);
      foreach (const Jet& j, jets) {
        const FourMomentum& pj = j.momentum();
        _hist_phi_jet->fill(mapAngleMPiToPi(pj.phi()), weight);
        _hist_eta_jet->fill(pj.eta(), weight);
        _hist_pt_jet->fill(pj.pT()/GeV, weight);
      }


      // Get bjets and fill bjet histos
      foreach (const Jet& j, jets) {
	    if (j.bTagged()) _hist_n_bjet->fill(jets.size(), weight);
        const FourMomentum& pj = j.momentum();
        _hist_phi_bjet->fill(mapAngleMPiToPi(pj.phi()), weight);
        _hist_eta_bjet->fill(pj.eta(), weight);
        _hist_pt_bjet->fill(pj.pT()/GeV, weight);
      }


      // Fill final state lepton histos
      _hist_n_l->fill(leptons.size(), weight);
      vector<FourMomentum> lpluses, lminuses;
      foreach (const Particle& l , leptons) {
        const FourMomentum& p = l.momentum();
        _hist_phi_l->fill(mapAngleMPiToPi(p.phi()), weight);
        _hist_eta_l->fill(p.eta(), weight);
        _hist_pt_l->fill(p.pT()/GeV, weight);
      }


      // lepton pair mass histo
      _hist_mll->fill(mll, weight);

	//HT histo
      _hist_ht->fill(HT, weight);
    }


    void finalize() {
      scale(_hist_met, crossSection()/sumOfWeights()); 
      normalize(_hist_met_incl);
      scale(_hist_mll, crossSection()/sumOfWeights());
      scale(_hist_ht, crossSection()/sumOfWeights());

      //normalize(_hist_deltaPhi);

      scale(_hist_deltaPhi_cross, crossSection()/sumOfWeights());
      normalize(_hist_n_trk);
      normalize(_hist_phi_trk);
      normalize(_hist_eta_trk);
      normalize(_hist_pt_trk);
      normalize(_hist_n_jet);
      scale(_hist_phi_jet, crossSection()/sumOfWeights());
      scale(_hist_eta_jet, crossSection()/sumOfWeights());
      scale(_hist_pt_jet, crossSection()/sumOfWeights());
      normalize(_hist_n_l);
      scale(_hist_phi_l, crossSection()/sumOfWeights());
      scale(_hist_eta_l, crossSection()/sumOfWeights());
      scale(_hist_pt_l, crossSection()/sumOfWeights());
      normalize(_hist_n_bjet);
      scale(_hist_phi_bjet, crossSection()/sumOfWeights());
      scale(_hist_eta_bjet, crossSection()/sumOfWeights());
      scale(_hist_pt_bjet, crossSection()/sumOfWeights());
      cout << "Vetoed " << nvetoed1 << " events" << endl;
      cout << "Vetoed " << nvetoed2 << " events" << endl;
      cout << "Vetoed " << nvetoed3 << " events" << endl;
      cout << "Vetoed " << nvetoed4 << " events" << endl;
      cout << "Vetoed " << nvetoed5 << " events" << endl;
      cout << "Vetoed " << nvetoed6 << " events" << endl;
      cout << "Vetoed " << nvetoed7 << " events" << endl;
      cout << "Total Vetoed " << nvetoed << " events" << endl;
      /// @todo Normalisations
    }

    //@}


  private:
    Histo1DPtr _hist_n_trk, _hist_phi_trk, _hist_eta_trk, _hist_pt_trk;
    Histo1DPtr _hist_n_jet, _hist_phi_jet, _hist_eta_jet, _hist_pt_jet;
    Histo1DPtr _hist_n_bjet, _hist_phi_bjet, _hist_eta_bjet, _hist_pt_bjet;
    Histo1DPtr _hist_n_l, _hist_phi_l, _hist_eta_l, _hist_pt_l;
    Histo1DPtr _hist_met;
    Histo1DPtr _hist_met_incl;
    Histo1DPtr _hist_mll;    
    Histo1DPtr _hist_ht;
    Histo1DPtr _hist_deltaPhi, _hist_deltaPhi_cross;    
    int nvetoed1;
    int nvetoed2;
    int nvetoed3;
    int nvetoed4;
    int nvetoed5;
    int nvetoed6;
    int nvetoed7;
    int nvetoed;
  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(SUSY_HH_ANALYSIS);

}
