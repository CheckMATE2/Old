#ifndef ATLAS_CONF_2017_019_H_
#define ATLAS_CONF_2017_019_H_
// AUTHOR: K. Rolbiecki
//  EMAIL: krolb@fuw.edu.pl
#include "AnalysisBase.h"

class Atlas_conf_2017_019 : public AnalysisBase {
  public:
    Atlas_conf_2017_019() : AnalysisBase()  {}               
    ~Atlas_conf_2017_019() {}
  
    void initialize();
    void analyze();        
    void finalize();

  private:
    
    bool pass_cuts3l(std::vector<Jet*> lightjets_in, std::vector<Jet*> bjets_in, double met_in, std::string sr, double jetpt, double bpt, int njets, double met, double ptll_in, double ptll_cut);
    
    bool pass_cuts1l(double mtrans_in, double ht_in, double met_in, std::string sr, double mtrans, double ht, double met);
    
    bool check_nTrack_jet(Jet* jet, std::vector<Track*> tracks, int nTracksMin) ;
    
    std::vector<Jet*> overlapRemoval_muon_jet_tracks(std::vector<Jet*> cand_jets, std::vector<Muon*> cand_muons, double deltaR, int nTracks);

    std::vector<Jet*> overlapRemoval_muon_jet_pt(std::vector<Jet*> cand_jets, std::vector<Muon*> cand_muons, double deltaR);    
};

#endif
