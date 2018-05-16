#include "atlas_conf_2017_019.h"
// AUTHOR: K. Rolbiecki
//  EMAIL: krolb@fuw.edu.pl
void Atlas_conf_2017_019::initialize() {
  setAnalysisName("atlas_conf_2017_019");          
  setInformation(""
    "# search for stops with Higgs or Z\n"
    "# 1 lepton and 4 b-jets or\n"
    "# 3 leptons and 1 b-jet\n"
  "");
  setLuminosity(36.1*units::INVFB);      
  bookSignalRegions("SR-3l1b-A;SR-3l1b-B;SR-3l1b-C;SR-1l4b-A;SR-1l4b-B;SR-1l4b-C");
  bookCutflowRegions("00_all;01_3leptons;02_SFOS;03_leadl>40;04_jet0>=60;05_>=1bjet;06_Zpair;07_MET>100;08_A_njets;08_B_njets;08_C_njets;09_A_LeadBjet;09_B_LeadBjet;09_C_LeadBjet;");
  bookCutflowRegions("10_A_LeadJetPT;10_B_LeadJetPT;10_C_LeadJetPT;11_A_ptll;11_B_ptll;11_C_ptll;12_A_met;12_B_met;12_C_met");
  bookCutflowRegions("21_1-2leptons;22_leadl>30;23_bjets>=4;24_A_nj60;24_B_nj60;24_C_nj30;");
  bookCutflowRegions("25_A_mbb;26_A_MET;27_A_HT;28_A_ptbb;");
  bookCutflowRegions("25_B_mT;26_B_MET;");
  bookCutflowRegions("25_C_MET;26_C_LeadB;27_C_MT;");
  // You can also book cutflow regions with bookCutflowRegions("CR1;CR2;..."). Note that the regions are
  //  always ordered alphabetically in the cutflow output files.

  // You should initialize any declared variables here
}

void Atlas_conf_2017_019::analyze() {

  missingET->addMuons(muonsCombined); 
  
  electronsLoose = filterPhaseSpace(electronsLoose, 10., -2.47, 2.47);
  electronsMedium = filterPhaseSpace(electronsMedium, 20., -2.47, 2.47);
  electronsTight = filterPhaseSpace(electronsTight, 20., -2.47, 2.47);
  muonsCombined = filterPhaseSpace(muonsCombined, 10., -2.5, 2.5);
  jets = filterPhaseSpace(jets, 20., -2.8, 2.8); 
  if (electronsMedium.size() + muonsCombined.size() >= 3 ) countCutflowEvent("000_all");
  std::vector<Jet*> bjets;
  std::vector<Jet*> lightjets;
  
  for (int i = 0; i < jets.size(); i++) 
    if ( fabs(jets[i]->Eta) < 2.5 && checkBTag(jets[i]) ) bjets.push_back(jets[i]);
    else lightjets.push_back(jets[i]);  
    
  lightjets = overlapRemoval(lightjets, electronsLoose, 0.2);  
  lightjets = overlapRemoval(lightjets, muonsCombined, 0.2);
  electronsLoose = overlapRemoval(electronsLoose, bjets, 0.4);  
  electronsLoose = overlapRemoval(electronsLoose, lightjets, 0.4); 
  electronsMedium = overlapRemoval(electronsMedium, bjets, 0.4);  
  electronsMedium = overlapRemoval(electronsMedium, lightjets, 0.4);  
  if (electronsMedium.size() + muonsCombined.size() >= 3 ) countCutflowEvent("001_all");
  electronsTight = overlapRemoval(electronsTight, bjets, 0.4);  
  electronsTight = overlapRemoval(electronsTight, lightjets, 0.4);    
  lightjets = overlapRemoval_muon_jet_tracks(lightjets, muonsCombined, 0.4, 2); 
  bjets = overlapRemoval_muon_jet_tracks(bjets, muonsCombined, 0.4, 2);  
  lightjets = overlapRemoval_muon_jet_pt(lightjets, muonsCombined, 0.4); 
  bjets = overlapRemoval_muon_jet_pt(bjets, muonsCombined, 0.4);    
  muonsCombined = overlapRemoval(muonsCombined, bjets, 0.4);  
  muonsCombined = overlapRemoval(muonsCombined, lightjets, 0.4);    
  
  muonsCombined = filterPhaseSpace(muonsCombined, 20., -2.5, 2.5);  
  bjets = filterPhaseSpace(bjets, 30., -2.5, 2.5);
  lightjets = filterPhaseSpace(lightjets, 30., -2.5, 2.5);
  
  electronsMedium = filterIsolation(electronsMedium);
  electronsTight = filterIsolation(electronsTight);
  muonsCombined = filterIsolation(muonsCombined);  
  if (electronsMedium.size() + muonsCombined.size() >= 3 ) countCutflowEvent("002_all");
  electronsMedium = Isolate_leptons_with_inverse_track_isolation_cone(electronsMedium, tracks, towers, 0.2, 10., 0.2, 0.06, 0.06, false);
  electronsTight = Isolate_leptons_with_inverse_track_isolation_cone(electronsTight, tracks, towers, 0.2, 10., 0.2, 0.06, 0.06, false);  
  muonsCombined = Isolate_leptons_with_inverse_track_isolation_cone(muonsCombined,  tracks, towers, 0.3, 10., 0.2, 0.06, 0.06, false);
  
  countCutflowEvent("00_all");
  
// 3l1b selection

  if ( electronsMedium.size() + muonsCombined.size() >= 3 ) {
    countCutflowEvent("01_3leptons");    
    
    bool Zpair = false;  
    bool sfos = false;
    double ptll = 0.;
    double mll = 0.;
    if ( electronsMedium.size() > 1 ) 
      for (int i = 0; i < electronsMedium.size()-1; i++)
        for (int j = i+1 ; j < electronsMedium.size(); j++)
          if ( electronsMedium[i]->Charge * electronsMedium[j]->Charge < 0 ) {
	    sfos = true;
	    double mll_temp = (electronsMedium[i]->P4() + electronsMedium[j]->P4()).M();
	    if ( fabs(mll_temp - 91.2) < 15. )
	      Zpair = true;
	      if ( fabs(mll_temp - 91.2) < fabs(mll - 91.2)  ) { //only keep the closest
		mll = mll_temp;
		ptll = (electronsMedium[i]->P4() + electronsMedium[j]->P4()).Perp();
	      }
	  }

    if ( muonsCombined.size() > 1 ) 
      for (int i = 0; i < muonsCombined.size()-1; i++)
        for (int j = i+1 ; j < muonsCombined.size(); j++)
          if ( muonsCombined[i]->Charge * muonsCombined[j]->Charge < 0 ) {
	    sfos = true;
	    double mll_temp = (muonsCombined[i]->P4() + muonsCombined[j]->P4()).M();
	    if ( fabs(mll_temp - 91.2) < 15. )
	      Zpair = true;
	      if ( fabs(mll_temp - 91.2) < fabs(mll - 91.2)  ) { //only keep the closest
		mll = mll_temp;
		ptll = (muonsCombined[i]->P4() + muonsCombined[j]->P4()).Perp();
	      }
	  } 

    if ( !sfos ) return;
    countCutflowEvent("02_SFOS");

    if ( (electronsMedium.size() > 0 and electronsMedium[0]->PT > 40.) or (muonsCombined.size() > 0 and muonsCombined[0]->PT > 40.) )
      countCutflowEvent("03_leadl>40");
    else return;      
	  
    if ( (bjets.size() > 0 and bjets[0]->PT > 60.) or (lightjets.size() > 0 and lightjets[0]->PT > 60.) ) countCutflowEvent("04_jet0>=60"); 
    else return;
    
    if ( bjets.size() == 0) return;
    countCutflowEvent("05_>=1bjet");    
    
    if ( !Zpair) return;
    countCutflowEvent("06_Zpair");
    
    double met = missingET->P4().Et();
    
    if (met < 100.) return;
    countCutflowEvent("07_MET>100");
    
    if (pass_cuts3l(lightjets, bjets, met, "A", 250., 40., 6, 100., ptll, 150.) )
        countSignalEvent("SR-3l1b-A");
    
    if (pass_cuts3l(lightjets, bjets, met, "B", 80., 40., 6, 180., ptll, 0.) ) 
        countSignalEvent("SR-3l1b-B");    
    
    if (pass_cuts3l(lightjets, bjets, met, "C", 60., 30., 5, 140., ptll, 80.) )
        countSignalEvent("SR-3l1b-C");
	  
    return; // orthogonality?
  }

// 1l4b selection  
  if ( electronsTight.size() + muonsCombined.size() >= 1 ) { // 1 or 2 leptons
    
    countCutflowEvent("21_1-2leptons");
    
    std::vector<FinalStateObject*> leptons;
    for(int e = 0; e < electronsTight.size(); e++) {
      FinalStateObject* lep = newFinalStateObject(electronsTight[e]);
      leptons.push_back(lep);
    }  
    for(int m = 0; m < muonsCombined.size(); m++) {
      FinalStateObject* lep = newFinalStateObject(muonsCombined[m]);
      leptons.push_back(lep);
    }
    std::sort(leptons.begin(), leptons.end(), FinalStateObject::sortByPT );
    
    if ( leptons[0]->PT > 30.) countCutflowEvent("22_leadl>30");
      else return;
    
//    if ( (electronsTight.size() > 0 and electronsTight[0]->PT > 30.) or (muonsCombined.size() > 0 and muonsCombined[0]->PT > 30.) )
//      countCutflowEvent("22_leadl>30");
//    else return;
    
    if (bjets.size() < 4 ) return;
    countCutflowEvent("23_bjets>=4");
    
    double met = missingET->P4().Et();
    
    double ht = 0.;
    int nj60 = 0;
    for (int i = 0; i < lightjets.size(); i++) {
      if ( lightjets[i]->PT > 60. ) nj60++;
      ht += lightjets[i]->PT;
    }  
      
    for (int i = 0; i < bjets.size(); i++) {
      if ( bjets[i]->PT > 60. ) nj60++;
      ht += bjets[i]->PT;
    }
    
//    TLorentzVector leading_lep;
//    if (electronsTight.size() == 0) leading_lep = muonsCombined[0]->P4();
//    else if (muonsCombined.size() == 0) leading_lep = electronsTight[0]->P4();
//    else if (muonsCombined[0]->PT > electronsTight[0]->PT) leading_lep = muonsCombined[0]->P4();
//    else leading_lep = electronsTight[0]->P4();
//    double mtrans = mT(leading_lep, missingET->P4()); 

    double mtrans_min = 1000.;
//    for (int i = 0; i < leptons.size(); i++) {
//       double mtrans_temp = mT(leptons[i]->P4(), missingET->P4());
//       if (mtrans_temp < mtrans_min) mtrans_min = mtrans_temp; 
//    }
    mtrans_min = mT(leptons[0]->P4(), missingET->P4());

    if (nj60 >= 5) {
      countCutflowEvent("24_B_nj60");
      if (mtrans_min > 150.) {
	countCutflowEvent("25_B_mT");
	if (met > 150.) {
          countCutflowEvent("26_B_MET");
          countSignalEvent("SR-1l4b-B");
        }
      }
    }
    
    if (lightjets.size() + bjets.size() >= 7) {
      countCutflowEvent("24_C_nj30");
      if ( met > 150. ) {
	countCutflowEvent("25_C_MET");
	if ( bjets[0]->PT < 140. ) {
	  countCutflowEvent("26_C_LeadB");
	  if ( mtrans_min > 125. ) {
	    countCutflowEvent("27_C_MT");
	    countSignalEvent("SR-1l4b-C");
	  }
	}
      }
    }

    
    double dr = 10.;
    double mbb = 0.;
    double ptbb = 0.;
    for ( int i = 0; i < bjets.size()-1; i++)
      for (int j = i+1; j < bjets.size(); j++) {
	double dr_temp = bjets[i]->P4().DeltaR(bjets[j]->P4());
	if ( dr_temp < dr ) {
	  dr = dr_temp;
	  mbb = (bjets[i]->P4() + bjets[j]->P4()).M(); 
	  ptbb = (bjets[i]->P4() + bjets[j]->P4()).Perp();
	}
      }
    
    if ( nj60 >= 6 ) {
      countCutflowEvent("24_A_nj60");
      if ( fabs(mbb - 125.) < 30. ) {
	countCutflowEvent("25_A_mbb");
	if ( met > 120 ) {
	  countCutflowEvent("26_A_MET");
	  if (ht > 1000. ) {
	    countCutflowEvent("27_A_HT");
	    if ( ptbb > 300.) {
	      countCutflowEvent("28_A_ptbb");
	      countSignalEvent("SR-1l4b-A");
	    }
	  }
	}
      }
    }
    
  }
  
    
  return;  
}

void Atlas_conf_2017_019::finalize() {
  // Whatever should be done after the run goes here
}       

bool Atlas_conf_2017_019::pass_cuts3l(std::vector<Jet*> lightjets_in, std::vector<Jet*> bjets_in, double met_in, std::string sr, double jetpt, double bpt, int njets, double met, double ptll_in, double ptll_cut) {
  
    double bpt_in = bjets_in[0]->PT;
    double lightpt_in = 0.;
    if ( lightjets_in.size() > 0) lightpt_in = lightjets_in[0]->PT;

    if (bjets_in.size() + lightjets_in.size() >= njets) countCutflowEvent("08_"+sr+"_njets");
    else return false;

    if (bpt_in > bpt ) countCutflowEvent("09_"+sr+"_LeadBjet");
    else return false;
    
    if ( bpt_in > jetpt or lightpt_in > jetpt ) countCutflowEvent("10_"+sr+"_LeadJetPT"); 
    else return false;

    if ( ( (sr == "A" or sr == "B") and ptll_in > ptll_cut) or (sr == "C" and ptll_in < ptll_cut) ) countCutflowEvent("11_"+sr+"_ptll");
    else return false;
    
    if (met_in > met) countCutflowEvent("12_"+sr+"_met");
    else return false;  
    
    return true;
}

bool Atlas_conf_2017_019::pass_cuts1l(double mtrans_in, double ht_in, double met_in, std::string sr, double mtrans, double ht, double met) {
  
  if (mtrans_in > mtrans ) countCutflowEvent("25_"+sr+"MT");
  else return false;
  
  if (ht_in > ht ) countCutflowEvent("26_"+sr+"HT");
  else return false;
  
  if (met_in > met ) countCutflowEvent("27_"+sr+"MET");
  else return false;
  
  return true;
}

std::vector<Jet*> Atlas_conf_2017_019::overlapRemoval_muon_jet_pt(std::vector<Jet*> cand_jets, std::vector<Muon*> cand_muons, double deltaR){
  
  std::vector<Jet*> passed;
  for (std::vector<Jet*>::iterator jet = cand_jets.begin(); jet != cand_jets.end(); jet++) {
    
    bool iso = true;
        
    for (std::vector<Muon*>::iterator mu = cand_muons.begin(); mu != cand_muons.end(); mu++) 
      if ((*jet)->P4().DeltaR((*mu)->P4()) < std::min(deltaR, 0.04 + ((*mu)->PT)/10.) and (*mu)->PT > 0.5*(*jet)->PT ) {
	iso = false;
	break;
      }
        		    
    if (iso) passed.push_back(*jet);
  }
  
  return passed;
}

bool Atlas_conf_2017_019::check_nTrack_jet(Jet* jet, std::vector<Track*> tracks, int nTracksMin) {
  
  int nTracks = 0;
  for (std::vector<Track*>::iterator it=tracks.begin(); it!=tracks.end(); it++) 
    for (int part = 0; part < jet->Particles.GetEntries(); part++)
      if (jet->Particles.At(part) == (*it)->Particle && (*it)->PT > 0.5) nTracks++;

    return nTracks > nTracksMin;
}

std::vector<Jet*> Atlas_conf_2017_019::overlapRemoval_muon_jet_tracks(std::vector<Jet*> cand_jets, std::vector<Muon*> cand_muons, double deltaR, int nTracks){
  
  std::vector<Jet*> passed;
  for (std::vector<Jet*>::iterator jet = cand_jets.begin(); jet != cand_jets.end(); jet++) {
  
    if (check_nTrack_jet(*jet, tracks, nTracks)) {
      passed.push_back(*jet);
      continue;
    }
    
    bool iso = true;
        
    for (std::vector<Muon*>::iterator mu = cand_muons.begin(); mu != cand_muons.end(); mu++) 
      if ((*jet)->P4().DeltaR((*mu)->P4()) < std::min(deltaR, 0.04 + (*mu)->PT/10.) ) {
	iso = false;
	break;
      }
        		    
    if (iso) passed.push_back(*jet);
  }
  
  return passed;
}
