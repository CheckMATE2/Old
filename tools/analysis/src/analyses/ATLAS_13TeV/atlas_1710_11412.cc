#include "atlas_1710_11412.h"
// AUTHOR: J.S.Kim
//  EMAIL: jsk@th.physik.uni-bonn.de
void Atlas_1710_11412::initialize() {
  setAnalysisName("atlas_1710_11412");          
  setInformation(""
    "# Search for dark matter produced in association with bottom or top quarks in sqrt{s}V pp collisions with the ATLAS detector\n"
  "");
  setLuminosity(36.1*units::INVFB);      
  bookSignalRegions("SRb1;SRb2-bin1;SRb2-bin2;SRb2-bin3;SRb2-bin4;SRt1;SRt2;SRt3");
  bookCutflowRegions("1SRb1;2SRb1;3SRb1;4SRb1;5SRb1;1SRb2;2SRb2;3SRb3;4SRb4;5SRb2;6SRb2;6SRb2");
  // You can also book cutflow regions with bookCutflowRegions("CR1;CR2;..."). Note that the regions are
  //  always ordered alphabetically in the cutflow output files.

  // You should initialize any declared variables here
}

void Atlas_1710_11412::analyze() {
  // Your eventwise analysis code goes here
  // The following objects are always defined unless they are 'ignored' above. They form std::vector objects of the respective Delphes class type (except for Etmiss which is a single object)
  // All std::vector members and etmiss have the common properties PT, Eta, Phi and P4() with the latter giving access to the full ROOT TLorentzVector.
  // Within a std::vector, all members are ordered with highest pt coming first.

  // electronsLoose, electronsMedium, electronsTight   are list of electrons that passed respective efficiency and reconstruction cuts
  // muonsCombinedPlus, muonsCombined                  as above for muons
  // photonsMedium                                     as above for photons
  // jets are all reconstructed jets                   as above for jets. Note that electrons are most certainly also reconstructed as a jet -> overlap removal do avoid double counting necessary!
  // tracks, towers                                    calorimeter and tracker information. Usually not needed.
  // missingET                                         rec missing ET EXCLUDING muons.

  
  // Here is a couple of useful functions and lines:  
  //------------Phase Space Cuts (defined for jets, electronsXYZ, muonsXYZ, photonsXYZ)
  // jets = filterPhaseSpace(jets, 20., -2.8, 2.8)  // The vector 'jets' only contains jets with pt >= 20 GeV and -2.8 < eta < 2.8. This function is applicable to other particles too (electronsMedium, ... ).
  // jets = overlapRemoval(jets, electronsLoose, 0.2) Removes all jets for which there exists any electron in 'electronsLoose' with deltaR < 0.2.
  // jets = overlapRemovel(jets, 0.2) If two jets overlap within deltaR < 0.2, only the harder jet is stored.
  
  //------------Isolation Checks (defined for electronsXYZ, muonsXYZ, photonsXYZ
  //------------        For each object, if the user entered N isolation conditions, they can be
  //------------        checked individually be the second argument (running from 0 to N-1).
  // electronsMedium = filterIsolation(electronsMedium, 0)            Removes electrons that do not pass the first isolation condition entered into the AnalysisManager by the user
  // std::vector<int> flags; flags.push_back(0); flags.push_back(2);
  // electronsMedium = filterIsolation(electronsMedium, flags)        Same as above, but both the first and the third condition have to be fulfilled
  // electronsMedium = filterIsolation(electronsMedium)               Same as above, but all conditions have to be fulfilled.
  
  //-----------Flavour Tag Checks (defined for jets only)
  //----------          Tau tags "loose", "medium" or "tight" can be individually checked if the user enabled tau tagging in the AM.
  //----------          For b-tags, if N working points have been defined, the ith condition can be tested by putting i-1 into the second argument (if there is only one, the argument can be omitted)
  // if checkTauTag(jets[0], "tight") leadingJetIsTagged = True;
  // if checkBTag(jets[0], 0) leadingJetIsBTagged = True;


  //-----------Auxiliary Information
  // - Always ensure that you don't access vectors out of bounds. E.g. 'if(jets[1]->PT > 150)' should rather be if (jets.size() > 1 && jets[1]->PT > 150). 
  // - Use rand()/(RAND_MAX+1.) for random numbers between 0 and 1. The random seed is determined from system time or by the RandomSeed parameter in CheckMATE.
  // - The 'return' statement will end this function for the current event and hence should be called whenever the current event is to be vetoed.
  // - Many advanced kinematical functions like mT2 are implemented. Check the manual for more information.
  // - If you need output to be stored in other files than the cutflow/signal files we provide, check the manual for how to do this conveniently.  

  missingET->addMuons(muonsCombined);  // Adds muons to missing ET. This should almost always be done which is why this line is not commented out.
  jets = filterPhaseSpace(jets, 20., -2.8, 2.8);
  muonsCombinedPlus = filterPhaseSpace(muonsCombinedPlus, 8., -2.5, 2.5);
  electronsLoose = filterPhaseSpace(electronsLoose, 7., -2.47, 2.47);
  electronsTight = filterPhaseSpace(electronsTight, 7., -2.47, 2.47);

  electronsLoose = filterIsolation(electronsLoose);
  electronsTight = filterIsolation(electronsTight);
  muonsCombinedPlus = filterIsolation(muonsCombinedPlus);

  jets=overlapRemoval(jets,muonsCombinedPlus,0.2);
  jets=overlapRemoval(jets,electronsLoose,0.2);
  muonsCombinedPlus=overlapRemoval(muonsCombinedPlus, jets, 0.4);
  electronsLoose=overlapRemoval(electronsLoose, jets, 0.4);
  electronsTight=overlapRemoval(electronsTight, jets, 0.4);

  muonsCombinedPlus = filterPhaseSpace(muonsCombinedPlus, 25., -2.5, 2.5);
  electronsLoose = filterPhaseSpace(electronsLoose, 25., -2.47, 2.47);
  electronsTight = filterPhaseSpace(electronsTight, 25., -2.47, 2.47);

  
  //JSK:overlap removel is missing
  //Photon definition
  //no distinction between baseline and signal leptons!
 
  std::vector<Jet*> bjetsL;
  for( int i=0; i<jets.size(); i++ ){
    if( checkBTag(jets[i],2) )  bjetsL.push_back(jets[i]);
  }

  std::vector<Jet*> bjetsM;
  for( int i=0; i<jets.size(); i++ ){
    if( checkBTag(jets[i],1) )  bjetsM.push_back(jets[i]);
  }

  std::vector<Jet*> bjetsT;
  for( int i=0; i<jets.size(); i++ ){
    if( checkBTag(jets[i],0) )  bjetsT.push_back(jets[i]);
  }

  std::vector<Jet*> nobjets;
  for( int i=0; i<jets.size(); i++ ){
    if( !checkBTag(jets[i],2) )  nobjets.push_back(jets[i]);
  }

  std::vector<Jet*> taujets;
  for(int i=0; i<jets.size(); i++) {
    if(checkTauTag(jets[i], "loose")) {
      taujets.push_back(jets[i]);
    }
  }

  
  //
  //SRb1
  //
    
  if( electronsLoose.size()+muonsCombinedPlus.size() == 0 && bjetsT.size() >= 1 && bjetsT[0]->PT > 160. ){
    countCutflowEvent("1SRb1");
    if( jets.size() > 1)
      {
	countCutflowEvent("2SRb1");
	if (jets[0]->PT > 160. )
	  {
	    countCutflowEvent("3SRb1");
	    if (jets[1]->PT > 160. ){
	      //calculation of scalar sum HT3
	      double HT3=0.;
	      for(int i=2; i<jets.size(); i++) {
		HT3 += jets[i]->PT;
	      }
	      if( HT3  < 100. ){
		countCutflowEvent("4SRb1");
		if( missingET->P4().Et() > 650. )
		  {
		    countCutflowEvent("5SRb1");
		    double deltaphi_min=10000.;
		    double deltaphi;
		    for(int i=0; i<jets.size(); i++){
		      deltaphi=fabs(jets[i]->P4().DeltaPhi(missingET->P4()));
		      if( deltaphi < deltaphi_min && i < 2)
			deltaphi_min=deltaphi;
		    } //loop
		    if( deltaphi_min > 0.6){
		      countSignalEvent("SRb1");
		    }//deltaphi cut
		    
		  }//met cut
	      } //HT3cut
	      
	    }//jet2cut
	  }//jet1cut
	
      }//numberjetscuts
    
  }//multiplicity cut


  //
  //SRb2
  //
    
  if( electronsLoose.size()+muonsCombinedPlus.size() == 0 && bjetsT.size() >= 2 && bjetsT[0]->PT > 150. ){
    countCutflowEvent("1SRb2");
    if( missingET->P4().Et() > 180. )
      {
	countCutflowEvent("2SRb2");
	if ( jets[0]->PT > 150. && ( jets.size() == 2 ) || ( jets.size() == 3 && jets[2]->PT < 60. ) )
	  {
	    countCutflowEvent("3SRb2");
	    double deltaphi;
	    double deltaphimin=1000.;
	    for(int i=0; i<jets.size(); i++){
	      deltaphi=jets[i]->P4().DeltaPhi(missingET->P4());
	      if( deltaphi < deltaphimin && i < 3)
		deltaphimin=deltaphi;
	    } //for loop
	    double deltaphimin_bb;
	    deltaphimin_bb= ( bjetsT[0]->P4().DeltaPhi(bjetsT[1]->P4()) );
	    double deltaplus,deltaminus;

	    deltaminus=deltaphimin - deltaphimin_bb;
	    deltaplus=fabs( - deltaphimin - deltaphimin_bb + 3.14159265359); 
	    
	    if ( deltaplus < 1000.5 ){
	      countCutflowEvent("4SRb2");
	      if( deltaminus < 1000.){
		countCutflowEvent("5SRb2");
		//calculation of scalar sum HT ratio
		double HT=0.;
		for(int i=0; i<jets.size(); i++) {
		  HT += jets[i]->PT;
		} //for loop
		double HTratio;
		HTratio=jets[0]->PT/HT;
		if( HTratio > 0.75){
		  countCutflowEvent("6SRb2");
		  //		  if ( deltaphimin > 0.0){
		  if ( 1 ){
		    countCutflowEvent("7SRb2");

		    double costheta_bb;
		    costheta_bb=fabs(tanh((bjetsT[0]->Eta-bjetsT[1]->Eta)/2.));

		    if( costheta_bb >= 0.00 && costheta_bb <= 0.25 )
		      countSignalEvent("SRb2-bin1");

		    if( costheta_bb = 0.25 && costheta_bb <= 0.5 )
		      countSignalEvent("SRb2-bin2");

		    if( costheta_bb = 0.5 && costheta_bb <= 0.75 )
		      countSignalEvent("SRb2-bin3");
		    
		    if( costheta_bb = 0.75 && costheta_bb <= 1.0 )
		      countSignalEvent("SRb2-bin4");
				    
		    
		    
		  }//deltaphi cut
		}//HTratio cut
	      }//deltaminus cut
	    }//deltaplus cut
	  }//jets multicplicity cut
      }//met cut
  }//lepton veto, bjets multiplicity cut
	      

  //
  //SRt1
  //

  if( electronsLoose.size()+muonsCombinedPlus.size() == 0 && bjetsM.size() > 1 && jets.size() > 3 && jets[1]->PT > 80 && jets[3]->PT > 40 ){
    countCutflowEvent("1SRt1");
    if( missingET->P4().Et()  > 250. ){
          countCutflowEvent("2SRt1");
      double deltaphi_min=10000.;
      for( int i=0; i<jets.size(); i++){
	if( deltaphi_min > fabs(jets[i]->P4().DeltaPhi(missingET->P4())) )
	  deltaphi_min=fabs(jets[i]->P4().DeltaPhi(missingET->P4()));
      }
      if(deltaphi_min > 0.4){
	countCutflowEvent("3SRt1");	
	countCutflowEvent("4SRt1");
	countCutflowEvent("5SRt1");
	
	if( taujets.size() == 0 ){
	  countCutflowEvent("6SRt1");

	  //now recluster jets
	  std::vector<fastjet::PseudoJet> particles;
	  for ( int i = 0; i < jets.size(); i++ )
	    particles.push_back( fastjet::PseudoJet(jets[i]->P4().Px(), jets[i]->P4().Py(), jets[i]->P4().Pz(), jets[i]->P4().E()) );

	  //	  fastjet::JetDefinition jet_def12(fastjet::antikt_algorithm, 1.2);
	  fastjet::JetDefinition jet_def08(fastjet::antikt_algorithm, 0.8);

	  //	  fastjet::ClusterSequence cs12(particles, jet_def12);
	  fastjet::ClusterSequence cs08(particles, jet_def08);
	  //	  std::vector<fastjet::PseudoJet> largeR12 = sorted_by_pt(cs12.inclusive_jets());
	  std::vector<fastjet::PseudoJet> largeR08 = sorted_by_pt(cs08.inclusive_jets());

	  if(largeR08.size() > 1 && largeR08[0].pt() > 80. ){
	    countCutflowEvent("7SRt1");
	    if(largeR08.size() > 1 && largeR08[1].pt() > 80. ){
	      countCutflowEvent("8SRt1");
	      double temp;
	      double deltaRmax=0.;
	      double deltaRmin=100000.;
	      int imin,imax;
	      for(int i=0;i < bjetsM.size(); i++){
		temp=bjetsM[i]->P4().DeltaR(missingET->P4());
		if(deltaRmin > temp){
		  deltaRmin = temp;
		  imin=i;
		}
		if(deltaRmax < temp){
		  deltaRmax = temp;
		  imax=i;
		}
	      }//loop
	      double mTmin,mTmax;
	      mTmin=mT(bjetsM[imin]->P4(), missingET->P4());
	      mTmax=mT(bjetsM[imax]->P4(), missingET->P4());
	      if( mTmax > 250 ){
		countCutflowEvent("9SRt1");
		if( mTmin > 150 ){
		  countCutflowEvent("10SRt1");
		  if( missingET->P4().Et() > 300. ){
		    countCutflowEvent("11SRt1");
		    
		    double deltaR_bb;
		    double deltaRmin_bb=1000.;
		    for(int i=0; i<bjetsM.size(); i++){
		      deltaR_bb=bjetsM[i]->P4().DeltaR(missingET->P4());
		      if( deltaR_bb < deltaRmin_bb )
			deltaRmin_bb=deltaR_bb;
		    } //for loop                                                                           
		    if( deltaRmin_bb > 1.5){
		      countCutflowEvent("12SRt1");
		      countSignalEvent("SRt1");
		    } // deltaRbb cut
		  } //met cut
		}//mTmin cut
	      }//mTmax cut
		 
		
	    }//large jet 2 PT cut
	  }//large jet 1 PT cut
	}//taujetscuts
      } //deltaphi_min cut
    } //met cut

  } // leptonveto, jets multiplicities


  //
  //SRt2
  //

  if( electronsLoose.size()+muonsCombinedPlus.size() == 0 && bjetsM.size() > 1 && jets.size() > 3 && jets[1]->PT > 80 && jets[3]->PT > 40 ){
    countCutflowEvent("1SRt2");
    if( missingET->P4().Et()  > 250. ){
      countCutflowEvent("2SRt2");
      double deltaphi_min=10000.;
      for( int i=0; i<jets.size(); i++){
	if( deltaphi_min > fabs(jets[i]->P4().DeltaPhi(missingET->P4())) )
	  deltaphi_min=fabs(jets[i]->P4().DeltaPhi(missingET->P4()));
      }
      if(deltaphi_min > 0.4){
	countCutflowEvent("3SRt2");
	//no emiss track
	countCutflowEvent("4SRt2");
	//no delta phi cut  on ptmiss, ptmisstrack 
	countCutflowEvent("5SRt2");
	
	if( taujets.size() == 0 ){
	  countCutflowEvent("6SRt2");

	  //now recluster jets
	  std::vector<fastjet::PseudoJet> particles;
	  for ( int i = 0; i < jets.size(); i++ )
	    particles.push_back( fastjet::PseudoJet(jets[i]->P4().Px(), jets[i]->P4().Py(), jets[i]->P4().Pz(), jets[i]->P4().E()) );

	  fastjet::JetDefinition jet_def12(fastjet::antikt_algorithm, 1.2);
	  //fastjet::JetDefinition jet_def08(fastjet::antikt_algorithm, 0.8);

	  fastjet::ClusterSequence cs12(particles, jet_def12);
	  //fastjet::ClusterSequence cs08(particles, jet_def08);
	  std::vector<fastjet::PseudoJet> largeR12 = sorted_by_pt(cs12.inclusive_jets());
	  //std::vector<fastjet::PseudoJet> largeR08 = sorted_by_pt(cs08.inclusive_jets());

	  if(largeR12.size() > 1 && largeR12[0].pt() > 140. ){
	    countCutflowEvent("7SRt2");
	    if(largeR12.size() > 1 && largeR12[1].pt() > 80. ){
	      countCutflowEvent("8SRt2");
	      double temp;
	      double deltaRmin=100000.;
	      int imin;
	      for(int i=0;i < bjetsM.size(); i++){
		temp=bjetsM[i]->P4().DeltaR(missingET->P4());
		if(deltaRmin > temp){
		  deltaRmin = temp;
		  imin=i;
		}
	      }//loop
	      double mTmin;
	      mTmin=mT(bjetsM[imin]->P4(), missingET->P4());

	      if( mTmin > 200 ){
		countCutflowEvent("9SRt2");
		double HT=0.;

		//calculation of scalar sum                                                                      
		for(int i=0; i<jets.size(); i++) {
		  HT += jets[i]->PT;
		}// llop

		if( missingET->P4().Et() / sqrt(HT) > 12. ){ 
		  countCutflowEvent("10SRt2");
		  if( missingET->P4().Et() > 300. ){
		    
		    
		    double deltaR_bb;
		    double deltaRmin_bb=1000.;
		    for(int i=0; i<bjetsM.size(); i++){
		      deltaR_bb=bjetsM[i]->P4().DeltaR(missingET->P4());
		      if( deltaR_bb < deltaRmin_bb )
			deltaRmin_bb=deltaR_bb;
		    } //for loop                                                                           
		    if( deltaRmin_bb > 1.5){
		      countCutflowEvent("11SRt2");
		      countSignalEvent("SRt2");
		    } // deltaRbb cut
		  } //met cut
		}//met,significance
	      }//mTmin cut
		 
		
	    }//large jet 2 PT cut
	  }//large jet 1 PT cut
	}//taujetscuts
      } //deltaphi_min cut
    } //met cut

  } // leptonveto, jets multiplicities
  

  //
  // SRt3
  //

  std::vector<FinalStateObject*> leptons;
  for(int e = 0; e < electronsTight.size(); e++) {
    FinalStateObject* lep = newFinalStateObject(electronsTight[e]);
    leptons.push_back(lep);
  }
  for(int m = 0; m < muonsCombinedPlus.size(); m++) {
    FinalStateObject* lep = newFinalStateObject(muonsCombinedPlus[m]);
    leptons.push_back(lep);
  }



  if( leptons.size() == 2 &&  leptons[0]->Charge * leptons[1]->Charge == -1 && leptons[0]->PT > 25. && leptons[1]->PT > 20. && fabs((leptons[0]->P4() + leptons[1]->P4()).M()) > 20. ){
    countCutflowEvent("1SRt3");

    if( fabs((leptons[0]->P4() + leptons[1]->P4()).M() - 91.) > 20. ){
      countCutflowEvent("2SRt3");
      if(bjetsM.size()>0 && bjetsM[0]->PT>30.){
	countCutflowEvent("3SRt3");
	double mT2_ll = AnalysisBase::mT2(leptons[0]->P4(), leptons[1]->P4(), 0., missingET->P4());
	if(mT2_ll > 100.){
	  countCutflowEvent("4SRt3");

	  
	  double m1=fabs((leptons[0]->P4() + bjetsM[0]->P4() ).M());
	  double m2=fabs((leptons[1]->P4() + bjetsM[0]->P4() ).M());

	  if( ( m1 < m2 ? m1 : m2 ) < 170.){
	    countCutflowEvent("5SRt3");
	    double xiplus=mT2_ll+0.2*missingET->P4().Et();
	    if( xiplus > 170.){
	      countCutflowEvent("6SRt3");
	      double dphi_boost = missingET->P4().DeltaPhi( missingET->P4() + leptons[0]->P4() + leptons[1]->P4() );
	      if(fabs(dphi_boost) < 0.8){
		countCutflowEvent("7SRt3");
		countSignalEvent("SRt3");
	      } //deltaphi_boost cut
	    }//xiplus cut
	  } //m2bl_min cut
	}//stransverse mass cut
      }//one b jet
    }//z veto
  }//2OS lepton
}

void Atlas_1710_11412::finalize() {
  // Whatever should be done after the run goes here
}       
