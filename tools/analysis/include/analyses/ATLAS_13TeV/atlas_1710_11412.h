#ifndef ATLAS_1710_11412_H_
#define ATLAS_1710_11412_H_
// AUTHOR: J.S.Kim
//  EMAIL: jsk@th.physik.uni-bonn.de
#include "AnalysisBase.h"

class Atlas_1710_11412 : public AnalysisBase {
  public:
    Atlas_1710_11412() : AnalysisBase()  {}               
    ~Atlas_1710_11412() {}
  
    void initialize();
    void analyze();        
    void finalize();

  private:
};

#endif
