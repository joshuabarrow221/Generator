//____________________________________________________________________________
/*!

\class    genie::BreitWignerI

\brief    Pure abstract base class. Defines the BreitWignerI interface to
          be implemented by any algorithmic class modeling a Breit Wigner
          function.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 04, 2004

*/
//____________________________________________________________________________


#ifndef _BREIT_WIGNER_I_H_
#define _BREIT_WIGNER_I_H_

#include "Algorithm/Algorithm.h"
#include "BaryonResonance/BaryonResonance.h"

namespace genie {

class BreitWignerI : public Algorithm {

public:
  virtual ~BreitWignerI();

  virtual double Eval(Resonance_t res, double W) const = 0;

protected:
  BreitWignerI();
  BreitWignerI(string name);
  BreitWignerI(string name, string config);
};

}         // genie namespace

#endif    // _BREIT_WIGNER_I_H_
