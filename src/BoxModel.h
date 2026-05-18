/*
  BoxModel.h

  Uniform-gamma variant of the StickySquare model for periodic-box simulations.
  Identical to NucleolusModel except that the coupling factor g is a single
  stored value (uniformGamma) rather than a position-dependent function.
  Set uniformGamma from the run loop to implement annealing or fixed coupling.

  Backbone bonds are NOT scaled by uniformGamma — they are always at full
  strength.  Hard-core repulsion (d < 1) is also never scaled.
*/

#ifndef _BOX_MODEL_H
#define _BOX_MODEL_H

#include "StickySquare.h"

class BoxModel : public StickySquare
{
public:
    double uniformGamma;  //!< Global coupling factor (0 = denatured, 1 = fully coupled).
    bool   singleChain = false;  //!< When true, disables same-chain suppression (single folding chain).

    BoxModel(Box&, std::vector<Particle>&, CellList&,
             unsigned int maxInteractions, double interactionEnergy,
             double interactionRange, Interactions&);

    //! Override: pair energy scaled by uniformGamma for all weak terms.
    double computePairEnergy(unsigned int, const double*, const double*,
                             unsigned int, const double*, const double*);

    //! Override: total particle energy (calls our computePairEnergy).
    double computeEnergy(unsigned int, const double*, const double*);

    //! Override: apply full periodic boundaries (including x) before cell update.
    void applyPostMoveUpdates(unsigned int, const double*, const double*);
};

#endif /* _BOX_MODEL_H */
