/*
  BoxModel.cpp

  Uniform-gamma variant of NucleolusModel for periodic-box simulations.
  The pair energy is identical to NucleolusModel except that all weak-coupling
  terms are multiplied by uniformGamma (a single global value) rather than
  gamma(x_i) * gamma(x_j).  Set BoxModel::uniformGamma from the run loop.

  Backbone bonds (magnitude ~1000) are left unscaled.
  Hard-core overlap (d < 1) is always forbidden (energy = INF).
*/

#include <algorithm>
#include <cmath>
#include "Box.h"
#include "CellList.h"
#include "Particle.h"
#include "BoxModel.h"

extern double INF;
extern double TOL;

BoxModel::BoxModel(
    Box& box_,
    std::vector<Particle>& particles_,
    CellList& cells_,
    unsigned int maxInteractions_,
    double interactionEnergy_,
    double interactionRange_,
    Interactions& interactions_)
    : StickySquare(box_, particles_, cells_,
                   maxInteractions_, interactionEnergy_, interactionRange_,
                   interactions_),
      uniformGamma(1.0)
{
}

double BoxModel::computePairEnergy(
    unsigned int particle1, const double* position1, const double* orientation1,
    unsigned int particle2, const double* position2, const double* orientation2)
{
    std::vector<double> sep(2);
    sep[0] = position1[0] - position2[0];
    sep[1] = position1[1] - position2[1];
    box.minimumImage(sep);

    double normSqd = sep[0]*sep[0] + sep[1]*sep[1];

    if (normSqd < 1.0 - TOL) return INF;

    if (interactions.backboneSpringK > 0.0 && !interactions.backbonePartners.empty()) {
        int n0 = interactions.n0_size > 0 ? interactions.n0_size : 1;
        int id1_ = (int)particle1 % n0;
        int id2_ = (int)particle2 % n0;
        bool isBackbone = false;
        for (int bp : interactions.backbonePartners[id1_])
            if (bp == id2_) { isBackbone = true; break; }
        if (isBackbone)
            return StickySquare::computePairEnergy(particle1, position1, orientation1,
                                                   particle2, position2, orientation2);
    }

    if (normSqd > 5.0 + TOL) return 0.0;

    int n0 = interactions.n0_size > 0 ? interactions.n0_size : 1;
    int id1 = (int)particle1 % n0;
    int id2 = (int)particle2 % n0;

    bool sameChain = ((int)particle1 / n0 == (int)particle2 / n0) &&
                     (id1 / 4 == id2 / 4);

    // Uniform coupling — no position dependence.
    double g = uniformGamma;

    double energy = 0.0;

    if (normSqd < 1.0 + TOL) {
        double backbone = Neighbours::cNone;
        if (mabs(sep[0] - 1.0) < TOL && mabs(sep[1]) < TOL)
            backbone = interactions.east[particle2].getVal(particle1);
        if (backbone == Neighbours::cNone &&
            mabs(sep[0] + 1.0) < TOL && mabs(sep[1]) < TOL)
            backbone = interactions.east[particle1].getVal(particle2);
        if (backbone == Neighbours::cNone &&
            mabs(sep[1] - 1.0) < TOL && mabs(sep[0]) < TOL)
            backbone = interactions.north[particle2].getVal(particle1);
        if (backbone == Neighbours::cNone &&
            mabs(sep[1] + 1.0) < TOL && mabs(sep[0]) < TOL)
            backbone = interactions.north[particle1].getVal(particle2);
        if (backbone != Neighbours::cNone)
            energy += backbone;
        else
            energy += interactions.crosstalk;

        bool isBackbone = (backbone != Neighbours::cNone);
        if (!interactions.weakD1.empty()) {
            bool patchOk = true;
            if (!isBackbone && interactions.patchesEnabled && !interactions.patchSlots.empty()) {
                int dx12 = (int)round(-sep[0]);
                int dy12 = (int)round(-sep[1]);
                int lx1 = (int)round( dx12 * orientation1[0] + dy12 * orientation1[1]);
                int ly1 = (int)round(-dx12 * orientation1[1] + dy12 * orientation1[0]);
                int lx2 = (int)round(-dx12 * orientation2[0] - dy12 * orientation2[1]);
                int ly2 = (int)round( dx12 * orientation2[1] - dy12 * orientation2[0]);
                auto toSlot = [](int lx, int ly) -> int {
                    if (lx== 1&&ly== 0) return 0;
                    if (lx== 0&&ly== 1) return 1;
                    if (lx==-1&&ly== 0) return 2;
                    if (lx== 0&&ly==-1) return 3;
                    return -1;
                };
                int s1 = toSlot(lx1, ly1);
                int s2 = toSlot(lx2, ly2);
                patchOk = (s1 >= 0 && s2 >= 0 &&
                           interactions.patchSlots[id1][s1] &&
                           interactions.patchSlots[id2][s2]);
            }
            if (isBackbone || (!sameChain && patchOk))
                energy += g * interactions.weakD1[id1][id2];
        }

    } else if (normSqd < 2.0 + TOL) {
        double backbone = interactions.east[particle1].getVal(particle2);
        if (backbone == Neighbours::cNone)
            backbone = interactions.east[particle2].getVal(particle1);
        if (backbone != Neighbours::cNone) {
            energy += backbone;
        } else {
            energy += interactions.crosstalk;
            if (!sameChain && !interactions.weakDsq2.empty())
                energy += g * interactions.weakDsq2[id1][id2];
        }

    } else if (normSqd < 4.0 + TOL) {
        if (!sameChain && !interactions.weakD2.empty())
            energy += g * interactions.weakD2[id1][id2];

    } else {
        if (!sameChain && !interactions.weakDsq5.empty())
            energy += g * interactions.weakDsq5[id1][id2];
    }

    return -energy;
}

double BoxModel::computeEnergy(
    unsigned int particle, const double* position, const double* orientation)
{
    double energy = 0.0;
    unsigned int nbrs[maxInteractions];
    unsigned int n = computeInteractions(particle, position, orientation, nbrs);
    for (unsigned int k = 0; k < n; k++) {
        unsigned int nbr = nbrs[k];
        energy += computePairEnergy(particle, position, orientation,
                     nbr, &particles[nbr].position[0], &particles[nbr].orientation[0]);
    }
    return energy;
}

void BoxModel::applyPostMoveUpdates(
    unsigned int particle, const double* position, const double* orientation)
{
    for (unsigned int i = 0; i < box.dimension; i++) {
        particles[particle].position[i]    = position[i];
        particles[particle].orientation[i] = orientation[i];
    }
    // VMMC only wraps y (i>=1) because NucleolusModel has a non-periodic x axis.
    // BoxModel is fully periodic, so wrap x too before computing the cell index.
    box.periodicBoundaries(particles[particle].position);

    unsigned int newCell = cells.getCell(particles[particle]);
    if (particles[particle].cell != newCell)
        cells.updateCell(newCell, particles[particle], particles);
}
