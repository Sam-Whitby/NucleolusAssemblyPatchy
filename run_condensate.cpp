/*
  run_condensate.cpp

  Circular condensate assembly simulation.

  Models assembly of nCopies copies of a Moore-curve target complex inside a
  2D circular condensate of radius R_c.

  Three simulation phases (all measured in outer iterations):
    --t-equil N   Equilibration: g=1 everywhere; complexes are stable.
    --t-denat N   Denaturation: g=0 everywhere; complexes fall apart.
    --steps   N   Main run: radial gradient γ(r)=γ0+(1-γ0)·r/R_c.

  Recycling (multi-chain):  exited polymer groups are placed on the 4 cardinal
    ring sites around the condensate centre.
  Recycling (single-chain): exited chains are threaded one particle at a time
    through the condensate centre.  The centre lattice point has a hard wall;
    no ring repulsion.  Energy is reported excluding backbone bonds (|e|<500)
    to avoid jumps as backbone bonds enter/leave interaction range.

  Options
    --steps       N        main-phase outer iterations                    [10000]
    --snapshots   N        total snapshots across all phases              [1000]
    --t-equil     N        equilibration iterations                       [0]
    --t-denat     N        denaturation iterations                        [0]
    --copies      N        number of target complex copies                [4]
    --radius      R        condensate radius in lattice units             [60]
    --gamma0      γ        minimum coupling at r=0 (gradient floor)      [0.0]
    --gradient             enable radial gradient (default: g=1 everywhere)
    --stokes               enable Stokes hydrodynamic drag (D ∝ 1/R)
    --coupling    MODE     gradient coupling: product or midpoint         [product]
    --phi-rot     φ        fraction of cluster rotation moves             [0.2]
    --phi-reorient φ       fraction of in-place reorientation moves       [0.2]
    --output      PREFIX   prefix for output files                        [condensate]
    --seed        S        RNG seed (0 = time-based)                     [1]
    --J           J        weak patch coupling strength (gradient-scaled) [8.0]
    --single               build a single self-folding Moore-curve chain
    --Moore-N     N        Moore curve order (4^N particles, N>=1)        [2]
*/

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <queue>
#include <algorithm>
#include <array>

#include "Demo.h"
#include "VMMC.h"
#include "CondensateModel.h"

using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

extern double INF;
extern double TOL;

// ============================================================
//  Multi-chain target (n=2 Moore curve, 4 polymers × 4 segments)
// ============================================================
static const int N0_MULTI   = 16;
static const int N_POLYMER  = 4;
static const int N_SEG      = 4;

static const int TARGET_X_MULTI[N0_MULTI] = { 0,0,1,2,  3,3,2,2,  0,0,1,1,  3,3,2,1 };
static const int TARGET_Y_MULTI[N0_MULTI] = { 2,3,3,3,  3,2,2,1,  0,1,1,2,  1,0,0,0 };

static const int BACKBONE_PAIRS[][2] = {
    {0,1},{1,2},{2,3},
    {4,5},{5,6},{6,7},
    {8,9},{9,10},{10,11},
    {12,13},{13,14},{14,15}
};
static const int N_BB_PAIRS = 12;

inline int polyType(int id) { return id / N_SEG; }

inline double targetDistSqd_multi(int id1, int id2) {
    double dx = TARGET_X_MULTI[id1] - TARGET_X_MULTI[id2];
    double dy = TARGET_Y_MULTI[id1] - TARGET_Y_MULTI[id2];
    return dx*dx + dy*dy;
}

// ============================================================
//  Moore curve L-system generator (same as run_box.cpp)
// ============================================================
static vector<pair<int,int>> generateMooreCurve(int N)
{
    string s = "LFL+F+LFL";
    for (int iter = 0; iter < N - 1; iter++) {
        string t;
        t.reserve(s.size() * 6);
        for (char c : s) {
            switch (c) {
                case 'L': t += "-RF+LFL+FR-"; break;
                case 'R': t += "+LF-RFR-FL+"; break;
                default:  t += c;             break;
            }
        }
        s = std::move(t);
    }
    int x = 0, y = 0, dx = 1, dy = 0;
    vector<pair<int,int>> path;
    path.reserve(1 << (2 * N));
    path.push_back({x, y});
    for (char c : s) {
        if (c == 'F') {
            x += dx; y += dy;
            path.push_back({x, y});
        } else if (c == '+') {
            int tmp = dx; dx = -dy; dy = tmp;
        } else if (c == '-') {
            int tmp = dx; dx = dy; dy = -tmp;
        }
    }
    int minX = path[0].first, minY = path[0].second;
    for (auto& p : path) { minX = min(minX, p.first); minY = min(minY, p.second); }
    for (auto& p : path) { p.first -= minX; p.second -= minY; }
    return path;
}

// ============================================================
//  Multi-chain build functions (unchanged from before)
// ============================================================
static void buildCouplingMatrices_multi(
    double J,
    vector<vector<double>>& wD1,
    vector<vector<double>>& wDsq2,
    vector<vector<double>>& wD2,
    vector<vector<double>>& wDsq5)
{
    wD1.assign(N0_MULTI, vector<double>(N0_MULTI, 0.0));
    wDsq2.assign(N0_MULTI, vector<double>(N0_MULTI, 0.0));
    wD2.assign(N0_MULTI, vector<double>(N0_MULTI, 0.0));
    wDsq5.assign(N0_MULTI, vector<double>(N0_MULTI, 0.0));
    for (int i = 0; i < N0_MULTI; i++)
        for (int j = 0; j < N0_MULTI; j++) {
            if (i == j || polyType(i) == polyType(j)) continue;
            if (targetDistSqd_multi(i, j) < 1.0 + TOL)
                wD1[i][j] = J;
        }
}

static vector<array<bool,4>> buildPatchSlots_multi()
{
    vector<array<bool,4>> slots(N0_MULTI);
    for (auto& s : slots) s.fill(false);
    for (int i = 0; i < N0_MULTI; i++)
        for (int j = 0; j < N0_MULTI; j++) {
            if (i == j || polyType(i) == polyType(j)) continue;
            if (targetDistSqd_multi(i, j) > 1.0 + 1e-6) continue;
            int dx = TARGET_X_MULTI[j] - TARGET_X_MULTI[i];
            int dy = TARGET_Y_MULTI[j] - TARGET_Y_MULTI[i];
            int slot = -1;
            if      (dx ==  1 && dy ==  0) slot = 0;
            else if (dx ==  0 && dy ==  1) slot = 1;
            else if (dx == -1 && dy ==  0) slot = 2;
            else if (dx ==  0 && dy == -1) slot = 3;
            if (slot >= 0) slots[i][slot] = true;
        }
    return slots;
}

static void buildBackboneTriples_multi(int nCopies, double bbEnergy,
                                        vector<Triple>& north, vector<Triple>& east)
{
    north.clear(); east.clear();
    for (int c = 0; c < nCopies; c++) {
        int base = c * N0_MULTI;
        for (int k = 0; k < N_BB_PAIRS; k++) {
            int gi = base + BACKBONE_PAIRS[k][0];
            int gj = base + BACKBONE_PAIRS[k][1];
            east.push_back({gi, gj, bbEnergy}); east.push_back({gj, gi, bbEnergy});
            north.push_back({gi, gj, bbEnergy}); north.push_back({gj, gi, bbEnergy});
        }
    }
}

static double referenceComplexEnergy_multi(double J, double bbEnergy)
{
    double E = N_BB_PAIRS * -bbEnergy;
    for (int i = 0; i < N0_MULTI; i++)
        for (int j = i+1; j < N0_MULTI; j++) {
            if (polyType(i) == polyType(j)) continue;
            if (targetDistSqd_multi(i, j) < 1.0 + 1e-6) E -= J;
        }
    return E;
}

// Returns energy of one perfect complex counting ONLY patch contacts (no backbone).
static double referenceComplexEnergyNoBB_multi(double J)
{
    double E = 0.0;
    for (int i = 0; i < N0_MULTI; i++)
        for (int j = i+1; j < N0_MULTI; j++) {
            if (polyType(i) == polyType(j)) continue;
            if (targetDistSqd_multi(i, j) < 1.0 + 1e-6) E -= J;
        }
    return E;
}

// ============================================================
//  Single-chain build functions (same as run_box.cpp)
// ============================================================
static void buildCouplingMatrices_single(
    int n0, const vector<pair<int,int>>& path, double J,
    vector<vector<double>>& wD1,
    vector<vector<double>>& wDsq2,
    vector<vector<double>>& wD2,
    vector<vector<double>>& wDsq5)
{
    wD1.assign(n0, vector<double>(n0, 0.0));
    wDsq2.assign(n0, vector<double>(n0, 0.0));
    wD2.assign(n0, vector<double>(n0, 0.0));
    wDsq5.assign(n0, vector<double>(n0, 0.0));
    for (int i = 0; i < n0; i++)
        for (int j = 0; j < n0; j++) {
            if (i == j || abs(i - j) == 1) continue;
            int dx = path[i].first  - path[j].first;
            int dy = path[i].second - path[j].second;
            if (dx*dx + dy*dy < 1.0 + 1e-6) wD1[i][j] = J;
        }
}

static vector<array<bool,4>> buildPatchSlots_single(
    int n0, const vector<pair<int,int>>& path)
{
    vector<array<bool,4>> slots(n0);
    for (auto& s : slots) s.fill(false);
    for (int i = 0; i < n0; i++)
        for (int j = 0; j < n0; j++) {
            if (i == j || abs(i - j) == 1) continue;
            int dx = path[j].first  - path[i].first;
            int dy = path[j].second - path[i].second;
            if (dx*dx + dy*dy > 1.0 + 1e-6) continue;
            int slot = -1;
            if      (dx ==  1 && dy ==  0) slot = 0;
            else if (dx ==  0 && dy ==  1) slot = 1;
            else if (dx == -1 && dy ==  0) slot = 2;
            else if (dx ==  0 && dy == -1) slot = 3;
            if (slot >= 0) slots[i][slot] = true;
        }
    return slots;
}

static void buildBackboneTriples_single(int n0, int nCopies, double bbEnergy,
                                         vector<Triple>& north, vector<Triple>& east)
{
    north.clear(); east.clear();
    for (int c = 0; c < nCopies; c++) {
        int base = c * n0;
        for (int k = 0; k < n0 - 1; k++) {
            int gi = base + k, gj = base + k + 1;
            east.push_back({gi, gj, bbEnergy}); east.push_back({gj, gi, bbEnergy});
            north.push_back({gi, gj, bbEnergy}); north.push_back({gj, gi, bbEnergy});
        }
    }
}

static double referenceComplexEnergy_single(
    int n0, const vector<pair<int,int>>& path, double J, double bbEnergy)
{
    double E = (n0 - 1) * -bbEnergy;
    for (int i = 0; i < n0; i++)
        for (int j = i + 2; j < n0; j++) {
            int dx = path[i].first - path[j].first;
            int dy = path[i].second - path[j].second;
            if (dx*dx + dy*dy < 1.0 + 1e-6) E -= J;
        }
    return E;
}

static double referenceComplexEnergyNoBB_single(
    int n0, const vector<pair<int,int>>& path, double J)
{
    double E = 0.0;
    for (int i = 0; i < n0; i++)
        for (int j = i + 2; j < n0; j++) {
            int dx = path[i].first - path[j].first;
            int dy = path[i].second - path[j].second;
            if (dx*dx + dy*dy < 1.0 + 1e-6) E -= J;
        }
    return E;
}

// ============================================================
//  Place particles as assembled target complexes at start
// ============================================================
static void placeParticlesTargetComplex(
    vector<Particle>& particles, CellList& cells, Box& box,
    int nCopies, int n0, const vector<pair<int,int>>& path,
    int spacing, double cx, double cy)
{
    cells.reset();
    int nParticles = nCopies * n0;
    for (int i = 0; i < nParticles; i++) {
        particles[i].index = i;
        particles[i].position.resize(2);
        particles[i].orientation.resize(2);
        particles[i].orientation[0] = 1.0;
        particles[i].orientation[1] = 0.0;
    }
    int nc_grid = (int)ceil(sqrt((double)nCopies));
    double ox0 = round(cx) - ((nc_grid - 1) * spacing) / 2.0;
    double oy0 = round(cy) - ((nc_grid - 1) * spacing) / 2.0;
    for (int c = 0; c < nCopies; c++) {
        int col = c % nc_grid, row = c / nc_grid;
        double ox = ox0 + col * spacing;
        double oy = oy0 + row * spacing;
        for (int lid = 0; lid < n0; lid++) {
            int gi = c * n0 + lid;
            particles[gi].position[0] = ox + path[lid].first;
            particles[gi].position[1] = oy + path[lid].second;
            box.periodicBoundaries(particles[gi].position);
            particles[gi].cell = cells.getCell(particles[gi]);
            cells.initCell(particles[gi].cell, particles[gi]);
        }
    }
}

// ============================================================
//  componentPairEnergy — sum all intra-component pair energies
// ============================================================
static double componentPairEnergy(CondensateModel& model,
                                   const vector<int>& comp,
                                   const vector<Particle>& particles)
{
    double E = 0.0;
    for (int ii = 0; ii < (int)comp.size(); ii++)
        for (int jj = ii+1; jj < (int)comp.size(); jj++) {
            int i = comp[ii], j = comp[jj];
            double e = model.computePairEnergy(
                i, &particles[i].position[0], &particles[i].orientation[0],
                j, &particles[j].position[0], &particles[j].orientation[0]);
            if (e < 1e5) E += e;
        }
    return E;
}

// ============================================================
//  Count copies with intra-copy energy ≈ refComplexEnergy
// ============================================================
static int countAssembled(CondensateModel& model, const vector<Particle>& particles,
                           int nCopies, int n0, double refComplexEnergy)
{
    int count = 0;
    for (int c = 0; c < nCopies; c++) {
        int base = c * n0;
        double E = 0.0;
        for (int i = base; i < base + n0; i++)
            for (int j = i+1; j < base + n0; j++) {
                double e = model.computePairEnergy(
                    i, &particles[i].position[0], &particles[i].orientation[0],
                    j, &particles[j].position[0], &particles[j].orientation[0]);
                if (e < 1e5) E += e;
            }
        if (fabs(E - refComplexEnergy) < 0.5) count++;
    }
    return count;
}

// ============================================================
//  handleExitsAndQueue
//  Detect components fully outside R_c, stage them, enqueue for recycling.
//  For multi-chain: groups are 4-particle polymers, pushed to injQueue
//                   for ring placement by tryInjectNext.
//  For single-chain: the whole chain is one group, pushed to injQueue
//                    for one-at-a-time centre threading by tryThreadNext.
// ============================================================
struct ExitCounts { int particles; int perfect; int fullComplex; };

static ExitCounts handleExitsAndQueue(
    CondensateModel& model,
    vector<Particle>& particles,
    int nParticles, int n0,
    CellList& cells, Box& box,
    double cx, double cy, double R_c,
    vmmc::VMMC& vmmc,
    queue<vector<int>>& injQueue,
    set<int>& staged,
    int stageX0,
    double refComplexEnergy,
    bool singleChain)
{
    const double Rc2 = R_c * R_c;
    vector<int> outsideSeeds;
    for (int i = 0; i < nParticles; i++) {
        if (staged.count(i)) continue;
        double dx = particles[i].position[0] - cx;
        double dy = particles[i].position[1] - cy;
        if (dx*dx + dy*dy > Rc2) outsideSeeds.push_back(i);
    }
    if (outsideSeeds.empty()) return {0, 0, 0};

    const int maxInt = 30;
    unsigned int nbrs[maxInt];
    vector<int> visited(nParticles, -1);
    vector<vector<int>> components;

    for (int seed : outsideSeeds) {
        if (visited[seed] != -1) continue;
        int compIdx = (int)components.size();
        vector<int> comp = {seed};
        visited[seed] = compIdx;
        for (int ci = 0; ci < (int)comp.size(); ci++) {
            int j = comp[ci];
            int nn = (int)model.computeInteractions(
                j, &particles[j].position[0], &particles[j].orientation[0], nbrs);
            for (int k = 0; k < nn; k++) {
                int nbr = (int)nbrs[k];
                if (visited[nbr] != -1 || staged.count(nbr)) continue;
                double e = model.computePairEnergy(
                    j,   &particles[j].position[0],   &particles[j].orientation[0],
                    nbr, &particles[nbr].position[0], &particles[nbr].orientation[0]);
                if (e != 0.0 && e < 1e5) {
                    visited[nbr] = compIdx;
                    comp.push_back(nbr);
                }
            }
        }
        components.push_back(comp);
    }

    int particlesExited = 0, perfectExited = 0, fullComplexExited = 0;
    int icy = (int)round(cy);

    for (auto& comp : components) {
        // Must be fully outside R_c
        bool allOut = true;
        for (int gi : comp) {
            double dx = particles[gi].position[0] - cx;
            double dy = particles[gi].position[1] - cy;
            if (dx*dx + dy*dy <= Rc2) { allOut = false; break; }
        }
        if (!allOut) continue;

        // Must be isolated (no non-staged bonded neighbours outside the component)
        set<int> compSet(comp.begin(), comp.end());
        bool isolated = true;
        for (int gi : comp) {
            int nn = (int)model.computeInteractions(
                gi, &particles[gi].position[0], &particles[gi].orientation[0], nbrs);
            for (int k = 0; k < nn; k++) {
                int nbr = (int)nbrs[k];
                if (compSet.count(nbr) || staged.count(nbr)) continue;
                double e = model.computePairEnergy(
                    gi,  &particles[gi].position[0],  &particles[gi].orientation[0],
                    nbr, &particles[nbr].position[0], &particles[nbr].orientation[0]);
                if (e != 0.0 && e < 1e5) { isolated = false; break; }
            }
            if (!isolated) break;
        }
        if (!isolated) continue;

        particlesExited += (int)comp.size();

        // Check for full/perfect complex exit
        if ((int)comp.size() == n0) {
            set<int> localIds;
            for (int gi : comp) localIds.insert(gi % n0);
            if ((int)localIds.size() == n0) {
                fullComplexExited++;
                double E = componentPairEnergy(model, comp, particles);
                if (fabs(E - refComplexEnergy) < 0.5) perfectExited++;
            }
        }

        // Group and stage
        if (singleChain) {
            // One group per copy: all n0 particles, sorted by local chain index
            map<int, vector<int>> groups;
            for (int gi : comp) groups[gi / n0].push_back(gi);

            for (auto& kv : groups) {
                auto& gids = kv.second;
                sort(gids.begin(), gids.end(),
                     [n0](int a, int b){ return (a % n0) < (b % n0); });
                int stageY = icy + kv.first * (n0 + 2);
                for (int s = 0; s < (int)gids.size(); s++) {
                    int gi = gids[s];
                    particles[gi].position[0] = stageX0 + s;
                    particles[gi].position[1] = stageY;
                    box.periodicBoundaries(particles[gi].position);
                    int newCell = cells.getCell(particles[gi]);
                    cells.updateCell(newCell, particles[gi], particles);
                    vmmc.syncPosition(gi, &particles[gi].position[0]);
                    staged.insert(gi);
                }
                injQueue.push(gids);
            }
        } else {
            // Multi-chain: group by polymer (copy × N_POLYMER + polymerWithinCopy)
            map<int, vector<int>> groups;
            for (int gi : comp) {
                int key = (gi / n0) * N_POLYMER + (gi % n0) / N_SEG;
                groups[key].push_back(gi);
            }
            for (auto& kv : groups) {
                auto& gids = kv.second;
                sort(gids.begin(), gids.end(),
                     [](int a, int b){ return (a % N_SEG) < (b % N_SEG); });
                int stageY = icy + kv.first * 2;
                for (int s = 0; s < (int)gids.size(); s++) {
                    int gi = gids[s];
                    particles[gi].position[0] = stageX0 + s;
                    particles[gi].position[1] = stageY;
                    box.periodicBoundaries(particles[gi].position);
                    int newCell = cells.getCell(particles[gi]);
                    cells.updateCell(newCell, particles[gi], particles);
                    vmmc.syncPosition(gi, &particles[gi].position[0]);
                    staged.insert(gi);
                }
                injQueue.push(gids);
            }
        }
    }

    return {particlesExited, perfectExited, fullComplexExited};
}

// ============================================================
//  tryInjectNext — multi-chain ring injection (unchanged logic)
//  Waits for all 4 ring sites to be clear, then places a 4-particle
//  polymer group at the ring positions CW or CCW.
// ============================================================
static void tryInjectNext(
    vector<Particle>& particles, int nParticles,
    CellList& cells, Box& box, vmmc::VMMC& vmmc,
    double cx, double cy,
    queue<vector<int>>& injQueue,
    set<int>& staged)
{
    if (injQueue.empty()) return;
    int icx = (int)round(cx), icy = (int)round(cy);
    pair<int,int> ring[4] = {{icx,icy+1},{icx+1,icy},{icx,icy-1},{icx-1,icy}};
    for (int ri = 0; ri < 4; ri++)
        for (int gi = 0; gi < nParticles; gi++) {
            if (staged.count(gi)) continue;
            double dx = particles[gi].position[0] - ring[ri].first;
            double dy = particles[gi].position[1] - ring[ri].second;
            if (fabs(dx) < 0.5 && fabs(dy) < 0.5) return;
        }
    vector<int>& gids = injQueue.front();
    sort(gids.begin(), gids.end(),
         [](int a, int b){ return (a % N_SEG) < (b % N_SEG); });
    int  start = (int)(vmmc.rng() * 4) % 4;
    bool cw    = vmmc.rng() < 0.5;
    for (int s = 0; s < (int)gids.size(); s++) {
        int ringIdx = cw ? (start + s) % 4 : (start - s + 4) % 4;
        int gi = gids[s];
        particles[gi].position[0] = ring[ringIdx].first;
        particles[gi].position[1] = ring[ringIdx].second;
        box.periodicBoundaries(particles[gi].position);
        int newCell = cells.getCell(particles[gi]);
        cells.updateCell(newCell, particles[gi], particles);
        vmmc.syncPosition(gi, &particles[gi].position[0]);
        staged.erase(gi);
    }
    injQueue.pop();
}

// ============================================================
//  tryThreadNext — single-chain threading
//  First particle is placed at the condensate centre.
//  Each subsequent particle is placed adjacent (cardinal) to the previous one
//  so backbone bonds are never broken during injection.
//  Starts a new chain from injQueue when threadingChain is empty.
// ============================================================
static void tryThreadNext(
    vector<Particle>& particles, int nParticles,
    CellList& cells, Box& box, vmmc::VMMC& vmmc,
    int icx, int icy,
    set<int>& staged,
    vector<int>& threadingChain, int& threadingIdx,
    queue<vector<int>>& injQueue)
{
    // Start a new chain if none is active
    if (threadingChain.empty()) {
        if (injQueue.empty()) return;
        threadingChain = injQueue.front();
        injQueue.pop();
        threadingIdx = 0;
    }
    if (threadingIdx >= (int)threadingChain.size()) {
        threadingChain.clear();
        return;
    }

    int tx, ty;  // target lattice position for the new particle

    if (threadingIdx == 0) {
        // First particle: enter at the condensate centre
        for (int gi = 0; gi < nParticles; gi++) {
            if (staged.count(gi)) continue;
            double dx = particles[gi].position[0] - icx;
            double dy = particles[gi].position[1] - icy;
            if (fabs(dx) < 0.5 && fabs(dy) < 0.5) return;  // centre occupied
        }
        tx = icx; ty = icy;
    } else {
        // Subsequent particles: place adjacent to the previous particle so the
        // backbone bond (distance 1) is immediately satisfied.
        int prev = threadingChain[threadingIdx - 1];
        int px = (int)round(particles[prev].position[0]);
        int py = (int)round(particles[prev].position[1]);

        const int ddx[4] = { 1, -1,  0,  0 };
        const int ddy[4] = { 0,  0,  1, -1 };
        int chosen = -1;
        for (int d = 0; d < 4; d++) {
            int cx_ = px + ddx[d];
            int cy_ = py + ddy[d];
            bool free = true;
            for (int gi = 0; gi < nParticles; gi++) {
                if (staged.count(gi)) continue;
                double ex = particles[gi].position[0] - cx_;
                double ey = particles[gi].position[1] - cy_;
                if (fabs(ex) < 0.5 && fabs(ey) < 0.5) { free = false; break; }
            }
            if (free) { chosen = d; break; }
        }
        if (chosen < 0) return;  // all cardinal neighbours occupied; wait
        tx = px + ddx[chosen];
        ty = py + ddy[chosen];
    }

    // Place particle at (tx, ty)
    int gi = threadingChain[threadingIdx++];
    particles[gi].position[0] = tx;
    particles[gi].position[1] = ty;
    int newCell = cells.getCell(particles[gi]);
    cells.updateCell(newCell, particles[gi], particles);
    vmmc.syncPosition(gi, &particles[gi].position[0]);
    staged.erase(gi);

    if (threadingIdx >= (int)threadingChain.size())
        threadingChain.clear();
}

// ============================================================
//  writeFrame
// ============================================================
static void writeFrame(FILE* fp, const vector<Particle>& particles,
                        int nCopies, int n0, bool singleMode,
                        double R_c, double cx, double cy,
                        long long step, double energy,
                        long long exitedParticles, long long exitedPerfect,
                        long long exitedFull, double exitQuality,
                        const string& couplingLabel, double gamma0,
                        const string& phase, double refEnergy, int nAssembled)
{
    int nParticles = (int)particles.size();
    fprintf(fp, "%d\n", nParticles);
    fprintf(fp,
        "step=%lld energy=%.6f exitedParticles=%lld exitedPerfect=%lld"
        " exitedFull=%lld exitQuality=%.4f"
        " R_c=%.1f cx=%.1f cy=%.1f nCopies=%d coupling=%s gamma0=%.4f phase=%s"
        " refEnergy=%.4f nAssembled=%d single=%d n0=%d\n",
        step, energy, exitedParticles, exitedPerfect,
        exitedFull, exitQuality,
        R_c, cx, cy, nCopies, couplingLabel.c_str(), gamma0, phase.c_str(),
        refEnergy, nAssembled, (int)singleMode, n0);
    for (int i = 0; i < nParticles; i++) {
        int copy  = i / n0;
        int lid   = i % n0;
        int ptype = singleMode ? 0 : polyType(lid);
        fprintf(fp, "%d %d %.4f %.4f %d %.4f %.4f\n",
                i, ptype,
                particles[i].position[0],
                particles[i].position[1],
                copy,
                particles[i].orientation[0],
                particles[i].orientation[1]);
    }
}

// ============================================================
//  main()
// ============================================================
int main(int argc, char** argv)
{
    long long nsteps      = 10000;
    long long nsnaps      = 1000;
    long long t_equil     = 0;
    long long t_denat     = 0;
    int       nCopies     = 4;
    double    R_c         = 60.0;
    double    gamma0      = 0.0;
    bool      useGradient = false;
    bool      useStokes   = false;
    string    couplingStr = "product";
    double    phi_rot     = 0.2;
    double    phi_reorient = 0.2;
    string    outPrefix   = "condensate";
    unsigned int seed     = 1;
    double    J           = 8.0;
    bool      singleChain = false;
    int       mooreN      = 2;

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i],"--steps")      && i+1<argc) { nsteps    = atoll(argv[++i]); }
        else if (!strcmp(argv[i],"--snapshots")  && i+1<argc) { nsnaps    = atoll(argv[++i]); }
        else if (!strcmp(argv[i],"--t-equil")    && i+1<argc) { t_equil   = atoll(argv[++i]); }
        else if (!strcmp(argv[i],"--t-denat")    && i+1<argc) { t_denat   = atoll(argv[++i]); }
        else if (!strcmp(argv[i],"--copies")     && i+1<argc) { nCopies   = atoi(argv[++i]); }
        else if (!strcmp(argv[i],"--radius")     && i+1<argc) { R_c       = atof(argv[++i]); }
        else if (!strcmp(argv[i],"--gamma0")     && i+1<argc) { gamma0    = atof(argv[++i]); }
        else if (!strcmp(argv[i],"--gradient"))               { useGradient = true; }
        else if (!strcmp(argv[i],"--stokes"))                 { useStokes   = true; }
        else if (!strcmp(argv[i],"--coupling")   && i+1<argc) { couplingStr = argv[++i]; }
        else if (!strcmp(argv[i],"--phi-rot")      && i+1<argc) { phi_rot     = atof(argv[++i]); }
        else if (!strcmp(argv[i],"--phi-reorient") && i+1<argc) { phi_reorient= atof(argv[++i]); }
        else if (!strcmp(argv[i],"--output")     && i+1<argc) { outPrefix  = argv[++i]; }
        else if (!strcmp(argv[i],"--seed")       && i+1<argc) { seed = (unsigned int)atoi(argv[++i]); }
        else if (!strcmp(argv[i],"--J")          && i+1<argc) { J         = atof(argv[++i]); }
        else if (!strcmp(argv[i],"--single"))                 { singleChain = true; }
        else if (!strcmp(argv[i],"--Moore-N")    && i+1<argc) { mooreN    = atoi(argv[++i]); }
        else { cerr << "Unknown argument: " << argv[i] << "\n"; }
    }
    if (mooreN < 1) { cerr << "--Moore-N must be >= 1\n"; return 1; }

    CouplingMode couplingMode = CouplingMode::Product;
    if (couplingStr == "midpoint") couplingMode = CouplingMode::Midpoint;
    else if (couplingStr != "product") {
        cerr << "[WARNING] Unknown --coupling '" << couplingStr << "'; defaulting to product.\n";
        couplingStr = "product";
    }

    // ============================================================
    //  Compute n0 and path
    // ============================================================
    int n0;
    vector<pair<int,int>> path;
    if (singleChain) {
        path = generateMooreCurve(mooreN);
        n0   = (int)path.size();
        int expected = 1 << (2 * mooreN);
        if (n0 != expected) {
            cerr << "Moore curve error: got " << n0 << " points, expected " << expected << "\n";
            return 1;
        }
    } else {
        n0 = N0_MULTI;
        for (int i = 0; i < N0_MULTI; i++)
            path.push_back({TARGET_X_MULTI[i], TARGET_Y_MULTI[i]});
        if (mooreN != 2)
            cerr << "Warning: --Moore-N ignored in multi-chain mode.\n";
    }

    // ============================================================
    //  Build interaction data
    // ============================================================
    const double bbEnergy = 1000.0;
    vector<vector<double>> wD1, wDsq2, wD2, wDsq5;
    vector<Triple> north0, east0;
    vector<array<bool,4>> patchSlots;
    double refComplexEnergy_full; // for exit detection
    double refComplexEnergy_nobb; // for energy reporting (no backbone)

    if (singleChain) {
        buildCouplingMatrices_single(n0, path, J, wD1, wDsq2, wD2, wDsq5);
        patchSlots            = buildPatchSlots_single(n0, path);
        buildBackboneTriples_single(n0, nCopies, bbEnergy, north0, east0);
        refComplexEnergy_full = referenceComplexEnergy_single(n0, path, J, bbEnergy);
        refComplexEnergy_nobb = referenceComplexEnergyNoBB_single(n0, path, J);
    } else {
        buildCouplingMatrices_multi(J, wD1, wDsq2, wD2, wDsq5);
        patchSlots            = buildPatchSlots_multi();
        buildBackboneTriples_multi(nCopies, bbEnergy, north0, east0);
        refComplexEnergy_full = referenceComplexEnergy_multi(J, bbEnergy);
        refComplexEnergy_nobb = referenceComplexEnergyNoBB_multi(J);
    }

    int gridSize = singleChain ? (1 << mooreN) : 4;
    int spacing  = gridSize + 3;

    int nFoldingContacts = 0;
    for (int i = 0; i < n0; i++)
        for (int j = i+1; j < n0; j++)
            if (wD1[i][j] != 0.0) nFoldingContacts++;

    const int nParticles = nCopies * n0;
    const double R_large = max(6.0 * R_c, 150.0);
    const double cx = R_large, cy = R_large;
    const double BOX = 2.0 * R_large;
    const int    icx = (int)round(cx), icy = (int)round(cy);

    const unsigned int dimension       = 2;
    const double interactionRange      = 2.5;
    const unsigned int maxInteractions = 30;
    const double interactionEnergy     = 0.0;
    const bool   isLattice             = true;

    // Snapshot distribution
    long long totalSteps  = t_equil + t_denat + nsteps;
    long long snaps_equil = (totalSteps > 0 && t_equil > 0)
                            ? max(1LL, nsnaps * t_equil / totalSteps) : 0;
    long long snaps_denat = (totalSteps > 0 && t_denat > 0)
                            ? max(1LL, nsnaps * t_denat / totalSteps) : 0;
    long long snaps_main  = nsnaps - snaps_equil - snaps_denat;
    if (snaps_main < 1) snaps_main = 1;
    auto saveEvery = [](long long steps, long long snaps) -> long long {
        return (steps <= 0 || snaps <= 0) ? steps : max(1LL, steps / snaps);
    };
    long long saveEvery_equil = saveEvery(t_equil, snaps_equil);
    long long saveEvery_denat = saveEvery(t_denat, snaps_denat);
    long long saveEvery_main  = saveEvery(nsteps,  snaps_main);

    cout << "=== Circular Condensate Assembly Simulation ===" << endl;
    cout << "  mode=" << (singleChain ? "single-chain" : "multi-chain")
         << "  N0=" << n0 << "  mooreN=" << mooreN
         << "  foldingContacts=" << nFoldingContacts << endl;
    cout << "  copies=" << nCopies
         << "  t_equil=" << t_equil << " t_denat=" << t_denat << " steps=" << nsteps
         << "  R_c=" << R_c << "  gamma0=" << gamma0
         << "  gradient=" << useGradient << " stokes=" << useStokes
         << "  coupling=" << couplingStr << "  J=" << J
         << "  phi_rot=" << phi_rot << "  phi_reorient=" << phi_reorient << endl;

    vector<vector<int>> bbPartners(n0);
    double springK = 0.0;
    Interactions interactions(nParticles, n0, north0, east0,
                               wD1, wDsq2, wD2, wDsq5, springK, bbPartners);
    interactions.patchesEnabled = true;
    interactions.patchSlots     = patchSlots;

    vector<double> boxSize    = { BOX, BOX };
    vector<bool>   isPeriodic = { false, false };
    Box box(boxSize, isPeriodic);
    box.isLattice = true;

    CellList cells;
    cells.setDimension(dimension);
    cells.initialise(box.boxSize, interactionRange);

    vector<Particle> particles(nParticles);

    CondensateModel model(box, particles, cells,
                          maxInteractions, interactionEnergy, interactionRange,
                          interactions, cx, cy, R_c, useGradient, couplingMode, gamma0);
    model.singleChain  = singleChain;

    cout << "Reference complex energy (full, g=1):    " << refComplexEnergy_full << endl;
    cout << "Reference complex energy (no-BB, g=1):   " << refComplexEnergy_nobb << endl;

    placeParticlesTargetComplex(particles, cells, box, nCopies, n0, path, spacing, cx, cy);

    // VMMC setup
    vector<double> coordinates(dimension * nParticles);
    vector<double> orientations(dimension * nParticles);
    unique_ptr<bool[]> isIsotropic(new bool[nParticles]());  // false = patchy

    for (int i = 0; i < nParticles; i++) {
        coordinates[2*i]   = particles[i].position[0];
        coordinates[2*i+1] = particles[i].position[1];
        orientations[2*i]   = 1.0;
        orientations[2*i+1] = 0.0;
    }

    double maxTrialTranslation = 1.5;
    double maxTrialRotation    = (phi_rot > 0.0 || phi_reorient > 0.0) ? M_PI : 0.0;
    double probTranslate       = 1.0 - phi_rot - phi_reorient;
    double referenceRadius     = 0.5;
    bool   isRepulsive         = true;
    int    nLatticeNeighbours  = 8;

    // Threading / injection state
    queue<vector<int>> injQueue;
    set<int>           staged;
    int stageX0 = (int)round(cx + BOX * 0.3);
    vector<int> threadingChain;
    int         threadingIdx = 0;

    using namespace std::placeholders;
    vmmc::CallbackFunctions callbacks;
    callbacks.energyCallback =
        std::bind(&CondensateModel::computeEnergy, &model, _1, _2, _3);
    callbacks.pairEnergyCallback =
        std::bind(&CondensateModel::computePairEnergy, &model, _1, _2, _3, _4, _5, _6);
    callbacks.interactionsCallback =
        std::bind(&CondensateModel::computeInteractions, &model, _1, _2, _3, _4);
    callbacks.postMoveCallback =
        std::bind(&CondensateModel::applyPostMoveUpdates, &model, _1, _2, _3);

    // Hard wall: centre lattice point (d²<0.5) is always forbidden.
    // Staged particles are frozen.
    callbacks.boundaryCallback =
        [cx, cy, &staged](unsigned int pid, const double* pos, const double*) -> bool {
            double dx = pos[0] - cx, dy = pos[1] - cy;
            if (dx*dx + dy*dy < 0.5) return true;
            if (staged.count((int)pid)) return true;
            return false;
        };

    if (!singleChain) {
        // Multi-chain: soft repulsion at the 4 cardinal ring sites keeps the
        // ring clear between injections.  Single-chain mode has no ring repulsion
        // — the chain simply diffuses outward from the centre.
        const double ringRepulsion = 50.0;
        callbacks.nonPairwiseCallback =
            [icx, icy, ringRepulsion](unsigned int, const double* pos, const double*) -> double {
                int px = (int)round(pos[0]) - icx;
                int py = (int)round(pos[1]) - icy;
                if ((px == 0 && (py == 1 || py == -1)) ||
                    (py == 0 && (px == 1 || px == -1)))
                    return ringRepulsion;
                return 0.0;
            };
        callbacks.isNonPairwise = true;
    }

    vmmc::VMMC vmmc(nParticles, dimension, coordinates.data(), orientations.data(),
                     maxTrialTranslation, maxTrialRotation,
                     probTranslate, referenceRadius,
                     maxInteractions, &boxSize[0], isIsotropic.get(), isRepulsive,
                     callbacks, isLattice, nLatticeNeighbours,
                     0.0, n0, phi_reorient);

    vmmc.hydrAlpha = useStokes ? 1.0 : 0.0;
    if (seed != 0) vmmc.rng.setSeed(seed);
    else { seed = std::random_device{}(); vmmc.rng.setSeed(seed); }
    fprintf(stderr, "[SEED] %u\n", seed);

    string trajFile = outPrefix + "_traj.txt";
    string statFile = outPrefix + "_stats.txt";
    FILE* fp_traj = fopen(trajFile.c_str(), "w");
    FILE* fp_stat = fopen(statFile.c_str(), "w");
    if (!fp_traj || !fp_stat) { cerr << "Cannot open output files.\n"; return 1; }
    fprintf(fp_stat, "# step  energy  exitedParticles  exitedPerfect  acceptRatio"
                     "  phase  perfectExited  exitedFull  exitQuality\n");

    string firstPhase   = (t_equil > 0) ? "equil" : (t_denat > 0) ? "denat" : "main";
    int    firstOverride = (t_equil > 0) ? 0 : (t_denat > 0) ? 1 : -1;
    model.phaseGOverride = firstOverride;

    double initEnergy = model.getEnergyExcludingBackbone();
    writeFrame(fp_traj, particles, nCopies, n0, singleChain,
               R_c, cx, cy, 0, initEnergy, 0, 0, 0, 0.0,
               couplingStr, gamma0, firstPhase, refComplexEnergy_nobb, 0);
    fprintf(fp_stat, "0  %.4f  0  0  0.0000  %s  0  0  0.0000\n",
            initEnergy, firstPhase.c_str());

    clock_t startTime = clock();
    long long globalStep             = 0;
    long long totalParticlesExited   = 0;
    long long totalPerfectExited     = 0;
    long long totalFullComplexExited = 0;
    long long lifetimeAccepts        = 0;
    long long lifetimeAttempts       = 0;

    auto runPhase = [&](long long phaseSteps, long long saveEveryN, const string& phaseName) {
        if (phaseSteps <= 0) return;
        cout << "--- Phase: " << phaseName << "  (" << phaseSteps << " steps) ---" << endl;
        for (long long s = 1; s <= phaseSteps; s++) {
            vmmc += nParticles;

            ExitCounts exitResult = handleExitsAndQueue(
                model, particles, nParticles, n0, cells, box,
                cx, cy, R_c, vmmc, injQueue, staged, stageX0,
                refComplexEnergy_full, singleChain);
            totalParticlesExited   += exitResult.particles;
            totalPerfectExited     += exitResult.perfect;
            totalFullComplexExited += exitResult.fullComplex;

            if (singleChain)
                tryThreadNext(particles, nParticles, cells, box, vmmc,
                              icx, icy, staged, threadingChain, threadingIdx, injQueue);
            else
                tryInjectNext(particles, nParticles, cells, box, vmmc,
                              cx, cy, injQueue, staged);

            globalStep++;
            bool doSave     = (s % saveEveryN == 0) || (s == phaseSteps);
            bool doProgress = (s % max(1LL, phaseSteps/10) == 0);

            if (doSave || doProgress) {
                double energy      = model.getEnergyExcludingBackbone();
                double acceptRatio = (vmmc.getAttempts() > 0)
                    ? (double)vmmc.getAccepts() / (double)vmmc.getAttempts() : 0.0;
                if (doSave) {
                    double exitQuality = (totalParticlesExited > 0)
                        ? (double)(totalPerfectExited * n0) / (double)totalParticlesExited : 0.0;
                    writeFrame(fp_traj, particles, nCopies, n0, singleChain,
                               R_c, cx, cy, globalStep, energy,
                               totalParticlesExited, totalPerfectExited,
                               totalFullComplexExited, exitQuality,
                               couplingStr, gamma0, phaseName, refComplexEnergy_nobb,
                               countAssembled(model, particles, nCopies, n0, refComplexEnergy_full));
                    fprintf(fp_stat, "%lld  %.4f  %lld  %lld  %.4f  %s  %lld  %lld  %.4f\n",
                            globalStep, energy,
                            totalParticlesExited, totalPerfectExited,
                            acceptRatio, phaseName.c_str(), totalPerfectExited,
                            totalFullComplexExited,
                            (totalParticlesExited > 0)
                                ? (double)(totalPerfectExited * n0) / totalParticlesExited : 0.0);
                    lifetimeAccepts  += vmmc.getAccepts();
                    lifetimeAttempts += vmmc.getAttempts();
                    vmmc.reset();
                }
                if (doProgress) {
                    cout << "  [" << phaseName << "] step " << s << "/" << phaseSteps
                         << "  E=" << energy
                         << "  exitParts=" << totalParticlesExited
                         << "  exitPerfect=" << totalPerfectExited
                         << "  accept=" << acceptRatio << "\n";
                }
            }
        }
    };

    model.phaseGOverride = 0;
    runPhase(t_equil, saveEvery_equil, "equil");
    model.phaseGOverride = 1;
    runPhase(t_denat, saveEvery_denat, "denat");
    model.phaseGOverride = -1;
    runPhase(nsteps, saveEvery_main, "main");

    double simTime = (clock() - startTime) / (double)CLOCKS_PER_SEC;
    cout << "Done! Time = " << simTime << " s (" << simTime/60.0 << " min)" << endl;
    cout << "Total particles exited: "      << totalParticlesExited   << endl;
    cout << "Full-size complexes exited: "  << totalFullComplexExited << endl;
    cout << "Perfect complexes exited: "    << totalPerfectExited     << endl;
    {
        double fq = (totalParticlesExited > 0)
            ? (double)(totalPerfectExited * n0) / totalParticlesExited : 0.0;
        cout << "Exit quality: " << fq << endl;
    }
    cout << "Acceptance ratio: "
         << (lifetimeAttempts > 0 ? (double)lifetimeAccepts / lifetimeAttempts : 0.0)
         << endl;

    fclose(fp_traj);
    fclose(fp_stat);
    cout << "Trajectory: " << trajFile << "\nStatistics: " << statFile << endl;
    cout << "\nTo visualize:\n  python3 visualize_condensate.py " << trajFile << endl;
    return EXIT_SUCCESS;
}
