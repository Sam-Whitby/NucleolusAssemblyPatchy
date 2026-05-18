/*
  run_box.cpp

  Periodic-box assembly simulation.

  All particles remain in the system at all times (no recycling).
  Starts with nCopies fully assembled target complexes.

  Two phases:
    --t-denat N   Denaturation: N steps at gamma=0 (breaks complexes apart)
    --steps   N   Main run:
                    without --anneal: gamma=1 throughout
                    with    --anneal: gamma ramps linearly from 0 to 1

  Options
    --steps        N    main-phase outer iterations                [10000]
    --snapshots    N    trajectory frames saved                    [1000]
    --t-denat      N    denaturation iterations at gamma=0         [0]
    --copies       N    number of target complex copies            [4]
    --box-size     L    side length of the square periodic box     [20]
    --t-after      N    steps at gamma=1 after main phase          [0]
    --anneal            ramp gamma 0→1 over --steps iterations
    --stokes            enable Stokes hydrodynamic drag (D ∝ 1/R)
    --phi-rot      φ    fraction of cluster rotation moves         [0.2]
    --phi-reorient φ    fraction of in-place reorientation moves   [0.2]
    --output   PREFIX   prefix for output files                    [box]
    --seed         S    RNG seed (0 = time-based)                  [1]
    --J            J    patch coupling strength                    [8]
    --single            build a single self-folding Moore-curve chain
    --Moore-N      N    Moore curve order (4^N particles, N>=1)    [2]
*/

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <array>

#include "Demo.h"
#include "VMMC.h"
#include "BoxModel.h"

using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

extern double INF;
extern double TOL;

// ============================================================
//  Multi-chain target complex T  (n=2 Moore curve)
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
//  Moore curve L-system generator
//  Axiom: LFL+F+LFL
//  L → -RF+LFL+FR-    R → +LF-RFR-FL+
//  F=forward, +=left 90°, -=right 90°
//  Apply N-1 expansions; path has 4^N points on a 2^N × 2^N grid.
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
            int tmp = dx; dx = -dy; dy = tmp;   // left 90°
        } else if (c == '-') {
            int tmp = dx; dx = dy; dy = -tmp;   // right 90°
        }
    }

    int minX = path[0].first, minY = path[0].second;
    for (auto& p : path) { minX = min(minX, p.first); minY = min(minY, p.second); }
    for (auto& p : path) { p.first -= minX; p.second -= minY; }
    return path;
}

// ============================================================
//  Multi-chain build functions (use static arrays, unchanged)
// ============================================================
static void buildCouplingMatrices(
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

    for (int i = 0; i < N0_MULTI; i++) {
        for (int j = 0; j < N0_MULTI; j++) {
            if (i == j) continue;
            if (polyType(i) == polyType(j)) continue;
            if (targetDistSqd_multi(i, j) < 1.0 + TOL)
                wD1[i][j] = J;
        }
    }
}

static vector<array<bool,4>> buildPatchSlots()
{
    vector<array<bool,4>> slots(N0_MULTI);
    for (auto& s : slots) s.fill(false);

    for (int i = 0; i < N0_MULTI; i++) {
        for (int j = 0; j < N0_MULTI; j++) {
            if (i == j) continue;
            if (polyType(i) == polyType(j)) continue;
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
    }
    return slots;
}

static void buildBackboneTriples(int nCopies, double bbEnergy,
                                  vector<Triple>& north, vector<Triple>& east)
{
    north.clear();
    east.clear();
    for (int c = 0; c < nCopies; c++) {
        int base = c * N0_MULTI;
        for (int k = 0; k < N_BB_PAIRS; k++) {
            int gi = base + BACKBONE_PAIRS[k][0];
            int gj = base + BACKBONE_PAIRS[k][1];
            east.push_back({gi, gj, bbEnergy});
            east.push_back({gj, gi, bbEnergy});
            north.push_back({gi, gj, bbEnergy});
            north.push_back({gj, gi, bbEnergy});
        }
    }
}

static double referenceComplexEnergy(double J, double bbEnergy)
{
    double E = (double)N_BB_PAIRS * -bbEnergy;
    for (int i = 0; i < N0_MULTI; i++) {
        for (int j = i+1; j < N0_MULTI; j++) {
            if (polyType(i) == polyType(j)) continue;
            if (targetDistSqd_multi(i, j) < 1.0 + 1e-6) E -= J;
        }
    }
    return E;
}

// ============================================================
//  Single-chain build functions (use generated Moore curve path)
// ============================================================

// coupling: d=1 non-backbone pairs only; backbone = |i-j|==1
static void buildCouplingMatricesSingle(
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

    for (int i = 0; i < n0; i++) {
        for (int j = 0; j < n0; j++) {
            if (i == j) continue;
            if (abs(i - j) == 1) continue;  // skip backbone neighbours
            int dx = path[i].first  - path[j].first;
            int dy = path[i].second - path[j].second;
            if (dx*dx + dy*dy < 1.0 + 1e-6)
                wD1[i][j] = J;
        }
    }
}

static vector<array<bool,4>> buildPatchSlotsSingle(
    int n0, const vector<pair<int,int>>& path)
{
    vector<array<bool,4>> slots(n0);
    for (auto& s : slots) s.fill(false);

    for (int i = 0; i < n0; i++) {
        for (int j = 0; j < n0; j++) {
            if (i == j) continue;
            if (abs(i - j) == 1) continue;  // backbone, no directed patch needed
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
    }
    return slots;
}

static void buildBackboneTriplesSingle(int n0, int nCopies, double bbEnergy,
                                        vector<Triple>& north, vector<Triple>& east)
{
    north.clear();
    east.clear();
    for (int c = 0; c < nCopies; c++) {
        int base = c * n0;
        for (int k = 0; k < n0 - 1; k++) {
            int gi = base + k;
            int gj = base + k + 1;
            east.push_back({gi, gj, bbEnergy});
            east.push_back({gj, gi, bbEnergy});
            north.push_back({gi, gj, bbEnergy});
            north.push_back({gj, gi, bbEnergy});
        }
    }
}

static double referenceComplexEnergySingle(
    int n0, const vector<pair<int,int>>& path, double J, double bbEnergy)
{
    double E = (double)(n0 - 1) * -bbEnergy;
    for (int i = 0; i < n0; i++) {
        for (int j = i + 2; j < n0; j++) {  // j>=i+2 skips backbone (|i-j|=1)
            int dx = path[i].first  - path[j].first;
            int dy = path[i].second - path[j].second;
            if (dx*dx + dy*dy < 1.0 + 1e-6) E -= J;
        }
    }
    return E;
}

// ============================================================
//  Count assembled copies
// ============================================================
static int countAssembled(BoxModel& model, const vector<Particle>& particles,
                           int nCopies, int n0, double refComplexEnergy)
{
    int count = 0;
    for (int c = 0; c < nCopies; c++) {
        int base = c * n0;
        double E = 0.0;
        for (int i = base; i < base + n0; i++) {
            for (int j = i+1; j < base + n0; j++) {
                double e = model.computePairEnergy(
                    i, &particles[i].position[0], &particles[i].orientation[0],
                    j, &particles[j].position[0], &particles[j].orientation[0]);
                if (e < 1e5) E += e;
            }
        }
        if (fabs(E - refComplexEnergy) < 0.5) count++;
    }
    return count;
}

// ============================================================
//  Place particles as assembled complexes
//  spacing = stride between copy origins (use gridSize + 3)
// ============================================================
static void placeParticlesAssembled(vector<Particle>& particles,
                                     CellList& cells, Box& box,
                                     int nCopies, int n0,
                                     const vector<pair<int,int>>& path,
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
        int col = c % nc_grid;
        int row = c / nc_grid;
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
//  Write one trajectory frame
// ============================================================
static void writeFrame(FILE* fp, const vector<Particle>& particles,
                        int nCopies, int n0, bool singleMode,
                        double L_box, long long step, double energy, double gamma,
                        double refEnergy, int nAssembled,
                        const string& phase = "main")
{
    int nParticles = (int)particles.size();
    fprintf(fp, "%d\n", nParticles);
    fprintf(fp, "step=%lld energy=%.6f gamma=%.6f L=%.1f nCopies=%d refEnergy=%.4f nAssembled=%d phase=%s single=%d n0=%d\n",
            step, energy, gamma, L_box, nCopies, refEnergy, nAssembled, phase.c_str(),
            (int)singleMode, n0);
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
//  Main
// ============================================================
int main(int argc, char* argv[])
{
    long long nsteps    = 10000;
    long long nsnaps    = 1000;
    long long t_equil   = 0;
    long long t_denat   = 0;
    long long t_after   = 0;
    int       nCopies   = 4;
    double    L_box     = 20.0;
    bool      useAnneal = false;
    bool      useStokes = false;
    double    phi_rot     = 0.2;
    double    phi_reorient = 0.2;
    string    outPrefix = "box";
    unsigned int seed   = 1;
    double    J         = 8.0;
    bool      singleChain = false;
    int       mooreN      = 2;

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i],"--steps")       && i+1<argc) { nsteps      = atoll(argv[++i]); }
        else if (!strcmp(argv[i],"--snapshots")   && i+1<argc) { nsnaps      = atoll(argv[++i]); }
        else if (!strcmp(argv[i],"--t-equil")     && i+1<argc) { t_equil     = atoll(argv[++i]); }
        else if (!strcmp(argv[i],"--t-denat")     && i+1<argc) { t_denat     = atoll(argv[++i]); }
        else if (!strcmp(argv[i],"--t-after")     && i+1<argc) { t_after     = atoll(argv[++i]); }
        else if (!strcmp(argv[i],"--copies")      && i+1<argc) { nCopies     = atoi(argv[++i]); }
        else if (!strcmp(argv[i],"--box-size")    && i+1<argc) { L_box       = atof(argv[++i]); }
        else if (!strcmp(argv[i],"--anneal"))                  { useAnneal   = true; }
        else if (!strcmp(argv[i],"--stokes"))                  { useStokes   = true; }
        else if (!strcmp(argv[i],"--phi-rot")      && i+1<argc) { phi_rot     = atof(argv[++i]); }
        else if (!strcmp(argv[i],"--phi-reorient") && i+1<argc) { phi_reorient= atof(argv[++i]); }
        else if (!strcmp(argv[i],"--output")      && i+1<argc) { outPrefix   = argv[++i]; }
        else if (!strcmp(argv[i],"--seed")        && i+1<argc) { seed = (unsigned int)atoi(argv[++i]); }
        else if (!strcmp(argv[i],"--J")           && i+1<argc) { J           = atof(argv[++i]); }
        else if (!strcmp(argv[i],"--single"))                  { singleChain = true; }
        else if (!strcmp(argv[i],"--Moore-N")     && i+1<argc) { mooreN      = atoi(argv[++i]); }
        else {
            cerr << "Unknown argument: " << argv[i] << "\n";
        }
    }

    if (mooreN < 1) { cerr << "--Moore-N must be >= 1\n"; return 1; }

    // ============================================================
    //  Compute N0 and path from mode selection
    // ============================================================
    int n0;
    vector<pair<int,int>> path;

    if (singleChain) {
        path = generateMooreCurve(mooreN);
        n0   = (int)path.size();
        int expected = 1 << (2 * mooreN);
        if (n0 != expected) {
            cerr << "Moore curve generation error: got " << n0
                 << " points, expected " << expected << "\n";
            return 1;
        }
    } else {
        n0 = N0_MULTI;
        for (int i = 0; i < N0_MULTI; i++)
            path.push_back({TARGET_X_MULTI[i], TARGET_Y_MULTI[i]});
        if (mooreN != 2)
            cerr << "Warning: --Moore-N is only used with --single; ignored.\n";
    }

    // ============================================================
    //  Build interaction data
    // ============================================================
    const double bbEnergy = 1000.0;

    vector<vector<double>> wD1, wDsq2, wD2, wDsq5;
    vector<Triple> north, east;
    vector<array<bool,4>> patchSlots;
    double refComplexEnergy;

    if (singleChain) {
        buildCouplingMatricesSingle(n0, path, J, wD1, wDsq2, wD2, wDsq5);
        patchSlots      = buildPatchSlotsSingle(n0, path);
        buildBackboneTriplesSingle(n0, nCopies, bbEnergy, north, east);
        refComplexEnergy = referenceComplexEnergySingle(n0, path, J, bbEnergy);
    } else {
        buildCouplingMatrices(J, wD1, wDsq2, wD2, wDsq5);
        patchSlots      = buildPatchSlots();
        buildBackboneTriples(nCopies, bbEnergy, north, east);
        refComplexEnergy = referenceComplexEnergy(J, bbEnergy);
    }

    // Spacing between copy origins: grid side length + 3-unit gap
    int gridSize = singleChain ? (1 << mooreN) : 4;
    int spacing  = gridSize + 3;

    // Count folding contacts for display
    int nFoldingContacts = 0;
    for (int i = 0; i < n0; i++)
        for (int j = i+1; j < n0; j++)
            if (wD1[i][j] != 0.0) nFoldingContacts++;

    const int    nParticles = nCopies * n0;
    const unsigned int dimension       = 2;
    const double interactionRange      = 2.5;
    const unsigned int maxInteractions = 100;
    const double interactionEnergy     = 0.0;
    const bool isLattice               = true;

    cout << "=== Box Assembly Simulation (patchy model) ===" << endl;
    cout << "  mode=" << (singleChain ? "single-chain" : "multi-chain")
         << "  N0=" << n0 << "  mooreN=" << mooreN
         << "  foldingContacts=" << nFoldingContacts << endl;
    cout << "  steps=" << nsteps << " snapshots=" << nsnaps
         << "  t_equil=" << t_equil << "  t_denat=" << t_denat << "  t_after=" << t_after
         << "  copies=" << nCopies << "  box=" << L_box << "x" << L_box
         << "  anneal=" << useAnneal << " stokes=" << useStokes
         << "  J=" << J
         << "  phi_rot=" << phi_rot << "  phi_reorient=" << phi_reorient << endl;

    vector<vector<int>> bbPartners(n0);
    double springK = 0.0;

    Interactions interactions(nParticles, n0, north, east,
                               wD1, wDsq2, wD2, wDsq5,
                               springK, bbPartners);

    interactions.patchesEnabled = true;
    interactions.patchSlots     = patchSlots;

    // --- Simulation box: fully periodic ---
    vector<double> boxSize    = { L_box, L_box };
    vector<bool>   isPeriodic = { true, true };
    Box box(boxSize, isPeriodic);
    box.isLattice = true;

    // --- Cell list ---
    CellList cells;
    cells.setDimension(dimension);
    cells.initialise(box.boxSize, interactionRange);

    // --- Particles ---
    vector<Particle> particles(nParticles);

    // --- BoxModel ---
    BoxModel model(box, particles, cells,
                   maxInteractions, interactionEnergy, interactionRange,
                   interactions);
    model.uniformGamma = 1.0;
    model.singleChain  = singleChain;

    cout << "Reference perfect-complex energy (g=1): " << refComplexEnergy << endl;

    // --- Place particles as assembled complexes ---
    double cx = L_box / 2.0;
    double cy = L_box / 2.0;
    placeParticlesAssembled(particles, cells, box, nCopies, n0, path, spacing, cx, cy);

    // --- VMMC setup ---
    vector<double> coordinates(dimension * nParticles);
    vector<double> orientations(dimension * nParticles);
    // vector<bool> has no .data(); use raw array (all false = not isotropic for patchy)
    unique_ptr<bool[]> isIsotropic(new bool[nParticles]());

    for (int i = 0; i < nParticles; i++) {
        coordinates[2*i]     = particles[i].position[0];
        coordinates[2*i + 1] = particles[i].position[1];
        orientations[2*i]     = 1.0;
        orientations[2*i + 1] = 0.0;
    }

    double maxTrialTranslation = 1.5;
    double maxTrialRotation    = (phi_rot > 0.0 || phi_reorient > 0.0) ? M_PI : 0.0;
    double probTranslate       = 1.0 - phi_rot - phi_reorient;
    double referenceRadius     = 0.5;
    bool   isRepulsive         = true;
    int    nLatticeNeighbours  = 8;

    using namespace std::placeholders;
    vmmc::CallbackFunctions callbacks;
    callbacks.energyCallback =
        std::bind(&BoxModel::computeEnergy, &model, _1, _2, _3);
    callbacks.pairEnergyCallback =
        std::bind(&BoxModel::computePairEnergy, &model, _1, _2, _3, _4, _5, _6);
    callbacks.interactionsCallback =
        std::bind(&BoxModel::computeInteractions, &model, _1, _2, _3, _4);
    callbacks.postMoveCallback =
        std::bind(&BoxModel::applyPostMoveUpdates, &model, _1, _2, _3);

    vmmc::VMMC vmmc(nParticles, dimension, &coordinates[0], &orientations[0],
                     maxTrialTranslation, maxTrialRotation,
                     probTranslate, referenceRadius,
                     maxInteractions, &boxSize[0], isIsotropic.get(), isRepulsive,
                     callbacks, isLattice, nLatticeNeighbours,
                     0.0 /*phi_sl: disabled*/, n0, phi_reorient);

    vmmc.hydrAlpha = useStokes ? 1.0 : 0.0;
    if (seed != 0) vmmc.rng.setSeed(seed);
    else { seed = std::random_device{}(); vmmc.rng.setSeed(seed); }
    fprintf(stderr, "[SEED] %u\n", seed);

    // Snapshot intervals distributed proportionally across all phases
    long long totalSteps    = t_equil + t_denat + nsteps + t_after;
    long long snaps_equil   = (totalSteps > 0 && t_equil > 0)
                              ? max(1LL, nsnaps * t_equil / totalSteps) : 0;
    long long snaps_denat   = (totalSteps > 0 && t_denat > 0)
                              ? max(1LL, nsnaps * t_denat / totalSteps) : 0;
    long long snaps_after   = (totalSteps > 0 && t_after > 0)
                              ? max(1LL, nsnaps * t_after / totalSteps) : 0;
    long long snaps_main    = nsnaps - snaps_equil - snaps_denat - snaps_after;
    if (snaps_main < 1) snaps_main = 1;

    auto saveEvery = [](long long steps, long long snaps) -> long long {
        if (steps <= 0 || snaps <= 0) return steps;
        return max(1LL, steps / snaps);
    };
    long long saveEvery_equil = saveEvery(t_equil, snaps_equil);
    long long saveEvery_denat = saveEvery(t_denat, snaps_denat);
    long long saveEvery_main  = saveEvery(nsteps,  snaps_main);
    long long saveEvery_after = saveEvery(t_after,  snaps_after);

    // --- Output files ---
    string trajFile = outPrefix + "_traj.txt";
    string statFile = outPrefix + "_stats.txt";
    FILE* fp_traj = fopen(trajFile.c_str(), "w");
    FILE* fp_stat = fopen(statFile.c_str(), "w");
    if (!fp_traj || !fp_stat) {
        cerr << "Cannot open output files.\n";
        return 1;
    }
    fprintf(fp_stat, "# step  energy  acceptRatio  gamma  nAssembled  phase\n");

    auto writeStep0 = [&](const string& phase, double gamma) {
        model.uniformGamma = gamma;
        double initEnergy  = model.getEnergy() * nParticles;
        int assembled = countAssembled(model, particles, nCopies, n0, refComplexEnergy);
        writeFrame(fp_traj, particles, nCopies, n0, singleChain,
                   L_box, 0, initEnergy, gamma, refComplexEnergy, assembled, phase);
        fprintf(fp_stat, "0  %.4f  0.0000  %.4f  %d  %s\n",
                initEnergy, gamma, assembled, phase.c_str());
    };

    {
        double gamma0_init = (t_equil > 0) ? 1.0 :
                             (t_denat > 0) ? 0.0 :
                             (useAnneal)   ? 0.0 : 1.0;
        string firstPhase  = (t_equil > 0) ? "equil" :
                             (t_denat > 0) ? "denat" : "main";
        writeStep0(firstPhase, gamma0_init);
    }

    long long stepOffset = 0;

    // ============================================================
    //  Phase 1: Equilibration (t_equil steps at gamma=1)
    // ============================================================
    if (t_equil > 0) {
        model.uniformGamma = 1.0;
        vmmc.reset();

        for (long long step = 1; step <= t_equil; step++) {
            for (int k = 0; k < nParticles; k++) vmmc.step();

            if (step % saveEvery_equil == 0 || step == t_equil) {
                double energy      = model.getEnergy() * nParticles;
                double acceptRatio = (vmmc.getAttempts() > 0)
                                     ? (double)vmmc.getAccepts() / vmmc.getAttempts()
                                     : 0.0;
                int assembled = countAssembled(model, particles, nCopies, n0, refComplexEnergy);
                writeFrame(fp_traj, particles, nCopies, n0, singleChain,
                           L_box, stepOffset + step, energy, 1.0,
                           refComplexEnergy, assembled, "equil");
                fprintf(fp_stat, "%lld  %.4f  %.4f  %.4f  %d  equil\n",
                        stepOffset + step, energy, acceptRatio, 1.0, assembled);
                fflush(fp_stat);
                vmmc.reset();
            }
        }
        stepOffset += t_equil;
    }

    // ============================================================
    //  Phase 2: Denaturation (t_denat steps at gamma=0)
    // ============================================================
    if (t_denat > 0) {
        model.uniformGamma = 0.0;
        vmmc.reset();

        for (long long step = 1; step <= t_denat; step++) {
            for (int k = 0; k < nParticles; k++) vmmc.step();

            if (step % saveEvery_denat == 0 || step == t_denat) {
                double energy      = model.getEnergy() * nParticles;
                double acceptRatio = (vmmc.getAttempts() > 0)
                                     ? (double)vmmc.getAccepts() / vmmc.getAttempts()
                                     : 0.0;
                int assembled = countAssembled(model, particles, nCopies, n0, refComplexEnergy);
                writeFrame(fp_traj, particles, nCopies, n0, singleChain,
                           L_box, stepOffset + step, energy, 0.0,
                           refComplexEnergy, assembled, "denat");
                fprintf(fp_stat, "%lld  %.4f  %.4f  %.4f  %d  denat\n",
                        stepOffset + step, energy, acceptRatio, 0.0, assembled);
                fflush(fp_stat);
                vmmc.reset();
            }
        }
        stepOffset += t_denat;
    }

    // ============================================================
    //  Phase 3: Main run
    // ============================================================
    {
        double gamma_main0 = useAnneal ? 0.0 : 1.0;
        model.uniformGamma = gamma_main0;
        vmmc.reset();

        for (long long step = 1; step <= nsteps; step++) {
            if (useAnneal && nsteps > 1)
                model.uniformGamma = (double)(step - 1) / (double)(nsteps - 1);
            else
                model.uniformGamma = 1.0;

            for (int k = 0; k < nParticles; k++) vmmc.step();

            if (step % saveEvery_main == 0 || step == nsteps) {
                double energy      = model.getEnergy() * nParticles;
                double acceptRatio = (vmmc.getAttempts() > 0)
                                     ? (double)vmmc.getAccepts() / vmmc.getAttempts()
                                     : 0.0;
                double g = model.uniformGamma;
                int assembled = countAssembled(model, particles, nCopies, n0, refComplexEnergy);
                writeFrame(fp_traj, particles, nCopies, n0, singleChain,
                           L_box, stepOffset + step, energy, g,
                           refComplexEnergy, assembled, "main");
                fprintf(fp_stat, "%lld  %.4f  %.4f  %.4f  %d  main\n",
                        stepOffset + step, energy, acceptRatio, g, assembled);
                fflush(fp_stat);
                vmmc.reset();
            }
        }
        stepOffset += nsteps;
    }

    // ============================================================
    //  Phase 4: After (t_after steps at gamma=1)
    // ============================================================
    if (t_after > 0) {
        model.uniformGamma = 1.0;
        vmmc.reset();

        for (long long step = 1; step <= t_after; step++) {
            for (int k = 0; k < nParticles; k++) vmmc.step();

            if (step % saveEvery_after == 0 || step == t_after) {
                double energy      = model.getEnergy() * nParticles;
                double acceptRatio = (vmmc.getAttempts() > 0)
                                     ? (double)vmmc.getAccepts() / vmmc.getAttempts()
                                     : 0.0;
                int assembled = countAssembled(model, particles, nCopies, n0, refComplexEnergy);
                writeFrame(fp_traj, particles, nCopies, n0, singleChain,
                           L_box, stepOffset + step, energy, 1.0,
                           refComplexEnergy, assembled, "after");
                fprintf(fp_stat, "%lld  %.4f  %.4f  %.4f  %d  after\n",
                        stepOffset + step, energy, acceptRatio, 1.0, assembled);
                fflush(fp_stat);
                vmmc.reset();
            }
        }
    }

    fclose(fp_traj);
    fclose(fp_stat);

    cout << "Trajectory written to: " << trajFile << endl;
    cout << "Statistics written to:  " << statFile << endl;

    return 0;
}
