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
    --t-after      N    steps at gamma=1 after main phase              [0]
    --anneal            ramp gamma 0→1 over --steps iterations
    --stokes            enable Stokes hydrodynamic drag (D ∝ 1/R)
    --phi-rot      φ    fraction of cluster rotation moves         [0.2]
    --phi-reorient φ    fraction of in-place reorientation moves   [0.2]
    --output   PREFIX   prefix for output files                    [box]
    --seed         S    RNG seed (0 = time-based)                  [1]
*/

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
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
//  Target complex T  (n=2 Moore curve, same as nucleolus/condensate)
// ============================================================

static const int N0        = 16;
static const int N_POLYMER = 4;
static const int N_SEG     = 4;

static const int TARGET_X[N0] = { 0,0,1,2,  3,3,2,2,  0,0,1,1,  3,3,2,1 };
static const int TARGET_Y[N0] = { 2,3,3,3,  3,2,2,1,  0,1,1,2,  1,0,0,0 };

static const int BACKBONE_PAIRS[][2] = {
    {0,1},{1,2},{2,3},
    {4,5},{5,6},{6,7},
    {8,9},{9,10},{10,11},
    {12,13},{13,14},{14,15}
};
static const int N_BB_PAIRS = 12;

inline int polyType(int id) { return id / N_SEG; }

inline double targetDistSqd(int id1, int id2) {
    double dx = TARGET_X[id1] - TARGET_X[id2];
    double dy = TARGET_Y[id1] - TARGET_Y[id2];
    return dx*dx + dy*dy;
}

// ============================================================
//  Build coupling matrices (patchy model — cross-type d=1 only)
// ============================================================
static void buildCouplingMatrices(
    double J,
    vector<vector<double>>& wD1,
    vector<vector<double>>& wDsq2,
    vector<vector<double>>& wD2,
    vector<vector<double>>& wDsq5)
{
    wD1.assign(N0, vector<double>(N0, 0.0));
    wDsq2.assign(N0, vector<double>(N0, 0.0));
    wD2.assign(N0, vector<double>(N0, 0.0));
    wDsq5.assign(N0, vector<double>(N0, 0.0));

    for (int i = 0; i < N0; i++) {
        for (int j = 0; j < N0; j++) {
            if (i == j) continue;
            if (polyType(i) == polyType(j)) continue;
            if (targetDistSqd(i, j) < 1.0 + TOL)
                wD1[i][j] = J;
        }
    }
}

// ============================================================
//  Build patch slot map
// ============================================================
static vector<array<bool,4>> buildPatchSlots()
{
    vector<array<bool,4>> slots(N0);
    for (auto& s : slots) s.fill(false);

    for (int i = 0; i < N0; i++) {
        for (int j = 0; j < N0; j++) {
            if (i == j) continue;
            if (polyType(i) == polyType(j)) continue;
            if (targetDistSqd(i, j) > 1.0 + 1e-6) continue;
            int dx = TARGET_X[j] - TARGET_X[i];
            int dy = TARGET_Y[j] - TARGET_Y[i];
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

// ============================================================
//  Build backbone Triples for all nCopies copies
// ============================================================
static void buildBackboneTriples(int nCopies, double bbEnergy,
                                  vector<Triple>& north, vector<Triple>& east)
{
    north.clear();
    east.clear();
    for (int c = 0; c < nCopies; c++) {
        int base = c * N0;
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

// ============================================================
//  Reference energy of one perfect complex (g=1, native geometry)
// ============================================================
static double referenceComplexEnergy(double J, double bbEnergy)
{
    double E = (double)N_BB_PAIRS * -bbEnergy;
    for (int i = 0; i < N0; i++) {
        for (int j = i+1; j < N0; j++) {
            if (polyType(i) == polyType(j)) continue;
            if (targetDistSqd(i, j) < 1.0 + 1e-6) E -= J;
        }
    }
    return E;
}

// ============================================================
//  Count copies whose intra-copy pairwise energy equals refComplexEnergy
//  within tolerance 0.5.  All complex energies are exact integer multiples
//  of J and bbEnergy, so 0.5 cleanly separates assembled from non-assembled.
// ============================================================
static int countAssembled(BoxModel& model, const vector<Particle>& particles,
                           int nCopies, double refComplexEnergy)
{
    int n = 0;
    for (int c = 0; c < nCopies; c++) {
        int base = c * N0;
        double E = 0.0;
        for (int i = base; i < base + N0; i++) {
            for (int j = i+1; j < base + N0; j++) {
                double e = model.computePairEnergy(
                    i, &particles[i].position[0], &particles[i].orientation[0],
                    j, &particles[j].position[0], &particles[j].orientation[0]);
                if (e < 1e5) E += e;
            }
        }
        if (fabs(E - refComplexEnergy) < 0.5) n++;
    }
    return n;
}

// ============================================================
//  Place particles as assembled target complexes in a square grid.
//  Copies are offset by 7 lattice units (3-unit gap > sqrt(5) range).
// ============================================================
static void placeParticlesAssembled(vector<Particle>& particles,
                                     CellList& cells, Box& box,
                                     int nCopies, double cx, double cy)
{
    cells.reset();
    int nParticles = nCopies * N0;
    for (int i = 0; i < nParticles; i++) {
        particles[i].index = i;
        particles[i].position.resize(2);
        particles[i].orientation.resize(2);
        particles[i].orientation[0] = 1.0;
        particles[i].orientation[1] = 0.0;
    }

    int nc_grid = (int)ceil(sqrt((double)nCopies));
    double span = nc_grid * 7 - 3;
    double ox0  = round(cx) - span / 2.0;
    double oy0  = round(cy) - span / 2.0;

    for (int c = 0; c < nCopies; c++) {
        int col = c % nc_grid;
        int row = c / nc_grid;
        double ox = ox0 + col * 7;
        double oy = oy0 + row * 7;
        for (int lid = 0; lid < N0; lid++) {
            int gi = c * N0 + lid;
            particles[gi].position[0] = ox + TARGET_X[lid];
            particles[gi].position[1] = oy + TARGET_Y[lid];
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
                        int nCopies, double L_box,
                        long long step, double energy, double gamma,
                        double refEnergy, int nAssembled,
                        const string& phase = "main")
{
    int nParticles = (int)particles.size();
    fprintf(fp, "%d\n", nParticles);
    fprintf(fp, "step=%lld energy=%.6f gamma=%.6f L=%.1f nCopies=%d refEnergy=%.4f nAssembled=%d phase=%s\n",
            step, energy, gamma, L_box, nCopies, refEnergy, nAssembled, phase.c_str());
    for (int i = 0; i < nParticles; i++) {
        int copy  = i / N0;
        int lid   = i % N0;
        int ptype = polyType(lid);
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
        else {
            cerr << "Unknown argument: " << argv[i] << "\n";
        }
    }

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

    cout << "=== Box Assembly Simulation (patchy model) ===" << endl;
    cout << "  steps=" << nsteps << " snapshots=" << nsnaps
         << "  t_equil=" << t_equil << "  t_denat=" << t_denat << "  t_after=" << t_after
         << "  copies=" << nCopies << "  box=" << L_box << "x" << L_box
         << "  anneal=" << useAnneal << " stokes=" << useStokes
         << "  phi_rot=" << phi_rot << "  phi_reorient=" << phi_reorient << endl;

    // --- Parameters ---
    const int    nParticles = nCopies * N0;
    const double J          = 8.0;
    const double bbEnergy   = 1000.0;

    const unsigned int dimension       = 2;
    const double interactionRange      = 2.5;
    const unsigned int maxInteractions = 100;
    const double interactionEnergy     = 0.0;
    const bool isLattice               = true;

    // Build coupling matrices and patch slots
    vector<vector<double>> wD1, wDsq2, wD2, wDsq5;
    buildCouplingMatrices(J, wD1, wDsq2, wD2, wDsq5);

    vector<Triple> north, east;
    buildBackboneTriples(nCopies, bbEnergy, north, east);

    vector<vector<int>> bbPartners(N0);
    double springK = 0.0;

    Interactions interactions(nParticles, N0, north, east,
                               wD1, wDsq2, wD2, wDsq5,
                               springK, bbPartners);

    auto patchSlots = buildPatchSlots();
    interactions.patchesEnabled = true;
    interactions.patchSlots     = patchSlots;

    // --- Simulation box: fully periodic ---
    vector<double> boxSize   = { L_box, L_box };
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

    const double refComplexEnergy = referenceComplexEnergy(J, bbEnergy);
    cout << "Reference perfect-complex energy (g=1): " << refComplexEnergy << endl;

    // --- Place particles as assembled complexes centred in the box ---
    double cx = L_box / 2.0;
    double cy = L_box / 2.0;
    placeParticlesAssembled(particles, cells, box, nCopies, cx, cy);

    // --- VMMC setup ---
    double coordinates[dimension * nParticles];
    double orientations[dimension * nParticles];
    bool   isIsotropic[nParticles];
    for (int i = 0; i < nParticles; i++) {
        coordinates[2*i]     = particles[i].position[0];
        coordinates[2*i + 1] = particles[i].position[1];
        orientations[2*i]     = 1.0;
        orientations[2*i + 1] = 0.0;
        isIsotropic[i]        = true;
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

    vmmc::VMMC vmmc(nParticles, dimension, coordinates, orientations,
                     maxTrialTranslation, maxTrialRotation,
                     probTranslate, referenceRadius,
                     maxInteractions, &boxSize[0], isIsotropic, isRepulsive,
                     callbacks, isLattice, nLatticeNeighbours,
                     0.0 /*phi_sl: disabled*/, N0, phi_reorient);

    vmmc.hydrAlpha = useStokes ? 1.0 : 0.0;
    if (seed != 0) vmmc.rng.setSeed(seed);
    else { seed = std::random_device{}(); vmmc.rng.setSeed(seed); }
    fprintf(stderr, "[SEED] %u\n", seed);

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

    // Write step-0 frame for the first active phase.
    auto writeStep0 = [&](const string& phase, double gamma) {
        model.uniformGamma = gamma;
        double initEnergy  = model.getEnergy() * nParticles;
        int assembled      = countAssembled(model, particles, nCopies, refComplexEnergy);
        writeFrame(fp_traj, particles, nCopies, L_box, 0, initEnergy, gamma,
                   refComplexEnergy, assembled, phase);
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

    // Cumulative step counter so trajectory steps increase monotonically across phases.
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
                int assembled = countAssembled(model, particles, nCopies, refComplexEnergy);
                writeFrame(fp_traj, particles, nCopies, L_box, stepOffset + step, energy, 1.0,
                           refComplexEnergy, assembled, "equil");
                fprintf(fp_stat, "%lld  %.4f  %.4f  %.4f  %d  equil\n",
                        stepOffset + step, energy, acceptRatio, 1.0, assembled);
                fflush(fp_stat);
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
                int assembled = countAssembled(model, particles, nCopies, refComplexEnergy);
                writeFrame(fp_traj, particles, nCopies, L_box, stepOffset + step, energy, 0.0,
                           refComplexEnergy, assembled, "denat");
                fprintf(fp_stat, "%lld  %.4f  %.4f  %.4f  %d  denat\n",
                        stepOffset + step, energy, acceptRatio, 0.0, assembled);
                fflush(fp_stat);
            }
        }
        stepOffset += t_denat;
    }

    // ============================================================
    //  Phase 3: Main run
    //    --anneal: gamma ramps linearly 0→1 over nsteps
    //    default:  gamma=1 throughout
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
                int assembled = countAssembled(model, particles, nCopies, refComplexEnergy);
                writeFrame(fp_traj, particles, nCopies, L_box, stepOffset + step, energy, g,
                           refComplexEnergy, assembled, "main");
                fprintf(fp_stat, "%lld  %.4f  %.4f  %.4f  %d  main\n",
                        stepOffset + step, energy, acceptRatio, g, assembled);
                fflush(fp_stat);
            }
        }
        stepOffset += nsteps;
    }

    // ============================================================
    //  Phase 4: After (t_after steps at gamma=1, post-anneal equilibration)
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
                int assembled = countAssembled(model, particles, nCopies, refComplexEnergy);
                writeFrame(fp_traj, particles, nCopies, L_box, stepOffset + step, energy, 1.0,
                           refComplexEnergy, assembled, "after");
                fprintf(fp_stat, "%lld  %.4f  %.4f  %.4f  %d  after\n",
                        stepOffset + step, energy, acceptRatio, 1.0, assembled);
                fflush(fp_stat);
            }
        }
    }

    fclose(fp_traj);
    fclose(fp_stat);

    cout << "Trajectory written to: " << trajFile << endl;
    cout << "Statistics written to:  " << statFile << endl;

    return 0;
}
