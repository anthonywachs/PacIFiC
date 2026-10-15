// =================================================================================================
/** @brief GJK Iteration Count vs. Distance Validation

    This program measures how the surface-to-surface distance between two convex
    particles affects the number of GJK iterations required for convergence.

    Four shape-pair combinations are tested:
      - cylinder vs cylinder
      - box vs box
      - cylinder vs box
      - superquadric vs box

    For each pair we sweep the nominal separation gap over several orders of
    magnitude (log-uniform from ~1e-4 to ~1e2 times the particle characteristic
    size).  At each nominal gap, many random orientations are sampled; only the
    non-overlapping trials are kept.

    Output: data/gjk_iterations.csv
      Columns: ShapePair, NominalGap, ActualDistance, NbIter

    @author A.Yazdani - 2026 - GJK Iterations Validation */
// =================================================================================================

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "Box.hh"
#include "Cylinder.hh"
#include "GJK.hh"
#include "Quaternion.hh"
#include "Superquadric.hh"
#include "Transform3.hh"
#include "Vector3.hh"

// -------------------------------------------------------------------------------------------------
/** @brief Generate a uniformly distributed random unit quaternion.

    Uses the Marsaglia method: draw four independent standard-normal samples and
    normalise the resulting 4-vector.
    @param rng Mersenne-Twister generator to draw from */
static Quaternion<double> randomUnitQuat(std::mt19937& rng)
{
    std::normal_distribution<double> dist(0.0, 1.0);
    const double                     x = dist(rng);
    const double                     y = dist(rng);
    const double                     z = dist(rng);
    const double                     w = dist(rng);
    const double                     n = std::sqrt(x * x + y * y + z * z + w * w);
    return Quaternion<double>(x / n, y / n, z / n, w / n);
}

// -------------------------------------------------------------------------------------------------
/** @brief Run the GJK iteration sweep for one shape pair.

    Shape A is kept at the origin with identity orientation.  Shape B is placed
    along the x-axis at distance
        d_center = rSumAB + nominalGap
    from the origin. Both shapes receive independent random orientations for each
    trial.  If the GJK-computed distance is positive (no overlap) and GJK
    converged (nbIter < 1000), the record is written to the CSV file.

    A small non-zero crust (1% of the shape's circumscribed radius) is used on
    each shape so that GJK is called under the same conditions as the actual
    simulation kernel, which always uses non-zero crust values.  The crust
    slightly shrinks each shape for the distance query and does not affect the
    nominal geometry.
    @param name          label for this shape pair
    @param shapeA        first convex shape (not modified)
    @param shapeB        second convex shape (not modified)
    @param rSumAB        sum of circumscribed radii (used to offset centres)
    @param crustA        crust thickness for shape A (metres)
    @param crustB        crust thickness for shape B (metres)
    @param gapValues     vector of nominal gap distances to test
    @param nTrials       number of random-orientation trials per gap
    @param rng           random-number generator
    @param out           open output file stream */
static void runPair(const std::string&         name,
                    const Convex<double>&      shapeA,
                    const Convex<double>&      shapeB,
                    double                     rSumAB,
                    double                     crustA,
                    double                     crustB,
                    const std::vector<double>& gapValues,
                    int                        nTrials,
                    std::mt19937&              rng,
                    std::ofstream&             out)
{
    Vector3<double> pa, pb;

    for(double gap : gapValues)
    {
        const double dCenter = rSumAB + gap;

        for(int t = 0; t < nTrials; ++t)
        {
            // Random orientations for A and B
            const Quaternion<double> qA = randomUnitQuat(rng);
            const Quaternion<double> qB = randomUnitQuat(rng);

            // World-frame transforms
            const Transform3<double> a2w(qA, Vector3<double>(0.0, 0.0, 0.0));
            const Transform3<double> b2w(qB, Vector3<double>(dCenter, 0.0, 0.0));

            // Build relative transform b2a = inv(a2w) * b2w
            const Transform3<double> b2a(a2w, b2w);

            uint         nbIter = 0u;
            const double dist   = computeClosestPoints_GJK<double, GJKType::SIGNEDVOLUME>(shapeA,
                                                                                        shapeB,
                                                                                        b2a,
                                                                                        crustA,
                                                                                        crustB,
                                                                                        pa,
                                                                                        pb,
                                                                                        nbIter);

            // Only record non-overlapping, converged configurations.
            // nbIter == 1000 indicates the iteration cap was hit (GJK did not
            // converge); these degenerate cases are excluded from statistics.
            if(dist > 0.0 && nbIter < 1000u)
            {
                out << name << "," << std::scientific << std::setprecision(6) << gap << ","
                    << std::scientific << std::setprecision(6) << dist << "," << nbIter << "\n";
            }
        }
    }
}

// -------------------------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    // ---- parameters ---------------------------------------------------------
    // Number of random-orientation trials per (pair, gap) combination
    const int nTrials = (argc > 1) ? std::stoi(argv[1]) : 100;
    // Random seed
    const unsigned int seed = (argc > 2) ? static_cast<unsigned int>(std::stoul(argv[2])) : 42u;
    // Output file
    const std::string csvPath = (argc > 3) ? argv[3] : "data/gjk_iterations.csv";

    std::cout << "=================================================================\n"
              << "GJK Iterations vs Distance Validation\n"
              << "  Trials per gap : " << nTrials << "\n"
              << "  Random seed    : " << seed << "\n"
              << "  Output CSV     : " << csvPath << "\n"
              << "=================================================================\n";

    // ---- shape definitions --------------------------------------------------
    // Characteristic size: all shapes fit inside a sphere of radius ~0.1 m
    //
    // Cylinder :  radius = 0.07 m,  halfHeight = 0.07 m
    //             circumscribed radius = sqrt(0.07^2 + 0.07^2) ~= 0.0990
    //
    // Box      :  half-extents (0.0577, 0.0577, 0.0577)
    //             circumscribed radius = sqrt(3 * 0.0577^2) ~= 0.0999
    //
    // Superquadric : semi-axes (0.0577, 0.0577, 0.0577), n1=n2=4
    //                circumscribed radius same as Box

    const double cyl_r  = 0.07;
    const double cyl_hh = 0.07;
    const double box_he = 0.0577;
    const double sq_he  = 0.0577;

    Cylinder<double>     cyl(cyl_r, cyl_hh);
    Box<double>          box(box_he, box_he, box_he);
    Superquadric<double> sqr(sq_he, sq_he, sq_he, 4.0, 4.0);

    const double rCyl = cyl.computeCircumscribedRadius();
    const double rBox = box.computeCircumscribedRadius();
    const double rSqr = sqr.computeCircumscribedRadius();

    // Zero crust: GJK operates on the full unshrunken shapes so that the
    // returned distance is the true surface-to-surface gap, matching the
    // Pacific implementation which has no crust parameter.
    const double crustCyl = 0.0;
    const double crustBox = 0.0;
    const double crustSqr = 0.0;

    std::cout << "Circumscribed radii:\n"
              << "  Cylinder      : " << rCyl << " m\n"
              << "  Box           : " << rBox << " m\n"
              << "  Superquadric  : " << rSqr << " m\n\n";

    // ---- gap sweep (log-uniform) -------------------------------------------
    // Reference length: rCyl (~ 0.1)
    // Gap range: [1e-4 * ref, 1e2 * ref]  => 7 decades
    const double refLen = rCyl;
    const int    nGaps  = 140;  // ~20 samples per decade
    const double gapMin = 1.0e-4 * refLen;
    const double gapMax = 1.0e+2 * refLen;
    const double logMin = std::log10(gapMin);
    const double logMax = std::log10(gapMax);

    std::vector<double> gapValues(static_cast<size_t>(nGaps));
    for(int i = 0; i < nGaps; ++i)
        gapValues[static_cast<size_t>(i)]
            = std::pow(10.0, logMin + (logMax - logMin) * i / (nGaps - 1));

    // ---- output -------------------------------------------------------------
    std::ofstream out(csvPath);
    if(!out.is_open())
    {
        std::cerr << "Error: cannot open output file '" << csvPath << "'\n";
        return 1;
    }
    out << "ShapePair,NominalGap,ActualDistance,NbIter\n";

    std::mt19937 rng(seed);

    // ---- shape pairs --------------------------------------------------------
    struct PairSpec
    {
        std::string           name;
        const Convex<double>* shapeA;
        const Convex<double>* shapeB;
        double                rSumAB;
        double                crustA;
        double                crustB;
    };

    const std::vector<PairSpec> pairs = {
        {"cylinder-cylinder", &cyl, &cyl, rCyl + rCyl, crustCyl, crustCyl},
        {"box-box", &box, &box, rBox + rBox, crustBox, crustBox},
        {"cylinder-box", &cyl, &box, rCyl + rBox, crustCyl, crustBox},
        {"superquadric-superquadric", &sqr, &sqr, rSqr + rSqr, crustSqr, crustSqr},
    };

    for(const auto& p : pairs)
    {
        std::cout << "Running pair: " << p.name << " ...\n";
        runPair(p.name,
                *p.shapeA,
                *p.shapeB,
                p.rSumAB,
                p.crustA,
                p.crustB,
                gapValues,
                nTrials,
                rng,
                out);
        out.flush();
    }

    out.close();
    std::cout << "\nDone. Results written to '" << csvPath << "'.\n";
    return 0;
}
