#ifndef _SUPERQUADRIC_CONTACT_HH_
#define _SUPERQUADRIC_CONTACT_HH_

#include "Basic.hh"
#include "Transform.hh"
#include "Matrix.hh"
#include "Superquadric.hh"
#include "PointContact.hh"

/** @brief Golden ratio inverse for line search */
constexpr double PHI_INV = 0.618033988749895;


/** @name Superquadric Contact Detection Functions */
//@{

/** @brief Detect contact between two superquadrics.

    Uses the contact point algorithm (Podlozhnyuk et al. 2017):
    1. calc_contact_point: 4D Newton-Raphson for a single point where
       the gradients are anti-parallel (gradA + mu^2*gradB = 0) and the
       implicit function values are equal (fA = fB).
    2. basic_overlap_algorithm: projects the contact point onto each
       surface and computes overlap from the projected distances along
       the contact normal.

    The 4D system variables are [point(3), mu] with Jacobian:
      [HA + mu^2*HB  | 2*mu*gradB ]
      [gradA - gradB |     0      ]

    @param sqA First superquadric shape
    @param transformA Transform for first superquadric
    @param radiusA Bounding radius
    @param sqB Second superquadric shape
    @param transformB Transform for second superquadric
    @param radiusB Bounding radius
    @return PointContact with contact point, overlap vector and distance
            (negative distance = penetration) */
PointContact detectSuperquadricContact(
    const Superquadric* sqA,
    const Transform& transformA,
    double radiusA,
    const Superquadric* sqB,
    const Transform& transformB,
    double radiusB,
    double* prevContactPt = nullptr
);

//@}

#endif // _SUPERQUADRIC_CONTACT_HH_
