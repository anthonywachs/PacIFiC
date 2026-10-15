#ifndef _INSERTIONWINDOW_HH_
#define _INSERTIONWINDOW_HH_

#include "ReaderXML.hh"
#include "Vector3.hh"
#include <random>

/** @name Enumerations */
//@{
// -------------------------------------------------------------------------------------------------
/** @brief Various supported insertion window types */
enum InsertionWindowType
{
    /** @brief Window with the shape of a box */
    BOXWINDOW,
    /** @brief Window with the shape of an annulus */
    ANNULUSWINDOW
};

// -------------------------------------------------------------------------------------------------
/** @brief Random generator seed */
enum RandomGeneratorSeed
{
    /** @brief initialized to default value (i.e., 1) */
    RGS_DEFAULT,
    /** @brief initialized to a user provided value */
    RGS_UDEF,
    /** @brief randomly initialized */
    RGS_RANDOM /**< randomly initialized */
};
//@}

// =================================================================================================
/** @brief The class InsertionWindow.

    This class defines a window of either box or annulus shape, and can return a
    random point within the window with a uniform distribution. This class is
    not optimized (performance-wise) for insertion over time.

    @author A.YAZDANI - 2024 - Construction */
// =================================================================================================
template <typename T>
class InsertionWindow
{
private:
    /** @name Parameters */
    //@{
    /** \brief First vector. In case of a box, it is the min point. In case
        of an annulus, it is the bottom point. */
    Vector3<T> m_v1;
    /** \brief Second vector. In case of a box, it is the max point. In case
        of an annulus, it is the height direction. */
    Vector3<T> m_v2;
    /** \brief The inner radius of the annulus. zero for the box. */
    T m_iRad;
    /** \brief The outer radius of the annulus. zero for the box. */
    T m_oRad;
    /** \brief Random generator engine */
    std::mt19937 m_randGenerator;
    /** \brief Uniform distribution [0, 1] */
    std::uniform_real_distribution<T> m_dist;
    /** \brief Type of the window */
    InsertionWindowType m_type;
    //@}

public:
    /**@name Contructors */
    //@{
    /** @brief Default constructor */
    __HOST__
    InsertionWindow();

    /** @brief Constructor for box window with min and max points
        @param minPoint Minimum corner of the box
        @param maxPoint Maximum corner of the box
        @param seed Random seed value (default: 0 for deterministic) */
    __HOST__
    InsertionWindow(const Vector3<T>& minPoint, const Vector3<T>& maxPoint, unsigned seed = 0);

    /** @brief Constructor with XML node and the type of the seed
        @param dn XML node
        @param seed random generator seed */
    __HOST__
    InsertionWindow(DOMNode* dn, RandomGeneratorSeed seed);

    /** @brief Destructor */
    __HOST__
    ~InsertionWindow();
    //@}

    /** @name Methods */
    //@{
    /** @brief Returns a random point with uniform distribution in window */
    __HOST__
    Vector3<T> generateRandomPoint();
    //@}
};

#endif
