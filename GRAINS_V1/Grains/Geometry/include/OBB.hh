#ifndef _OBB_HH_
#define _OBB_HH_

#include "BVolume.hh"
#include "Matrix.hh"
#include "Transform.hh"

using namespace solid;

/** @brief The class OBB.

    Oriented Bounding box used prior to GJK.

    @author A.YAZDANI - 2023 - Creation */
// ============================================================================
class OBB : public BVolume
{
  public:
    /**@name Constructors */
    //@{
    /** @brief Default constructor */
    OBB();

    /** @brief Constructor with half-lengths extent and initial orientation ori
    @param extent half-lengths
    @param ori initial orientation */
    OBB( Vector3 const& extent, Matrix const& ori );

    /** @brief Copy constructor
    @param obc_ reference OBB */
    OBB( OBB const& obb_ );

    /** @brief Destructur */
    ~OBB();
    //@}

    /**@name Methods */
    //@{
    /** @brief Returns the type of the bounding volume */
    BVolumeType getBVolumeType() const;

    // Returns a clone of the OBB
    BVolume* clone() const;

    /** @brief Returns the OBB half-lengths */
    Vector3 const& getExtent() const;

    /** @brief Returns the OBB initial orientation */
    Matrix const& getInitOrientation() const;

    /** @brief Sets the OBB half-lenghts
    @param extent new half-lengths */
    void setExtent( Vector3 const& extent );

    /** @brief Sets the OBB initial orientation
    @param ori new initial orientation */
    void setInitOrientation( Matrix const& ori );

    /** @brief Output operator (is called by <<)
    @param fileOut output stream */
    void writeShape( ostream& fileOut ) const;
    //@}


  private:
    /** @name Parameters */
    //@{
    Vector3 m_extent; /**< OBB half-lenghts */
    Matrix m_initOrientation; /**< OBB initial orientation */
    //@}
};



/**@name OBB : External methods */
//@{
/** @brief Returns whether two OBBs are in contact
    @param obbA 1st OBB
    @param obbB 2nd OBB
    @param a2w transformation of the 1st OBB
    @param b2w transformation of the 2nd OBB */
    bool isContactBVolume( OBB const& obbA,
                           OBB const& obbB,
                           Transform const& a2w, 
                           Transform const& b2w );
//@}

#endif
