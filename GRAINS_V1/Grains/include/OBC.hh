#ifndef _OBC_HH_
#define _OBC_HH_

#include "BVolume.hh"
#include "Matrix.hh"
#include "Transform.hh"

using namespace solid;

/** @brief The class OBC.

    Oriented Bounding cylinder used prior to GJK.

    @author A.YAZDANI - 2022 - Creation */
// ============================================================================
class OBC : public BVolume
{
  public:
    /**@name Constructors */
    //@{
    /** @brief Default constructor */
    OBC();

    /** @brief Constructor with radius, height, and orientation
    @param r radius
    @param h height
    @param ori initial orientation */
    OBC( double r, double h, Vector3 const& ori );

    /** @brief Copy constructor
    @param obc_ reference OBC */
    OBC( OBC const& obc_ );

    /** @brief Destructur */
    ~OBC();
    //@}

    /**@name Methods */
    //@{
    /** @brief Returns the type of the bounding volume */
    BVolumeType getBVolumeType() const;

    // Returns a clone of the OBC
    BVolume* clone() const;

    /** @brief Returns the OBC radius */
    double getRadius() const;

    /** @brief Returns the OBC height */
    double getHeight() const;

    /** @brief Returns the OBC initial orientation */
    Vector3 const& getInitOrientation() const;

    /** @brief Sets the OBC radius
    @param r new radius */
    void setRadius( double r );

    /** @brief Sets the OBC height
    @param h new height */
    void setHeight( double h );

    /** @brief Sets the OBC initial orientation
    @param ori new initial orientation */
    void setInitOrientation( Vector3 const& ori );

    /** @brief Output operator (is called by <<)
    @param fileOut output stream */
    void writeShape( ostream& fileOut ) const;
    //@}


  private:
    /** @name Parameters */
    //@{
    double m_radius; /**< OBC radius */
    double m_height; /**< OBC height */
    Vector3 m_initOrientation; /**< OBC initial orientation */
    //@}
};



/**@name OBC : External methods */
//@{
/** @brief Returns whether two OBCs are in contact
    @param obcA 1st OBC
    @param obcB 2nd OBC
    @param a2w transformation of the 1st OBC
    @param b2w transformation of the 2nd OBC */
    bool isContactBVolume( OBC const& obcA,
                           OBC const& obcB,
                           Transform const& a2w, 
                           Transform const& b2w );
//@}
#endif
