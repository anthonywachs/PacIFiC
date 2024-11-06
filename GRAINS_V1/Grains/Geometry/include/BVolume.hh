#ifndef _BVOLUME_HH_
#define _BVOLUME_HH_

#include "Transform.hh"
#include "ReaderXML.hh"

enum BVolumeType {
  typeAABB,
  typeOBB,
  typeOBC
};

/** @brief The class BVolume.

    Bounding volume used prior to GJK.

    @author A.YAZDANI - 2023 - Creation */
// ============================================================================
class BVolume
{
  public:
    /**@name Constructors */
    //@{
    /** @brief Default constructor */
    BVolume();

    /** @brief Destructeur */
    virtual ~BVolume();
    //@}

    /**@name Virtual methods */
    //@{
    /** @brief Returns the type of the bounding volume */
    virtual BVolumeType getBVolumeType() const = 0;

    // Returns a clone of the OBB
    virtual BVolume* clone() const = 0;

    /** @brief Output operator ( is called by << )
    @param fileOut output stream */
    virtual void writeShape( ostream &fileOut ) const = 0;
    //@}


    /** @name Friend methods */
    //@{
    /** @brief Output operator
    @param f output stream
    @param bvol_ BVolume object */
    friend ostream& operator << ( ostream& f, BVolume const& bvol_ );
    //@}
};


/**@name BVolume : External methods */
//@{
/** @brief Returns whether two bounding volumes are in contact
    @param a 1st bounding volume
    @param b 2nd bounding volume
    @param a2w transformation of the 1st bounding volume
    @param b2w transformation of the 2nd bounding volume */
    bool isContact( BVolume const& a,
                    BVolume const& b,
                    Transform const& a2w, 
                    Transform const& b2w );
//@}
    
#endif
