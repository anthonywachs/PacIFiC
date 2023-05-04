#include "BVolume.hh"
#include "OBB.hh"
#include "OBC.hh"

// --------------------------------------------------------------------
// Default constructor
BVolume::BVolume()
{}




// --------------------------------------------------------------------
// Destructor
BVolume::~BVolume()
{}




// ----------------------------------------------------------------------------
// Output operator
ostream& operator << ( ostream& fileOut, BVolume const& bvol_ )
{
  bvol_.writeShape( fileOut );

  return ( fileOut );
}




// ----------------------------------------------------------------------------
// Returns whether the bounding volumes are in contact
bool isContact( BVolume const& bvolA, 
                BVolume const& bvolB,
                Transform const& a2w,
                Transform const& b2w )
{
  if ( bvolA.getBVolumeType() == typeOBB )
    return( isContactBVolume( (OBB const&) bvolA, 
                              (OBB const&) bvolB,
                              a2w, 
                              b2w ) );

  if ( bvolA.getBVolumeType() == typeOBC )
    return( isContactBVolume( (OBC const&) bvolA, 
                              (OBC const&) bvolB,
                              a2w, 
                              b2w ) );

  return ( true );
}




