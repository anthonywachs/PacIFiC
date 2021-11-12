#include <FS_RigidBody_BuilderFactory.hh>
#include <MAC.hh>
#include <FS_RigidBody.hh>
#include <FS_Sphere.hh>
#include <FS_3Dcylinder.hh>
using std::endl;


//---------------------------------------------------------------------------
FS_RigidBody* FS_RigidBody_BuilderFactory:: create( size_t& dimens, 
	istream& in, size_t& ncorners_, size_t& id_ )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_RigidBody_BuilderFactory:: create" ) ;

  FS_RigidBody* prb = NULL;
  
  // Build the rigid body
  if ( dimens == 2 ) 
  {
  
  }
  else
  {
    switch ( ncorners_ )
    {
      case 1:
        prb = new FS_Sphere( in, id_ );
	break;
	
      case 777:
        prb = new FS_3Dcylinder( in, id_ );
	break;        

      default:
        MAC::out() << "Unknown rigid body shape" << endl;
    }
  }
      
  return ( prb );
  
}
