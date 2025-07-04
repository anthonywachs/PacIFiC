#include <DS_RigidBody_BuilderFactory.hh>
#include <MAC.hh>
#include <DS_RigidBody.hh>
#include <FS_RigidBody.hh>
#include <DS_Sphere.hh>
#include <DS_3Dcylinder.hh>
#include <DS_2Dcylinder.hh>
#include <DS_3Dbox.hh>
#include <DS_GeneralPolyhedron.hh>
#include <DS_STL.hh>
using std::endl;


//---------------------------------------------------------------------------
DS_RigidBody* DS_RigidBody_BuilderFactory:: create( FS_RigidBody* pgrb )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_RigidBody_BuilderFactory:: create" ) ;

  DS_RigidBody* dsrb = NULL;

  // Build the Direction Splitting rigid body
  switch ( pgrb->get_shape_type() )
  {
    case GEOM_SPHERE:
      dsrb = new DS_Sphere( pgrb );
      break;

    case GEOM_2DCYLINDER:
      dsrb = new DS_2Dcylinder( pgrb );
      break;

    case GEOM_3DCYLINDER:
      dsrb = new DS_3Dcylinder( pgrb );
      break;

    case GEOM_3DBOX:
      dsrb = new DS_3Dbox( pgrb );
      break;

    case GEOM_GENERAL_POLYHEDRON:
      dsrb = new DS_GeneralPolyhedron( pgrb );
      break;

    default:
      MAC::out() << "Unknown geometric shape in "
      	"DS_RigidBody_BuilderFactory::create" << endl;
  }

  return ( dsrb );

}




//---------------------------------------------------------------------------
DS_RigidBody* DS_RigidBody_BuilderFactory:: create( FV_Mesh const* MESH,
                                                    istream& STL_input)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_RigidBody_BuilderFactory:: create" ) ;

  DS_RigidBody* dsrb = NULL;

  // Build the Direction Splitting STL
  dsrb = new DS_STL( MESH, STL_input );

  return ( dsrb );

}
