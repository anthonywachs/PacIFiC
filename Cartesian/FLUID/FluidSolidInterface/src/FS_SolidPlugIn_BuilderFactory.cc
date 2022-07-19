#include <FS_SolidPlugIn_BuilderFactory.hh>
#include <FS_SolidPlugIn.hh>
#include <FS_Grains3DPlugIn.hh>
#include <MAC.hh>
#include <iostream>
using std::endl;


//---------------------------------------------------------------------------
FS_SolidPlugIn* FS_SolidPlugIn_BuilderFactory:: create( string const& type,
	string const& insertion_file_,
        string const& simulation_file_,
        double const& fluid_density,
        bool const& b_restart,
        bool const& b_initializeClonePer,
        double const& grid_size,
        bool const& is_solidsolver_parallel,
	int& error )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_SolidPlugIn_BuilderFactory:: create" ) ;

  FS_SolidPlugIn* psolid = NULL;
  
  if ( type == "Grains3D" )
    	psolid = new FS_Grains3DPlugIn( insertion_file_, simulation_file_, 
		fluid_density, b_restart, b_initializeClonePer, grid_size, 
		is_solidsolver_parallel, error );
  else 
  {
    MAC::out() << "Warning: unknown solid solver type" << endl;
    error = 1;
  }
      
  return ( psolid );
  
}
