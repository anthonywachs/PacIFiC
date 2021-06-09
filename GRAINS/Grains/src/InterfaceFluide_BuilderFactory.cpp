#include "Grains_BuilderFactory.H"
#include "InterfaceFluide_BuilderFactory.hh"
#include "InterfaceFluide2D.hh"
#include "InterfaceFluide3D.hh"


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Creation de l'interface de couplage avec le Fluide
InterfaceFluide* InterfaceFluide_BuilderFactory::create()
{
  InterfaceFluide* appFluide = NULL;

  EAPPLI context = Grains_BuilderFactory::getContext();

  switch (context)
  {
    case DIM_2:
      appFluide = new InterfaceFluide2D();
      break;
    case DIM_3:
      appFluide = new InterfaceFluide3D();
      break;
    case UNDEFINED:
      appFluide = NULL;
      break;    
  }

  return appFluide;
}
