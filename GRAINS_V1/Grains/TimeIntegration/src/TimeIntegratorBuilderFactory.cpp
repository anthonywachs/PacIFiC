#include "GrainsExec.hh"
#include "TimeIntegratorBuilderFactory.hh"
#include "TimeIntegrator.hh"
#include "SecondOrderLeapFrog.hh"
#include "FirstOrderExplicit.hh"
#include "SecondOrderExplicit.hh"
#include "SecondOrderAdamsBashforth.hh"
#include <assert.h>


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Creation de de l'integrateur en time
TimeIntegrator* TimeIntegratorBuilderFactory::create()
{
  TimeIntegrator* TI = NULL;

  if ( GrainsExec::m_TIScheme == "SecondOrderLeapFrog" )
    TI = new SecondOrderLeapFrog();
  else if ( GrainsExec::m_TIScheme == "FirstOrderExplicit" )
    TI = new FirstOrderExplicit();
  else if ( GrainsExec::m_TIScheme == "SecondOrderExplicit" )
    TI = new SecondOrderExplicit();
  else if ( GrainsExec::m_TIScheme == "SecondOrderAdamsBashforth" )
    TI = new SecondOrderAdamsBashforth();  
  else {
    cout << "Wrong type of time integration schema";
    cout << " in <TimeIntegration Type=\"xxx\"/>" << endl; 
    cout << "Allowed entries for xxx are:" << endl;
    cout << "  * SecondOrderLeapFrog" << endl;
    cout << "  * FirstOrderExplicit" << endl;    
    cout << "  * SecondOrderExplicit" << endl;
    cout << "  * SecondOrderAdamsBashforth" << endl;    
  }       

  assert( TI != NULL ); 

  return ( TI );
}
