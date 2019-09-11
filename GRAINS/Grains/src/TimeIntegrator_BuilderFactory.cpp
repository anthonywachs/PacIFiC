#include "Grains_Exec.hh"
#include "TimeIntegrator_BuilderFactory.hh"
#include "TimeIntegrator.hh"
#include "SecondOrderLeapFrog.hh"
#include "FirstOrderExplicit.hh"
#include "SecondOrderExplicit.hh"
#include "SecondOrderAdamsBashforth.hh"
#include "SecondOrderRungeKutta.hh"
#include "ThirdOrderAdamsBashforth.hh"
#include <assert.h>

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Creation de de l'integrateur en temps
TimeIntegrator* TimeIntegrator_BuilderFactory::create()
{
  TimeIntegrator* TI = NULL;

  if (Grains_Exec::m_TIScheme == "SecondOrderLeapFrog")
    TI = new SecondOrderLeapFrog();
  else if (Grains_Exec::m_TIScheme == "FirstOrderExplicit")
    TI = new FirstOrderExplicit();
  else if (Grains_Exec::m_TIScheme == "SecondOrderExplicit")
    TI = new SecondOrderExplicit();
  else if (Grains_Exec::m_TIScheme == "SecondOrderAdamsBashforth")
    TI = new SecondOrderAdamsBashforth();
else if (Grains_Exec::m_TIScheme == "SecondOrderRungeKutta")
    TI = new SecondOrderRungeKutta();
else if (Grains_Exec::m_TIScheme == "ThirdOrderAdamsBashforth")
    TI = new ThirdOrderAdamsBashforth();
  else {
    cout << "Wrong type of time integration schema";
    cout << " in <TimeIntegration Type=\"xxx\"/>" << endl;
    cout << "Allowed entries for xxx are:" << endl;
    cout << "  * SecondOrderLeapFrog" << endl;
    cout << "  * FirstOrderExplicit" << endl;
    cout << "  * SecondOrderExplicit" << endl;
    cout << "  * SecondOrderAdamsBashforth" << endl;
    cout << "  * SecondOrderRungeKutta" << endl;
    cout << "  * ThirdOrderAdamsBashforth" << endl;
  }

  assert(TI != NULL);

  return TI;
}
