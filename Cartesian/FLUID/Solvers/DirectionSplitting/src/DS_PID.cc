#include <DS_PID.hh>

//---------------------------------------------------------------------------
DS_PID*
DS_PID:: create( double const& Kp, double const& Kd, double const& Ki )
//---------------------------------------------------------------------------
{

 DS_PID* result = new DS_PID( Kp, Kd, Ki ) ;

 return( result ) ;

}




//---------------------------------------------------------------------------
DS_PID:: DS_PID( double const& Kp, double const& Kd, double const& Ki )
//---------------------------------------------------------------------------
: _Kp (Kp)
, _Kd (Kd)
, _Ki (Ki)
{

   _integral = 0.;
   _pre_error = 0.;

}




//---------------------------------------------------------------------------
DS_PID:: ~DS_PID( void )
//---------------------------------------------------------------------------
{


}




//---------------------------------------------------------------------------
double DS_PID:: calculate( double const& setpoint,
                           double const& pv,
                           double const& dt )
//---------------------------------------------------------------------------
{

   // Calculate error
   double error = (setpoint - pv)/setpoint;

   // Proportional term
   double Pout = _Kp * error;

   // Integral term
   _integral += error * dt;
   double Iout = _Ki * _integral;

   // Derivative term
   double derivative = (error - _pre_error) / dt;
   double Dout = _Kd * derivative;

   // Calculate total output
   double output = Pout + Iout + Dout;

   // Save error to previous error
   _pre_error = error;

   return output;
}
