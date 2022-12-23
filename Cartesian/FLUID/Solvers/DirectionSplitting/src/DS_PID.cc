#include <DS_PID.hh>

//---------------------------------------------------------------------------
DS_PID*
DS_PID:: create( double const& max, double const& min
               , double const& Kp, double const& Kd, double const& Ki )
//---------------------------------------------------------------------------
{

 DS_PID* result = new DS_PID( max, min, Kp, Kd, Ki ) ;

 return( result ) ;

}




//---------------------------------------------------------------------------
DS_PID:: DS_PID( double const& max, double const& min
               , double const& Kp, double const& Kd, double const& Ki )
//---------------------------------------------------------------------------
: _max (max)
, _min (min)
, _Kp (Kp)
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
   double error = setpoint - pv;

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

   // Restrict to max/min
   if( output > _max )
     output = _max;
   else if( output < _min )
     output = _min;

   // Save error to previous error
   _pre_error = error;

   return output;
}
