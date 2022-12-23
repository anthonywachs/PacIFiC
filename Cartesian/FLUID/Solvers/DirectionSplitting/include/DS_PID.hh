#ifndef _DS_PID_
#define _DS_PID_

#include <iostream>
#include <cmath>

/** @brief The Class controller to maintain constant flow rate

@author A. Goyal - Pacific project 2022 */

class DS_PID
{
   private: //-----------------------------------------------------------------

      // Kp -  proportional gain
      // Ki -  Integral gain
      // Kd -  derivative gain
      // dt -  loop interval time
      // max - maximum value of manipulated variable
      // min - minimum value of manipulated variable
      DS_PID( double const& max,
              double const& min,
              double const& Kp,
              double const& Kd,
              double const& Ki );


      ~DS_PID(void);

      double _max, _min;
      double _Kp;
      double _Kd;
      double _Ki;
      double _integral;
      double _pre_error;

   public:
      static DS_PID* create(double const& max,
                         double const& min,
                         double const& Kp,
                         double const& Kd,
                         double const& Ki);

      // Returns the manipulated variable given a setpoint and current process value
      double calculate( double const& setpoint,
                        double const& pv,
                        double const& dt );
};

#endif
