#ifndef _BASIC_HH_
#define _BASIC_HH_

#include <math.h>
#include <stdlib.h>


using namespace std;

/** @brief Various constants and type definitions.

    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================   

    /**@name Constants */
    //@{
    const int long_string = 255; /**< long string size */
    const double DEGS_PER_RAD = 57.29577951308232286465; /**< degree per 
    	radian */
    const double RADS_PER_DEG =  0.01745329251994329547; /**< radian per 
    	degree */
    const double PI = 3.14159265358979323846; /**< pi number */
    const double TWO_PI = 6.28318530717958623200; /**< 2 times pi number */
    const double LOWEPS = 1.0e-6; /**< very low precision approximation 
    	constant */
    const double EPSILON = 1.0e-10; /**< low precision approximation constant */
    const double EPSILON2 = 1.0e-15; /**< high precision approximation 
    	constant */
    const double EPSILON3 = 1.0e-20; /**< very high precision approximation 
    	constant */
    const int TIMEFORMAT = 10; /**< number of significant digits after the
    	decimal point to write time in high precision format */
    const int POSITIONFORMAT = 16; /**< number of significant digits after the
    	decimal point to write component position in high precision format */
    //@}


    /** @name Basic methods */
    //@{
    /** @brief Returns whether a real number is approximately 0 with respect to
    EPSILON 
    @param x the real number */
    inline bool eqz( double x ) { return ( fabs(x) <= EPSILON ); }
    
    /** @brief Returns the minimum of 2 real numbers defined as double
    @param x 1st real number 
    @param y 2nd real number */    
    inline double min( double x, double y ) { return ( x > y ? y : x ); }
    
    /** @brief Returns the maximum of 2 real numbers defined as double
    @param x 1st real number 
    @param y 2nd real number */      
    inline double max( double x, double y ) { return ( x < y ? y : x ); }
    
    /** @brief Sests the minimum of 2 real numbers defined as double to these 2
    numbers
    @param x 1st real number 
    @param y 2nd real number */     
    inline void set_min( double& x, double y ) { if (x > y) x = y; }
    
    /** @brief Sests the maximum of 2 real numbers defined as double to these 2
    numbers
    @param x 1st real number 
    @param y 2nd real number */ 
    inline void set_max( double& x, double y ) { if (x < y) x = y; }
    
    /** @brief Returns an angle in radians given an angle in degrees
    @param x angle in degrees */
    inline double rads( double x ) { return ( x * RADS_PER_DEG ); }

    /** @brief Returns an angle in degrees given an angle in radians
    @param x angle in radians */
    inline double degs( double x ) { return ( x * DEGS_PER_RAD ); }
    //@}


    /** @name Enumerations */
    //@{    
    /** @brief Space dimensions */
    enum Direction 
    {
      X, // x direction
      Y, // y direction
      Z, // z direction
      W, // scalar component of quaternions
      NONE // no direction
    };

#endif
