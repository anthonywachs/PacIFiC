#include <MAC.hh>

#include <MAC_assertions.hh>
#include <MAC_Error.hh>
#include <MAC_Exec.hh>
#include <MAC_Randomizer.hh>
#include <MAC_Root.hh>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <limits>
#include <sstream>
#include <unistd.h>

#include <math.h>
// #ifdef __SUNPRO_CC
#if defined(SOLARIS)
   #include <sunmath.h>
#endif

#ifdef _WIN32
extern "C" { 
	double erf( double x ) ;
	double erfc( double x ) ;
} 
#endif

struct MAC_ERROR
{
   static void n0( std::string const& mac_name,
                   std::string const& c_name ) ;
} ;


std::string MAC::undefined_string = "___undef_string___";


//-------------------------------------------------------------------------
double 
MAC:: pi( void )
//-------------------------------------------------------------------------
{
   static double result = 
      3.14159265358979323846264338327950288419716939937510 ;
   return( result ) ;
}




//-------------------------------------------------------------------------
double 
MAC:: e( void )
//-------------------------------------------------------------------------
{
   static double result =
      2.71828182845904523536028747135266249775724709369995 ;
   return( result ) ;
}




//-------------------------------------------------------------------------
double 
MAC:: euler( void )
//-------------------------------------------------------------------------
{
   static double result = 
      0.57721566490153286060651209008240243104215933593992 ;
   return( result ) ;
}




//-------------------------------------------------------------------------
double 
MAC:: epsilon_double( void )
//-------------------------------------------------------------------------
{
   return( std::numeric_limits<double>::epsilon() ) ;
}




//-------------------------------------------------------------------------
double 
MAC::min_double( void )     
//-------------------------------------------------------------------------
{
   return( std::numeric_limits<double>::min() ) ;
}




//-------------------------------------------------------------------------
double 
MAC:: max_double( void )     
//-------------------------------------------------------------------------
{
   return( std::numeric_limits<double>::max() ) ;
}




//-------------------------------------------------------------------------
int 
MAC:: max_int( void )
//-------------------------------------------------------------------------
{
   return( std::numeric_limits<int>::max() ) ;
}




//-------------------------------------------------------------------------
double 
MAC:: bad_double( void )
//-------------------------------------------------------------------------
{
   return( std::numeric_limits<double>::max() ) ;
}




//-------------------------------------------------------------------------
int 
MAC::bad_int( void )
//-------------------------------------------------------------------------
{
   return( std::numeric_limits<int>::max() ) ;
}




//-------------------------------------------------------------------------
double
MAC:: sqrt( double const x )
//-------------------------------------------------------------------------
{
   using namespace std ;
   return( ::sqrt( x ) ) ;  
}




//-------------------------------------------------------------------------
double
MAC:: pow( double const x, double const exp )
//-------------------------------------------------------------------------
{
   using namespace std ;
   return( ::pow( x, exp ) ) ;  
}




//-------------------------------------------------------------------------
double
MAC:: exp( double const x )
//-------------------------------------------------------------------------
{
   using namespace std ;
   return( ::exp( x ) ) ;  
}




//-------------------------------------------------------------------------
double
MAC:: log( double const x )
//-------------------------------------------------------------------------
{
   using namespace std ;
   return( ::log( x ) ) ;  
}




//-------------------------------------------------------------------------
double
MAC:: log10( double const x )
//-------------------------------------------------------------------------
{
   using namespace std ;
   return( ::log10( x ) ) ;  
}




//-------------------------------------------------------------------------
double
MAC:: sin( double const x )
//-------------------------------------------------------------------------
{
   using namespace std ;
   return( ::sin( x ) ) ;  
}




//------------------------------------------------------------------------
double
MAC:: cos( double const x )
//-------------------------------------------------------------------------
{
   using namespace std ;
   return( ::cos( x ) ) ;  
}




//-------------------------------------------------------------------------
double
MAC:: tan( double const x )
//-------------------------------------------------------------------------
{
   using namespace std ;
   return( ::tan( x ) ) ;  
}




//-------------------------------------------------------------------------
double
MAC:: asin( double const x )
//-------------------------------------------------------------------------
{
   using namespace std ;
   return( ::asin( x ) ) ;  
}




//------------------------------------------------------------------------
double
MAC:: acos( double const x )
//-------------------------------------------------------------------------
{
   using namespace std ;
   return( ::acos( x ) ) ;  
}




//-------------------------------------------------------------------------
double
MAC:: atan( double const x )
//-------------------------------------------------------------------------
{
   using namespace std ;
   return( ::atan( x ) ) ;  
}




//-------------------------------------------------------------------------
double
MAC:: atan2( double const x, double const y )
//-------------------------------------------------------------------------
{
   using namespace std ;
   return( ::atan2( x, y ) ) ;  
}




//-------------------------------------------------------------------------
double
MAC:: sinh( double const x )
//-------------------------------------------------------------------------
{
   using namespace std ;
   return( ::sinh( x ) ) ;  
}




//-------------------------------------------------------------------------
double
MAC:: cosh( double const x )
//-------------------------------------------------------------------------
{
   using namespace std ;
   return( ::cosh( x ) ) ;  
}




//-------------------------------------------------------------------------
double
MAC:: tanh( double const x )
//-------------------------------------------------------------------------
{
   using namespace std ;
   return( ::tanh( x ) ) ;  
}




//------------------------------------------------------------------------
bool
MAC:: double_equality( double const x, double const y, 
                       double const a_dbl_eps, double const a_dbl_min )
//------------------------------------------------------------------------
{
   MAC_LABEL( "MAC:: double_equality" ) ;
   MAC_CHECK_PRE( a_dbl_eps >= 0.0 ) ;
   MAC_CHECK_PRE( a_dbl_min >= 0.0 ) ;

   double status = MAC::bad_double() ;
   return( MAC::double_equality( x, y, a_dbl_eps, a_dbl_min, status ) ) ;
}




//------------------------------------------------------------------------
bool
MAC:: double_equality( double const x, double const y, 
                       double const a_dbl_eps, double const a_dbl_min,
                       double& status )
//------------------------------------------------------------------------
{
   MAC_LABEL( "MAC:: double_equality" ) ;
   MAC_CHECK_PRE( a_dbl_eps >= 0.0 ) ;
   MAC_CHECK_PRE( a_dbl_min >= 0.0 ) ;

   status = MAC::bad_double() ;

   if( x == y )
   {
      status = 0.0 ;
      return true ; // <----  easy way out
   }

   bool result = false ;

   double abs_x = MAC::abs( x ) ;
   double abs_y = MAC::abs( y ) ;

   if( abs_y < a_dbl_min ) // test for closeness to 0.0
   {
      // |y| is undistinguishable from 0.0
      if( abs_x < a_dbl_min )
      {
         // |x| is undistinguishable from 0.0
         status = 0.0 ;
         result = true ;
      }
      else
      {
         result = ratio_close_to_one( x, a_dbl_min, a_dbl_eps, status ) ;
      }
   }
   else if( abs_x < a_dbl_min )
   {
      result = ratio_close_to_one( a_dbl_min, y, a_dbl_eps, status ) ;
   }
   else  // test for closeness to each other
   {
      result = ratio_close_to_one( x, y, a_dbl_eps, status ) ;
   }

   MAC_CHECK_POST( status != MAC::bad_double() ) ;
   return( result ) ;
}




//-------------------------------------------------------------------------
double
MAC:: ceil( double const x )
//-------------------------------------------------------------------------
{
   using namespace std ;
   return( ::ceil( x ) ) ;  
}




//-------------------------------------------------------------------------
double
MAC:: floor( double const x )
//-------------------------------------------------------------------------
{
   using namespace std ;
   return( ::floor( x ) ) ;  
}




//-------------------------------------------------------------------------
double
MAC:: asinh( double const x )
//-------------------------------------------------------------------------
{
#ifdef EXTENDED_MATH
   return( ::asinh( x ) ) ;
#else
   return( log(x + sqrt(x * x + 1.0)) ) ;
#endif
}




//-------------------------------------------------------------------------
double
MAC:: acosh( double const x )
//-------------------------------------------------------------------------
{
#ifdef EXTENDED_MATH
   return( ::acosh( x ) ) ;  
#else
   //return log( x + sqrt(x * x – 1.0 ) ) ;
	return log( x + sqrt( x*x - 1.0 ) ) ;
#endif
}




//-------------------------------------------------------------------------
double
MAC:: atanh( double const x )
//-------------------------------------------------------------------------
{
#ifdef EXTENDED_MATH
   return( ::atanh( x ) ) ;
#else
	return log((1.0+x)/(1.0-x))/2.0 ;
#endif   
}




//-------------------------------------------------------------------------
double
MAC:: j0( double const x )
//-------------------------------------------------------------------------
{
#ifdef EXTENDED_MATH
   return( ::j0( x ) ) ;
#else
#ifdef _WIN32
	return _j0(x) ;
#else
   MAC_ERROR:: n0( "j0", "j0" ) ;
   return( MAC::bad_double() ) ;
#endif
#endif     
}




//-------------------------------------------------------------------------
double
MAC:: j1( double const x )
//-------------------------------------------------------------------------
{
#ifdef EXTENDED_MATH
   return( ::j1( x ) ) ;
#else
#ifdef _WIN32
	return _j1(x) ;
#else
   MAC_ERROR:: n0( "j1", "j1" ) ;
   return( MAC::bad_double() ) ;
#endif
#endif    
}




//-------------------------------------------------------------------------
double
MAC:: jn( int const n, double const x )
//-------------------------------------------------------------------------
{
#ifdef EXTENDED_MATH
   return( ::jn( n, x ) ) ;
#else
#ifdef _WIN32
	return _jn(n,x) ;
#else
   MAC_ERROR:: n0( "jn", "jn" ) ;
   return( MAC::bad_double() ) ;
#endif
#endif    
}




//-------------------------------------------------------------------------
double
MAC:: y0( double const x )
//-------------------------------------------------------------------------
{
#ifdef EXTENDED_MATH
   return( ::y0( x ) ) ;
#else
#ifdef _WIN32
	return _y0(x) ;
#else
   MAC_ERROR:: n0( "y0", "y0" ) ;
   return( MAC::bad_double() ) ;
#endif
#endif    
}




//-------------------------------------------------------------------------
double
MAC:: y1( double const x )
//-------------------------------------------------------------------------
{
#ifdef EXTENDED_MATH
   return( ::y1( x ) ) ;
#else
#ifdef _WIN32
	return _y1(x) ;
#else
   MAC_ERROR:: n0( "y1", "y1" ) ;
   return( MAC::bad_double() ) ;
#endif
#endif    
}




//-------------------------------------------------------------------------
double
MAC:: yn( int const n, double const x )
//-------------------------------------------------------------------------
{
#ifdef EXTENDED_MATH
   return( ::yn( n, x ) ) ;
#else
#ifdef _WIN32
	return _yn(n,x) ;
#else
   MAC_ERROR:: n0( "yn", "yn" ) ;
   return( MAC::bad_double() ) ;
#endif
#endif     
}




//-------------------------------------------------------------------------
double
MAC:: gamma( double const x )
//-------------------------------------------------------------------------
{
    int n = x < 1.5 ? -((int) (2.5 - x)) : (int) (x - 1.5);
    double w = x - (n + 2);
    double y = ((((((((((((-1.99542863674e-7 * w + 1.337767384067e-6) * w - 
        2.591225267689e-6) * w - 1.7545539395205e-5) * w + 
        1.45596568617526e-4) * w - 3.60837876648255e-4) * w - 
        8.04329819255744e-4) * w + 0.008023273027855346) * w - 
        0.017645244547851414) * w - 0.024552490005641278) * w + 
        0.19109110138763841) * w - 0.233093736421782878) * w - 
        0.422784335098466784) * w + 0.99999999999999999;
    if (n > 0) {
        w = x - 1;
        for (int k = 2; k <= n; k++) {
            w *= x - k;
        }
    } else {
        w = 1;
        for (int k = 0; k > n; k--) {
            y *= x - k;
        }
    }
    double result =  w / y ;
    return result ;
}




//-------------------------------------------------------------------------
double
MAC:: lgamma( double const x )
//-------------------------------------------------------------------------
{
    int k;
    double w, t, result, v;
    static double a[22] = {
        9.967270908702825e-5, -1.9831672170162227e-4, 
        -0.00117085315349625822, 0.00722012810948319552, 
        -0.0096221300936780297, -0.04219772092994235254, 
        0.16653861065243609743, -0.04200263501129018037, 
        -0.65587807152061930091, 0.57721566490153514421, 
        0.99999999999999999764, 
        4.67209725901142e-5, -6.812300803992063e-5, 
        -0.00132531159076610073, 0.0073352117810720277, 
        -0.00968095666383935949, -0.0421764281187354028, 
        0.16653313644244428256, -0.04200165481709274859, 
        -0.65587818792782740945, 0.57721567315209190522, 
        0.99999999973565236061
    };
    static double b[98] = {
        -4.587497028e-11, 1.902363396e-10, 
        8.6377323367e-10, 1.15513678861e-8, 
        -2.556403058605e-8, -1.5236723372486e-7, 
        -3.1680510638574e-6, 1.22903704923381e-6, 
        2.334372474572637e-5, 0.00111544038088797696, 
        0.00344717051723468982, 0.03198287045148788384, 
        -0.32705333652955399526, 0.40120442440953927615, 
        -5.184290387e-11, -8.3355121068e-10, 
        -2.56167239813e-9, 1.455875381397e-8, 
        1.3512178394703e-7, 2.9898826810905e-7, 
        -3.58107254612779e-6, -2.445260816156224e-5, 
        -4.417127762011821e-5, 0.00112859455189416567, 
        0.00804694454346728197, 0.04919775747126691372, 
        -0.24818372840948854178, 0.11071780856646862561, 
        3.0279161576e-10, 1.60742167357e-9, 
        -4.05596009522e-9, -5.089259920266e-8, 
        -2.029496209743e-8, 1.35130272477793e-6, 
        3.91430041115376e-6, -2.871505678061895e-5, 
        -2.3052137536922035e-4, 4.5534656385400747e-4, 
        0.01153444585593040046, 0.07924014651650476036, 
        -0.12152192626936502982, -0.07916438300260539592, 
        -5.091914958e-10, -1.15274986907e-9, 
        1.237873512188e-8, 2.937383549209e-8, 
        -3.0621450667958e-7, -7.7409414949954e-7, 
        8.16753874325579e-6, 2.412433382517375e-5, 
        -2.60612176060637e-4, -9.1000087658659231e-4, 
        0.01068093850598380797, 0.11395654404408482305, 
        0.07209569059984075595, -0.10971041451764266684, 
        4.0119897187e-10, -1.3224526679e-10, 
        -1.002723190355e-8, 2.569249716518e-8, 
        2.0336011868466e-7, -1.1809768272606e-6, 
        -3.00660303810663e-6, 4.402212897757763e-5, 
        -1.462405876235375e-5, -0.0016487379559600128, 
        0.00513927520866443706, 0.13843580753590579416, 
        0.32730190978254056722, 0.08588339725978624973, 
        -1.5413428348e-10, 6.4905779353e-10, 
        1.60702811151e-9, -2.655645793815e-8, 
        7.619544277956e-8, 4.7604380765353e-7, 
        -4.90748870866195e-6, 8.21513040821212e-6, 
        1.4804944070262948e-4, -0.00122152255762163238, 
        -8.7425289205498532e-4, 0.1443870369965796831, 
        0.61315889733595543766, 0.55513708159976477557, 
        1.049740243e-11, -2.5832017855e-10, 
        1.39591845075e-9, -2.1177278325e-10, 
        -5.082950464905e-8, 3.7801785193343e-7, 
        -7.3982266659145e-7, -1.088918441519888e-5, 
        1.2491810452478905e-4, -4.9171790705139895e-4, 
        -0.0042570708944826646, 0.13595080378472757216, 
        0.89518356003149514744, 1.31073912535196238583
    };
    static double c[65] = {
        1.16333640008e-8, -8.33156123568e-8, 
        3.832869977018e-7, -1.5814047847688e-6, 
        6.50106723241e-6, -2.74514060128677e-5, 
        1.209015360925566e-4, -5.666333178228163e-4, 
        0.0029294103665559733, -0.0180340086069185819, 
        0.1651788780501166204, 1.1031566406452431944, 
        1.2009736023470742248, 
        1.3842760642e-9, -6.9417501176e-9, 
        3.42976459827e-8, -1.785317236779e-7, 
        9.525947257118e-7, -5.2483007560905e-6, 
        3.02364659535708e-5, -1.858396115473822e-4, 
        0.0012634378559425382, -0.0102594702201954322, 
        0.1243625515195050218, 1.3888709263595291174, 
        2.4537365708424422209, 
        1.298977078e-10, -8.02957489e-10, 
        4.945484615e-9, -3.17563534834e-8, 
        2.092136698089e-7, -1.4252023958462e-6, 
        1.01652510114008e-5, -7.74550502862323e-5, 
        6.537746948291078e-4, -0.006601491253552183, 
        0.0996711934948138193, 1.6110931485817511402, 
        3.9578139676187162939, 
        1.83995642e-11, -1.353537034e-10, 
        9.984676809e-10, -7.6346363974e-9, 
        5.99311464148e-8, -4.868554120177e-7, 
        4.1441957716669e-6, -3.77160856623282e-5, 
        3.805693126824884e-4, -0.0045979851178130194, 
        0.0831422678749791178, 1.7929113303999329439, 
        5.6625620598571415285, 
        3.4858778e-12, -2.97587783e-11, 
        2.557677575e-10, -2.2705728282e-9, 
        2.0702499245e-8, -1.954426390917e-7, 
        1.9343161886722e-6, -2.0479024910257e-5, 
        2.405181940241215e-4, -0.0033842087561074799, 
        0.0713079483483518997, 1.9467574842460867884, 
        7.5343642367587329552
    };
    static double d[7] = {
        -0.00163312359200500807, 8.3644533703385956e-4, 
        -5.9518947575728181e-4, 7.9365057505415415e-4, 
        -0.00277777777735463043, 0.08333333333333309869, 
        0.91893853320467274178
    };

    w = x;
    if (x < 0) {
        w = 1 - x;
    }
    if (w < 0.5) {
        k = w < 0.25 ? 0 : 11;
        result = ((((((((((a[k] * w + a[k + 1]) * w + 
            a[k + 2]) * w + a[k + 3]) * w + a[k + 4]) * w + 
            a[k + 5]) * w + a[k + 6]) * w + a[k + 7]) * w + 
            a[k + 8]) * w + a[k + 9]) * w + a[k + 10]) * w;
        result = -log(result);
    } else if (w < 3.5) {
        t = w - 4.5 / (w + 0.5);
        k = (int) (t + 4);
        t -= k - 3.5;
        k *= 14;
        result = ((((((((((((b[k] * t + b[k + 1]) * t + 
            b[k + 2]) * t + b[k + 3]) * t + b[k + 4]) * t + 
            b[k + 5]) * t + b[k + 6]) * t + b[k + 7]) * t + 
            b[k + 8]) * t + b[k + 9]) * t + b[k + 10]) * t + 
            b[k + 11]) * t + b[k + 12]) * t + b[k + 13];
    } else if (w < 8) {
        k = ((int) w) - 3;
        t = w - (k + 3.5);
        k *= 13;
        result = (((((((((((c[k] * t + c[k + 1]) * t + 
            c[k + 2]) * t + c[k + 3]) * t + c[k + 4]) * t + 
            c[k + 5]) * t + c[k + 6]) * t + c[k + 7]) * t + 
            c[k + 8]) * t + c[k + 9]) * t + c[k + 10]) * t + 
            c[k + 11]) * t + c[k + 12];
    } else {
        v = 1 / w;
        t = v * v;
        result = (((((d[0] * t + d[1]) * t + d[2]) * t + 
            d[3]) * t + d[4]) * t + d[5]) * v + d[6];
        result += (w - 0.5) * log(w) - w;
    }
    if (x < 0) {
        result = log(3.141592653589793238 / sin(x * 3.141592653589793238)) - result;
    }
    return result;
  
}




//-------------------------------------------------------------------------
double
MAC:: erf( double const x )
//-------------------------------------------------------------------------
{
#if defined( EXTENDED_MATH ) || defined( _WIN32 )
   return( ::erf( x ) ) ;   
#else
   MAC_ERROR:: n0( "erf", "erf" ) ;
   return( MAC::bad_double() ) ;
#endif      
}




//-------------------------------------------------------------------------
double
MAC:: erfc( double const x )
//-------------------------------------------------------------------------
{
#if defined( EXTENDED_MATH ) || defined( _WIN32 )
   return( ::erfc( x ) ) ;
#else
   MAC_ERROR:: n0( "erfc", "erfc" ) ;
   return( MAC::bad_double() ) ;
#endif     
}




//-------------------------------------------------------------------------
double
MAC:: incomplete_gamma( double const a,
                        double const x,
                        double const epsilon_min,
                        double const error_bound,
                        size_t const iteration_max )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC:: incomplete_gamma" ) ;
   MAC_CHECK_PRE( a>0. ) ;
   MAC_CHECK_PRE( x>0. ) ;
   MAC_CHECK_PRE( epsilon_min>0. ) ;
   MAC_CHECK_PRE( error_bound>0. ) ;
   
   double gammaa = MAC::gamma(a) ;
   double result = gammaa ;
   if( x<100.0*a*sqrt(a) )
   {
      result = exp(-x+a*log(x)) ;
   
      bool conv = false ;
      if( x<a+1 )
      {
         double s=0.0 ;
         double den = gammaa ;
         double xi = 1.0 ;
         for( size_t i=0 ; i<iteration_max ; i++ )
         {
            den = (a+i)*den ;
            double ds = gammaa*xi/den ;
            xi = xi*x ;
            s += ds ;
            conv = MAC::abs(ds)<error_bound ;
            if( conv ) break ;
         }
         result *= s ;
      }
      else
      {
         double b = x+1-a ;
         double d = 1.0/b ;
         double h=d ;
         double c = 1.0/epsilon_min ;
         for( size_t i=1 ; i<iteration_max ; i++ )
         {
            double an = -( i * ( i - a ) ) ;
            b += 2.0 ;
            d = an*d + b ;
            if( MAC::abs(d)<epsilon_min ) d = epsilon_min ;
            c = b + an/c ;
            if( MAC::abs(c)<epsilon_min ) c = epsilon_min ;
            d = 1.0/d ;
            double del = d*c ;
            h *= del ;
            conv = MAC::abs(del-1.0)<error_bound ;
            if( conv ) break ;
         }
         result = gammaa - result*h ;
      }
      if( !conv )
      {
         std::ostringstream msg ;
         msg << "*** MAC::incomplete_gamma failed" << std::endl ;
         msg << "    a = " << a << " x = " << x << std::endl ;
         msg << "    result = " << result << std::endl ;
         MAC_Error::object()->raise_plain( msg.str() ) ;
      }
   }
   return( result ) ;
}




//-------------------------------------------------------------------------
double
MAC:: En( size_t const n, double const x,
          double const epsilon_min,
          double const error_bound,
          size_t const iteration_max )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC:: En" ) ;
   MAC_CHECK_PRE( x>=0. ) ;
   MAC_CHECK_PRE( epsilon_min>0. ) ;
   MAC_CHECK_PRE( error_bound>0. ) ;
   MAC_CHECK_PRE( IMPLIES( x<epsilon_min, n != 0 && n != 1 ) ) ;
   
   double result = MAC::bad_double() ;
   
   if( n == 0 )
   {
      result = MAC::exp(-x)/x ;
   }
   else if( x<epsilon_min )
   {
      result = 1./(n-1.) ;
   }
   else if( x > 1. )
   {
      size_t const nm1 = n-1 ;
      double b = x+n ;
      double c = 1./epsilon_min ;
      double d = 1./b ;
      double h = d ;
      bool conv = false ;
      double a,del ;
      for( size_t i=1 ; !conv && i<= iteration_max ; ++i )
      {
         a = -( (double) i*(nm1+i) ) ;
         b += 2. ;
//         d = 1./(a*d+b) ;
         d= a*d+b ;
         if( MAC::abs(d)<epsilon_min ) d = epsilon_min ;
         d = 1./d ;
         c = b+a/c ;
         if( MAC::abs(c)<epsilon_min ) c = epsilon_min ;
         del = c*d ;
         h *= del ;
         conv = ( MAC::abs(del-1.)<error_bound ) ;
      }
      if( !conv )
      {
         MAC_Error::object()->raise_plain(
            "*** MAC::En computation failed\n"
            "    continued fraction failed " ) ;
      }
      else
      {
         result = h*MAC::exp(-x) ;
      }
   }
   else
   {
      int const nm1 = (int) (n-1) ;
      result = ( nm1 != 0 ) ? 1.0/nm1 : (-MAC::log( x ) - MAC::euler() ) ;
      double fact = 1. ;
      bool conv = false ;
      double del ;
      for( int i=1 ; !conv && i<=(int) iteration_max ; ++i )
      {
         fact *= -x/i ;
         if( i != nm1 )
         {
            del = -fact/(i-nm1) ;
         }
         else
         {
            double psi = -MAC::euler() ;
            for( int ii=1 ; ii<=nm1 ; ++ii ) psi += 1./ii ;
            del = fact*(-MAC::log(x)+psi) ;
         }
         result += del ;
         conv = ( MAC::abs(del)<error_bound*MAC::abs(result) ) ;
      }
      if( !conv )
      {
         MAC_Error::object()->raise_plain(
            "*** MAC::En computation failed\n"
            "    series failed " ) ;
      }
   }
   return( result ) ;
}




//-------------------------------------------------------------------------
double
MAC:: Ei( double const x,
          double const epsilon_min,
          double const error_bound,
          size_t const iteration_max )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC:: Ei" ) ;
   MAC_CHECK_PRE( x>=0. ) ;
   MAC_CHECK_PRE( epsilon_min>0. ) ;
   MAC_CHECK_PRE( error_bound>0. ) ;
   
   double result = MAC::bad_double() ;
   
   if( x < epsilon_min )
   {
      result = MAC::log(x)+MAC::euler() ;
   }
   else if( x<=-MAC::log(error_bound) )
   {
      double sum = 0. ;
      double fact = 1. ;
      bool conv = false ;
      double term ;
      for( size_t k=1 ; !conv && k<=iteration_max ; ++k )
      {
         fact *= x/k ;
         term = fact/k ;
         sum += term ;
         conv = ( term < error_bound*sum ) ;
      }
      if( !conv )
      {
         MAC_Error::object()->raise_plain(
            "*** MAC::Ei computation failed\n"
            "    series failed " ) ;
      }
      else
      {
         result = sum+MAC::log(x)+MAC::euler() ;
      }
   }
   else
   {
      double sum = 0. ;
      double term = 1. ;
      double prev ;
      bool conv = false ;
      for( size_t k=1 ; !conv && k<=iteration_max ; ++k )
      {
         prev = term ;
         term *= k/x ;
         if( term<error_bound )
         {
            conv = true ;
         }
         else if( term<prev )
         {

            sum += term ;
         }
         else
         {
            sum -= prev ;
            conv = true ;
         }  
      }
      result = MAC::exp(x)*(1.+sum)/x ;
   }
   return( result ) ;
}




//-------------------------------------------------------------------------
double
MAC:: random_double( void )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC:: random_double" ) ;
   
   static MAC_Randomizer* randomizer = 0 ;
   if( randomizer == 0 ) 
   {
      randomizer = MAC_Randomizer::create( MAC_Root::object(), 24091971 ) ;
      randomizer->start() ;
   }
   double result = randomizer->item() ;
   randomizer->go_next() ;
   
   MAC_CHECK_POST( result>=0.0 && result<=1.0 ) ;
   return( result ) ;
}




//-------------------------------------------------------------------------
int
MAC:: rand( void )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC:: rand" ) ;

   int result = (int) ( random_double()*max_int() ) ;

   MAC_CHECK_POST( result>0 ) ;
   return( result ) ;
}




//------------------------------------------------------------------------
std::ostream&
MAC::out( void )
//------------------------------------------------------------------------
{
   return( MAC_Exec::out() ) ;
}




//------------------------------------------------------------------------
std::ostream&
MAC::err( void )
//------------------------------------------------------------------------
{
   return( MAC_Exec::err() ) ;
}




//------------------------------------------------------------------------
std::istream&
MAC::in( void )
//------------------------------------------------------------------------
{
   return( MAC_Exec::in() ) ;
}




//----------------------------------------------------------------------------
void
MAC:: print_double( std::ostream& os, double data  )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC:: print_double" ) ;
   
   std::ostringstream stream ;
   stream.precision( os.precision() ) ;
   stream << std::setiosflags( std::ios::scientific )
          << data ;
   std::string s = stream.str() ;
   size_t i = MAC::min( s.find_first_of( 'E' ), s.find_first_of( 'e' ) );
   if( i != std::string::npos )
   {
      while( i>0 && s[i-1]=='0' )
      {
         s.erase( i-1, 1 ) ;
         --i ;
      }
   }
   os << s ;
}




//-------------------------------------------------------------------------
double
MAC:: relative( double const v1, double const v2 )
//-------------------------------------------------------------------------
{
   double r = MAC::abs( v1-v2 ) ;
   double s = MAC::abs( v1 ) + MAC::abs( v2 ) ;
   if( s>1e-10 ) r/=s ;
   return( r ) ;  
}




//-------------------------------------------------------------------------
bool
MAC:: toler( double const val, double const eps )
//-------------------------------------------------------------------------
{
   return( MAC::abs(val) < eps ) ;  
}




//-------------------------------------------------------------------------
bool
MAC:: equal( double const val1, double const val2, double const eps )
//-------------------------------------------------------------------------
{
   return( toler( relative( val1, val2 ), eps ) ) ;  
}




//------------------------------------------------------------------------
bool
MAC:: ratio_close_to_one( double const x, double const y, 
                          double const a_dbl_eps,
                          double& status )
//------------------------------------------------------------------------
{
   MAC_LABEL( "MAC:: ratio_close_to_one" ) ;

   static double dbl_max = MAC::max_double() ;
   static double dbl_min = MAC::min_double() ;

   bool result = false ;

   double abs_x = MAC::abs( x ) ;
   double abs_y = MAC::abs( y ) ;

   // avoid overflow if |y| is very small and |x| is very large
   // (the part (|y|<1.0) is required to avoid overflow in (|y|*dbl_max))
   if( (abs_y < 1.0) && (abs_x > abs_y*dbl_max) )
   {
      result = false ;
      status = -1.0 ;
   }
   // avoid underflow if |y| is very large and |x| is very small
   // (the part (|y|>1.0) is required to avoid underflow in (|y|*dbl_min))
   else if( (abs_y > 1.0) && (abs_x < abs_y*dbl_min) )
   {
      result = false ;
      status = -1.0 ;
   }
   else
   {
      double ratio = x/y ;
      if( ( ratio > (1.0-a_dbl_eps) ) && ( ratio < (1.0+a_dbl_eps) ) ) 
      {
         result = true ;
         status = 0.0 ;
      }
      else
      {
         result = false ;
         status = MAC::abs( x/y-1.0 ) ;
      }
   }
   return( result ) ;
}




//---------------------------------------------------------------------------
void
MAC:: replace( std::string& str,
               std::string const& src,
               std::string const& dest )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC:: replace" ) ;
   MAC_CHECK_PRE( dest.find(src) >= dest.length() ) ;
   
   size_t idx ;
   while( ( idx = str.find(src) ) < str.length() )
   {
      if( dest.length() > 0 )
         str.replace( idx, src.length(), dest ) ;
      else
         str.erase( idx, src.length() ) ;
   }
   MAC_CHECK_POST( str.find(src) >= str.length() ) ;
}




//---------------------------------------------------------------------------
void
MAC:: remove_enclosing_characters( std::string& str,
                                  char const car )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "MAC:: remove_enclosing_characters" ) ;
   
   if( str.length()>2 && str[0]==car && str[str.length()-1]==car)
      str = str.substr(1,str.length()-2) ;
}




//---------------------------------------------------------------------------
std::string 
MAC::doubleToString( std::ios_base::fmtflags format, int digits,
      	double const& number )
//--------------------------------------------------------------------------- 	
{
   MAC_LABEL( "MAC:: doubleToString" ) ;
   
   std::ostringstream oss;
   oss.setf( format, std::ios::floatfield );
   oss.precision( digits );
   oss << number; 
   
   return oss.str();  
}




//---------------------------------------------------------------------------
std::string 
MAC::intToString( int const& int_number )
//--------------------------------------------------------------------------- 	
{
   MAC_LABEL( "MAC:: intToString" ) ;
   
   std::ostringstream oss;
   oss << int_number; 
   
   return oss.str();    
}




//---------------------------------------------------------------------------
std::string 
MAC::extract_root_directory( std::string const& name )
//--------------------------------------------------------------------------- 	
{
   MAC_LABEL( "MAC:: extract_root_directory" ) ;
   
   std::string result;
   size_t pos = name.find_last_of( "/" );  
   if ( pos == std::string::npos ) result = "."; 
   else result = name.substr( 0, pos );   
   
   return( result ) ;    
}




//---------------------------------------------------------------------------
std::string 
MAC::extract_root_file_name( std::string const& name )
//--------------------------------------------------------------------------- 	
{
   MAC_LABEL( "MAC:: extract_root_file_name" ) ;
   
   std::string result;
   size_t pos = name.find_last_of( "/" );  
   if ( pos == std::string::npos ) result = name; 
   else result = name.substr( pos + 1 );
   pos = result.find_first_of( "." );
   if ( pos != std::string::npos ) result = result.substr( 0, pos );   
   
   return( result ) ;    
}




//---------------------------------------------------------------------------
std::string 
MAC::extract_file_name( std::string const& name )
//--------------------------------------------------------------------------- 	
{
   MAC_LABEL( "MAC:: extract_file_name" ) ;
   
   std::string result;
   size_t pos = name.find_last_of( "/" );  
   if ( pos == std::string::npos ) result = name; 
   else result = name.substr( pos + 1 );
   
   return( result ) ;    
}




//-------------------------------------------------------------------------
unsigned long long int
MAC:: used_memory( void )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC:: used_memory" ) ;
   unsigned long long int result = 0 ;
 
   // For Linux systems only !!!  
   std::ostringstream os ;
   std::string word ;
   
   os << "/proc/" << getpid() << "/status" ;
      
   std::ifstream in( os.str().c_str() ) ;
   if( !in )
     MAC::out() << "MAC_System : Unable to open " << os.str() << std::endl ;
   else
   {
     while(!in.eof())
     {
       in >> word ;
       if( word == "VmSize:" )
       {
         in >> result ;
         in >> word ;
         if( !( word == "kB" ) )
           MAC_Error::object()->raise_plain( "Unit is "+word ) ;
         result *= 1000 ;
	 break ;
       }
     }
     in.close() ;
   }

   return( result ) ;
}




//----------------------------------------------------------------------
void
MAC:: display_memory( std::ostream& os, 
	unsigned long long int const& memory )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC:: display_memory" ) ;
   MAC_CHECK_PRE( os ) ;

   static size_t const mo = 1024*1024 ;
   static size_t const go = 1024*1024*1024 ;

   if( memory > go )
   {
      os << ( (double) memory )/go << " Go" ;
   }
   else if( memory > mo )
   {
      os << ( (double) memory )/mo << " Mo" ;
   }
   else
   {
      os << memory << " octets" ;
   }
}




//internal--------------------------------------------------------------
void 
MAC_ERROR:: n0( std::string const& mac_name, std::string const& c_name )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "*** MAC::" << mac_name << " not supported" << std::endl ;
   mesg << "*** it uses the c extended mathematical function ::" << c_name
        << std::endl ;
   mesg << "*** that is not supported by this system." << std::endl ;
   MAC_Error::object()->raise_plain( mesg.str() ) ;
}
