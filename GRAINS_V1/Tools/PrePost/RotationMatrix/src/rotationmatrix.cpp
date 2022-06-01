#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <cmath>
#include <Basic.hh>
#include <Matrix.hh>
using namespace std;


string doubleToString( ios_base::fmtflags format, int digits,
      	double const& number )
{ 
  ostringstream oss;
  if ( number != 0. )
  {
    oss.setf( format, ios::floatfield );
    oss.precision( digits );
  }
  oss << number; 
   
  return ( oss.str() ); 
}


int main(int argc, char *argv[])
{
  double angleX,angleY,angleZ;
  string filename="rotationmatrix.txt"; 
  
  cout << "Rotation angle around X : ";
  cin >> angleX;
  cout << "Rotation angle around Y : ";
  cin >> angleY;  
  cout << "Rotation angle around Z : ";
  cin >> angleZ;

  angleX *= PI / 180.;
  angleY *= PI / 180.;  
  angleZ *= PI / 180.; 
  
  Matrix rZ(cos(angleZ),-sin(angleZ),0.,sin(angleZ),cos(angleZ),0.,0.,0.,1.);
  Matrix rX(1.,0.,0.,0.,cos(angleX),-sin(angleX),0.,sin(angleX),cos(angleX));
  Matrix rY(cos(angleY),0.,sin(angleY),0.,1.,0.,-sin(angleY),0.,cos(angleY));
  Matrix tmp = rY * rZ;
  Matrix rotation = rX * tmp; 
  
  cout << "The computed rotation matrix is M = Rx * Ry * Rz" << endl;
  cout << "Output is written in local file \"rotationmatrix.txt\" in XML and"
  	<< " can be copied & pasted as is in the Grains input file" << endl;  

  ofstream fileOUT(filename.c_str(),ios::out);
  fileOUT << "        <AngularPosition Type=\"Matrix\">" << endl;
  fileOUT << "        " << 
  	doubleToString( ios::scientific, 12, rotation[X][X] ) << "   " <<
	doubleToString( ios::scientific, 12, rotation[X][Y] ) << "   " <<	 
	doubleToString( ios::scientific, 12, rotation[X][Z] ) << endl;
  fileOUT << "        " << 
  	doubleToString( ios::scientific, 12, rotation[Y][X] ) << "   " <<
	doubleToString( ios::scientific, 12, rotation[Y][Y] ) << "   " <<	 
	doubleToString( ios::scientific, 12, rotation[Y][Z] ) << endl;
  fileOUT << "        " << 
  	doubleToString( ios::scientific, 12, rotation[Z][X] ) << "   " <<
	doubleToString( ios::scientific, 12, rotation[Z][Y] ) << "   " <<	 
	doubleToString( ios::scientific, 12, rotation[Z][Z] ) << endl;		 
  fileOUT << "        </AngularPosition>" << endl;    
          
  return(0);
  
}
