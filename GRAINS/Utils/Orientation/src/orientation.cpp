#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <cmath>
#include <Basic.H>
#include <Matrix.H>
using namespace std;


int main(int argc, char *argv[])
{
  double angleX,angleY,angleZ;
  string filename="orientation.txt"; 
  
  cout << "Angle de rotation axe X : ";
  cin >> angleX;
  cout << "Angle de rotation axe Y : ";
  cin >> angleY;  
  cout << "Angle de rotation axe Z : ";
  cin >> angleZ;

  angleX *= PI / 180.;
  angleY *= PI / 180.;  
  angleZ *= PI / 180.; 
  
  Matrix rZ(cos(angleZ),-sin(angleZ),0.,sin(angleZ),cos(angleZ),0.,0.,0.,1.);
  Matrix rX(1.,0.,0.,0.,cos(angleX),-sin(angleX),0.,sin(angleX),cos(angleX));
  Matrix rY(cos(angleY),0.,sin(angleY),0.,1.,0.,-sin(angleY),0.,cos(angleY));
  Matrix tmp = rY * rZ;
  Matrix rotation = rX * tmp;   

  ofstream fileOUT(filename.c_str(),ios::out);
  fileOUT << "          <Orientation Type=\"Matrice\">" << endl;
  fileOUT << "          " << rotation[X][X] << "   " << rotation[X][Y] 
  	<< "   " << rotation[X][Z] << endl;
  fileOUT << "          " << rotation[Y][X] << "   " << rotation[Y][Y] 
  	<< "   " << rotation[Y][Z] << endl;	
  fileOUT << "          " << rotation[Z][X] << "   " << rotation[Z][Y] 
  	<< "   " << rotation[Z][Z] << endl;	 
  fileOUT << "          </Orientation>" << endl;    
          
  return(0);
  
}
