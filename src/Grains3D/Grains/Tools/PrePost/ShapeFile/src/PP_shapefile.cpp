#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <cmath>
using namespace std;


int main(int argc, char *argv[])
{
  int choice, i, np ;
  bool good_choice = false;
  string filename; 
  ofstream fileOUT;
  double pi = fabs(acos(-1.));
  double x, y, z ;
  double edge, sphere_radius, half_edge, length, width, height, 
  	coord, vol ;  
  double dtheta, theta, theta0, rc, rd, surface, zzz ; 
  double rcirc, h, l ;
  double gold = 0.5 * ( 1. + sqrt(5.) ), sqrt2o2 = sqrt( 2. ) / 2. ;
     
  while(!good_choice)
  {
    cout << "Polygon/polyhedron type : " << endl;
    cout << "  (1) Isometric polygon" << endl;
    cout << "  (2) 3D box" << endl;
    cout << "  (3) Cube" << endl;    
    cout << "  (4) Regular tetrahedron" << endl;
    cout << "  (5) Regular octahedron" << endl;
    cout << "  (6) Regular dodecahedron" << endl;        
    cout << "  (7) Regular icosahedron" << endl;    
    cout << "Choice : ";
    cin >> choice;
  
    if ( choice < 0 || choice > 7 )
      cout << "Wrong choice, select again" << endl;
    else
      good_choice = true;
  }
  
  switch(choice)
  {
    case 1:
      good_choice = false;
      while(!good_choice)
      { 
        cout << "Number of edges = ";
        cin >> np;
  
        if (np<3)
          cout << "Wrong choice, number of edges > 2, enter again" << endl;
        else
          good_choice = true;
      } 
      good_choice = false;
      while(!good_choice)
      { 
        cout << "Equivalent disc radius = ";
        cin >> rd;
  
        if (rd<0.)
          cout << "Wrong choice, radius > 0., enter again" << endl;
        else
          good_choice = true;
      }          
      dtheta = 2. * pi / np;
      theta0 = dtheta / 2.;
      zzz = pi / np;
      rc = sqrt( zzz / (sin(zzz)*cos(zzz)) ) * rd; 
      filename="polygon";
      if (np==3) filename="triangle";
      else if (np==4) filename="square";
      else if (np==5) filename="pentagon";      
      else if (np==6) filename="hexagon";      
      else if (np==7) filename="heptagon";      
      else if (np==8) filename="octogon"; 
      else if (np==9) filename="nonagon";      
      else if (np==10) filename="decagon";       
      filename+=".insert";     
      fileOUT.open(filename.c_str(),ios::out);
      fileOUT << "3" << endl;
      fileOUT << np << endl;
      for (i=0;i<np;++i)
      {
        theta = theta0+i*dtheta;
        x = rc * cos(theta);
	y = rc * sin(theta);
	x = fabs(x) > 1e-12 ? x : 0.;
	y = fabs(y) > 1e-12 ? y : 0.;	
	fileOUT << x << " " << y << " 0.0" << endl;
	if (i==0) surface = x*y;
      } 
      fileOUT << "1" << endl;
      for (i=0;i<np;++i)  fileOUT << i << " ";
      fileOUT << "0" << endl; 
      fileOUT.close(); 
      surface *= np;
      cout << "Polygon surface area = " << surface << endl; 
      cout << "Circumscribed radius = " << rc << endl;
      cout << "Circularity = " << sqrt( zzz / tan(zzz) ) << endl;
      break;

    case 2:
      good_choice = false;
      while(!good_choice)
      { 
        cout << "Length = ";
        cin >> length;
  
        if ( length < 0. )
          cout << "Wrong choice, length > 0, enter again" << endl;
        else
          good_choice = true;
      }
      good_choice = false;
      while(!good_choice)
      { 
        cout << "Depth = ";
        cin >> width;
  
        if ( width < 0. )
          cout << "Wrong choice, depth > 0, enter again" << endl;
        else
          good_choice = true;
      }
      good_choice = false;
      while(!good_choice)
      { 
        cout << "Height = ";
        cin >> height;
  
        if ( height < 0. )
          cout << "Wrong choice, height > 0, enter again" << endl;
        else
          good_choice = true;
      }
      
      cout << "3D Box" << endl;
      cout << "   Length = " << length << endl;
      cout << "   Depth = " << width << endl;      
      cout << "   Height = " << height << endl;      
      
      if ((length==width)&&(length==height)) filename="cube.insert";
      else filename="3Dbox.insert";
      
      x=length/2.;
      y=width/2.;
      z=height/2.;
      fileOUT.open(filename.c_str(),ios::out);
      
      fileOUT << "3" << endl << "8" << endl;
      fileOUT << -x << " " << -y << " " << -z << endl;
      fileOUT << x << " " << -y << " " << -z << endl;      
      fileOUT << x << " " << -y << " " << z << endl;      
      fileOUT << -x << " " << -y << " " << z << endl;      
      fileOUT << -x << " " << y << " " << -z << endl;
      fileOUT << x << " " << y << " " << -z << endl;      
      fileOUT << x << " " << y << " " << z << endl;      
      fileOUT << -x << " " << y << " " << z << endl;
      fileOUT << endl << "6" << endl;
      fileOUT << "4 5 6 7" << endl;
      fileOUT << "5 1 2 6" << endl;              
      fileOUT << "1 0 3 2" << endl;      
      fileOUT << "0 4 7 3" << endl;      
      fileOUT << "4 0 1 5" << endl;      
      fileOUT << "3 7 6 2" << endl;      
      
      fileOUT.close(); 
      
      cout << "   Volume-equivalent sphere radius = " << 
      	pow( 3. * length * width * height / ( 4. * pi ) , 1./3. ) << endl;
      break;
      
    case 3:
      good_choice = false; 
      while(!good_choice)
      { 
        cout << "Volume-equivalent sphere radius = ";
        cin >> sphere_radius;
  
        if ( sphere_radius < 0. )
          cout << "Wrong choice, radius > 0, enter again" << endl;
        else
          good_choice = true;
      } 
      edge = sphere_radius * pow( 4. * pi / 3., 1./3. );      

      cout << "Cube" << endl;
      cout << "   edge length = " << edge << endl;

      filename="cube.insert";
      
      x=edge/2.;
      fileOUT.open(filename.c_str(),ios::out);
      
      fileOUT << "3" << endl << "8" << endl;
      fileOUT << -x << " " << -x << " " << -x << endl;
      fileOUT << x << " " << -x << " " << -x << endl;      
      fileOUT << x << " " << -x << " " << x << endl;      
      fileOUT << -x << " " << -x << " " << x << endl;      
      fileOUT << -x << " " << x << " " << -x << endl;
      fileOUT << x << " " << x << " " << -x << endl;      
      fileOUT << x << " " << x << " " << x << endl;      
      fileOUT << -x << " " << x << " " << x << endl;
      fileOUT << endl << "6" << endl;
      fileOUT << "4 5 6 7" << endl;
      fileOUT << "5 1 2 6" << endl;              
      fileOUT << "1 0 3 2" << endl;      
      fileOUT << "0 4 7 3" << endl;      
      fileOUT << "4 0 1 5" << endl;      
      fileOUT << "3 7 6 2" << endl;      
      
      fileOUT.close();  
      break;
                  
    case 4:
      good_choice = false; 
      while(!good_choice)
      { 
        cout << "Volume-equivalent sphere radius = ";
        cin >> sphere_radius;
  
        if ( sphere_radius < 0. )
          cout << "Wrong choice, radius > 0, enter again" << endl;
        else
          good_choice = true;
      } 
      edge = sphere_radius * pow( 8. * sqrt(2.) * pi, 1./3. );      

      cout << "Tetrahedron" << endl;
      cout << "   edge length = " << edge  << endl;

      filename="regulartetrahedron.insert";  
      
      theta = 2. * atan( sqrt(2.) ) - pi / 2.;
      rcirc = sqrt( 3. / 8. ) * edge;
      h = rcirc * sin (theta);
      l = rcirc * cos (theta);      

      fileOUT.open(filename.c_str(),ios::out);
      
      fileOUT << "3" << endl << "4" << endl;
      fileOUT << l << " " << "0." << " " << -h << endl;
      fileOUT << l*cos(2.*pi/3.) << " " << l*sin(2.*pi/3.) << " " << -h << endl;
      fileOUT << l*cos(4.*pi/3.) << " " << l*sin(4.*pi/3.) << " " << -h << endl;
      fileOUT << "0. 0. " << rcirc << endl;      
      fileOUT << endl << "4" << endl;
      fileOUT << "0 1 2" << endl;
      fileOUT << "0 3 1" << endl;              
      fileOUT << "2 1 3" << endl;      
      fileOUT << "2 3 0" << endl;      
      
      fileOUT.close();          
      break;

    case 5:
      good_choice = false; 
      while(!good_choice)
      { 
        cout << "Volume-equivalent sphere radius = ";
        cin >> sphere_radius;
  
        if ( sphere_radius < 0. )
          cout << "Wrong choice, radius > 0, enter again" << endl;
        else
          good_choice = true;
      } 
      edge = pow( 4. * pi / sqrt(2.) , 0.3333333333 ) * sphere_radius ;
      coord = edge / sqrt(2.);

      cout << "Octahedron" << endl;
      cout << "   edge length = " << edge  << endl;
      
      filename="regularoctahedron.insert";
      
      fileOUT.open(filename.c_str(),ios::out);
      fileOUT.setf(std::ios::scientific,std::ios::floatfield);
      fileOUT.precision(8);
      
      fileOUT << "3" << endl << "6" << endl;
      fileOUT << sqrt2o2 * coord << " " << sqrt2o2 * coord << " 0." << endl;
      fileOUT << - sqrt2o2 * coord << " " << - sqrt2o2 * coord << " 0." << endl;
      fileOUT << - sqrt2o2 * coord << " " << sqrt2o2 * coord << " 0." << endl;
      fileOUT << sqrt2o2 * coord << " " << - sqrt2o2 * coord << " 0." << endl;
      fileOUT << "0. 0. " << coord << endl;
      fileOUT << "0. 0. " << - coord << endl;  
      fileOUT << endl << "8" << endl; 
      fileOUT << "4 1 3" << endl;
      fileOUT << "4 2 1" << endl;      
      fileOUT << "4 0 2" << endl;      
      fileOUT << "4 3 0" << endl;
      fileOUT << "5 3 1" << endl;
      fileOUT << "5 1 2" << endl;
      fileOUT << "5 2 0" << endl;      
      fileOUT << "5 0 3" << endl;                                   
                  
      fileOUT.close();        
      
      break;
      
    case 6:
      good_choice = false; 
      while(!good_choice)
      { 
        cout << "Volume-equivalent sphere radius = ";
        cin >> sphere_radius;
  
        if ( sphere_radius < 0. )
          cout << "Wrong choice, radius > 0, enter again" << endl;
        else
          good_choice = true;
      } 
      vol = ( 4. / 3. ) * pi * pow( sphere_radius, 3. );
      edge = pow( 2. * vol / ( sqrt( 5. ) * pow( gold, 4. ) ), 
      	0.3333333333 ) ;
      zzz = 0.5 * gold * edge ;

      cout << "Dodecahedron" << endl;
      cout << "   edge length = " << edge  << endl;
      
      filename="regulardodecahedron.insert";
      
      fileOUT.open(filename.c_str(),ios::out);
      fileOUT.setf(std::ios::scientific,std::ios::floatfield);
      fileOUT.precision(8);      
      fileOUT << "3" << endl << "20" << endl;
      fileOUT << zzz << " " << - zzz << " " << zzz << endl;
      fileOUT << zzz << " " << zzz << " " << zzz << endl;      
      fileOUT << - zzz << " " << zzz << " " << zzz << endl;
      fileOUT << - zzz << " " << - zzz << " " << zzz << endl;        
      fileOUT << zzz << " " << - zzz << " " << - zzz << endl;
      fileOUT << zzz << " " << zzz << " " << - zzz << endl;      
      fileOUT << - zzz << " " << zzz << " " << - zzz << endl;
      fileOUT << - zzz << " " << - zzz << " " << - zzz << endl;      
      fileOUT << "0. " << - zzz * gold << " " << zzz / gold << endl;
      fileOUT << "0. " << - zzz * gold << " " << - zzz / gold << endl;      
      fileOUT << "0. " << zzz * gold << " " << - zzz / gold << endl;      
      fileOUT << "0. " << zzz * gold << " " << zzz / gold << endl;      
      fileOUT << zzz / gold << " 0. " << " " << zzz * gold << endl;      
      fileOUT << zzz / gold << " 0. " << " " << - zzz * gold << endl;
      fileOUT << - zzz / gold << " 0. " << " " << - zzz * gold << endl;   
      fileOUT << - zzz / gold << " 0. " << " " << zzz * gold << endl;  
      fileOUT << zzz * gold << " " << zzz / gold << " 0." << endl;
      fileOUT << zzz * gold << " " << - zzz / gold << " 0." << endl;
      fileOUT << - zzz * gold << " " << - zzz / gold << " 0." << endl;
      fileOUT << - zzz * gold << " " << zzz / gold << " 0." << endl;
      fileOUT << endl << "12" << endl; 
      fileOUT << "15 2 19 18 3" << endl;
      fileOUT << "15 12 1 11 2" << endl;      
      fileOUT << "12 0 17 16 1" << endl;      
      fileOUT << "0 12 15 3 8" << endl;
      fileOUT << "2 11 10 6 19" << endl;
      fileOUT << "1 16 5 10 11" << endl;
      fileOUT << "18 19 6 14 7" << endl;      
      fileOUT << "10 5 13 14 6" << endl;      
      fileOUT << "17 4 13 5 16" << endl;
      fileOUT << "9 7 14 13 4" << endl;
      fileOUT << "0 8 9 4 17" << endl;
      fileOUT << "3 18 7 9 8" << endl;                                    
                  
      fileOUT.close();        
      
      break; 
      
    case 7:
      good_choice = false; 
      while(!good_choice)
      { 
        cout << "Volume-equivalent sphere radius = ";
        cin >> sphere_radius;
  
        if ( sphere_radius < 0. )
          cout << "Wrong choice, radius > 0, enter again" << endl;
        else
          good_choice = true;
      } 
      edge = pow( 48. * pi / ( 15. * ( 3. + sqrt(5.) ) ) , 0.3333333333 ) 
      		* sphere_radius ;
      half_edge = 0.5 * edge ;

      cout << "Icosahedron" << endl;
      cout << "   edge length = " << edge  << endl;
      
      filename="regularicosahedron.insert";
      
      fileOUT.open(filename.c_str(),ios::out);
      fileOUT.setf(std::ios::scientific,std::ios::floatfield);
      fileOUT.precision(8);
      
      fileOUT << "3" << endl << "12" << endl;
      fileOUT << - gold * half_edge << " " << - half_edge << " 0." << endl;
      fileOUT << gold * half_edge << " " << - half_edge << " 0." << endl;
      fileOUT << gold * half_edge << " " << half_edge << " 0." << endl;	
      fileOUT << - gold * half_edge << " " << half_edge << " 0." << endl;
      fileOUT << "0. " << - gold * half_edge << " " << - half_edge << endl;
      fileOUT << "0. " << gold * half_edge << " " << - half_edge << endl;
      fileOUT << "0. "  << gold * half_edge << " " << half_edge << endl;
      fileOUT << "0. "  << - gold * half_edge << " " << half_edge << endl;
      fileOUT << - half_edge << " 0. " << - gold * half_edge << endl;
      fileOUT << - half_edge << " 0. " << gold * half_edge << endl;
      fileOUT << half_edge << " 0. " << gold * half_edge << endl;	
      fileOUT << half_edge << " 0. " << - gold * half_edge << endl;
      fileOUT << endl << "20" << endl; 
      fileOUT << "10 2 1" << endl;
      fileOUT << "2 11 1" << endl;      
      fileOUT << "3 9 0" << endl;      
      fileOUT << "8 3 0" << endl;
      fileOUT << "7 1 4" << endl;
      fileOUT << "0 7 4" << endl;
      fileOUT << "2 6 5" << endl;      
      fileOUT << "6 3 5" << endl;      
      fileOUT << "6 10 9" << endl;
      fileOUT << "10 7 9" << endl;
      fileOUT << "4 11 8" << endl;
      fileOUT << "11 5 8" << endl;      
      fileOUT << "6 2 10" << endl;      
      fileOUT << "3 6 9" << endl;
      fileOUT << "0 9 7" << endl;
      fileOUT << "1 7 10" << endl;
      fileOUT << "8 0 4" << endl;      
      fileOUT << "11 4 1" << endl;      
      fileOUT << "8 5 3" << endl;
      fileOUT << "11 2 5" << endl;                              
                  
      fileOUT.close();        
      
      break;                                

  }  
        
  return(0);
  
}
