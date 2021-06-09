#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <cmath>
#include <cstdlib>
using namespace std;

double random_value(const double &min,const double &max)
{  
  return(min+(double(rand())/double(RAND_MAX))*(max-min));  
}

int main(int argc, char *argv[])
{
  string tline;
  unsigned int nx = 0, ny = 1, nz = 1, nxmax, nymax, nzmax = 1,
  	n, ncalc, nmax, i, j, k, ipart, buf ;
  double lx, ly, lz = 1., radius, ry, rz, hx, hy, hz, posx, posy, posz,
  	origx, origy, origz, deltax, deltay, deltaz, 
	magnitude = 1.1, particle_measure = 0., equivol_radius = 0. ;
  istringstream iss;
  bool is_3D = false, random = false ;
  string output_file, random_keyword, filling_mode ;
  size_t ndigits = 16;

    
  ifstream file_IN(argv[1],ios::in);
  getline(file_IN,tline,'\n');  
  getline(file_IN,tline,'\n');
  iss.str(tline);
  iss >> origx >> origy;
  if (iss >> origz)
  {
    cout << "3D geometry" << endl;
    is_3D=true;
  }
  else
    cout << "2D geometry" << endl;    
  iss.clear();

  getline(file_IN,tline,'\n');  
  getline(file_IN,tline,'\n');
  iss.str(tline);
  iss >> lx >> ly;
  if (is_3D) iss >> lz;   
  iss.clear();

  cout << "Origin = " << origx << " " << origy;
  if (is_3D) cout << " " << origz;
  cout << endl;  
  cout << "Box geometry = " << lx << " " << ly;
  if (is_3D) cout << " " << lz;
  cout << endl;
  
  getline(file_IN,tline,'\n');  
  getline(file_IN,tline,'\n');  
  iss.str(tline);
  iss >> radius >> equivol_radius;
  iss.clear();  

  getline(file_IN,tline,'\n');  
  getline(file_IN,tline,'\n');  
  iss.str(tline);
  iss >> filling_mode >> buf;
  if ( filling_mode != "STRUCTURED" ) filling_mode = "TOTAL";
  if ( filling_mode == "TOTAL" ) n = buf ;
  else 
  {
    nx = buf;
    iss >> ny;
    if (is_3D) iss >> nz; 
  }  
  iss.clear(); 

  getline(file_IN,tline,'\n');  
  getline(file_IN,tline,'\n');  
  iss.str(tline);
  iss >> random_keyword;
  if ( random_keyword == "true" || random_keyword == "truetotal" )
  { 
    random = true;
    if ( random_keyword == "true" ) srand(0);
    else srand( time(NULL) );
  }   
  iss.clear();     
  
  getline(file_IN,tline,'\n');  
  getline(file_IN,tline,'\n');  
  iss.str(tline);
  iss >> output_file;
  iss.clear();
  
  getline(file_IN,tline,'\n'); 
  if ( !file_IN.eof() )
  {
    getline(file_IN,tline,'\n');  
    iss.str(tline);
    iss >> ndigits;
  }   
      
  file_IN.close(); 
  
  // Max number of particles to be inserted
  nxmax=int((1.+1e-12)*lx/(2.*radius));
  nymax=int((1.+1e-12)*ly/(2.*radius));  
  if (is_3D) nzmax=int((1.+1e-12)*lz/(2.*radius));   
  cout << "Maximum number of particles to be inserted as an array = ";
  cout << nxmax << " x " << nymax;
  if (is_3D) cout << " x " << nzmax;
  nmax=nxmax*nymax;
  if (is_3D) nmax*=nzmax;
  cout << " = " << nmax << endl; 
  cout << "Number of digits to write positions = " << ndigits << endl;

  if ( filling_mode == "TOTAL" )
  {  
    // Number of particles to be inserted
    ry=ly/lx;
    if (is_3D) rz=lz/lx;  
    cout << "Number of particles to be inserted as an array = " << n << endl;
    ncalc=0;  
    while (ncalc<n)
    {
      nx++;
      ny=int(nx*ry);
      if (!ny) ny=1;
      ncalc=nx*ny;
      if (is_3D) 
      {
        nz=int(nx*rz);
        if (!nz) nz=1;
        ncalc*=nz;
      }
      cout << nx << " " << ny;
      if (is_3D) cout << " " << nz;
      cout << endl;
    } 

    if (is_3D) 
    {
      // Correction for homogeneity
      unsigned int nplan=nx*ny;
      if ( double(n)/double(nplan)-int(n/nplan) !=0. )
        nz=int(n/nplan)+1;
    }  

    // Check that ncalc < nmax
    // If yes, write position in an output file
    if ( ncalc>nmax && n<nmax )
    {
      nx=nxmax;
      ny=nymax;      
      nz=nzmax;
      ncalc = nmax;
    }

    cout << "Box meshing = " << nx << " x " << ny;
    if (is_3D) cout << " x " << nz;
    cout << endl; 
     
    // Check that ncalc < nmax
    // If yes, write position in an output file
    if (ncalc<=nmax)
    {
      hx=lx/nx;
      hy=ly/ny;
      if (random) deltax=hx/2.-magnitude*radius;
      else deltax=hx/1e8;
      if (random) deltay=hy/2.-magnitude*radius; 
      else deltay=hy/1e8;       
      if (is_3D) 
      {
        hz=lz/nz;
        if (random) deltaz=hz/2.-magnitude*radius;
        else deltaz=hz/1e8; 
      } 
    
      ofstream file_OUT(output_file.c_str(),ios::out);
      file_OUT.setf(ios::scientific,ios::floatfield);
      file_OUT.precision(ndigits);
      ipart=0;
      if (is_3D) 
      {
        k=0;
        while((k<nz)&&(ipart<n))      
        {
	  i=0;
	  while((i<nx)&&(ipart<n))
          {
            j=0;
	    while((j<ny)&&(ipart<n))      
            { 
              posx=origx+(i+0.5)*hx;
	      posx+=random_value(-deltax,deltax); 
	      posy=origy+(j+0.5)*hy;
	      posy+=random_value(-deltay,deltay); 
	      posz=origz+(k+0.5)*hz;
	      posz+=random_value(-deltaz,deltaz);	    
	      file_OUT << posx << "\t" << posy << "\t" << posz << endl;
	      ++j;
	      ++ipart;
	    }
	    ++i;	  
	  }
	  ++k;
        }   
      }
      else
      {
        j=0;
        while((j<ny)&&(ipart<n))      
        {         
	  i=0;
	  while((i<nx)&&(ipart<n))
          {
            posx=origx+(i+0.5)*hx;
	    posx+=random_value(-deltax,deltax); 
	    posy=origy+(j+0.5)*hy;
	    posy+=random_value(-deltay,deltay); 
	    file_OUT << posx << "\t" << posy << "\t 0." << endl;
	    ++ipart;
	    ++i;
	  }
	  ++j;
        }    
      }
      file_OUT.close();    
    }
    else
    {
      cout << "!!! Number of particles to be inserted as a structured array"
    	<< " is too high !!!" << endl;
    }
  }
  else
  {
    bool goon = true ;
    if ( nx > nxmax )
    {
      cout << "!!! Number of particles to be inserted in the x direction"
    	<< " is too high !!!" << endl; 
      goon = false ;
    }
    else if ( ny > nymax )
    {
      cout << "!!! Number of particles to be inserted in the y direction"
    	<< " is too high !!!" << endl;
      goon = false ;
    }	
    else if ( nz > nzmax )
    {
      cout << "!!! Number of particles to be inserted in the z direction"
    	<< " is too high !!!" << endl;
      goon = false ;
    }
    
    if ( !goon ) 
    {
      ofstream file_OUT(output_file.c_str(),ios::out);
      file_OUT.close();          
    }	
    else
    {
      magnitude = 1.00000001;
      hx = lx / nx;
      hy = ly / ny;
      if ( random ) 
      { 
        deltax = hx / 2. - magnitude * radius > 0. ? 
		hx / 2. - magnitude * radius : 0. ;
        deltay = hy / 2. - magnitude * radius > 0. ? 
		hy / 2. - magnitude * radius : 0. ; 
      }
      else deltax = deltay = 0. ;      
      if ( is_3D ) 
      {
        hz = lz / nz;
        if ( random ) deltaz = hz / 2. - magnitude * radius > 0. ? 
		hz / 2. - magnitude * radius : 0. ;
	particle_measure = ( 4. /3. ) * acos( -1. ) 
		* pow( equivol_radius, 3. ) ;	
      }
      else particle_measure = acos( -1. ) * pow( equivol_radius, 2. ) ;
      
      ofstream file_OUT(output_file.c_str(),ios::out);
      file_OUT.setf(ios::scientific,ios::floatfield);
      file_OUT.precision(ndigits);
      if ( is_3D )
        for (i=0;i<nx;++i)
          for (j=0;j<ny;++j)
	    for (k=0;k<nz;++k)
	    {
              posx = origx + ( i + 0.5 ) * hx + random_value( -deltax, deltax );
	      posy = origy + ( j + 0.5 ) * hy + random_value( -deltay, deltay );
	      posz = origz + ( k + 0.5 ) * hz + random_value( -deltaz, deltaz );
	      file_OUT << posx << "\t" << posy << "\t" << posz << endl;	    
	    }
      else
        for (i=0;i<nx;++i)
          for (j=0;j<ny;++j)
	  {
            posx = origx + ( i + 0.5 ) * hx + random_value( -deltax, deltax );
	    posy = origy + ( j + 0.5 ) * hy + random_value( -deltay, deltay );
	    file_OUT << posx << "\t" << posy << "\t 0." << endl;
	  }              
      
      file_OUT.close();
      n = nx * ny * nz ;
      cout << "Total number of particles inserted = " << n << endl;
      cout << "Solid fraction = " << n * particle_measure / ( lx * ly * lz )
      	<< endl;
    }          
  }	
        
  return(0);
  
}
