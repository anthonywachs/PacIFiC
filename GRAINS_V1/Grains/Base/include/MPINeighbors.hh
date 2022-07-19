#ifndef _MPINEIGHBORS_HH_
#define _MPINEIGHBORS_HH_

#include <mpi.h>
#include <Cell.hh>
#include <vector>
#include <iostream>
#include <list>
using std::list;
using std::vector;
using std::ostream;
using std::endl;


struct MPICartInfos
{
  int rank; /**< process rank in the MPI communicator */
  int *coords; /**< process coordinates in the Cartesian MPI topology */
};


/** @brief The class MPINeighbors.

    Manages neighbor processes in a Cartesian MPI topology. 
    
    @author A.WACHS - Institut Francais du Petrole - 2009 - Creation 
    @author A.WACHS - 2019 - Major cleaning & refactoring */
//=============================================================================
class MPINeighbors
{
  public:
    /**@name Constructors */
    //@{
    /** @brief Constructor with input parameters 
    @param commgrainsMPI_3D Cartesian MPI communicator
    @param center_coords MPI coordinates of this process in the Cartesian MPI 
    	topology
    @param nprocsdir numbers of processes in each direction
    @param period periodicity of the Cartesian MPI topology in each 
    	direction */
    MPINeighbors( MPI_Comm& commgrainsMPI_3D, int const* center_coords,
  	int const* nprocsdir, int const* period );

    /** @brief Destructor */
    ~MPINeighbors();
    //@}

  
    /**@name Access */
    //@{
    /** @brief Returns the MPI rank of a neighbor based on a relative 
    position (-1,0 ou +1) 
    @param i relative position in the x direction (-1,0 ou 1) 
    @param j relative position in the y direction (-1,0 ou 1)   
    @param k relative position in the z direction (-1,0 ou 1) */
    int rank( int i, int j, int k ) const;
  
    /** @brief Returns the coordinates in the Cartesian MPI topology of a 
    neighbor based on a relative position (-1,0 ou +1) 
    @param i relative position in the x direction (-1,0 ou 1) 
    @param j relative position in the y direction (-1,0 ou 1)   
    @param k relative position in the z direction (-1,0 ou 1) */
    int const* coordinates( int i, int j, int k ) const;
  
    /** @brief Returns the number of neighbors including the process itself */
    int nbMPINeighbors() const;
  
    /** @brief Returns the MPI rank of all neighbors (including the process
    itself) in a array of integers */
    int* rankMPINeighbors() const; 
  
    /** @brief Returns the MPI rank of all neighbors (excluding the process
    itself) in a list of integers */
    list<int> const* rankMPINeighborsOnly() const;
  
    /** @brief Returns the MPI geolocalisation of all neighbors (excluding the 
    process itself) in a list of GeoPosition */
    list<GeoPosition> const* geolocMPINeighborsOnly() const;     
    //@} 

  
    /**@name I/O methods */
    //@{
    /** @brief Output operator 
    @param os output stream
    @param M the MPINeighbors object */
    friend ostream& operator <<( ostream& os, MPINeighbors const& M );
    //@}    


  private:
    /** @name Parameters */
    //@{ 
    vector<struct MPICartInfos> m_data; /**< rank and coordinates of processes
    	in the Cartesian MPI topology */
    int m_nneighbors; /**< number of MPI neighbors */
    list<int> m_rankNeighborsOnly; /**< MPI rank of neighbors (excluding the
    	process itself) */
    list<GeoPosition> m_geolocNeighborsOnly; /**< MPI geolocalisation of
    	neighbors (excluding the process itself) */	
    //@}
    

    /**@name Constructors */
    //@{
    /** @brief Default constructor (forbidden) */
    MPINeighbors();
    //@}
};

#endif
