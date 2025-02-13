#include "AllInsertionWindows.hh"


// ----------------------------------------------------------------------------
// Default constructor
AllInsertionWindows::AllInsertionWindows()
{}




// ----------------------------------------------------------------------------
// Destructor
AllInsertionWindows::~AllInsertionWindows()
{
  m_insertion_windows.clear();
  // Note: the loadings are destroyed by the destructors of the class
  // ObstacleKinematicsVelocity 
  // Hence we free the lists here but do not destroy the pointed objects
  m_AllImposedVelocitiesOnWindows.clear();
}




// ----------------------------------------------------------------------------
// Adds a new insertion window
void AllInsertionWindows::addWindow( Window const& iwindow )
{
  m_insertion_windows.insert( m_insertion_windows.begin(), iwindow );
} 




// ----------------------------------------------------------------------------
// Converts the vector of insertion windows to a list of rigid bodies
void AllInsertionWindows::convertToRigidBodies( list<RigidBody*>& iwlist ) const
{  
  for (vector<Window>::const_iterator iv=m_insertion_windows.cbegin();
  	iv!=m_insertion_windows.cend();iv++)
    iv->addAsRigidBody( iwlist );
} 




// ----------------------------------------------------------------------------
// Returns a pointer to a randomly selected window
Window const* AllInsertionWindows::getRandomWindow() const
{
  size_t nWindow = 0, nbreWindows = m_insertion_windows.size();

  // Random selection of an insertion window
  if ( nbreWindows != 1 )
  {
    double n = double(random()) / double(INT_MAX);
    nWindow = size_t( n * double(nbreWindows) );
    if ( nWindow == nbreWindows ) nWindow--;
  }
  
  return ( &(m_insertion_windows[nWindow]) );
}




// ----------------------------------------------------------------------------
// Associates the imposed velocity to an insertion window
void AllInsertionWindows::LinkImposedMotion( ObstacleImposedVelocity* impvel )
{
  bool status = false;
  for (vector<Window>::iterator iv=m_insertion_windows.begin();
  	iv!=m_insertion_windows.end() && !status;iv++)
    status = iv->LinkImposedMotion( impvel );
  m_AllImposedVelocitiesOnWindows.push_back( impvel );    
}




// ----------------------------------------------------------------------------
// Moves the insertion windows if any has an imposed motion
void AllInsertionWindows::Move( double time, double dt )
{
  static bool anyactive_previousdt = false;

  if ( !m_AllImposedVelocitiesOnWindows.empty() )
  {
    // Check if any imposed velocity/ is active over this time interval
    bool anyactive = false;
    double subinterval = 0.;  
    list<ObstacleImposedVelocity*>::iterator il;
    for (il=m_AllImposedVelocitiesOnWindows.begin();
	il!=m_AllImposedVelocitiesOnWindows.end();il++)
      if ( (*il)->isActif( time - dt, time, dt, subinterval ) )
        anyactive = true;

    // Move insertion windows
    if ( anyactive || anyactive_previousdt )
      for (vector<Window>::iterator iv=m_insertion_windows.begin();
  	iv!=m_insertion_windows.end();iv++)
        iv->resetKinematics(); 
    anyactive_previousdt = anyactive;
    if ( anyactive )
      for (vector<Window>::iterator iv=m_insertion_windows.begin();
  	iv!=m_insertion_windows.end();iv++)
        iv->Move( time, dt );    
      	
    // Update imposed velocity kinematics
    for (il=m_AllImposedVelocitiesOnWindows.begin();
	il!=m_AllImposedVelocitiesOnWindows.end(); )
      if ( (*il)->isCompleted( time, dt ) )
        il = m_AllImposedVelocitiesOnWindows.erase( il );
      else il++;	
  }
} 
