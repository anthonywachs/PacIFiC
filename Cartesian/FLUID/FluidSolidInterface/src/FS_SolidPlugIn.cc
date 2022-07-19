#include <FS_SolidPlugIn.hh>


//---------------------------------------------------------------------------
FS_SolidPlugIn:: FS_SolidPlugIn()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_SolidPlugIn:: FS_SolidPlugIn" ) ;

  m_insertion_file = "undefined";
  m_simulation_file = "undefined";
}




//---------------------------------------------------------------------------
FS_SolidPlugIn::FS_SolidPlugIn( string const& insertion_file_,
        string const& simulation_file_ ) 
//---------------------------------------------------------------------------
  : m_insertion_file( insertion_file_ )
  , m_simulation_file( simulation_file_ )
{
  MAC_LABEL( "FS_SolidPlugIn:: FS_SolidPlugIn(x,x,x,x)" ) ;

}




//---------------------------------------------------------------------------
FS_SolidPlugIn:: ~FS_SolidPlugIn()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_SolidPlugIn:: ~FS_SolidPlugIn" ) ;

}
