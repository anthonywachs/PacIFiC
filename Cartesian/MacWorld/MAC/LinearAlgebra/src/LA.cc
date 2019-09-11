#include <LA.hh>

#include <string>
#include <iostream>

//----------------------------------------------------------------------
std::ostream& operator<<( std::ostream& out, LA::SyncState state )
//----------------------------------------------------------------------
{
   if( state == LA::Sync )
   {
      out << "Sync" ;
   }
   else if( state == LA::NotSync_add )
   {
      out << "NotSync_add" ;
   }
   else if( state == LA::NotSync_set )
   {
      out << "NotSync_set" ;
   }
   else if( state == LA::NotSync_undef )
   {
      out << "NotSync_undef" ;
   }
   return( out ) ;
}

//----------------------------------------------------------------------
std::ostream& operator<<( std::ostream& out, LA::DistributionStrategy dist )
//----------------------------------------------------------------------
{
   if( dist == LA::FromLocalSize )
   {
      out << "from_local_size" ;
   }
   else if( dist == LA::FromGlobalSize )
   {
      out << "from_global_size" ;
   }
   else if( dist == LA::NoDistribution )
   {
      out << "no_distribution" ;
   }
   else if( dist == LA::InvalidDistribution )
   {
      out << "invalid_distribution" ;
   }
   return( out ) ;
}
