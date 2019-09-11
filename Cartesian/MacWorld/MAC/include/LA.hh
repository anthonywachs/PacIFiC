#ifndef LA_HH
#define LA_HH

#include <iosfwd>

/*
PUBLISHED
*/

class LA
{
   public: //-----------------------------------------------------------

      /*
      Synchronization state for distributed matrices and vectors
         Sync          : is synchronized
         NotSync_add   : is not synchronized because of a previous
                         call to "add" methods
         NotSync_set   : is not synchronized because of a previous
                         call to "set" methods
         NotSync_undef : is in an undefined not synchronized state
                         (eg just after the creational process) 
      */
      enum SyncState
      {
         Sync,
         NotSync_add,
         NotSync_set,
         NotSync_undef
      } ;

      /*
      Strategy used for the distribution of rows (or columns)
      of matrices and vectors
         NoDistribution : all rows (or columns) are handled
                          by the current process
         FromLocalSize  : the distribution over all processes is set
                          from the number of local rows (or columns),
                          i.e. handled by each process
         FromGlobalSize : the global number of rows (or columns) is
                          equi-distributed over all processes
      */
      enum DistributionStrategy
      {
         NoDistribution,
         FromLocalSize,
         FromGlobalSize,
         InvalidDistribution
      } ;

   //-- Input - Output

      friend std::ostream& operator<<(
                              std::ostream& out, LA::SyncState state ) ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      LA( void ) ;
     ~LA( void ) ;
      LA( LA const& other ) ;
      LA& operator=( LA const& other ) ;
} ;

#endif



