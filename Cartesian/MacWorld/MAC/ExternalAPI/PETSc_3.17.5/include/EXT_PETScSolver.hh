#ifndef EXT_PETSc_SOLVER_HH
#define EXT_PETSc_SOLVER_HH

#include <LA_Solver.hh>

#include <EXT_PETScAPI.hh>
#include <EXT_PETScMatrix.hh>


/*
Wrappers around PETSc iterative solvers.
*/

class EXT_PETScSolver : public LA_Solver
{
   public: //-----------------------------------------------------------------
      
   //-- Instance delivery and initialization

      virtual EXT_PETScSolver* create_clone( MAC_Object* a_owner ) const ;
      
   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

     ~EXT_PETScSolver( void ) ;
      EXT_PETScSolver( EXT_PETScSolver const& other ) ;
      EXT_PETScSolver& operator=( EXT_PETScSolver const& other ) ;

      EXT_PETScSolver( MAC_Object* a_owner, MAC_ModuleExplorer const* exp ) ;
      EXT_PETScSolver( MAC_Object* a_owner, EXT_PETScSolver const* other ) ;

   //-- Plug in

      EXT_PETScSolver( void ) ;

      virtual EXT_PETScSolver* create_replica( 
                                       MAC_Object* a_owner,
                                       MAC_ModuleExplorer const* exp ) const ;

   //-- Linear system resolution
      
      virtual void set_matrix_self( LA_Matrix const* mat, bool &ok, bool same_pattern ) ;

      virtual void solve_self( LA_Vector const* b, LA_Vector* x,
                               size_t &nb_iter, bool &ok ) ;
      
      virtual void unset_matrix_self( void ) ;

      void build_ksp( KSP &ksp, MAC_ModuleExplorer const* A_EXP ) ;
      
      void build_pc( KSP & ksp, MAC_ModuleExplorer const* A_EXP ) ;
      
      std::string mat_ordering( MAC_ModuleExplorer const* exp,
                                std::string const& name ) const ;
      
   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   //-- Class attributes

      static EXT_PETScSolver const* PROTOTYPE ;
      
    //-- Attributes

      MAC_ModuleExplorer const* const EXP ;
      MAC_ModuleExplorer* SUBEXP ;
      KSP MY_KSP ;
      bool HAS_TO_DESTROY_KSP ;
      EXT_PETScMatrix* MATRIX ;
      bool VERB ;
} ;

#endif
