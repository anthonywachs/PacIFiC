#include <EXT_PETScSolver.hh>
#include <EXT_PETScImplementation.hh>

#include <MAC_assertions.hh>
#include <MAC_Exec.hh>
#include <MAC_Communicator.hh>
#include <MAC_Error.hh>
#include <MAC_ModuleExplorer.hh>
#include <stringVector.hh>

#include <LA_SeqMatrix.hh>
#include <LA_SeqVector.hh>

#include <EXT_PETScMatrix.hh>

#include <iostream>
#include <sstream>

using std::endl ;
using std::string ;

EXT_PETScSolver const* EXT_PETScSolver::PROTOTYPE = new EXT_PETScSolver() ;


//----------------------------------------------------------------------------
EXT_PETScSolver:: EXT_PETScSolver( void ) 
//----------------------------------------------------------------------------
   : LA_Solver( "EXT_PETScSolver" )
   , EXP( 0) 
   , SUBEXP( 0 )
   , MY_KSP( 0 )
   , HAS_TO_DESTROY_KSP( true )
   , MATRIX( 0 )
   , VERB( false )
{   
}




//----------------------------------------------------------------------------
EXT_PETScSolver:: ~EXT_PETScSolver( void )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScSolver:: ~EXT_PETScSolver" ) ;
   if( HAS_TO_DESTROY_KSP && MY_KSP!=0  ) PETSc_do( KSPDestroy( &MY_KSP ) ) ;
   MY_KSP = 0 ;  
   if( is_a_prototype() )
   {
      PROTOTYPE = 0 ;
   }
}




//----------------------------------------------------------------------------
EXT_PETScSolver*
EXT_PETScSolver:: create_replica( MAC_Object* a_owner, 
                                  MAC_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScSolver:: create_replica" ) ;
   MAC_CHECK( create_replica_PRE(  a_owner, exp ) ) ;
   
   EXT_PETScSolver* result = new EXT_PETScSolver( a_owner, exp ) ;

   MAC_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------------
EXT_PETScSolver:: EXT_PETScSolver( MAC_Object* a_owner,
                                   MAC_ModuleExplorer const* exp ) 
//----------------------------------------------------------------------------
   : LA_Solver( a_owner )
   , EXP( exp->create_clone( this ) )
   , SUBEXP( 0 )
   , MY_KSP( 0 )
   , HAS_TO_DESTROY_KSP( true )
   , MATRIX( 0 )
   , VERB( false )
{
   MAC_LABEL( "EXT_PETScSolver:: EXT_PETScSolver" ) ;

   if( exp->has_module( "LA_Matrix" ) )
   {
      MAC_ModuleExplorer const* sexp = exp->create_subexplorer( 0, 
      	"LA_Matrix" ) ;
      LA_Matrix* mat = LA_Matrix::make( this, sexp ) ;
      if( mat->implementation() != EXT_PETScImplementation::object() )
      {
         MAC_Error::object()->raise_bad_data_value( sexp, 
	 	"concrete_name", "PETSc kind" ) ;
      }
      MATRIX = static_cast<EXT_PETScMatrix*>( mat ) ;
      
      sexp->destroy() ;
   }

   bool iter = ! ( EXP->has_module( "PETSc_preconditioner" ) &&
                    EXP->string_data( "PETSc_preconditioner/type" )=="lu" ) ;
   set_iterative( iter ) ;
}




//----------------------------------------------------------------------------
EXT_PETScSolver*
EXT_PETScSolver:: create_clone( MAC_Object* a_owner ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScSolver:: create_clone" ) ;
   
   EXT_PETScSolver* result = new EXT_PETScSolver( a_owner, this ) ;

   MAC_CHECK( create_clone_POST( result, a_owner ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------------
EXT_PETScSolver:: EXT_PETScSolver( MAC_Object* a_owner,
                                   EXT_PETScSolver const* other ) 
//----------------------------------------------------------------------------
   : LA_Solver( a_owner, other )
   , EXP( other->EXP->create_clone( this ) )
   , SUBEXP( 0 )
   , MY_KSP( 0 )
   , HAS_TO_DESTROY_KSP( true )
   , MATRIX( 0 )
{
   MAC_LABEL( "EXT_PETScSolver:: EXT_PETScSolver" ) ;

   if( other->MATRIX!=0 ) MATRIX=other->MATRIX->create_matrix( this ) ; 
}




//----------------------------------------------------------------------------
void
EXT_PETScSolver:: set_matrix_self( LA_Matrix const* mat,
                                   bool &ok, bool same_pattern ) 
//----------------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScSolver:: set_matrix_self" ) ;
   MAC_CHECK( set_matrix_self_PRE( mat, same_pattern ) ) ;

   MAC_Communicator const* com = MAC_Exec::communicator() ;
   Mat PETSc_mat = 0 ;
   EXT_PETScMatrix const* pmat = dynamic_cast<EXT_PETScMatrix const*>( mat ) ;
   if( pmat!=0 )
   {
      PETSc_mat = pmat->matrix() ;
   }
   else
   {
      LA_SeqMatrix const* smat = dynamic_cast<LA_SeqMatrix const*>( mat ) ;
      if( MATRIX==0 ) 
      {
         MAC_Error::object()->raise_plain(
            "You should specify a valid PETSc matrix (LA_Matrix module) \n"\
            "to use PETSc algorithm with MAC matrices" ) ;
      }
      MAC_ASSERT( MATRIX!=0 ) ;
      if( smat!=0 )
      {
         MATRIX->set( smat ) ;
         PETSc_mat = MATRIX->matrix() ;
      }
   
      else
      {
         MAC_Error::object()->raise_plain( "Not yet implemented" ) ;
      }
   }
   MAC_ASSERT( PETSc_mat!=0 ) ;

   if( !same_pattern )
   {
      if( HAS_TO_DESTROY_KSP && MY_KSP!=0  ) PETSc_do( KSPDestroy( &MY_KSP ) ) ;
      MY_KSP = 0 ;  
      PETSc_do( KSPCreate( PETSC_COMM_WORLD, &MY_KSP ) ) ;
      build_ksp( MY_KSP, EXP ) ;
      build_pc( MY_KSP, EXP ) ;
   }
   PETSc_do( KSPSetOperators( MY_KSP, PETSc_mat, PETSc_mat,
	( same_pattern ? SAME_NONZERO_PATTERN : DIFFERENT_NONZERO_PATTERN ) ) );

   if( SUBEXP != 0 )
   {         
      const PCType pctype ;
      PC pc ;
      PETSc_do( KSPGetPC( MY_KSP, &pc ) ) ;
      PETSc_do( PCGetType(pc, &pctype) ) ;
      
      MAC_ASSERT( std::string(pctype) == std::string(PCBJACOBI) ) ;

      int nlocal, mlocal ;
      
      PETSc_do( MatGetLocalSize(PETSc_mat,&nlocal,&mlocal) ) ;

      size_t nb = com->nb_ranks() ;
      
      size_t_vector vec( nb ) ;
      if( com->rank() < nb-1 )
      {
         com->send( nb-1, (size_t) nlocal ) ;
         com->receive( nb-1, vec ) ;
      }
      else 
      {
         for( size_t i=0 ; i< nb-1 ; i++ )
            com->receive( i, vec(i) ) ;
         vec( nb-1 ) = (size_t) nlocal ;
         for( size_t i=0 ; i< nb-1 ; i++ )
            com->send( i, vec ) ;
      }
      int * blks = new int [nb ] ;
      for( size_t i=0 ; i< nb ; i++ )
         blks[i] = (int) vec(i) ;
      PETSc_do( PCBJacobiSetTotalBlocks(pc,nb,blks) ) ;
      delete [] blks ;

      PETSc_do( KSPSetUp(MY_KSP) ) ;
      KSP *subksp;
      int first ;
      PCBJacobiGetSubKSP(pc,&nlocal,&first,&subksp);
      for( size_t i=0 ; i<nb ; i++ )
      {
         MAC_ASSERT( nlocal==1 ) ;
         build_ksp( subksp[0], SUBEXP ) ;
         build_pc( subksp[0], SUBEXP ) ;
      }
      SUBEXP = 0 ;
   }
   ok = true ;

   MAC_CHECK_POST( set_matrix_self_POST( mat, ok ) ) ;
}




//----------------------------------------------------------------------------
void
EXT_PETScSolver:: unset_matrix_self( void ) 
//----------------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScSolver:: unset_matrix_self" ) ;
   MAC_CHECK_PRE( unset_matrix_self_PRE() ) ;
}




//----------------------------------------------------------------------------
void
EXT_PETScSolver:: solve_self( LA_Vector const* b, LA_Vector * x,
                              size_t &nb_iter, bool &ok ) 
//----------------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScSolver:: solve_self" ) ;
   MAC_CHECK( solve_self_PRE( b, x ) ) ;
   
   EXT_PETScVector* px = 0 ;
   EXT_PETScVector const* prhs = 0 ;
   EXT_PETScVector* my_prhs = 0 ;
   LA_Scatter* scatter = 0 ;
   LA_SeqVector const* brhs = 0 ;
   LA_SeqVector* bx = 0 ;
   
   prhs = dynamic_cast< EXT_PETScVector const* >( b ) ;
   if( prhs!=0 )
   {
      px = static_cast< EXT_PETScVector* >( x ) ;
      MAC_CHECK( dynamic_cast< EXT_PETScVector const* >( x ) !=0 ) ;
   }
   else
   {
      MAC_ASSERT( MATRIX!=0 ) ;
      prhs = my_prhs = MATRIX->create_vector( 0 ) ;
      MAC_ASSERT( prhs->nb_rows() == MATRIX->nb_rows() ) ;
      
      size_t_vector id( prhs->nb_rows() ) ;
      for( size_t i=0 ; i<prhs->nb_rows() ; i++ ) id(i)=i ;
      
      scatter = prhs->create_scatter( my_prhs, id, id ) ;
      
      brhs = dynamic_cast<LA_SeqVector const* >( b ) ;
      if( brhs!=0 )
      {
         scatter->set( brhs, my_prhs ) ;
         
         bx = static_cast<LA_SeqVector* >( x ) ;
         MAC_CHECK( dynamic_cast<LA_SeqVector const* >( x ) != 0 ) ;
         px = prhs->create_vector( my_prhs ) ;
         if( !zero_initial_guess() )
         {
            scatter->set( bx, px ) ;
         }
      }
      else
      {
         MAC_Error::object()->raise_plain( "Not yet implemented" ) ;
      }
   }
   Vec PETSc_rhs = prhs->vector() ;
   Vec PETSc_x = px->vector() ;
   
   // Initial guess :
   if( !zero_initial_guess() )
   {
      const KSPType type ;
      PETSc_do( KSPGetType( MY_KSP, &type ) ) ;
      
      if( string(type) != string(KSPPREONLY) )
         PETSc_do( KSPSetInitialGuessNonzero( MY_KSP, PETSC_TRUE ) );
   }

   // Solve :
   PETSc_do( KSPSolve( MY_KSP, PETSc_rhs, PETSc_x ) ) ;
   if( VERB ) PETSc_do( KSPView( MY_KSP, PETSC_VIEWER_STDOUT_WORLD ) ) ;

   PetscInt nb ;
   PETSc_do( KSPGetIterationNumber( MY_KSP, &nb ) ) ;
   nb_iter = nb ;
      
   KSPConvergedReason reason ;
   PETSc_do( KSPGetConvergedReason( MY_KSP, &reason ) ) ;

   ok = ( reason>0 ) ;

   if( !ok )
   {
      std::ostringstream mesg ;
      mesg << "*** EXT_PETScSolver:" << endl ;
      mesg << "    the converged reason: " << reason << endl ;
      mesg << "    has been returned by: \"KSPGetConvergedReason\"" << endl ;
      mesg << "    which means failure to converge" << endl ;
      mesg << "    (see PETSc documentation for details)." ;
      MAC_Error::object()->display_info( mesg.str() ) ;
   }
   
   if( scatter!=0 )
   {
      MAC_ASSERT( bx!=0 ) ;
      MAC_ASSERT( my_prhs!=0 ) ;
      scatter->get( px, bx ) ;
      my_prhs->destroy() ;   
   }

   MAC_CHECK_POST( solve_self_POST( b, x, nb_iter, ok ) ) ;
}




//----------------------------------------------------------------------------
void
EXT_PETScSolver:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScSolver:: print" ) ;
   MAC_CHECK_INV( invariant() ) ;

   std::string s( indent_width, ' ' ) ;

   // Solver :
   MAC_ModuleExplorer* exp =
      EXP->create_subexplorer( 0, "PETSc_Krylov_subspace_method" ) ;
   std::string const& algo = exp->string_data( "type" ) ;
   os << s << "linear solver: \"PETSc_" << algo << "\"" << std::endl ;
   if( algo != "preonly" )
   {
      os << s << "   ksp_rtol : "
         << exp->double_data( "ksp_rtol" ) << std::endl ;
      os << s << "   ksp_atol : "
         << exp->double_data( "ksp_atol" ) << std::endl ;
      os << s << "   ksp_max_it : "
         << exp->int_data( "ksp_max_it" ) << std::endl ;
   }
   if( algo == "gmres" )
   {
      os << s << "   ksp_gmres_restart : "
         << exp->int_data( "ksp_gmres_restart" ) << std::endl ;
   }
   else if( algo == "richardson" )
   {
      os << s << "   ksp_richardson_scale : "
         << exp->double_data( "ksp_richardson_scale" ) << std::endl ;
      
   }
   exp->destroy() ; exp = 0 ;
   
   // Preconditioner :
   if( EXP->has_module( "PETSc_preconditioner" ) )
   {
      exp = EXP->create_subexplorer( 0, "PETSc_preconditioner" ) ;
      std::string const& prec = exp->string_data( "type" ) ;
      os << s << "preconditioner: \"PETSc_" << prec << "\"" << std::endl ;
      exp->print( os, indent_width+3 ) ;
      exp->destroy() ; exp = 0 ;
   }


   else
   {
      os << s << "preconditioner: PETSc default" << std::endl ;
   }
}




//----------------------------------------------------------------------------
void
EXT_PETScSolver:: build_ksp( KSP &ksp, MAC_ModuleExplorer const* A_EXP ) 
//----------------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScSolver:: build_ksp" ) ;
   
   MAC_ModuleExplorer* exp =
      A_EXP->create_subexplorer( 0, "PETSc_Krylov_subspace_method" ) ;
   if( SUBEXP==0 )
   {
      VERB = exp->has_entry( "verbose" ) && exp->bool_data( "verbose" ) ;
   }

   bool has_opt = EXT_PETScAPI::parse_options( exp, VERB ) ;
   
   HAS_TO_DESTROY_KSP = 
   	! ( exp->has_entry( "destroy") && ! exp->bool_data( "destroy") ) ;
   
   std::string const& type = exp->string_data( "type" ) ;
   if( type=="richardson" || type=="chebychev" || type=="cg" || type=="gmres"
       || type=="tcqmr" || type=="bcgs" || type=="cgs" || type=="tfqmr"
       || type=="cr" || type=="lsqr" || type=="preonly" || type=="qcg"
       || type=="bicg" || type=="fgmres" || type=="minres"|| type=="symmlq" )
   {
      
      PETSc_do( KSPSetType( ksp, const_cast<KSPType>( type.c_str() ) ) ) ;
      
      if( type != "preonly" )
      {
         double const rtol = exp->double_data( "ksp_rtol" ) ;
         double const atol = exp->has_entry( "ksp_atol" ) ?
	 	exp->double_data( "ksp_atol" ) : 1.e-50 ;
         int const max_it = exp->int_data( "ksp_max_it" ) ;
   
         PETSc_do( KSPSetTolerances( ksp, rtol, atol, PETSC_DEFAULT, max_it ) );
      }
      
      if( type == "gmres" )
      {
         int ierr = KSPGMRESSetRestart( ksp, 
                                        exp->int_data( "ksp_gmres_restart" ) ) ;
         PETSc_do( ierr ) ;
         if( exp->has_entry( "ksp_gmres_orthogonalization" ) )
         {
            std::string const& nn = 
                        exp->string_data( "ksp_gmres_orthogonalization" ) ;
            exp->test_data_in( "ksp_gmres_orthogonalization",
            "ksp_gmres_classicalgramschmidt,ksp_gmres_modifiedgramschmidt") ;
            if( nn == "ksp_gmres_classicalgramschmidt" )
            {
               PETSc_do( KSPGMRESSetOrthogonalization( ksp, 
                         KSPGMRESClassicalGramSchmidtOrthogonalization ) ) ;
            }
            else if( nn == "ksp_gmres_modifiedgramschmidt" )
            {
               PETSc_do( KSPGMRESSetOrthogonalization( ksp, 
                         KSPGMRESModifiedGramSchmidtOrthogonalization ) ) ;
            }
         }
      }
      else if( type == "richardson" )
      {
         PETSc_do( KSPRichardsonSetScale( ksp,
		exp->double_data( "ksp_richardson_scale" ) ) ) ;
      }
      
      if( type == "richardson" || type == "chebychev" || type == "cg" )
      {
         if( exp->has_entry( "ksp_norm_type" ) )
         {
            std::string const& nn = exp->string_data( "ksp_norm_type" ) ;
            exp->test_data_in( "ksp_norm_type", 
                               "preconditioned,unpreconditioned" ) ;
            if( nn == "preconditioned" )
            {
               PETSc_do( KSPSetNormType( ksp, KSP_NORM_PRECONDITIONED ) ) ;
            }
            else if( nn == "unpreconditioned" )
            {
               PETSc_do( KSPSetNormType( ksp, KSP_NORM_UNPRECONDITIONED ) ) ;
            }                    
         }
      }
      
      if( exp->has_entry( "preconditioner_side" ) )
      {
         string nn = exp->string_data( "preconditioner_side" ) ;
         PCSide pcs = PC_LEFT ;
         if( nn == "ksp_right_pc" )
         {
            pcs = PC_RIGHT ;  
         }
         else if( nn == "ksp_left_pc" )
         {
            pcs = PC_LEFT ;  
         }
         else if( nn == "ksp_symmetric_pc" )
         {
            pcs = PC_SYMMETRIC ;
         }
         else MAC_Error::object()->raise_bad_data_value(
                                   exp, "preconditioner_side",
                                   "  - \"ksp_right_pc\"\n"
                                   "  - \"ksp_left_pc\"\n"
                                   "  - \"ksp_symmetric_pc\"\n" ) ;
         PETSc_do( KSPSetPCSide( ksp, pcs ) ) ;
      }
   }
   else
   {
      MAC_Error::object()->raise_bad_data_value(
         exp, "type",
         "  - \"richardson\"\n"
         "  - \"chebychev\"\n"
         "  - \"cg\"\n"
         "  - \"gmres\"\n"
         "  - \"tcqmr\"\n"
         "  - \"bcgs\"\n"
         "  - \"cgs\"\n"
         "  - \"tfqmr\"\n"
         "  - \"cr\"\n"
         "  - \"lsqr\"\n"
         "  - \"preonly\"\n"
         "  - \"qcg\"\n"
         "  - \"bicg\"\n"
         "  - \"fgmres\"\n"
         "  - \"minres\"\n"
         "  - \"symmlq" ) ;
   }
   if( VERB && SUBEXP==0 )
   {
//     PETSc_do( KSPMonitorSet( ksp, KSPMonitorDefault, 
//                              PETSC_NULL, PETSC_NULL ) ) ;
//for a more extensive monitoring
      PETSc_do( KSPMonitorSet( ksp,  KSPMonitorTrueResidualNorm, 
                               PETSC_NULL, PETSC_NULL ) ) ;
//
   }   
   if( has_opt ) PETSc_do( KSPSetFromOptions( ksp ) ) ;

   exp->destroy() ; exp = 0 ;
}




//----------------------------------------------------------------------------
void
EXT_PETScSolver:: build_pc( KSP & ksp, MAC_ModuleExplorer const* A_EXP ) 
//----------------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScSolver:: build_pc" ) ;
   
   PC pc ;
   PETSc_do( KSPGetPC( ksp, &pc ) ) ;
   
   if( A_EXP->has_module( "PETSc_preconditioner" ) )
   {
      MAC_ModuleExplorer* exp =
         A_EXP->create_subexplorer( 0, "PETSc_preconditioner" ) ;
      bool has_opt = EXT_PETScAPI::parse_options( exp, VERB ) ;
      
      std::string const& type = exp->string_data( "type" ) ;
      
      if( type=="none" || type=="jacobi" || type=="sor" || type=="lu"
          || type=="cholesky" || type=="eisenstat" || type=="icc" 
          || type=="ilu" || type=="hypre" || type=="prometheus" 
	  || type ==PCBJACOBI )
      {
         PETSc_do( PCSetType( pc, const_cast<PCType>( type.c_str() ) ) ) ;
      
         if( type=="sor" )
         {
            PETSc_do( PCSORSetOmega( pc, exp->double_data( "pc_sor_omega" ) ) );
            exp->test_data( "pc_sor_omega",
                            "pc_sor_omega>0. && pc_sor_omega<2." ) ;
            
            
            std::string const& s = exp->string_data( "pc_sor_mat_type" ) ;
            exp->test_data_in(
               "pc_sor_mat_type",
               "SOR_FORWARD_SWEEP,SOR_BACKWARD_SWEEP,SOR_SYMMETRIC_SWEEP,"
               "SOR_LOCAL_FORWARD_SWEEP,SOR_LOCAL_BACKWARD_SWEEP,"
               "SOR_LOCAL_SYMMETRIC_SWEEP" ) ;
            if( s=="SOR_FORWARD_SWEEP" )
            {
               PETSc_do( PCSORSetSymmetric( pc, SOR_FORWARD_SWEEP ) ) ;
            }
            else if( s=="SOR_BACKWARD_SWEEP" )
            {
               PETSc_do( PCSORSetSymmetric( pc, SOR_BACKWARD_SWEEP ) ) ;
            }
            else if( s=="SOR_SYMMETRIC_SWEEP" )
            {
               PETSc_do( PCSORSetSymmetric( pc, SOR_SYMMETRIC_SWEEP ) ) ;
            }
            else if( s=="SOR_LOCAL_FORWARD_SWEEP" )
            {
               PETSc_do( PCSORSetSymmetric( pc, SOR_LOCAL_FORWARD_SWEEP ) ) ;
            }
            else if( s=="SOR_LOCAL_BACKWARD_SWEEP" )
            {
               PETSc_do( PCSORSetSymmetric( pc, SOR_LOCAL_BACKWARD_SWEEP ) ) ;
            }
            else if( s=="SOR_LOCAL_SYMMETRIC_SWEEP" )
            {
               PETSc_do( PCSORSetSymmetric( pc, SOR_LOCAL_SYMMETRIC_SWEEP ) ) ;
            }
            else
            {
               MAC_Error::object()->raise_bad_data_value(
                  exp, "pc_sor_mat_type",
                  "  - \"SOR_FORWARD_SWEEP\"\n"
                  "  - \"SOR_BACKWARD_SWEEP\"\n"
                  "  - \"SOR_SYMMETRIC_SWEEP\"\n"
                  "  - \"SOR_LOCAL_FORWARD_SWEEP\"\n"
                  "  - \"SOR_LOCAL_BACKWARD_SWEEP\"\n"
                  "  - \"SOR_LOCAL_SYMMETRIC_SWEEP\"" ) ;
            }
         
            int its = 1 ;
            if( exp->has_entry( "pc_sor_its" ) )
            {
               its = exp->int_data( "pc_sor_its" ) ;
               exp->test_data( "pc_sor_its", "pc_sor_its>0" ) ;
               exp->set_default( "pc_sor_its", "1" ) ;
            }
            int lits = 1 ;
            if( exp->has_entry( "pc_sor_lits" ) )
            {
               lits = exp->int_data( "pc_sor_lits" ) ;
               exp->test_data( "pc_sor_lits", "pc_sor_lits>0" ) ;
               exp->set_default( "pc_sor_lits", "1" ) ;
            }
            PETSc_do( PCSORSetIterations( pc, its, lits ) ) ;
         }
         else if( type=="eisenstat" )
         {
            PETSc_do( PCEisenstatSetOmega(
                         pc, exp->double_data( "pc_eisenstat_omega") ) );
         }
         else if( type==PCBJACOBI || type==PCASM )
         {
            SUBEXP = exp->create_subexplorer( this, "sub_linear_solvers" ) ;
         }
         else if( type=="ilu" || type=="icc" )
         {
            int level = exp->int_data( "pc_factor_levels" ) ;
            
            PETSc_do( PCFactorSetLevels( pc, level ) ) ;
            
            if( exp->has_entry( "pc_factor_fill" )  )
               PETSc_do( PCFactorSetFill(
                            pc, exp->double_data( "pc_factor_fill" ) ) ) ;
         
            // PCFactorSetUseDropTolerance has disappeared in Petsc-3.2.0 
	    // (and already in Petsc-3.1.xxx)
	    // and this is not documented in "Changes", not very clear why ...
// 	    if( exp->has_module( "pc_factor_use_drop_tolerance" ) )
//             {
//                MAC_ModuleExplorer* subexp = exp->create_subexplorer(
//                   0, "pc_factor_use_drop_tolerance" ) ;
//                PETSc_do( PCFactorSetUseDropTolerance(
//                             pc,
//                             subexp->double_data( "tolerance" ),
//                             subexp->double_data( "pivot_tolerance" ),
//                             subexp->int_data( "max_non_zeros_in_row" ) ) ) ;
//                subexp->destroy() ;
//             }
            
         }
         
         if( exp->has_entry( "pc_factor_ordering_type" ) )
         {
            char* ordering =
               const_cast<char*>(
                  mat_ordering( exp, "pc_factor_ordering_type" ).c_str() ) ;
            PETSc_do( PCFactorSetMatOrderingType( pc, ordering ) ) ;
         }
         
         if( exp->has_entry( "pc_factor_zero_pivot" ) )
            PETSc_do( PCFactorSetZeroPivot(
                         pc, exp->double_data( "pc_factor_zero_pivot" ) ) ) ;

         if( exp->has_entry( "pc_factor_pivoting" ) )
            PETSc_do( PCFactorSetColumnPivot(
                         pc, exp->double_data( "pc_factor_pivoting" ) ) ) ;

         if( exp->has_entry( "pc_factor_reorder_for_nonzero_diagonal" ) )
            PETSc_do( PCFactorReorderForNonzeroDiagonal(
                         pc, exp->double_data( 
			 "pc_factor_reorder_for_nonzero_diagonal" ) ) ) ;

         // Modified to do the same as with Petsc-3.0.0
	 // In Petsc-3.2.0, values of ShiftType are:
	 // MAT_SHIFT_NONE, MAT_SHIFT_NONZERO, MAT_SHIFT_POSITIVE_DEFINITE  
	 // and MAT_SHIFT_INBLOCKS 
	 if( exp->has_entry( "pc_factor_shift_nonzero" ) )
	 {
            PETSc_do( PCFactorSetShiftType(
                         pc, MAT_SHIFT_NONZERO ) ) ;
	    PETSc_do( PCFactorSetShiftAmount(
                         pc, exp->double_data( "pc_factor_shift_nonzero" ) ) ) ;
	 }

         if( exp->has_entry( "pc_factor_allow_diagonal_fill" ) &&
             exp->bool_data( "pc_factor_allow_diagonal_fill" ) )
            PETSc_do( PCFactorSetAllowDiagonalFill( pc ) ) ;

         if( exp->has_entry( "pc_factor_pivot_in_blocks" ) &&
             exp->bool_data( "pc_factor_pivot_in_blocks" ) )
            PETSc_do( PCFactorSetPivotInBlocks( pc, PETSC_TRUE ) ) ;

      }
      else if ( type == "mumps" )
      {
        PETSc_do( PCSetType( pc, "lu" ) ) ;
	PETSc_do( PCFactorSetMatSolverPackage( pc, MATSOLVERMUMPS ) );
      }
      else
      {
         MAC_Error::object()->raise_bad_data_value(
            exp, "type",
            "  - \"none\"\n"
            "  - \"jacobi\"\n"
            "  - \"sor\"\n"
            "  - \"lu\"\n"
            "  - \"ilu\"\n"
            "  - \"icc\"\n"
            "  - \"cholesky\"\n"
            "  - \"prometheus\"\n"
            "  - \"eisenstat\""
            "  - \"mumps\"" ) ;	    
      }
      exp->destroy() ; exp = 0 ;
      if(has_opt) PETSc_do( PCSetFromOptions( pc ) ) ;
   }
}




//---------------------------------------------------------------------------
std::string
EXT_PETScSolver:: mat_ordering( MAC_ModuleExplorer const* exp,
                                    std::string const& name ) const
//---------------------------------------------------------------------------
{
   static std::string const natural = MATORDERINGNATURAL ;
   static std::string const nd = MATORDERINGND ;
   static std::string const wd = MATORDERING1WD ;
   static std::string const rcm = MATORDERINGRCM ;
   static std::string const qmd = MATORDERINGQMD ;
   
   std::string result = natural ;
   std::string const& ordering = exp->string_data( name ) ;
   if( ordering == "natural" )
   {
      result = natural ;
   }
   else if( ordering == "nd" )
   {
      result = nd ;
   }
   else if( ordering == "1wd" )
   {
      result = wd ;
   }
   else if( ordering == "rcm" )
   {
      result = rcm ;
   }
   else if( ordering == "qmd" )
   {
      result = qmd ;
   }
   else
   {
      MAC_Error::object()->raise_bad_data_value(
         exp, name,
         "  - \"natural\"\n"
         "  - \"nd\"\n"
         "  - \"1wd\"\n"
         "  - \"rcm\"\n"
         "  - \"qmd\"" ) ;
   }
   return( result ) ;
}
