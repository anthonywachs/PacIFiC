#include <PostProcessing.hh>
#include <MAC_ModuleExplorer.hh>
#include <math.h>
#include <MAC.hh>
#include <MAC_Error.hh>
#include <FV_Mesh.hh>
using std::endl;



//---------------------------------------------------------------------------
PostProcessing:: PostProcessing()
//---------------------------------------------------------------------------
  : m_is_solids( false )
  , m_allrigidbodies( NULL )
{
  MAC_LABEL( "DS_AllRigidBodies:: DS_AllRigidBodies" ) ;

}




//---------------------------------------------------------------------------
PostProcessing:: PostProcessing( bool is_solids_
                               , DS_AllRigidBodies const* allrigidbodies_
                               , size_t const& dim_
                               , MAC_Communicator const* macCOMM_)
//---------------------------------------------------------------------------
{
   MAC_LABEL( "PostProcessing:: PostProcessing" ) ;

   m_is_solids = is_solids_ ;
   m_allrigidbodies = allrigidbodies_ ;
   m_dim = dim_ ;
   m_macCOMM = macCOMM_ ;

}




//---------------------------------------------------------------------------
PostProcessing:: ~PostProcessing()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "PostProcessing:: ~PostProcessing" ) ;



}




//---------------------------------------------------------------------------
double
PostProcessing:: compute_mean( FV_DiscreteField const* FF
                             , size_t const& level
                             , size_t const& comp )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "PostProcessing:: compute_mean" ) ;

   size_t_vector min_unknown_index(3,0);
   size_t_vector max_unknown_index(3,0);

   for (size_t l = 0; l < 3; ++l) {
      min_unknown_index(l) =
                   FF->get_min_index_unknown_handled_by_proc( comp, l ) ;
      max_unknown_index(l) =
                   FF->get_max_index_unknown_handled_by_proc( comp, l ) ;
   }

   double value = 0.;
   double volume = 0.;

   for (size_t i = min_unknown_index(0); i <= max_unknown_index(0); ++i) {
      double dx = FF->get_cell_size( i, comp, 0 );
      for (size_t j = min_unknown_index(1); j <= max_unknown_index(1); ++j) {
         double dy = FF->get_cell_size( j, comp, 1 );
         for (size_t k = min_unknown_index(2); k <= max_unknown_index(2); ++k) {
            double dz = (m_dim == 3) ? FF->get_cell_size( k, comp, 2 ) : 1;
            value += FF->DOF_value(i, j, k, comp, level)*dx*dy*dz;
            volume += dx*dy*dz;
         }
      }
   }

   value = m_macCOMM->sum( value ) ;
   volume = m_macCOMM->sum( volume ) ;
   std::cout << "Mean value: " << value/volume << endl;

   return (value/volume);

}




//---------------------------------------------------------------------------
void
PostProcessing::prepare_fieldVolumeAverageInBox( MAC_ModuleExplorer const* exp
                                               , FV_DomainAndFields const* dom)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "PostProcessing::prepare_fieldVolumeAverageInBox" ) ;

  // Read the module for box averaging
  MAC_ModuleExplorer* se =
                  exp->create_subexplorer( 0,"fieldVolumeAverageInBox" ) ;

  for (se->start_module_iterator(); se->is_valid_module();
                                    se->go_next_module() ) {
     struct fieldVolumeAverageInBox fva ;
     MAC_ModuleExplorer* sse =
                     se->create_subexplorer( 0 ) ;
     string field_name = "none";

     // Read the field name
     if ( sse->has_entry( "field_name" ) ) {
        field_name = sse->string_data( "field_name" );
        fva.field_name = field_name;
     } else {
        MAC_Error::object()->raise_missing_keyword( sse, "field_name" ) ;
     }

     if ( !dom->has_discrete_field( field_name ) ) {
        MAC_Error::object()->raise_bad_data_value( sse,
                            "field_name", field_name+" does not exist" );
     } else {
        fva.FF = dom->discrete_field( field_name );
     }

     // Read the box name
     if ( sse->has_entry( "box_name" ) ) {
        string box_name = sse->string_data( "box_name" );
        fva.box_name = box_name;
     } else {
        MAC_Error::object()->raise_missing_keyword( sse, "box_name" ) ;
     }

     // Read box parameters
     geomVector ctemp(3), ltemp(3);
     if ( sse->has_entry( "GravityCenterX" ) ) {
        ctemp(0) = sse->double_data( "GravityCenterX" ) ;
     } else {
        MAC_Error::object()->raise_missing_keyword( sse,
                                          "GravityCenterX" ) ;
     }
     if ( sse->has_entry( "GravityCenterY" ) ) {
        ctemp(1) = sse->double_data( "GravityCenterY" ) ;
     } else {
        MAC_Error::object()->raise_missing_keyword( sse,
                                          "GravityCenterY" ) ;
     }
     if ( sse->has_entry( "GravityCenterZ" ) ) {
        ctemp(2) = sse->double_data( "GravityCenterZ" ) ;
     } else {
        MAC_Error::object()->raise_missing_keyword( sse,
                                          "GravityCenterZ" ) ;
     }

     if ( sse->has_entry( "LX" ) ) {
        ltemp(0) = sse->double_data( "LX" ) ;
     } else {
        MAC_Error::object()->raise_missing_keyword( sse,
                                                      "LX" ) ;
     }
     if ( sse->has_entry( "LY" ) ) {
        ltemp(1) = sse->double_data( "LY" ) ;
     } else {
        MAC_Error::object()->raise_missing_keyword( sse,
                                                      "LY" ) ;
     }
     if ( sse->has_entry( "LZ" ) ) {
        ltemp(2) = sse->double_data( "LZ" ) ;
     } else {
        MAC_Error::object()->raise_missing_keyword( sse,
                                                      "LZ" ) ;
     }

     fva.center = ctemp;
     fva.length = ltemp;

     for (size_t i = 0; i < m_dim; i++) {
        double global_min = dom->primary_grid()
                               ->get_main_domain_min_coordinate(i);
        double global_max = dom->primary_grid()
                               ->get_main_domain_max_coordinate(i);

        if (fva.length(i) > (global_max - global_min)) {
           std::ostringstream msg;
           msg << "In MAC_Utils:: prepare_fieldVolumeAverageInBox()" <<endl;
           msg << "--> Control Volume larger than domain in " << i
                << " direction" << endl;
           MAC_Error::object()->raise_plain( msg.str() ) ;
        }
     }

     m_fieldVolumeAverageInBox_list.push_back(fva);

     sse->destroy(); sse = 0;
  }
  se->destroy(); se = 0;
}




//---------------------------------------------------------------------------
void
PostProcessing::compute_fieldVolumeAverageInBox( )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "PostProcessing::compute_fieldVolumeAverageInBox" ) ;



}




//---------------------------------------------------------------------------
double
PostProcessing::kernel( double distance, double radius, double boxWidth,
                        size_t kernel_type)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "PostProcessing::kernel" ) ;
  double result = 1.;
  double maxL = 0.5*pow(3.,0.5)*boxWidth ;

  switch ( kernel_type ) {
    case 0 :
    {
     result=1.;
     break;
    }
    case 1 :
    {
     result=1.-distance/(2.*radius*maxL);
     break;
    }
    case 2 :
    {
     result=1.-pow(distance/(2.*radius*maxL),2.);
     break;
    }
    case 3 :
    {
     result=pow(1.-pow(distance/(2.*radius*maxL),2.),2.);
     break;
    }
    case 4 :
    {
     result=pow(1.-pow(distance/(2.*radius*maxL),2.),3.);
     break;
    }
    case 5 :
    {
     result=1./pow(2.*MAC::pi(),0.5)*exp(-1./2.*pow(distance/radius,2.));
     break;
    }
    case 6 :
    {
      result=exp(-distance/radius);
      break;
    }
    case 7 :
    {
      if ( distance >= radius )
        result=exp(-(distance-radius)/radius);
      else result=0.;
      break;
    }
    case 8 :
    {
      result=exp(-2.*distance/radius);
      break;
    }
    case 9 :
    {
      result=exp(-distance/(2.*radius));
      break;
    }
  }

  return (result);
}
