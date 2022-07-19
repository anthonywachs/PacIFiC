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
  MAC_LABEL( "PostProcessing:: PostProcessing" ) ;

}




//---------------------------------------------------------------------------
PostProcessing:: PostProcessing( bool is_solids_
                               , FV_DomainAndFields * dom
                               , MAC_Communicator const* macCOMM_)
//---------------------------------------------------------------------------
{
   MAC_LABEL( "PostProcessing:: PostProcessing" ) ;

   m_is_solids = is_solids_ ;
   m_allrigidbodies = NULL ;
   m_dim = dom->primary_grid()->nb_space_dimensions();
   m_macCOMM = macCOMM_ ;
   if ( !dom->has_discrete_field( "PP_epsilon" ) ) {
      MAC_Error::object()->raise_plain( "field_name PP_epsilon does not exist");
   } else {
      EPSILON = dom->discrete_field( "PP_epsilon" );
   }
   // EPSILON = porosity_ ;
   MAC_ASSERT( EPSILON->discretization_type() == "centered" ) ;
}




//---------------------------------------------------------------------------
PostProcessing:: PostProcessing( bool is_solids_
                               , FS_AllRigidBodies const* allrigidbodies_
                               , FV_DomainAndFields * dom
                               , MAC_Communicator const* macCOMM_)
//---------------------------------------------------------------------------
{
   MAC_LABEL( "PostProcessing:: PostProcessing" ) ;

   m_is_solids = is_solids_ ;
   m_allrigidbodies = allrigidbodies_ ;
   m_dim = dom->primary_grid()->nb_space_dimensions();
   m_macCOMM = macCOMM_ ;
   if ( !dom->has_discrete_field( "PP_epsilon" ) ) {
      MAC_Error::object()->raise_plain( "field_name PP_epsilon does not exist");
   } else {
      EPSILON = dom->discrete_field( "PP_epsilon" );
   }
   // EPSILON = porosity_ ;
   MAC_ASSERT( EPSILON->discretization_type() == "centered" ) ;
}




//---------------------------------------------------------------------------
PostProcessing:: ~PostProcessing()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "PostProcessing:: ~PostProcessing" ) ;



}




//---------------------------------------------------------------------------
void
PostProcessing::prepare_fieldVolumeAverageAroundRB( MAC_Object* a_owner
                                               , MAC_ModuleExplorer const* exp
                                               , FV_DomainAndFields const* dom)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "PostProcessing::prepare_fieldVolumeAverageAroundRB" ) ;

  // Read the module for box averaging
  MAC_ModuleExplorer* se =
            exp->create_subexplorer( 0,"fieldVolumeAverageAroundRB" ) ;

  MESH = dom->primary_grid();

  if (!m_is_solids) {
     std::ostringstream msg;
     msg << "Average around RB not possible in single phase flow " << endl;
     MAC_Error::object()->raise_plain( msg.str() ) ;
  }

  for (se->start_module_iterator(); se->is_valid_module();
                                    se->go_next_module() ) {
     struct fieldVolumeAverageAroundRB fva ;
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

     // Read if mean is with or without porosity
     if ( sse->has_entry( "withPorosity" ) ) {
        fva.withPorosity = sse->bool_data( "withPorosity" );
     } else {
        MAC_Error::object()->raise_missing_keyword( sse, "withPorosity" ) ;
     }

     // Read kernel
     if ( sse->has_entry( "kernel_type" ) )
        fva.kernelType = sse->int_data( "kernel_type" ) ;
     else fva.kernelType = 0 ;

     // Read volume width
     if ( sse->has_entry( "volume_width" ) )
      fva.volumeWidth = sse->double_data( "volume_width" ) ;
     else MAC_Error::object()->raise_missing_keyword( sse,"volume_width (in Dp)" ) ;

     m_fieldVolumeAverageAroundRB_list.push_back(fva);

     sse->destroy(); sse = 0;
  }
  se->destroy(); se = 0;
}




//---------------------------------------------------------------------------
void
PostProcessing::prepare_fieldVolumeAverageInBox( MAC_Object* a_owner
                                               , MAC_ModuleExplorer const* exp
                                               , FV_DomainAndFields const* dom)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "PostProcessing::prepare_fieldVolumeAverageInBox" ) ;

  // Read the module for box averaging
  MAC_ModuleExplorer* se =
                  exp->create_subexplorer( 0,"fieldVolumeAverageInBox" ) ;

  MESH = dom->primary_grid();

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

     // Read if mean is with or without porosity
     if ( sse->has_entry( "withPorosity" ) ) {
        fva.withPorosity = sse->bool_data( "withPorosity" );
     } else {
        MAC_Error::object()->raise_missing_keyword( sse, "withPorosity" ) ;
     }
     if (!m_is_solids) fva.withPorosity = false;

     // Read the gravity center of box
     doubleVector gc( m_dim, 0 );
     if ( sse->has_entry( "GravityCenter" ) ) {
        gc = sse->doubleVector_data( "GravityCenter" );
        fva.center = MAC_DoubleVector::create(a_owner,gc);
     } else {
        MAC_Error::object()->raise_missing_keyword( sse, "GravityCenter" ) ;
     }

     // Read the box length
     doubleVector bl( m_dim, 0 );
     if ( sse->has_entry( "BoxLength" ) ) {
        bl = sse->doubleVector_data( "BoxLength" ) ;
        fva.length = MAC_DoubleVector::create(a_owner,bl);
     } else {
        MAC_Error::object()->raise_missing_keyword( sse, "BoxLength" ) ;
     }

     doubleVector const& box_length = fva.length->to_double_vector();

     for (size_t i = 0; i < m_dim; i++) {
        double global_min = MESH->get_main_domain_min_coordinate(i);
        double global_max = MESH->get_main_domain_max_coordinate(i);

        if (box_length(i) > (global_max - global_min)) {
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
PostProcessing::compute_fieldVolumeAverageInBox()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "PostProcessing::compute_fieldVolumeAverageInBox" ) ;

  string fileName = "./Res/fieldVolumeAverageInBox.res";
  std::ofstream MyFile;
  if (m_macCOMM->rank() == 0) {
     MyFile.open(fileName.c_str(), std::ios::out);
     MyFile.setf(std::ios::scientific,std::ios::floatfield);
     MyFile.precision(6);
  }

  list<struct fieldVolumeAverageInBox>::const_iterator it;

  for ( it = m_fieldVolumeAverageInBox_list.begin()
      ; it != m_fieldVolumeAverageInBox_list.end()
      ; it++ ) {

      doubleVector const& box_length = it->length->to_double_vector();
      doubleVector const& box_center = it->center->to_double_vector();
      size_t ncomps = it->FF->nb_components();

      if (m_macCOMM->rank() == 0) {
         MyFile << "# Average "
                << it->field_name
                << " in "
                << it->box_name
                << " with porosity "
                << it->withPorosity << endl;
      }

      intVector min_local_index(m_dim,0);
      intVector max_local_index(m_dim,0);

      for (size_t dir = 0; dir < m_dim; dir++) {
         // Control volume min and max, if any
         doubleVector box_extents(2,0.);
         box_extents(0) = box_center(dir) - box_length(dir)/2.;
         box_extents(1) = box_center(dir) + box_length(dir)/2.;

         intVector temp =
         get_local_index_of_extents(box_extents, EPSILON, dir, 0);

         min_local_index(dir) = temp(0);
         max_local_index(dir) = temp(1);
      }

      for (size_t comp = 0; comp < ncomps; comp++) {
         double value = 0.;
         double volume = 0.;

         if (min_local_index(0) != -1) {
            for (int i = min_local_index(0);
                     i <= max_local_index(0); i++) {
               double dx = it->FF->get_cell_size( i, comp, 0 );
               if (min_local_index(1) != -1) {
                  for (int j = min_local_index(1);
                           j <= max_local_index(1); j++) {
                     double dy = it->FF->get_cell_size( j, comp, 1 );
                     if (m_dim == 2 || min_local_index(2) != -1) {
                        for (int k = min_local_index(2);
                                 k <= max_local_index(2); k++) {
                           double dz = (m_dim == 3) ?
                                       it->FF->get_cell_size( k, comp, 2 ) : 1;
                           double epsilon = 1;
                           if (it->withPorosity) {
                              epsilon = EPSILON->DOF_value(i, j, k, 0, 0);
                           }
                           value += field_value(it->FF, i , j, k, comp, 0)
                                  * epsilon * dx * dy * dz;
                           volume += dx * dy * dz * epsilon;
                        }
                     }
                  }
               }
            }
         }

         value = m_macCOMM->sum( value ) ;
         volume = m_macCOMM->sum( volume ) ;
         if (m_macCOMM->rank() == 0) {
            if (volume != 0.) {
               MyFile << value/volume << '\t';
            } else {
               MyFile << 0. << '\t';
            }
         }

     }  // comp
     if (m_macCOMM->rank() == 0) MyFile << endl;
  }  // box
}




//---------------------------------------------------------------------------
double
PostProcessing::field_value(FV_DiscreteField const* FF, size_t const& i
                                                      , size_t const& j
                                                      , size_t const& k
                                                      , size_t const& comp
                                                      , size_t const& level)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "PostProcessing::field_value" ) ;

  FV_SHIFT_TRIPLET shift = EPSILON->shift_staggeredToCentered() ;

  double field = 0.;

  size_t ncomps = FF->nb_components();

  if (ncomps == 1) {   //centered
     field = FF->DOF_value(i, j, k, comp, level);
  } else { //staggered
    if (comp == 0) {
       field = (FF->DOF_value(i+shift.i-1, j, k, comp, level)
              + FF->DOF_value(i+shift.i, j, k, comp, level)) * 0.5;
    } else if (comp == 1) {
       field = (FF->DOF_value(i, j+shift.j-1, k, comp, level)
              + FF->DOF_value(i, j+shift.j, k, comp, level)) * 0.5;
    } else if (comp == 2) {
       field = (FF->DOF_value(i, j, k+shift.k-1, comp, level)
              + FF->DOF_value(i, j, k+shift.k, comp, level)) * 0.5;
    }
  }

  return(field);

}




//---------------------------------------------------------------------------
void
PostProcessing::compute_fieldVolumeAverageAroundRB()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "PostProcessing::compute_fieldVolumeAverageAroundRB" ) ;

  string fileName = "./Res/fieldVolumeAverageAroundRB.res";
  std::ofstream MyFile;
  if (m_macCOMM->rank() == 0) {
     MyFile.open(fileName.c_str(), std::ios::out);
     MyFile.setf(std::ios::scientific,std::ios::floatfield);
     MyFile.precision(6);
  }
  size_t m_nrb = m_allrigidbodies->get_number_rigid_bodies();
  list<struct fieldVolumeAverageAroundRB>::const_iterator it;

  for (size_t parID = 0; parID < m_nrb; parID++) {
     for ( it = m_fieldVolumeAverageAroundRB_list.begin()
         ; it != m_fieldVolumeAverageAroundRB_list.end()
         ; it++ ) {

         if (m_macCOMM->rank() == 0) {
            MyFile << "# Average "
                   << it->field_name
                   << " around rigid body "
                   << parID
                   << " with kernel "
                   << it->kernelType
                   << ", the box size "
                   << it->volumeWidth
                   << ", and porosity "
                   << it->withPorosity
                   << endl;
         }

         FS_RigidBody const* rigidBody = m_allrigidbodies
                                             ->get_ptr_rigid_body(parID);

         geomVector const* ptgc = rigidBody->get_ptr_to_gravity_centre();
         double cr = rigidBody->get_circumscribed_radius();
         size_t ncomps = it->FF->nb_components();

         intVector min_local_index(m_dim,0);
         intVector max_local_index(m_dim,0);

         for (size_t dir = 0; dir < m_dim; dir++) {
            // Control volume min and max, if any
            doubleVector box_extents(2,0.);
            box_extents(0) = ptgc->operator()(dir)
                           - (0.5*it->volumeWidth*(2.*cr));
            box_extents(1) = ptgc->operator()(dir)
                           + (0.5*it->volumeWidth*(2.*cr));

            intVector temp =
            get_local_index_of_extents(box_extents, EPSILON, dir, 0);

            min_local_index(dir) = temp(0);
            max_local_index(dir) = temp(1);
         }

         for (size_t comp = 0; comp < ncomps; comp++) {
            double value = 0.;
            double volume = 0.;

            if (min_local_index(0) != -1) {
               for (int i = min_local_index(0);
                        i <= max_local_index(0); i++) {
                  double dx = it->FF->get_cell_size( i, comp, 0 );
                  double xC = it->FF->get_DOF_coordinate( i, comp, 0) ;
                  xC = delta_periodic_transformation(
                                 fabs(xC - ptgc->operator()(0)), 0);
                  if (min_local_index(1) != -1) {
                     for (int j = min_local_index(1);
                              j <= max_local_index(1); j++) {
                        double dy = it->FF->get_cell_size( j, comp, 1 );
                        double yC = it->FF->get_DOF_coordinate( j, comp, 1) ;
                        yC = delta_periodic_transformation(
                                       fabs(yC - ptgc->operator()(1)), 1);
                        if (m_dim == 2 || min_local_index(2) != -1) {
                           for (int k = min_local_index(2);
                                    k <= max_local_index(2); k++) {
                              double dz = (m_dim == 3) ?
                                       it->FF->get_cell_size( k, comp, 2 ) : 1;
                              double zC = (m_dim == 3) ?
                                  it->FF->get_DOF_coordinate( k, comp, 2) : 0.;
                              zC = delta_periodic_transformation(
                                           fabs(zC - ptgc->operator()(2)), 2);

                              double distance = MAC::sqrt(xC*xC + yC*yC + zC*zC);
                              double weight = kernel(distance
                                                   , cr
                                                   , it->volumeWidth
                                                   , it->kernelType);
                              double epsilon = 1;
                              if (it->withPorosity) {
                                 epsilon = EPSILON->DOF_value(i, j, k, 0, 0);
                              }
                              value += field_value(it->FF, i , j, k, comp, 0)
                                       * epsilon * weight * dx * dy * dz;
                              volume += epsilon * weight * dx * dy * dz;
                           }
                        }
                     }
                  }
               }
            }

            value = m_macCOMM->sum( value ) ;
            volume = m_macCOMM->sum( volume ) ;
            if (m_macCOMM->rank() == 0) {
               if (volume != 0.) {
                  MyFile << value/volume << '\t';
               } else {
                  MyFile << 0. << '\t';
               }
            }

         }
         if (m_macCOMM->rank() == 0) MyFile << endl;
      }
   }
}




//---------------------------------------------------------------------------
intVector
PostProcessing::get_local_index_of_extents( class doubleVector& bounds
                                          , FV_DiscreteField const* FF
                                          , size_t const& dir
                                          , size_t const& comp)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "PostProcessing::get_grid_index_of_extents" ) ;

  intVector value(2,0);

  bounds(0) = periodic_transformation(bounds(0),dir);
  bounds(1) = periodic_transformation(bounds(1),dir);

  double global_min = MESH->get_main_domain_min_coordinate(dir);
  double global_max = MESH->get_main_domain_max_coordinate(dir);
  double local_min = MESH->get_min_coordinate_on_current_processor(dir);
  double local_max = MESH->get_max_coordinate_on_current_processor(dir);

  boolVector const* is_periodic = MESH->get_periodic_directions();

  // Warnings
  if (m_macCOMM->rank() == 0) {
     if (!(*is_periodic)(dir)) {
        if ((bounds(0) < global_min)
        || (bounds(1) > global_max))
           std::cout << endl <<
             " WARNING : Box Averaging Control Volume overlaps a " <<
             " non-periodic BC" << endl <<
             " Control volume will be reduced" << endl << endl;
     }
  }

  // Getting the minimum grid index in control volume (CV)
  if (bounds(0) < bounds(1)) {// Non-periodic CV
     if (bounds(0) > local_min) {
        if (bounds(0) < local_max) {
           size_t i0_temp = 0;
           bool found = FV_Mesh::between(
                          FF->get_DOF_coordinates_vector(comp,dir)
                        , bounds(0)
                        , i0_temp) ;
           value(0) = (found) ? (int)i0_temp : 0;
        } else {
           value(0) = -1;
        }
     } else {
        if (bounds(1) > local_min) {
           value(0) = (int)FF->get_min_index_unknown_handled_by_proc(comp,dir);
        } else {
           value(0) = -1;
        }
     }
  } else {// periodic control volume
     if (bounds(0) > local_min) {
        if (bounds(0) < local_max) {
           size_t i0_temp = 0;
           bool found = FV_Mesh::between(
                          FF->get_DOF_coordinates_vector(comp,dir)
                        , bounds(0)
                        , i0_temp) ;
           value(0) = (found) ? (int)i0_temp : 0;
        } else if (bounds(1) > local_min) {
           value(0) = (int)FF->get_min_index_unknown_handled_by_proc(comp,dir);
        } else {
           value(0) = -1;
        }
     } else {
        value(0) = (int)FF->get_min_index_unknown_handled_by_proc(comp,dir);
     }
  }

  // Getting the maximum grid index in control volume (CV)
  if (bounds(0) < bounds(1)) {// Non-periodic CV
     if (bounds(1) < local_max) {
        if (bounds(1) > local_min) {
           size_t i0_temp = 0;
           bool found = FV_Mesh::between(
                          FF->get_DOF_coordinates_vector(comp,dir)
                        , bounds(1)
                        , i0_temp) ;
           value(1) = (found) ? (int)i0_temp : 0;
        } else {
           value(1) = -1;
        }
     } else {
        if (bounds(0) < local_max) {
           value(1) = (int)FF->get_max_index_unknown_handled_by_proc(comp,dir);
        } else {
           value(1) = -1;
        }
     }
  } else {// periodic control volume
     if (bounds(1) < local_max) {
        if (bounds(1) > local_min) {
           size_t i0_temp = 0;
           bool found = FV_Mesh::between(
                          FF->get_DOF_coordinates_vector(comp,dir)
                        , bounds(1)
                        , i0_temp) ;
           value(1) = (found) ? (int)i0_temp : 0;
        } else if (bounds(0) < local_max) {
           value(1) = (int)FF->get_max_index_unknown_handled_by_proc(comp,dir);
        } else {
           value(1) = -1;
        }
     } else {
        value(1) = (int)FF->get_max_index_unknown_handled_by_proc(comp,dir);
     }
  }

  return (value);

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




//---------------------------------------------------------------------------
double PostProcessing:: periodic_transformation( double const& x
                                                , size_t const& dir)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "PostProcessing:: periodic_transformation" ) ;

  double value = x;

  boolVector const* is_periodic = MESH->get_periodic_directions();
  // Control volume periodic treatment, if any
  if ((*is_periodic)(dir)) {
     double isize = MESH->get_main_domain_max_coordinate(dir)
                  - MESH->get_main_domain_min_coordinate(dir);
     double imin = MESH->get_main_domain_min_coordinate(dir);
     value = x - MAC::floor((x-imin)/isize)*isize;
  }

  return (value);
}




//---------------------------------------------------------------------------
double PostProcessing:: delta_periodic_transformation( double const& delta
                                                     , size_t const& dir)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "PostProcessing:: delta_periodic_transformation" ) ;

  double value = delta;

  boolVector const* is_periodic = MESH->get_periodic_directions();
  // Control volume periodic treatment, if any
  if ((*is_periodic)(dir)) {
     double isize = MESH->get_main_domain_max_coordinate(dir)
                  - MESH->get_main_domain_min_coordinate(dir);
     value = delta - round(delta/isize) * isize;
  }

  return (value);
}
