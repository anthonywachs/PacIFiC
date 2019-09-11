#include <MAC_NumberedDoubleVectors.hh>

#include <MAC_BalancedBinaryTree.hh>
#include <MAC_BalancedBinaryTreeIterator.hh>
#include <MAC_DoubleComparator.hh>
#include <MAC_Error.hh>

#include <iostream>
#include <sstream>

class MAC_NumberedDoubleVectorsItem : public MAC_Object
{
   public:
      
      MAC_NumberedDoubleVectorsItem( MAC_NumberedDoubleVectors* a_owner,
                                     doubleVector const& vector,
                                     size_t n ) ;
     ~MAC_NumberedDoubleVectorsItem( void ) ;

      void set( doubleVector const& vector, size_t n ) ;
      virtual bool is_equal( MAC_Object const* other ) const ;
      virtual int three_way_comparison( MAC_Object const* other ) const ;
      
      size_t N ;
      doubleVector VECTOR ;
} ;

//-----------------------------------------------------------------------------
MAC_NumberedDoubleVectors*
MAC_NumberedDoubleVectors:: create( MAC_Object* a_owner,
                                    MAC_DoubleComparator const* dbl_comp,
                                    size_t a_size_of_items )
//-----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_NumberedDoubleVectors:: create" ) ;
   MAC_CHECK_PRE( dbl_comp != 0 ) ;
   
   MAC_NumberedDoubleVectors* result =
      new MAC_NumberedDoubleVectors( a_owner, dbl_comp, a_size_of_items ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( result->size_of_items() == a_size_of_items ) ;
   MAC_CHECK_POST( result->nb_items() == 0 ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
MAC_NumberedDoubleVectors:: MAC_NumberedDoubleVectors(
                                    MAC_Object* a_owner,
                                    MAC_DoubleComparator const* dbl_comp,
                                    size_t a_size_of_items ) 
//-----------------------------------------------------------------------------
   : MAC_Object( a_owner )
   , DBL_COMP( dbl_comp )
   , DIM( a_size_of_items )
   , NB_PTS( 0 )
   , PT_TREE( 0 )
   , PT_TREE_IT( 0 )
   , ALL_ITEMS( 0, 0 )
   , ORDER( 0 )
{
   MAC_LABEL( "MAC_NumberedDoubleVectors:: MAC_NumberedDoubleVectors" ) ;
   PT_TREE = MAC_BalancedBinaryTree::create( this ) ;
   PT_TREE_IT = PT_TREE->create_iterator( PT_TREE ) ;
   TMP = new MAC_NumberedDoubleVectorsItem( this, doubleVector( DIM ), 0 ) ;
}

//-----------------------------------------------------------------------------
MAC_NumberedDoubleVectors:: ~MAC_NumberedDoubleVectors( void ) 
//-----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_NumberedDoubleVectors:: ~MAC_NumberedDoubleVectors" ) ;
}

//-----------------------------------------------------------------------------
size_t
MAC_NumberedDoubleVectors:: nb_items( void ) const
//-----------------------------------------------------------------------------
{
   return( NB_PTS ) ;
}

//-----------------------------------------------------------------------------
size_t
MAC_NumberedDoubleVectors:: size_of_items( void ) const
//-----------------------------------------------------------------------------
{
   return( DIM ) ;
}

//-----------------------------------------------------------------------------
bool
MAC_NumberedDoubleVectors:: has( doubleVector const& item ) const
//-----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_NumberedDoubleVectors:: has" ) ;
   MAC_CHECK_PRE( item.size()==size_of_items() ) ;

   TMP->set( item, 0 ) ;
   bool result = PT_TREE->has( TMP ) ;
   
   return( result ) ;
}

//-----------------------------------------------------------------------------
size_t
MAC_NumberedDoubleVectors:: index( doubleVector const& item ) const
//-----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_NumberedDoubleVectors:: index" ) ;
   MAC_CHECK_PRE( item.size()==size_of_items() ) ;
   MAC_CHECK_PRE( has( item ) ) ;

   TMP->set( item, 0 ) ;
   MAC_NumberedDoubleVectorsItem const* vec =
      static_cast<MAC_NumberedDoubleVectorsItem*>( PT_TREE->item( TMP ) ) ;
   
   size_t result = vec->N ;

   MAC_CHECK_POST( result < nb_items() ) ;

   return( result ) ;
}

//-----------------------------------------------------------------------------
doubleArray2D const& 
MAC_NumberedDoubleVectors:: ordered_items( void ) const 
//-----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_NumberedDoubleVectors:: ordered_items" ) ;

   ALL_ITEMS.re_initialize( size_of_items(), nb_items() ) ;
   size_t cpt = 0 ;
   for( PT_TREE_IT->start() ; PT_TREE_IT->is_valid() ; PT_TREE_IT->go_next() )
   {
      MAC_NumberedDoubleVectorsItem const* item =
         static_cast<MAC_NumberedDoubleVectorsItem*>( PT_TREE_IT->item() ) ;
      for( size_t j=0 ; j<DIM ; j++ ) ALL_ITEMS( j, cpt ) = item->VECTOR( j ) ;
      cpt++ ;
   }
   MAC_CHECK( cpt==nb_items() ) ;
   doubleArray2D const& result = ALL_ITEMS ;
   
   MAC_CHECK_POST( result.index_bound(1) == nb_items() ) ;
   MAC_CHECK_POST( result.index_bound(0) == size_of_items() ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
size_t_vector const& 
MAC_NumberedDoubleVectors:: order( void ) const 
//-----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_NumberedDoubleVectors:: order" ) ;

   ORDER.re_initialize( nb_items() ) ;
   size_t cpt=0 ;
   for( PT_TREE_IT->start() ; PT_TREE_IT->is_valid() ; PT_TREE_IT->go_next() )
   {
      MAC_NumberedDoubleVectorsItem const* item =
         static_cast<MAC_NumberedDoubleVectorsItem*>( PT_TREE_IT->item() ) ;
      size_t i = item->N ;
      ORDER( i ) = cpt++ ;
   }
   size_t_vector const& result = ORDER ;
   
   MAC_CHECK_POST( result.size() == nb_items() ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
void
MAC_NumberedDoubleVectors:: extend( doubleVector const& item )
//-----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_NumberedDoubleVectors:: extend" ) ;
   MAC_CHECK_PRE( item.size() == size_of_items() ) ;
   if( !has( item ) )
   {
      MAC_NumberedDoubleVectorsItem* new_one =
         new MAC_NumberedDoubleVectorsItem( this, item, NB_PTS++ ) ;
      PT_TREE->extend( new_one ) ;
      MAC_CHECK( PT_TREE->count()==NB_PTS ) ;
   }
   
   MAC_CHECK_POST( has( item ) ) ;
}

//internal--------------------------------------------------------------
MAC_NumberedDoubleVectorsItem:: MAC_NumberedDoubleVectorsItem(
                                       MAC_NumberedDoubleVectors* a_owner,
                                       doubleVector const& vector,
                                       size_t n )
//internal--------------------------------------------------------------
   : MAC_Object( a_owner )
   , N( n )
   , VECTOR( vector )
{
}

//internal--------------------------------------------------------------
MAC_NumberedDoubleVectorsItem:: ~MAC_NumberedDoubleVectorsItem( void )
//internal--------------------------------------------------------------
{
}

//internal--------------------------------------------------------------
void
MAC_NumberedDoubleVectorsItem:: set( doubleVector const& vector, size_t n )
//internal--------------------------------------------------------------
{
   VECTOR = vector ;
   N = n ;
}

//internal--------------------------------------------------------------
bool
MAC_NumberedDoubleVectorsItem:: is_equal( MAC_Object const* other ) const
//internal--------------------------------------------------------------
{
   MAC_CHECK( is_equal_PRE( other ) ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   bool result = three_way_comparison( other ) == 0 ;

   return( result ) ;
}

//internal--------------------------------------------------------------
int
MAC_NumberedDoubleVectorsItem:: three_way_comparison( 
                                              MAC_Object const* other ) const
//internal--------------------------------------------------------------
{
   MAC_CHECK( three_way_comparison_PRE( other ) ) ;
   MAC_CHECK( 
    dynamic_cast<MAC_NumberedDoubleVectorsItem const*>( other ) != 0 ) ;
   MAC_CHECK( VECTOR.size() == 
    static_cast<MAC_NumberedDoubleVectorsItem const*>(other)->VECTOR.size() ) ;
   MAC_CHECK( 
      dynamic_cast<MAC_NumberedDoubleVectors const*>( owner() ) != 0 ) ;

   MAC_DoubleComparator const* DBL_COMP =
      static_cast<MAC_NumberedDoubleVectors const*>( owner() )->DBL_COMP ;
   
   MAC_NumberedDoubleVectorsItem const* dbl_vec =
                   static_cast<MAC_NumberedDoubleVectorsItem const*>( other ) ;
   
   int result = 0 ;
   size_t const nb_coords = VECTOR.size() ;
   for( size_t i = 0 ; i<nb_coords && result==0 ; ++i )
   {
      result = DBL_COMP->three_way_comparison( VECTOR(i), dbl_vec->VECTOR(i) ) ;
   }

   MAC_CHECK_POST( three_way_comparison_POST( result, other ) ) ;
   return( result ) ;
}


