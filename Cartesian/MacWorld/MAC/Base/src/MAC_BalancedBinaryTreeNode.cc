#include <MAC_BalancedBinaryTreeNode.hh>

#include <MAC.hh>
#include <MAC_assertions.hh>
#include <MAC_Error.hh>

//-------------------------------------------------------------------------
MAC_BalancedBinaryTreeNode:: MAC_BalancedBinaryTreeNode( MAC_Object * value )
//-------------------------------------------------------------------------
   :   MAC_Object(0), height(0), val(value), less(0), more(0)
{
   MAC_CHECK_INV( invariant() ) ;
}



//-------------------------------------------------------------------------
MAC_BalancedBinaryTreeNode:: ~MAC_BalancedBinaryTreeNode( void )
//-------------------------------------------------------------------------
{
   MAC_CHECK_INV( invariant() ) ;
   clear() ;
}



//-------------------------------------------------------------------------
MAC_Object *
MAC_BalancedBinaryTreeNode:: search( const MAC_Object * value ) const
//-------------------------------------------------------------------------
{
   MAC_CHECK_PRE( value!=0 ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_Object * ret = 0 ;
   int cmp = val->three_way_comparison(value) ;
   
   if( cmp < 0 )
   {
      if( more!=0 )
      {
         ret = more->search( value );
      }
   }
   else
   {
      if( cmp == 0 )
      {
         ret = val ;
      }
      else
      {
         if( less!=0 )
         {
            ret = less->search( value ) ;
         }
      }
   }
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( ret==0 || ret->is_equal( value ) ) ;
   return( ret ) ;
}



//-------------------------------------------------------------------------
MAC_Object *
MAC_BalancedBinaryTreeNode:: addChild(
   MAC_Object * value,
   MAC_BalancedBinaryTreeNode * & top ) 
//-------------------------------------------------------------------------
{
   MAC_CHECK_PRE( value!=0 ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_Object * ret = 0 ;
   if( val == 0 )
   {
      val = value ;
      ret = val ;
   }
   else
   {
      int cmp = val->three_way_comparison( value ) ;
   
      if( cmp<0 )
      {
         if( more!=0 )
         {
            ret = more->addChild( value, more ) ;
         }
         else
         {
            more = new MAC_BalancedBinaryTreeNode( value ) ;
            ret = value ;
         }
      }
      else
      {
         if( cmp==0 )
         {
            ret=val ;
         }
         else
         {
            if( less!=0 )
            {
               ret = less->addChild( value, less )  ;
            }
            else
            {
               less =  new MAC_BalancedBinaryTreeNode( value ) ;
               ret = value ;
            }
         }
      }
      balance( top ) ;
   }
   MAC_CHECK_INV( invariant() ) ;
   return( ret ) ;
}



//-------------------------------------------------------------------------
MAC_BalancedBinaryTreeNode *
MAC_BalancedBinaryTreeNode:: remove(
   const MAC_Object * value,
   MAC_BalancedBinaryTreeNode * & top ) 
//-------------------------------------------------------------------------
{
   MAC_CHECK_PRE( value!=0 ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_BalancedBinaryTreeNode * ret = 0 ;
   int cmp = val->three_way_comparison( value ) ;
   
   if( cmp<0 )
   {
      if( more!=0 )
      {
         ret = more->remove( value, more ) ;
         balance( top ) ;
      }
      else
      {
         MAC_Error::object()->raise_plain( "Remove : error : item doesn't exist" ) ;
      }
   }
   else if( cmp>0 )
   {
      if( less!=0 )
      {
         ret = less->remove( value, less )  ;
         balance( top ) ;
      }
      else
      {
         MAC_Error::object()->raise_plain( "Remove : error : item doesn't exist" ) ;
      }
   }
   else // cmp == 0
   {
      top = 0 ;
      if( less!=0 )
      {
         top = less ;
         if( more!=0 )
         {
            less->copySubTree( more, top ) ;
         }
      }
      else if( more!=0 )
      {
         top = more ;
      }
      ret = this ;
      less = 0 ;
      more = 0 ;
   }
   MAC_CHECK_INV( invariant() ) ;   
   return( ret ) ;
}



//-------------------------------------------------------------------------
void
MAC_BalancedBinaryTreeNode:: balance( MAC_BalancedBinaryTreeNode * & top )
//-------------------------------------------------------------------------
{
   MAC_CHECK_INV( invariant() ) ;
   // Balancing
   int lessH = -1 ;
   if( less!=0 )
   {
      lessH = less->height ;
   }
   int moreH = -1 ;
   if( more!=0 )
   {
      moreH = more->height ;
   }
   //calculHeight() ;
   if( lessH > moreH+1 )
   {
      //MAC::out() << "Before less balancing : " << *top ;
      MAC_BalancedBinaryTreeNode * lessMore = less->more ;
      less->more = this ;
      if( this==top )
      {
         top = less ;
      }
      MAC_BalancedBinaryTreeNode * oldLess = less ;
      less = lessMore ;
      calculHeight() ;
      oldLess->calculHeight() ;
      //MAC::out() << "After balancing : " << *top ;
   } 
   else if( moreH > lessH+1 )
   {
      //MAC::out() << "Before more balancing : " << *top ;
      MAC_BalancedBinaryTreeNode * moreLess = more->less ;
      more->less = this ;
      if( this==top )
      {
         top = more ;
      }
      MAC_BalancedBinaryTreeNode * oldMore = more ;
      more = moreLess ;
      calculHeight() ;
      oldMore->calculHeight() ;
      //MAC::out() << "After balancing : " << *top ;
   }
   else
   {
      calculHeight() ;
   }
   MAC_CHECK_INV( invariant() ) ;
}



//-------------------------------------------------------------------------
void
MAC_BalancedBinaryTreeNode:: copySubTree( MAC_BalancedBinaryTreeNode * source,
                                          MAC_BalancedBinaryTreeNode *& top)
//-------------------------------------------------------------------------
{
   if( source->less )
   {
      copySubTree( source->less, top ) ;
   }
   addChild( source->val, top ) ;
   if( source->more )
   {
      copySubTree( source->more, top ) ;
   }
}

      
//-------------------------------------------------------------------------
void
MAC_BalancedBinaryTreeNode:: calculHeight( void )
//-------------------------------------------------------------------------
{
   int lessH = -1 ;
   if( less!=0 )
   {
      lessH = less->height ;
   }
   int moreH = -1 ;
   if( more!=0 )
   {
      moreH = more->height ;
   }
   height=(size_t)MAC::max( lessH, moreH )+1 ;
}



//-------------------------------------------------------------------------
void
MAC_BalancedBinaryTreeNode:: clear( void )
//-------------------------------------------------------------------------
{
   val = 0 ;
   
   if( more!=0 )
   {
      more->destroy() ;
      more = 0 ;
   }
   if( less!=0 )
   {
      less->destroy() ;
      less = 0 ;
   }
}



//-------------------------------------------------------------------------
size_t
MAC_BalancedBinaryTreeNode:: getNbElem( void ) const
//-------------------------------------------------------------------------
{
   // These methodis used for invariant verification purpose
   size_t nbElem = 1 ;
   MAC_ASSERT( val!=0 ) ;
   
   if( less!=0 )
   {
      nbElem += less->getNbElem() ;
   }
   if( more!=0 )
   {
      nbElem += more->getNbElem() ;
   }
   return nbElem ;
}



//-------------------------------------------------------------------------
bool
MAC_BalancedBinaryTreeNode:: invariant( void ) const
//-------------------------------------------------------------------------
{
   MAC_ASSERT( MAC_Object::invariant() ) ;
   MAC_ASSERT( val!=0 ) ;

   return( true ) ;
}
