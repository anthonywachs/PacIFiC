#include "CollisionHistory.hh"


// -----------------------------------------------------------------------------
// Default constructor
CollisionHistory::CollisionHistory()
{}




// -----------------------------------------------------------------------------
// Destructor
CollisionHistory::~CollisionHistory()
{}




// -----------------------------------------------------------------------------
// Updates (or adds) a pair in (to) the map
void CollisionHistory::updateCollision( int id, Vector3 contactVec )
{
    m_collisionMap.insert_or_assign( id, contactVec );
}




// -----------------------------------------------------------------------------
// Removes a pair from the map 
void CollisionHistory::removeCollision( int id )
{
    m_collisionMap.erase( id );
}




// -----------------------------------------------------------------------------
// Returns the overlap vector if the data for the id is available
Vector3 CollisionHistory::lookupCollision( int id ) const
{
	auto it = m_collisionMap.find( id );
	if ( it == m_collisionMap.end() ) // the id is not in the map
		return ( Vector3Null );
	else
		return ( it->second );
}