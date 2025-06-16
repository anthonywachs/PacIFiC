#ifndef _COLLISIONHISTORY_HH_
#define _COLLISIONHISTORY_HH_


#include <map>
#include "Vector3.hh"

using namespace solid;


// =============================================================================
/** @brief The class CollisionHistory.

    A class to track the collision direction between two rigid bodies using 
    their IDs. We use map here, so lookups are constant in time. 
    Insertion/Removal is expensive, but we hope that the benefits we get from 
    lookups is worth it.
    Note that there is no need to create a comparison class for pairs as it
    already implemented in the std library.

    @author A.YAZDANI - 2024 - Construction */
// =============================================================================
class CollisionHistory
{
    public:
		/**@name Constructors */
		//@{
		/** @brief Default constructor */
		CollisionHistory();

		/** @brief Destructor */
		~CollisionHistory();
		//@}


		/**@name Methods */
		//@{
		/** @brief Updates (or adds) a pair in (to) the map
		@param id the id
		@param contactVec contact vector in the local coordinate system */
		void updateCollision( int id, Vector3 contactVec );
	
		/** @brief Removes a pair from the map 
		@param id the id */
		void removeCollision( int id ); 
	
		/** @brief Returns the overlap vector if the data for the ids is 
		available
		@param id the id */
		Vector3 lookupCollision( int id ) const;
		//@}


	private:
		/**@name Parameter */
		//@{
		/** @brief */
		map<int, Vector3> m_collisionMap;
		//@}
};

#endif
