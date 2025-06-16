#ifndef _TruncatedConicalShell_HH_
#define _TruncatedConicalShell_HH_

#include "CompositeObstacle.hh"
#include "SimpleObstacle.hh"

#include <list>
using namespace std;

#include "ReaderXML.hh"


/** @brief The class TruncatedConicalShell

    Truncated conical shell made of trapezoidal prisms.

    @author A.WACHS - 2024 - Creation */
// ============================================================================
class TruncatedConicalShell : public CompositeObstacle
{
  public:
    /** @name Constructors */
    //@{
    /** @brief Constructor with an XML node as an input parameter
    @param root XML node */
    TruncatedConicalShell( DOMNode* root );

    /** @brief Constructor with name as input parameter
    @param s obstacle name */
    TruncatedConicalShell( string const& s );

    /** @brief Destructor */
    ~TruncatedConicalShell();
    //@}


    /** @name I/O methods */
    //@{
    /** @brief Reloads the truncated conical shell and links it to the higher 
    level obstacle in the obstacle tree
    @param mother higher level obstacle
    @param file input stream */
    virtual void reload( Obstacle& mother, istream& file ) ;    

    /** @brief Outputs the truncated conical shell for reload
    @param fileSave output stream */
    virtual void write( ostream& fileSave ) const;
    //@}


  private:
    /** @name Parameters */
    //@{  
    double m_height; /**< height */
    double m_innerLargeRadius; /**< inner large radius */
    double m_innerSmallRadius; /**< inner small radius */    
    double m_shellWidth; /**< shell radial width */
    size_t m_nbBoxes; /**< number of boxes in the angular direction */   
    //@}


    /** @name Constructors */
    //@{
    /** @brief Copy constructor
    @param copy copied TruncatedConicalShell */
    TruncatedConicalShell( TruncatedConicalShell const& copy );    
    //@}   
};

#endif
