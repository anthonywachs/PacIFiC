#ifndef _CYLINDRICALSHELL_HH_
#define _CYLINDRICALSHELL_HH_

#include "CompositeObstacle.hh"
#include "SimpleObstacle.hh"

#include <list>
using namespace std;

#include "ReaderXML.hh"


/** @brief The class CylindricalShell

    Cylindrical shell made of boxes.

    @author A.WACHS - 2024 - Creation */
// ============================================================================
class CylindricalShell : public CompositeObstacle
{
  public:
    /** @name Constructors */
    //@{
    /** @brief Constructor with an XML node as an input parameter
    @param root XML node */
    CylindricalShell( DOMNode* root );

    /** @brief Constructor with name as input parameter
    @param s obstacle name */
    CylindricalShell( string const& s );

    /** @brief Destructor */
    ~CylindricalShell();
    //@}


    /** @name I/O methods */
    //@{
    /** @brief Reloads the cylindrical shell and links it to the higher level 
    obstacle in the obstacle tree
    @param mother higher level obstacle
    @param file input stream */
    virtual void reload( Obstacle& mother, istream& file ) ;    

    /** @brief Outputs the cylindrical shell for reload
    @param fileSave output stream */
    virtual void write( ostream& fileSave ) const;
    //@}


  private:
    /** @name Parameters */
    //@{  
    double m_height; /**< height */
    double m_innerRadius; /**< inner radius */
    double m_shellWidth; /**< shell radial width */
    size_t m_nbBoxes; /**< number of boxes in the angular direction */   
    //@}


    /** @name Constructors */
    //@{
    /** @brief Copy constructor
    @param copy copied CylindricalShell */
    CylindricalShell( CylindricalShell const& copy );    
    //@}   
};

#endif
