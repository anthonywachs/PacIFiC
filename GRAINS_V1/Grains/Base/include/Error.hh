#ifndef _ERROR_HH_
#define _ERROR_HH_

#include <string>
#include <iostream>
#include <list>
using namespace std;

class Component;


/** @brief The class ContactError.

    To manage contact errors as exceptions.

    @author G.FERRER - Institut Francais du Petrole - 2000 - Creation
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class ContactError
{
  public:
    /**@name Constructors */
    //@{
    /** @brief Default contructor */
    ContactError();

    /** @brief Destructor */
    ~ContactError();
    //@}


    /**@name Methods */
    //@{
    /** @brief Outputs message when exception is caught
    @param fileOut output stream */
    void Message( ostream& fileOut ) const;

    /** @brief Sets the pointers to the 2 components involved in the contact
    error and the time of the contact error
    @param id0_ 1st component
    @param id1_ 2nd component
    @param time_ physical time */
    void setComponents( Component* id0_, Component* id1_, double time_ );

    /** @brief Sets header message
    @param mes header message of exception */
    void setMessage( string const& mes );

    /** @brief Returns the pointers to the 2 components involved in the contact
    error in a list for further post-processing */
    list<Component*> getComponents();
    //@}


  private:
    /**@name Parameters */
    //@{
    string m_message; /**< Header message of the exception */
    Component *m_id0; /**< 1st component involved in the contact error */
    Component *m_id1; /**< 2nd component involved in the contact error */
    double m_time; /**< physical time */
    //@}
};




/** @brief The class DisplacementError.

    To manage displacement errors as exceptions.

    @author GRAINS Project - IFP - 2007
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class DisplacementError
{
  public:
    /** @name Constructors */
    //@{
    /** @brief Constructor with input parameters
    @param id0_ component
    @param depl component displacement
    @param deplMax maximum displacement allowed
    @param time_ physical time */
    DisplacementError( Component* id0_, double depl, double deplMax,
  	double time_ );

    /** @name Contructors */
    //@{
    /** @brief Default constructor (forbidden) */
    DisplacementError();
    //@}

    /** @brief Destructor */
    ~DisplacementError();
    //@}


    /**@name Methods */
    //@{
    /** @brief Outputs message when exception is caught
    @param fileOut output stream */
    void Message( ostream& fileOut ) const;

    /** @brief Returns the pointer to the component involved in the displacement
    error in a list for further post-processing */
    list<Component*> getComponent();
    //@}


  private:
    /** @name Parameters */
    //@{
    Component *m_id0; /**< Component */
    double m_depl; /**< Component displacement */
    double m_deplMax; /**< Maximum displacement allowed */
    double m_time; /**< physical time */
    //@}
};




/** @brief The class SimulationError.

    To manage simulation errors as exceptions.

    @author GRAINS Project - IFP - 2007
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class SimulationError
{
  public:
    /**@name Constructors */
    //@{
    /** @brief Default constructor */
    SimulationError();

    /** @brief Constructor with the name of the method involved in the error
    as an input parameter
    @param str method name */
    SimulationError( string const& str );

    /** @brief Destructor */
    ~SimulationError();
    //@}

    /**@name Methods */
    //@{
    /** @brief Outputs message when exception is caught
    @param fileOut output stream */
    void Message( ostream& fileOut ) const;
    //@}


  private:
    /**@name Parameters */
    //@{
    string m_method; /**< name of the method involved in the error  */
    //@}
};

#endif
