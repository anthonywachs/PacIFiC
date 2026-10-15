#ifndef _RAWDATAPOSTPROCESSINGWRITER_HH_
#define _RAWDATAPOSTPROCESSINGWRITER_HH_

#include <fstream>

#include "PostProcessingWriter.hh"

// =================================================================================================
/** @brief The class RawPostProcessingWriter.

    Writes particle data in raw format for post-processing.

    @author A.Yazdani - 2024 - Construction */
// =================================================================================================
template <typename T>
class RawDataPostProcessingWriter : public PostProcessingWriter<T>
{
private:
    /**@name Parameters */
    //@{
    /** \brief center of mass x-coordinate output stream */
    std::ofstream m_gc_coordinates_x;
    /** \brief center of mass y-coordinate output stream */
    std::ofstream m_gc_coordinates_y;
    /** \brief center of mass z-coordinate output stream */
    std::ofstream m_gc_coordinates_z;
    /** \brief x-component of translational velocity along output stream */
    std::ofstream m_translational_velocity_x;
    /** \brief y-component of translational velocity along output stream */
    std::ofstream m_translational_velocity_y;
    /** \brief z-component of translational velocity along output stream */
    std::ofstream m_translational_velocity_z;
    /** \brief x-component of angular velocity along output stream */
    std::ofstream m_angular_velocity_x;
    /** \brief y-component of angular velocity along output stream */
    std::ofstream m_angular_velocity_y;
    /** \brief z-component of angular velocity along output stream */
    std::ofstream m_angular_velocity_z;
    /** \brief coordination number output stream */
    std::ofstream m_coordination_number;
    /** \brief particle class output stream */
    std::ofstream m_particle_class;
    /** \brief files root name */
    std::string m_rootName;
    /** \brief output directory */
    std::string m_directory;
    /** \brief No. digits after the decimal in the scientific format */
    int m_ndigits;
    //@}

public:
    /** @name Constructors */
    //@{
    /** @brief Default constructor */
    RawDataPostProcessingWriter();

    /** @brief Constructor with XML node
        @param dn XML node */
    RawDataPostProcessingWriter(DOMNode* dn);

    /** @brief Destructor */
    ~RawDataPostProcessingWriter();
    //@}

    /** @name Get methods */
    //@{

    PostProcessingWriterType getPostProcessingWriterType() const;
    //@}

    /** @name Methods */
    //@{
    /** @brief Removes post-processing files already in the directory */
    void clearPostProcessingFiles() const;

    /** @brief Initializes the post-processing writer */
    void PostProcessing_start();

    /** @brief Writes post-processing data
        @param rb Arrays of rigid bodies
        @param cm Component manager
        @param currentTime Current simulation time */
    void PostProcessing(const GrainsMemBuffer<RigidBody<T>*>&       rb,
                        const std::unique_ptr<ComponentManager<T>>& cm,
                        const T                                     currentTime) final;

    /** @brief Finalizes writing data */
    void PostProcessing_end();

    /** @brief Creates output files and open streams */
    void prepareResultFiles(ios_base::openmode mode);
    //@}
};

#endif
