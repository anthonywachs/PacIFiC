#include "ParaviewPostProcessingWriter.hh"
#include "GrainsUtils.hh"
#include "VectorMath.hh"

/* ============================================================================================== */
/* Low-Level Methods                                                                              */
/* ============================================================================================== */
// Writes obstacles data
template <typename T>
void writeObstacles_Paraview(const GrainsMemBuffer<RigidBody<T>*>&       rb,
                             const std::unique_ptr<ComponentManager<T>>& cm,
                             const std::string&                          obsFileName)
{
    std::ofstream f(obsFileName, std::ios::out);
    if(!f.is_open())
    {
        std::cerr << "Failed to open file: " << obsFileName << std::endl;
        throw std::runtime_error("Cannot open file for writing: " + obsFileName);
    }
    const uint                            numObstacles = cm->getNumberOfObstacles();
    const GrainsMemBuffer<Vector3<T>>&    position     = cm->getPosition();
    const GrainsMemBuffer<Quaternion<T>>& quaternion   = cm->getQuaternion();
    const GrainsMemBuffer<Kinematics<T>>& kin          = cm->getVelocity();
    GrainsMemBuffer<Transform3<T>>        tr(numObstacles);
    for(uint i = 0; i < numObstacles; ++i)
        tr[i] = Transform3<T>(quaternion[i], position[i]);

    f << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" " << "byte_order=\"LittleEndian\" ";
    f << ">" << endl;
    f << "<UnstructuredGrid>" << endl;
    uint nbpts = 0, nbcells = 0;
    for(uint i = 0; i < numObstacles; ++i)
    {
        nbpts += rb[i]->getConvex()->numberOfPoints_PARAVIEW();
        nbcells += rb[i]->getConvex()->numberOfCells_PARAVIEW();
    }
    f << "<Piece NumberOfPoints=\"" << nbpts << "\"" << " NumberOfCells=\"" << nbcells << "\">"
      << endl;
    f << "<Points>" << endl;
    f << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ";
    f << "format=\"ascii\">";
    f << endl;
    for(uint i = 0; i < numObstacles; ++i)
        rb[i]->getConvex()->writePoints_PARAVIEW(f, tr[i]);
    f << "</DataArray>" << endl;
    f << "</Points>" << endl;

    list<uint>           connectivity, offsets, cellstype;
    list<uint>::iterator ii;
    uint                 firstpoint_globalnumber = 0, last_offset = 0;
    for(uint i = 0; i < numObstacles; ++i)
        rb[i]->getConvex()->writeConnection_PARAVIEW(connectivity,
                                                     offsets,
                                                     cellstype,
                                                     firstpoint_globalnumber,
                                                     last_offset);
    f << "<Cells>" << endl;
    f << "<DataArray type=\"Int32\" Name=\"connectivity\" ";
    f << "format=\"ascii\">";
    f << endl;
    for(ii = connectivity.begin(); ii != connectivity.end(); ii++)
        f << *ii << " ";
    f << endl;
    f << "</DataArray>" << endl;
    f << "<DataArray type=\"Int32\" Name=\"offsets\" ";
    f << "format=\"ascii\">";
    f << endl;
    for(ii = offsets.begin(); ii != offsets.end(); ii++)
        f << *ii << " ";
    f << endl;
    f << "</DataArray>" << endl;
    f << "<DataArray type=\"Int32\" Name=\"types\" ";
    f << "format=\"ascii\">";
    f << endl;
    for(ii = cellstype.begin(); ii != cellstype.end(); ii++)
        f << *ii << " ";
    f << endl;
    f << "</DataArray>" << endl;
    f << "</Cells>" << endl;

    f << "<CellData Scalars=\"Indicator\">" << endl;
    f << "<DataArray type=\"Float32\" Name=\"Indicator\" ";
    f << "format=\"ascii\">";
    f << endl;
    for(uint i = 0; i < numObstacles; ++i)
    {
        // double indic = obstacleRB[i]->getIndicator();
        double indic = 0;
        int    nc    = rb[i]->getConvex()->numberOfCells_PARAVIEW();
        for(uint j = 0; j < nc; ++j)
            f << indic << " ";
    }
    f << endl;
    f << "</DataArray>" << endl;
    f << "</CellData>" << endl;

    f << "</Piece>" << endl;
    f << "</UnstructuredGrid>" << endl;
    f << "</VTKFile>" << endl;
    f.close();
}

// -------------------------------------------------------------------------------------------------
// Writes particles data
template <typename T>
void writeParticles_Paraview(const GrainsMemBuffer<RigidBody<T>*>&       rb,
                             const std::unique_ptr<ComponentManager<T>>& cm,
                             const std::string&                          parFileName)
{
    std::ofstream f(parFileName, std::ios::out);
    if(!f.is_open())
    {
        std::cerr << "Failed to open file: " << parFileName << std::endl;
        throw std::runtime_error("Cannot open file for writing: " + parFileName);
    }
    const uint                            numObstacles = cm->getNumberOfObstacles();
    const uint                            numParticles = cm->getNumberOfParticles();
    const GrainsMemBuffer<Vector3<T>>&    position     = cm->getPosition();
    const GrainsMemBuffer<Quaternion<T>>& quaternion   = cm->getQuaternion();
    const GrainsMemBuffer<Kinematics<T>>& kin          = cm->getVelocity();
    GrainsMemBuffer<Transform3<T>>        tr(numParticles);
    for(uint i = 0; i < numParticles; ++i)
        tr[i] = Transform3<T>(quaternion[numObstacles + i], position[numObstacles + i]);

    f << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" " << "byte_order=\"LittleEndian\" ";
    f << ">" << endl;
    f << "<UnstructuredGrid>" << endl;
    uint nbpts = 0, nbcells = 0;
    for(uint i = numObstacles; i < numObstacles + numParticles; ++i)
    {
        nbpts += rb[i]->getConvex()->numberOfPoints_PARAVIEW();
        nbcells += rb[i]->getConvex()->numberOfCells_PARAVIEW();
    }

    f << "<Piece NumberOfPoints=\"" << nbpts << "\"" << " NumberOfCells=\"" << nbcells << "\">"
      << endl;

    f << "<Points>" << endl;
    f << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ";
    f << "format=\"ascii\">" << endl;

    for(uint i = 0; i < numParticles; ++i)
    {
        rb[numObstacles + i]->getConvex()->writePoints_PARAVIEW(f, tr[i]);
    }
    f << "</DataArray>" << endl;
    f << "</Points>" << endl;

    list<uint>           connectivity, offsets, cellstype;
    list<uint>::iterator ii;
    uint                 firstpoint_globalnumber = 0, last_offset = 0;
    for(uint i = 0; i < numParticles; ++i)
        rb[numObstacles + i]->getConvex()->writeConnection_PARAVIEW(connectivity,
                                                                    offsets,
                                                                    cellstype,
                                                                    firstpoint_globalnumber,
                                                                    last_offset);
    f << "<Cells>" << endl;
    f << "<DataArray type=\"Int32\" Name=\"connectivity\" ";
    f << "format=\"ascii\">" << endl;
    for(ii = connectivity.begin(); ii != connectivity.end(); ii++)
        f << *ii << " ";
    f << endl;
    f << "</DataArray>" << endl;
    f << "<DataArray type=\"Int32\" Name=\"offsets\" ";
    f << "format=\"ascii\">" << endl;

    for(ii = offsets.begin(); ii != offsets.end(); ii++)
        f << *ii << " ";
    f << endl;

    f << "</DataArray>" << endl;
    f << "<DataArray type=\"Int32\" Name=\"types\" ";
    f << "format=\"ascii\">" << endl;

    for(ii = cellstype.begin(); ii != cellstype.end(); ii++)
        f << *ii << " ";
    f << endl;

    f << "</DataArray>" << endl;
    f << "</Cells>" << endl;

    f << "<CellData Scalars=\"NormU,NormOm,CoordNumb\">" << endl;

    f << "<DataArray type=\"Float32\" Name=\"NormU\" ";
    f << "format=\"ascii\">" << endl;

    for(uint i = numObstacles; i < numObstacles + numParticles; ++i)
    {
        T    normU = norm(kin[i].getTranslationalComponent());
        uint nc    = rb[i]->getConvex()->numberOfCells_PARAVIEW();
        for(uint j = 0; j < nc; ++j)
            f << normU << " ";
    }
    f << endl;
    f << "</DataArray>" << endl;

    f << "<DataArray type=\"Float32\" Name=\"NormOm\" ";
    f << "format=\"ascii\">" << endl;

    for(uint i = numObstacles; i < numObstacles + numParticles; ++i)
    {
        T    normOm = norm(kin[i].getAngularComponent());
        uint nc     = rb[i]->getConvex()->numberOfCells_PARAVIEW();
        for(uint j = 0; j < nc; ++j)
            f << normOm << " ";
    }
    f << endl;
    f << "</DataArray>" << endl;

    f << "<DataArray type=\"Float32\" Name=\"CoordNumb\" ";
    f << "format=\"ascii\">" << endl;

    for(uint i = numObstacles; i < numObstacles + numParticles; ++i)
    {
        T    coordNum = 0;
        uint nc       = rb[i]->getConvex()->numberOfCells_PARAVIEW();
        for(uint j = 0; j < nc; ++j)
            f << coordNum << " ";
    }
    f << endl;
    f << "</DataArray>" << endl;
    f << "</CellData>" << endl;
    f << "</Piece>" << endl;

    f << "</UnstructuredGrid>" << endl;
    f << "</VTKFile>" << endl;
    f.close();
}

/* ============================================================================================== */
/* High-Level Methods                                                                             */
/* ============================================================================================== */
// Default constructor
template <typename T>
ParaviewPostProcessingWriter<T>::ParaviewPostProcessingWriter()
{
}

// -------------------------------------------------------------------------------------------------
// Constructor with XML node, rank and number of processes as input parameters
template <typename T>

ParaviewPostProcessingWriter<T>::ParaviewPostProcessingWriter(DOMNode* dn)
{
    m_rootName  = ReaderXML::getNodeAttr_String(dn, "RootName");
    m_directory = ReaderXML::getNodeAttr_String(dn, "Directory");

    GoutWI(9, "Type = Paraview");
    GoutWI(12, "Output file directory name =", m_directory);
    GoutWI(12, "Output file root name =", m_rootName);
    // GoutWI(12, "Writing mode =", (m_binary ? "Binary" : "Text"));
}

// -------------------------------------------------------------------------------------------------
// Destructor
template <typename T>
ParaviewPostProcessingWriter<T>::~ParaviewPostProcessingWriter()
{
}

// ------------------------------------------------------------------------------------------------
// Gets the post-processing writer type
template <typename T>
PostProcessingWriterType ParaviewPostProcessingWriter<T>::getPostProcessingWriterType() const
{
    return (PARAVIEW);
}

// ------------------------------------------------------------------------------------------------
// Removes post-processing files already in the directory
template <typename T>
void ParaviewPostProcessingWriter<T>::clearPostProcessingFiles() const
{
    std::string              directory   = m_directory;
    std::vector<std::string> patternsStr = {"^" + m_rootName + R"(_.*\.pvd$)",
                                            "^" + m_rootName + R"(_.*\.vtu$)",
                                            "^" + m_rootName + R"(_.*\.pvtu$)",
                                            "^" + m_rootName + R"(_.*\.vtp$)",
                                            "^" + m_rootName + R"(_.*\.pvtp$)"};
    std::vector<std::regex>  patternsReg;
    for(const auto& pattern : patternsStr)
        patternsReg.push_back(std::regex(pattern));
    Gout("Removing Paraview post-processing files in", directory);
    PostProcessingWriter<T>::clearPostProcessingFiles(directory, patternsReg);
}

// -------------------------------------------------------------------------------------------------
template <typename T>
void ParaviewPostProcessingWriter<T>::PostProcessing_start()
{
    clearPostProcessingFiles();
    // Obstacles
    m_Paraview_saveObstacles_pvd << "<?xml version=\"1.0\"?>" << endl;
    m_Paraview_saveObstacles_pvd << "<VTKFile type=\"Collection\" version=\"0.1\""
                                 << " byte_order=\"LittleEndian\"";
    m_Paraview_saveObstacles_pvd << ">" << endl;
    m_Paraview_saveObstacles_pvd << "<Collection>" << endl;

    // Particles
    ostringstream* ossNULL = NULL;
    m_Paraview_saveParticles_pvd.reserve(1);
    m_Paraview_saveParticles_pvd.push_back(ossNULL);
    m_Paraview_saveParticles_pvd[0] = new ostringstream;
    *m_Paraview_saveParticles_pvd[0] << "<?xml version=\"1.0\"?>" << endl;
    *m_Paraview_saveParticles_pvd[0] << "<VTKFile type=\"Collection\" version=\"0.1\""
                                     << " byte_order=\"LittleEndian\"";
    *m_Paraview_saveParticles_pvd[0] << ">" << endl;
    *m_Paraview_saveParticles_pvd[0] << "<Collection>" << endl;
}

// -------------------------------------------------------------------------------------------------
// Writes data
template <typename T>
void ParaviewPostProcessingWriter<T>::PostProcessing(const GrainsMemBuffer<RigidBody<T>*>&       rb,
                                                     const std::unique_ptr<ComponentManager<T>>& cm,
                                                     const T currentTime)
{
    // list<string> Scalars;
    // Scalars.push_back("NormU");
    // Scalars.push_back("NormOm");
    // Scalars.push_back("CoordNumb");
    std::ostringstream ossCN;
    ossCN << m_ParaviewCycleNumber;

    // Obstacles
    std::string obsFileName     = m_rootName + "_Obstacles_T" + ossCN.str() + ".vtu";
    std::string obsFileNamePath = m_directory + "/" + obsFileName;
    m_Paraview_saveObstacles_pvd << "<DataSet timestep=\"" << currentTime << "\" "
                                 << "group=\"\" part=\"0\" file=\"" << obsFileName << "\"/>\n";

    std::string   obstacleFile = m_directory + "/" + m_rootName + "_Obstacles.pvd";
    std::ofstream f(obstacleFile, std::ios::out);
    if(!f.is_open())
    {
        std::cerr << "Failed to open file: " << obstacleFile << std::endl;
        throw std::runtime_error("Cannot open file for writing: " + obstacleFile);
    }
    f << m_Paraview_saveObstacles_pvd.str();
    f << "</Collection>" << endl;
    f << "</VTKFile>" << endl;
    f.close();
    writeObstacles_Paraview(rb, cm, obsFileNamePath);

    // Particles
    std::string parFileName     = m_rootName + "_Particles_T" + ossCN.str() + ".vtu";
    std::string parFileNamePath = m_directory + "/" + parFileName;
    *m_Paraview_saveParticles_pvd[0] << "<DataSet timestep=\"" << currentTime << "\" "
                                     << "group=\"\" part=\"0\" file=\"" << parFileName << "\"/>\n";

    std::string   particlesPvdFile = m_directory + "/" + m_rootName + "_Particles.pvd";
    std::ofstream g(particlesPvdFile, std::ios::out);
    if(!g.is_open())
    {
        std::cerr << "Failed to open file: " << particlesPvdFile << std::endl;
        throw std::runtime_error("Cannot open file for writing: " + particlesPvdFile);
    }
    g << m_Paraview_saveParticles_pvd[0]->str();
    g << "</Collection>" << endl;
    g << "</VTKFile>" << endl;
    g.close();

    writeParticles_Paraview(rb, cm, parFileNamePath);
    m_ParaviewCycleNumber++;
}

// ------------------------------------------------------------------------------------------------
// Finalizes writing data
template <typename T>
void ParaviewPostProcessingWriter<T>::PostProcessing_end()
{
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class ParaviewPostProcessingWriter<float>;
template class ParaviewPostProcessingWriter<double>;