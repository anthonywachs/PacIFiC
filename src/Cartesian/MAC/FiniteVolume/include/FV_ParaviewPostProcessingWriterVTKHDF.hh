#ifndef FV_PARAVIEW_POSTPROCESSING_WRITER_VTKHDF_HH
#define FV_PARAVIEW_POSTPROCESSING_WRITER_VTKHDF_HH

#include <FV_PostProcessingWriter.hh>

#include <list>
#include <string>
#include <vector>
using std::list;
using std::string;

class MAC_Communicator;
class MAC_ModuleExplorer;
class FV_DiscreteField;
class FV_Mesh;
class FV_TimeIterator;

/**
 * @brief The class ParaviewPostProcessingWriterVTKHDF
 *
 * @author C. Olive 2026
 */
class FV_ParaviewPostProcessingWriterVTKHDF : public FV_PostProcessingWriter {
public: //-----------------------------------------------------------
  void write_cycle(FV_TimeIterator const *t_it, size_t cycle_number);

  //-- Data clearing

  void clearResultFiles(void);

  size_t getPreviousCycleNumber(void);

  void readTimeFile(FV_TimeIterator const *t_it, size_t &cycle_number);

protected: //--------------------------------------------------------
private:   //----------------------------------------------------------
  FV_ParaviewPostProcessingWriterVTKHDF(void);

  FV_ParaviewPostProcessingWriterVTKHDF(MAC_Object *a_owner,
                                        MAC_ModuleExplorer const *exp,
                                        MAC_Communicator const *com,
                                        list<FV_DiscreteField const *> a_fields,
                                        FV_Mesh const *a_primary_mesh,
                                        bool a_binary);

  ~FV_ParaviewPostProcessingWriterVTKHDF(void);

  FV_ParaviewPostProcessingWriterVTKHDF(
      FV_ParaviewPostProcessingWriterVTKHDF const &other);

  FV_ParaviewPostProcessingWriterVTKHDF &
  operator=(FV_ParaviewPostProcessingWriterVTKHDF const &other);

  FV_ParaviewPostProcessingWriterVTKHDF(MAC_Object *a_owner);

  void write_imagedata_vtkhdf(std::string file, int compression_level, bool overwrite);
  void write_vtkhdf_series();
  void write_vtkhdf_series_idx();
  void read_vtkhdf_series_idx();

  //-- Plug in

  FV_ParaviewPostProcessingWriterVTKHDF(std::string const &a_name);

  FV_PostProcessingWriter *
  create_replica(MAC_Object *a_owner, MAC_ModuleExplorer const *exp,
                 MAC_Communicator const *com,
                 list<FV_DiscreteField const *> a_fields,
                 FV_Mesh const *a_primary_mesh, bool a_binary) const;

  //-- Class attributes

  static FV_ParaviewPostProcessingWriterVTKHDF const *PROTOTYPE;

  //-- Attributes

  MAC_ModuleExplorer const *EXP;
  MAC_Communicator const *COM;
  list<FV_DiscreteField const *> FVFIELDS;
  FV_Mesh const *PRIMARY_GRID;

  std::string RES_DIRECTORY;
  std::string BASE_FILENAME;
  std::string VTKHDF_SERIES_FILENAME;
  std::string VTKHDF_SERIES_IDX_FILENAME;

  struct TimeSeriesRecord {
    double time;
    size_t cycle;
    std::string file;
  };
  std::vector<TimeSeriesRecord> TIME_SERIES_RECORDS;

  bool TIME_SERIES_LOADED;
  bool TIME_SERIES_RESTART_MODE;
  size_t TIME_SERIES_PREVIOUS_CYCLE;

  size_t CYCLE_NUMBER;
  bool BINARY;
};

#endif
