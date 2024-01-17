/*=========================================================================
// .NAME vtkSalvusHDF5Reader - read unstructured grid data 
// .SECTION Description
// vtkSalvusHDF5Reader is a source object that reads HDF5 
// unstructured grid data files.
*/
#ifndef __vtkSalvusHDF5Reader_h
#define __vtkSalvusHDF5Reader_h

#include "SalvusHDF5ReaderModule.h" // for export macro
#include "vtkUnstructuredGridAlgorithm.h"

#include <vector> // for vector

class vtkDataArraySelection;

class SALVUSHDF5READER_EXPORT vtkSalvusHDF5Reader : public vtkUnstructuredGridAlgorithm
{
public:
  static vtkSalvusHDF5Reader *New();
  vtkTypeMacro(vtkSalvusHDF5Reader, vtkUnstructuredGridAlgorithm);

  vtkSetMacro(ModelName, int);
  vtkGetMacro(ModelName, int);

  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  vtkGetMacro(NbCells,int);
  vtkGetMacro(NbNodes,int);

  vtkGetObjectMacro(ELASTIC_PointDataArraySelection, vtkDataArraySelection);
  vtkGetObjectMacro(ACOUSTIC_PointDataArraySelection, vtkDataArraySelection);

  int GetNumberOfPointArrays();
  const char* GetPointArrayName(int index);
  int GetPointArrayStatus(const char* name);
  void SetPointArrayStatus(const char* name, int status);
  void DisableAllPointArrays();  
  void EnableAllPointArrays();
  void DisablePointArray(const char* name);  
  void EnablePointArray(const char* name);
  
  int GetNumberOf_Acoustic_PointArrays();
  const char* Get_Acoustic_PointArrayName(int index);
  int Get_Acoustic_PointArrayStatus(const char* name);
  void Set_Acoustic_PointArrayStatus(const char* name, int status);
  void DisableAll_Acoustic_PointArrays();  
  void EnableAll_Acoustic_PointArrays();
  void Disable_Acoustic_PointArray(const char* name);  
  void Enable_Acoustic_PointArray(const char* name);

protected:
  vtkSalvusHDF5Reader();
  ~vtkSalvusHDF5Reader();

  int RequestData(vtkInformation*,
                  vtkInformationVector**,
                  vtkInformationVector*);
  int RequestInformation(vtkInformation*,
                         vtkInformationVector**,
                         vtkInformationVector*);

  int CanReadFile(const char* fname);
  
  char *FileName;
  vtkDataArraySelection* ELASTIC_PointDataArraySelection;
  vtkDataArraySelection* ACOUSTIC_PointDataArraySelection;
  int NbNodes;
  int NbCells;
 private:
  vtkSalvusHDF5Reader(const vtkSalvusHDF5Reader&) = delete;
  void operator=(const vtkSalvusHDF5Reader&) = delete;

  int ModelName; // 0 = ELASTIC, 1 = ACOUSTIC
  std::vector<double> TimeStepValues;
  int NumberOfTimeSteps;
  int TimeStep;
  int ActualTimeStep;
  double TimeStepTolerance;
};

#endif
