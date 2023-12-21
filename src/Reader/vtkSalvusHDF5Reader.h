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

  vtkGetObjectMacro(PointDataArraySelection, vtkDataArraySelection);

  int GetNumberOfPointArrays();
  const char* GetPointArrayName(int index);
  int GetPointArrayStatus(const char* name);
  void SetPointArrayStatus(const char* name, int status);

  void DisableAllPointArrays();  
  void EnableAllPointArrays();
  void DisablePointArray(const char* name);  
  void EnablePointArray(const char* name);

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
  vtkDataArraySelection* PointDataArraySelection;
  int NbNodes;
  int NbCells;
 private:
  vtkSalvusHDF5Reader(const vtkSalvusHDF5Reader&) = delete;
  void operator=(const vtkSalvusHDF5Reader&) = delete;

  int ModelName; // 0 = ELASTIC, 1 = ACOUSTIC
  double *LoadStepScalingFactors;
  int NumLoadSteps;
};

#endif
