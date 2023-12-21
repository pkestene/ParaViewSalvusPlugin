#include "vtkAutoInit.h"

//#include "vtkDataSet.h"
//#include "vtkDataSetWriter.h"
#include "vtkSalvusHDF5Reader.h"
#include "vtkInformation.h"
#include "vtkNew.h"
#include "vtkPointData.h"
#include "vtkUnstructuredGrid.h"

#include <vtksys/SystemTools.hxx>
#include <vtksys/CommandLineArguments.hxx>

int
vtkIOSalvusCxxTests(int argc, char **argv)
{
  std::string filein;
  std::string varname;

  double TimeStep = 4.2898e-05;

  vtksys::CommandLineArguments args;
  args.Initialize(argc, argv);
  args.AddArgument(
    "-f", vtksys::CommandLineArguments::SPACE_ARGUMENT, &filein, "(the names of the Gadget (HDF5) files to read)");
  args.AddArgument(
    "-var", vtksys::CommandLineArguments::SPACE_ARGUMENT, &varname, "(the name of the SCALAR variable to display)");

  if ( !args.Parse() || argc == 1 || filein.empty())
    {
    cerr << "\nTestSalvusHDF5Reader: Written by Jean M. Favre\n"
         << "options are:\n";
    cerr << args.GetHelp() << "\n";
    return EXIT_FAILURE;
    }

  if(!vtksys::SystemTools::FileExists(filein.c_str()))
    {
    cerr << "\nFile " << filein.c_str() << " does not exist\n\n";
    return EXIT_FAILURE;
    }

  vtkNew<vtkSalvusHDF5Reader> reader;
  reader->DebugOff();
  reader->SetFileName(filein.c_str());
  reader->UpdateInformation();
  
  for(auto i=0; i < reader->GetNumberOfPointArrays(); i++)
    cout << "found array (" << i << ") = " << reader->GetPointArrayName(i) << endl;
  reader->DisableAllPointArrays();
  reader->SetPointArrayStatus("stress_xx", 1);

  reader->UpdateTimeStep(TimeStep); // time value
  reader->Update();

  double range[2];
      
  if(varname.size())
    {
    reader->GetOutput()->GetPointData()->GetArray(0)->GetRange(range);
    cerr << varname.c_str() << ": scalar range = [" << range[0] << ", " << range[1] << "]\n";
    }
  //cout << *reader;
  /*
  VTK_CREATE(vtkDataSetWriter, writer);
  writer->SetInputData(FirstBlock);
  writer->SetFileTypeToBinary();
  writer->SetFileName("/tmp/foo.vtk");
  writer->Write();
  */
  
  return EXIT_SUCCESS;
}

int
main(int argc, char **argv)
{
  vtkIOSalvusCxxTests(argc, argv);
}

