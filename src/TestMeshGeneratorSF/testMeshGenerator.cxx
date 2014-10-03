#include <Framework/ShockFittingAPI.hh>
#include <Framework/IOFunctions.hh>
#include <Framework/MeshGenerator.hh>

using namespace std;
using namespace ShockFitting;

int main (int argc, char** argv)
{
  int federateID = 0;
  cout << "### Creating coupling tool in federate [" << federateID << "]\n";
  SF_create_(&federateID);

  system("cp ../../../src/TestMeshGeneratorSF/input.case .");

  // link the I/O files 
  system("ln -sf ../../../src/TestMeshGeneratorSF/na00.1.node .");

  system("ln -sf ../../../src/TestMeshGeneratorSF/na00.1.poly .");

  system("ln -sf ../../../src/TestMeshGeneratorSF/na00.1.ele .");

  system("ln -sf ../../../src/TestMeshGeneratorSF/na00.1.neigh .");

  system("ln -sf ../../../src/TestMeshGeneratorSF/na00.1.edge .");

  string inputFile = "input.case";
  if (!fileExists(inputFile.c_str())) {
    cout << "ERROR: file <" << inputFile << "> does not exist in current directory!\n"; abort();
  }

  SF_configure_(&federateID, (char*)inputFile.c_str(), &argc, &argv);

  SF_process_(&federateID);


  cout << "### Destroy coupling tools in federate [" << federateID << "]\n";
  SF_destroy_(&federateID);

  return 1;
}

