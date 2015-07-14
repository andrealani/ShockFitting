#include <Framework/ShockFittingAPI.hh>
#include <Framework/IOFunctions.hh>

using namespace std;
using namespace ShockFitting;

int main (int argc, char** argv)
{
  int federateID = 0;

  cout << "### Creating coupling tool in federate [" << federateID << "]\n";
  /// create the shock fitting for the given federate
  /// @param federateID   ID of the federate (< NFEDERATES)
  SF_create_(&federateID);

  if (!fileExists(argv[1])) {
    cout << "ERROR: file <" << argv[1] << "> does not exist in current directory!\n"; abort();
  }

  /// configure the shock fitting
  /// @param federateID      ID of the federate (< NFEDERATES)
  /// @param inputFile       configuration filename or string
  /// @param argc        number of command line options
  /// @param argv        command line options
  SF_configure_(&federateID, (char*)argv[1], &argc, &argv);

  /// run the shock fitting 
  /// @param federateID   ID of the federate (< NFEDERATES)
  SF_process_(&federateID);


  cout << "### Destroy coupling tools in federate [" << federateID << "]\n";
  /// finalize the shock fitting for the given federate
  /// @param federateID    ID of the federate (< NFEDERATES)
  SF_destroy_(&federateID);

  return 1;
}


