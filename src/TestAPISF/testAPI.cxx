#include <Framework/ShockFittingAPI.hh>
#include <Framework/IOFunctions.hh>

using namespace std;
using namespace ShockFitting;

int main(int argc, char** argv)
{
  int federateID = 0;
  cout << "### Creating coupling tool in federate [" << federateID << "]\n";
  SF_create_(&federateID);
  
  system("cp ../../../src/TestAPISF/input.case .");
 
  string inputFile = "input.case";
  if (!fileExists(inputFile.c_str())) {
    cout << "ERROR: file <" << inputFile << "> does not exist in current directory!\n"; abort(); 
  }
  
  SF_configure_(&federateID, (char*)inputFile.c_str(), &argc, &argv);
  
  int inNbElems = 0;
  int inNbStates = 0;
  int inStateStride = 0;
  int outNbElems = 0;
  int outNbStates = 0;
  int outStateStride = 0;
  vector<int> inElementState;  
  vector<int> inElementStatePtr;
  vector<int> outElementState; 
  vector<int> outElementStatePtr;
  vector<double> inStateField;
  vector<double> outStateField;
  
  cout << "### Interpolate field in federate [" << federateID << "]\n";
  SF_interpolate_field_(&federateID,
			(char*)"DummyFieldInterpolator",
			&inNbElems,
			&inNbStates,
			&inStateStride,
			&inElementState[0], 
			&inElementStatePtr[0],
			&inStateField[0],
			&outNbElems,
			&outNbStates,
			&outStateStride,
			&outElementState[0], 
			&outElementStatePtr[0],
			&outStateField[0]);
  
  cout << "### Process file in federate [" << federateID << "]\n";
  SF_process_file_(&federateID, (char*)"DummyFileProcessing");
  
  cout << "### Destroy coupling tools in federate [" << federateID << "]\n";
  SF_destroy_(&federateID);
  
  return 1;
}
