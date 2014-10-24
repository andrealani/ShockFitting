#include <Framework/ShockFittingAPI.hh>
#include <Framework/IOFunctions.hh>

using namespace std;
using namespace ShockFitting;

int main (int argc, char** argv)
{
  int federateID = 0;
  cout << "### Creating coupling tool in federate [" << federateID << "]\n";
  SF_create_(&federateID);

  // (0) CircularCylinder_UnDiFiv_1.2
  // (1) CircularCylinder_VKI_LRD_UnDiFiv_2.1 
  // (2) CircularCylinder_VKI_UnDiFiv_2.1 
  const unsigned nbTest = 4;
  vector<string> testDir(nbTest);
  testDir.at(0) = "CircularCylinder_1.2";
  testDir.at(1) = "CircularCylinder_VKI_2.1";
  testDir.at(2) = "CircularCylinder_VKI_LRD_2.1";
  testDir.at(3) = "CircularCylinder_VKI_LDA_2.1";

  const unsigned i = 2; // number of executing test
  string pwdTestDir = "../../../src/TestStandardSF/"+testDir.at(i);

  string commandcp = "cp "+pwdTestDir+"/input.case .";
  system(commandcp.c_str());

  // link the I/O files 
  string commandln = "ln -sf "+pwdTestDir+"/na00.1.node .";
  system(commandln.c_str());
  commandln = "ln -sf "+pwdTestDir+"/na00.1.poly .";
  system(commandln.c_str());
  commandln = "ln -sf "+pwdTestDir+"/na00.1.ele .";
  system(commandln.c_str());
  commandln = "ln -sf "+pwdTestDir+"/na00.1.neigh .";
  system(commandln.c_str());
  commandln = "ln -sf "+pwdTestDir+"/na00.1.edge .";
  system(commandln.c_str());

  if(i==1 || i==2 || i==3) {
     commandln = "ln -sf "+pwdTestDir+"/nitrogen2.dat .";
     system(commandln.c_str());
  }
  commandln = "ln -sf "+pwdTestDir+"/sh00.dat .";
  system(commandln.c_str());

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


