#include <Framework/ShockFittingAPI.hh>
#include <Framework/IOFunctions.hh>

using namespace std;
using namespace ShockFitting;

int main (int argc, char** argv)
{
  int federateID = 0;
  cout << "### Creating coupling tool in federate [" << federateID << "]\n";
  SF_create_(&federateID);

  // (0) CircularCylinder_Unibas_N_inv_M20_2.0 : perfect gas
  // (1) CircularCylinder_VKI_LRD_2.1 : TCneq
  const unsigned nbTest = 4;
  vector<string> testDir(nbTest);
  testDir.at(0) = "CircularCylinder_Unibas_N_inv_M20_2.0";
  testDir.at(1) = "CircularCylinder_VKI_LRD_2.1";

  const unsigned i = 0; // number of executing test
  string pwdTestDir = "../../../src/TestStandardSF/"+testDir.at(i);

  string commandcp = "cp "+pwdTestDir+"/input.case .";
  system(commandcp.c_str());

  if(i==1) {
  string commandcp = "cp "+pwdTestDir+"/cylRDS.inter .";
  system(commandcp.c_str());
  }

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

  if(i==1) {
     commandln = "ln -sf "+pwdTestDir+"/nitrogen2.dat .";
     system(commandln.c_str());
  }
  commandln = "ln -sf "+pwdTestDir+"/sh00.dat .";
  system(commandln.c_str());

  // link the coolfluid files
  commandln = "ln -sf "+pwdTestDir+"/coolfluid-solver.xml .";
  system(commandln.c_str());
  commandln = "ln -sf "+pwdTestDir+"/cf00.CFcase .";
  system(commandln.c_str());

  string inputFile = "input.case";
  if (!fileExists(inputFile.c_str())) {
    cout << "ERROR: file <" << inputFile << "> does not exist in current directory!\n"; abort();
  }

  SF_configure_(&federateID, (char*)inputFile.c_str(), &argc, &argv);

cout << "-------------------------------------------" << endl;
cout << "--------------------------------------------" << endl;
cout << endl <<"Processing test " << testDir.at(i) << endl << endl;
cout << "-------------------------------------------" << endl;
cout << "-------------------------------------------" << endl;


  SF_process_(&federateID);


  cout << "### Destroy coupling tools in federate [" << federateID << "]\n";
  SF_destroy_(&federateID);

  return 1;
}


