#include <Framework/ShockFittingAPI.hh>
#include <Framework/IOFunctions.hh>

using namespace std;
using namespace ShockFitting;

int main (int argc, char** argv)
{
  int federateID = 0;
  cout << "### Creating coupling tool in federate [" << federateID << "]\n";
  SF_create_(&federateID);

  // (0) CircularCylinder_Pg_inv_N_M15  : perfect gas M=15 inviscid
  // (1) CircularCylinder_Pg_inv_N_M20  : perfect gas M=20 inviscid
  // (2) CircularCylinder_Pg_inv_N_M25  : perfect gas M=25 inviscid
  // (3) CircularCylinder_VKI_LRD_2.1   : TCneq       M6   inviscid
  // (4) CircularCylinder_Pg_vis_Bx_M17 : perfect gas M=17 viscous
  //                                      from the folder vki
  // (5) CircularCylinder_Pg_vis_Bx_M17 : perfect gas M=6 inviscid
  const unsigned nbTest = 10;
  vector<string> testDir(nbTest);
  testDir.at(0) = "CircularCylinder_Pg_inv_N_M15";
  testDir.at(1) = "CircularCylinder_Pg_inv_N_M20";
  testDir.at(2) = "CircularCylinder_Pg_inv_N_M25";
  testDir.at(3) = "CircularCylinder_VKI_LRD_2.1";
  testDir.at(4) = "CircularCylinder_Pg_inv_N_M6";
  testDir.at(5) = "CircularCylinder_Pg_vis_Bx_M17";

  // number of executing test
  const unsigned i = 4;

  string pwdTestDir = "../../../src/TestStandardSF/"+testDir.at(i);

  string commandcp = "cp "+pwdTestDir+"/input.case .";
  system(commandcp.c_str());

  if(i==3 || i==4) {
  string commandcp = "cp "+pwdTestDir+"/cylRDS.inter .";
  system(commandcp.c_str());
  }

  if(i==5) {
  string commandcp = "cp "+pwdTestDir+"/cyl.inter .";
  system(commandcp.c_str());
  }

  // link the I/O files 
  string commandln;
  commandln = "cp -rf "+pwdTestDir+"/na00.1.node .";
  system(commandln.c_str());
  commandln = "cp -rf "+pwdTestDir+"/na00.1.poly .";
  system(commandln.c_str());
  commandln = "cp -rf "+pwdTestDir+"/na00.1.ele .";
  system(commandln.c_str());
  commandln = "cp -rf "+pwdTestDir+"/na00.1.neigh .";
  system(commandln.c_str());
  commandln = "cp -rf "+pwdTestDir+"/na00.1.edge .";
  system(commandln.c_str());
  commandln = "cp -rf "+pwdTestDir+"/sh00.dat .";
  system(commandln.c_str());

  // link the starting captured solution
/*  commandln = "cp -rf "+pwdTestDir+"/StartCapturedSolution/CFresults/cyl_TCneq_noChem-P0.CFmesh .";
  system(commandln.c_str());
  commandln = "cp -rf "+pwdTestDir+"/StartCapturedSolution/CFresults/shock.dat .";
  system(commandln.c_str());
*/
  // link the chemical info file
  if(i==3) {
     commandln = "ln -sf "+pwdTestDir+"/nitrogen2.dat .";
     system(commandln.c_str());
  }

  // link the coolfluid files
  commandln = "cp -rf " + pwdTestDir + "/coolfluid-solver.xml .";
  system(commandln.c_str());
  commandln = "cp -rf "+pwdTestDir+"/cf00.CFcase .";
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


