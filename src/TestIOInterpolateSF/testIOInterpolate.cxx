#include <ctime>
#include <Framework/ShockFittingAPI.hh>
#include <Framework/IOFunctions.hh>

using namespace std;
using namespace ShockFitting;

int main(int argc, char** argv)
{  
  int federateID = 0;
  cout << "### Create coupling tool in federate [" << federateID << "]\n";
  SF_create_(&federateID);

  // copy the input file for the coupling tools   
  system("cp ../../../src/TestIOInterpolateSF/input.case .");
  // link the I/O (iPIC3D in this example) file 
  system("ln -sf ../../../src/TestIOInterpolateSF/ipic3d.out.0 .");
  // unpack the donor solution file 
  system("bunzip2 ../../../src/TestIOInterpolateSF/output-P7-iter_0.CFmesh.bz2"); 
  // link the donor file
  system("ln -sf ../../../src/TestIOInterpolateSF/output-P7-iter_0.CFmesh .");
  // some link-based  trickery to reuse existing I/O files
  system("ln -sf ipic3d.out.0 ipic3d.out.100");
  system("ln -sf output-P7-iter_0.CFmesh output-P7-iter_100.CFmesh"); 
 
  string inputFile = "input.case";
  if (!fileExists(inputFile.c_str())) {
    cout << "ERROR: file <" << inputFile << "> does not exist in current directory!\n"; abort(); 
  }
  
  SF_configure_(&federateID, (char*)inputFile.c_str(), &argc, &argv);
  
  clock_t start = clock();  
  cout << "### Process coupling tools in federate [" << federateID << "] coupling iter 0\n";
  SF_process_(&federateID);
  clock_t finish = clock();
  cout << "### Process coupling tools took " << float(finish-start)/CLOCKS_PER_SEC << " seconds\n";
  
  start = clock();  
  cout << "### Process coupling tools in federate [" << federateID << "] coupling iter 1\n";
  SF_process_(&federateID);
  finish = clock();
  cout << "### Process coupling tools took " << float(finish-start)/CLOCKS_PER_SEC << " seconds\n";
  
  cout << "### Destroy coupling tools in federate [" << federateID << "]\n";
  SF_destroy_(&federateID);
  
  return 1;
}
