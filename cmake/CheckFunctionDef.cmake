# try to compile a program
CHECK_CXX_SOURCE_COMPILES (
" #include <iostream>
  int main(int argc, char* argv[])
  {
    std::cout << __FUNCTION__ << std::endl;
    return 0;
  }"
  SF_HAVE_FUNCTION_DEF )
