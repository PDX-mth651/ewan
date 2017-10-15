
#include "CSRclass.h"
using namespace std;


int main()
{
  string filePath,destination,name,answer;
  bool keepGoing = true;
  while (keepGoing)
  {
    bool dataInFile = false;
    bool useData = false;
    cout << "File must either consist of lines of the form 'i j', with i and j positive integers, or 'i j d', with data d a numeric type with 0,1,multiplication and addition. \n\r";
    cout << "Enter file path:\n\r";
    cin >> filePath;
    cout << "Does you file include data?\n\r";
    cin >> answer;
    if (answer=="yes")
    {
      dataInFile = true;
    }
    if (dataInFile)
    {
      cout << "Should I use this data?  If not, then each non zero entry will be assumed to be 1.\n\r";
      cin >> answer;
      if (answer == "yes")
      {
        useData = true;
      }
    } 
    cout << "importing\n\r";
    const CSRMat<double> bob(filePath,dataInFile, useData);
  
    print(bob);
    destroy(bob);
    cout << "load another matrix?\n\r";
    cin >> answer;
    if (answer!="yes")
    {
      keepGoing=false;
    }  
      
  }
  

  cout<<"Finished\n\r";

}