#include<iostream> 
#include <string>  
#include"argparser/argparser.h"
using namespace std;
int main(int argc, const char ** argv)
{
	ArgParser a;
	const std::string optString = "stringTest";
	std::string stringVariable;
	std::string dftString = "default";
	const std::string helpString = "help String..........";

	const std::string optInt = "intTest";
	int intVariable;
	int dftInt = 0;
	const std::string helpInt = "help Int..........";

	const std::string optDouble = "doubleTest";
	double doubleVariable;
	double dftDouble = 0.0;
	const std::string helpDouble = "help Double..........";

	const  char* optBool = "boolTest";
	bool boolVariable;
	bool dftBool = false;
	const char* helpBool = "help Bool..........";

	a.refOption(optString, helpString, stringVariable, dftString, true);
	a.refOption(optInt, helpInt, intVariable, dftInt, true);
	a.refOption(optDouble, helpDouble, doubleVariable, dftDouble, true);
	a.refOption(optBool, helpBool, boolVariable, dftBool, true);
	
	a.showOptions();
	a.showUsages();
	a.run(argc, argv);
	cout << "---------------" << endl;
	cout << stringVariable << endl;
	cout << intVariable << endl;
	cout << doubleVariable << endl;
	cout << boolVariable << endl;
	cout << "---------------"<< endl;
	for (int i = 0; i < argc; i++)
	{
		cout << argv[i] << endl;
	}

	int count = 1;
	cin.get();
}