/****************************************************************/
#include "..\include\VCLInput.h"

#include "..\include\Database.h"
#include "..\include\CalcEquilibrium.h"
#include "..\include\LineSolver.h"

using namespace std;
using namespace VCLab;

int main(int argc, char *argv[])
{
	ProjectInput				VCLInput;
	Database                    VCLDatabase;
	Condition					VCLCondition;
	CalcEquilibrium				VCLCalEq;
	MinEnergy					MinE;

	string quit;
	string ProjectFileName;
	clock_t start, finish;
	double totaltime;

	PrintTitle();

	// Read the contral file, "VCLInput.txt"
	VCLInput.ReadVCLInput();

    //******************       Read TDB   ***************************************//
	//FilTDB = "ALZNMG.TDB";
	VCLDatabase.readDatabase(VCLInput.FilTDB);
	VCLDatabase.showDatabase();
	//******************    Define System    ***********************************//
	//ElementsInputInput = "AL ZN MG";
	VCLDatabase.defineSystem(VCLInput.EleList, VCLInput.FilTDB);
	VCLDatabase.showSystem();
	//*****************     Calculate Energy    **************************//

	//*****************     Calculate Equilibrium     **************************//


	start = clock();
	VCLCondition = VCLInput.Cond;
	if (VCLInput.VCLmode == "POINT")
	{
		// Calc the equilibrium
		VCLCalEq.Solver(1, VCLDatabase.SysPhases, VCLCondition);
		// Print the equilibrium
		VCLCalEq.ShowEquilibrium();
	}
	else if(VCLInput.VCLmode == "LINE")
	{
		LineSolver Line;
		Line.Solver(VCLInput, VCLDatabase.SysPhases);
		
	}
	else if (VCLInput.VCLmode == "MAP")
	{
		cout << "Map module is not inplemented now!" << endl;
	}
	finish = clock();
	totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
	cout << "\n Running times:" << totaltime << "s." << endl;

	cout << "Enter any key to stop" << endl;
	cin >> quit;
	return 0;
}
