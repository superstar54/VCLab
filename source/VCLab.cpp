<<<<<<< HEAD
/****************************************************************/
#include "..\include\VCLInput.h"

#include "..\include\Database.h"
#include "..\include\CalcEquilibrium.h"
#include "..\include\LineSolver.h"
#include "..\include\SinglePoint.h"
#include "..\include\EnergySurface.h"


using namespace std;
using namespace VCLab;

int main(int argc, char *argv[])
{
	ProjectInput				VCLInput;
	Database                    VCLDatabase;
	Condition					VCLCondition;
	CalcEquilibrium				VCLCalEq;
	MinEnergy					MinE;
	GlobalGrid					VCLGGP;

	string quit;
	string ProjectFileName;
	clock_t start, finish;
	double totaltime;

	PrintTitle();

	// Read the contral file, "VCLInput.txt"
	VCLInput.ReadVCLInput();
	VCLCondition = VCLInput.Cond;

    //******************       Read TDB   ***************************************//
	//FilTDB = "ALZNMG.TDB";
	VCLDatabase.readDatabase(VCLInput.FilTDB);
	VCLDatabase.showDatabase();
	//******************    Define System    ***********************************//
	//ElementsInputInput = "AL ZN MG";
	VCLDatabase.defineSystem(VCLInput.EleList, VCLInput.SelPhaList, 
		VCLInput.RejPhaList, VCLInput.FilTDB);
	VCLDatabase.showSystem();
	//*****************     Calculate Energy    **************************//
	
	//*****************     Calculate Equilibrium     **************************//


	start = clock();
	if (VCLInput.Mode == "EQUILIBRIUM")
	{
		if (VCLInput.Dimension == "0D")
		{
			SinglePoint SingleP;
			SingleP.Solver(VCLInput, VCLDatabase.SysPhases);
		}
		else if (VCLInput.Dimension == "1D")
		{
			LineSolver Line;
			Line.Solver(VCLInput, VCLDatabase.SysPhases);

		}
		else if (VCLInput.Dimension == "2D")
		{
			cout << "Map module is not inplemented now!" << endl;
		}
	}
	else if (VCLInput.Mode == "ENERGY_SURFACE")
	{
		EnergySurface VCLESurf;
		VCLESurf.Solver(1, VCLDatabase.SysPhases, VCLCondition);
	}
	else
	{
		cout << "This module: "<< VCLInput.Mode << " is not inplemented now!" << endl;
	}

	finish = clock();
	totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
	cout << "\n Running times:" << totaltime << "s." << endl;

	cout << "\n Enter any word to stop" << endl;
	cin >> quit;
	return 0;
}
=======
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
>>>>>>> b34a88e731ed6049f6c364b37b469ae02b913009
