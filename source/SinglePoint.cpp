#include "..\include\SinglePoint.h"
#include "..\include\VCLInput.h"
namespace VCLab
{
	void SinglePoint::Solver(ProjectInput VCLInput, vector<Phase> SysPhases)
	{
		int i;
		double temp, tempT;
		SPCond = VCLInput.Cond;
		Phn = 0;
		nele = SPCond.nele - 2;

		SPPhases = SysPhases;
		Phn = SPPhases.size();
		//Grid Point
		SPGGP.GeneGPy(SPPhases, SPCond);
		// Calc the equilibrium
		SPCalEq.Solver(1, SPPhases, SPCond);
		// Print the equilibrium

		WriteResults();

		/*
		for (i = 0;i < Phn;i++)
		{
			if(SPPhases[i].GP!=NULL)
				delete[] SPPhases[i].GP;
		}
		if (SPGGP.GP != NULL)
		delete[] SPGGP.GP;
		*/

	}

	void SinglePoint::WriteResults()
	{
		int i, j, kk;
		fstream VCLOut("VCLOutput.txt", ios::out);

		cout << endl << endl;
		PrintLine("=");
		PrintLineFile("=", VCLOut);

		cout << "Conditions:";
		VCLOut << "Conditions:";

		cout << "   T = " << SPCond.T
			<< "   P = " << SPCond.P
			<< "   N = " << SPCond.N << endl;
		VCLOut << "   T = " << SPCond.T
			<< "   P = " << SPCond.P
			<< "   N = " << SPCond.N << endl;
		//cout << "   Compositions:";

		for (i = 2; i < SPCond.nele; i++)
		{
			cout << "   X(" << SPCond.ele[i] << ") = " << SPCond.x[i] << endl;
			VCLOut << "   X(" << SPCond.ele[i] << ") = " << SPCond.x[i] << endl;
		}
		PrintLine("=");
		PrintLineFile("=", VCLOut);

		cout << "\nEquilibrium: " << endl;
		VCLOut << "\nEquilibrium: " << endl;
		cout << "\nChemical potential: " << endl;
		VCLOut << "\nChemical potential: " << endl;

		for (i = 2; i < SPCond.nele; i++)
		{
			cout << SPCond.ele[i] << "   " << SPCalEq.Chp[i] << endl;
			VCLOut << SPCond.ele[i] << "   " << SPCalEq.Chp[i] << endl;

		}
		cout << "\nPhase: " << endl;
		VCLOut << "\nPhase: " << endl;

		for (i = 0; i < SPCalEq.SPhases.size(); i++)
		{
			cout << SPCalEq.SPhases[i].name << endl;
			VCLOut << SPCalEq.SPhases[i].name << endl;

			cout << "Moles: " << SPCalEq.SPhases[i].vm << endl;
			VCLOut << "Moles: " << SPCalEq.SPhases[i].vm << endl;

			for (j = 0; j < nele; j++)
			{
				if (SPCalEq.SPhases[i].x[j]>0)
				{
					cout << "X(" << SPCond.ele[j + 2] << ") = " << SPCalEq.SPhases[i].x[j] << endl;
					VCLOut << "X(" << SPCond.ele[j + 2] << ") = " << SPCalEq.SPhases[i].x[j] << endl;

				}
			}
			if (debug >= 1)
			{
				cout << endl;
				kk = SPCalEq.SPhases[i].yids[0];
				for (j = 0; j < SPCalEq.SPhases[i].yn; j++)
				{
					if (SPCalEq.SPhases[i].yids[j] != kk) cout << "   :   ";
					cout << "  " << SPCalEq.SPhases[i].y[j];
					kk = SPCalEq.SPhases[i].yids[j];
				}
			}
			cout << endl;
			VCLOut << endl;
		}
		VCLOut.close();
	}
}