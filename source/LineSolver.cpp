
#include "..\include\LineSolver.h"
#include "..\include\VCLInput.h"
namespace VCLab
{
	void LineSolver::Solver(ProjectInput VCLInput, vector<Phase> SysPhases)
	{
		int i;
		double temp, tempT;
		LCond = VCLInput.Cond;
		Phn = 0;
		Lineid = 0;
		//
		
		//double lvar = VCLInput.V_s;
		if (VCLInput.Var == "TEMPERATURE")
		{
			LCond.T = VCLInput.Var_s;
			while (LCond.T <= VCLInput.Var_e)
			{
				cout << LCond.T << endl;
				// Calc the equilibrium
				VCLCalEq.Solver(1, SysPhases, LCond);
				//Save
				SaveEq(VCLCalEq.SPhases);
				//
				if (Lineid>0)
				{
					if (PhasesChange(VCLCalEq.SPhases) == 10)
					{
						oldSPhases.clear();
						for (i = 0;i < VCLCalEq.SPhases.size();i++)
							oldSPhases.push_back(VCLCalEq.SPhases[i]);
						Lineid++;
						LCond.Tb = 1;
						tempT = LCond.T;
						VCLCalEq.Solver(0, LPhasesChage, LCond);
						LCond.T = VCLCalEq.CEConditions.T;
						cout << "Phases set change: " << LCond.T << endl;;

						SaveEq(VCLCalEq.SPhases);
						for (i = 0;i < 10;i++)
						{
							temp = LPhfrac[Lineid][i];
							LPhfrac[Lineid][i] = LPhfrac[Lineid - 1][i];
							LPhfrac[Lineid - 1][i] = temp;
						}
						temp = LV[Lineid];
						LV[Lineid] = LV[Lineid - 1];
						LV[Lineid - 1] = temp;
						LCond.Tb = 0;
						LCond.T = tempT;
					}
				}
				LCond.T += VCLInput.Var_d;
				oldSPhases.clear();
				for (i = 0;i < VCLCalEq.SPhases.size();i++)
					oldSPhases.push_back(VCLCalEq.SPhases[i]);
				Lineid++;
			}
			// Print the equilibrium
			WriteResults();
		}
		else
		{
			cout << "Other type of variblies is not inplemented now!" << endl;
		}
	}
	void LineSolver::SaveEq(vector<Phase> SPhases)
	{
		int i, j, phaseid;
		//
		for (i = 0;i < 10;i++)
		{
			LPhb[Lineid][i] = 0;
		}
		// Find phase id in Lphases
		LPhn[Lineid] = SPhases.size();
		for (i = 0; i < LPhn[Lineid]; i++)
		{
			phaseid = FindPhase(SPhases[i].name, LPhases);
			if (phaseid == -1)  // new Phases
			{
				LPhases.push_back(SPhases[i]);
				phaseid = Phn;
				Phn++;
			}
			LPhfrac[Lineid][phaseid] = SPhases[i].vm;
			LPhid[Lineid][phaseid] = phaseid;
			LPhb[Lineid][phaseid] = 1;
			LV[Lineid] = LCond.T;
		}
	}
	//
	int LineSolver::PhasesChange(vector<Phase> Phases)
	{
		int i, j, phaseid;
		double pro, sum;
		int PhCh = 0;

		LPhasesChage.clear();
		// Find phase id in Lphases
		for (i = 0; i < Phases.size(); i++)
		{
			phaseid = FindPhase(Phases[i].name, LPhasesChage);
			if (phaseid == -1)
			{
				LPhasesChage.push_back(Phases[i]);
			}
		}
		// Find phase id in Lphases
		for (i = 0; i < oldSPhases.size(); i++)
		{
			phaseid = FindPhase(oldSPhases[i].name, LPhasesChage);
			if (phaseid == -1)
			{
				LPhasesChage.push_back(oldSPhases[i]);
			}
		}
		//
		for (i = 0; i < Phn; i++)
		{
			pro = LPhfrac[Lineid - 1][i] * LPhfrac[Lineid][i];
			sum = LPhfrac[Lineid - 1][i] + LPhfrac[Lineid][i];
			if (pro<eps&&sum>eps)  // new Phases
			{
				phaseid = FindPhase(LPhases[i].name, LPhasesChage);
				if (phaseid == -1)
				{
					PrintError("error in finding the changed phase.");
				}
				else
				{
					LPhasesChage[phaseid].vmb = 0;
					LPhasesChage[phaseid].vm = 0;
					return PhCh = 1;
					break;
				}
			}
		}
		return PhCh;
	}
	//
	void LineSolver::SortResults(ProjectInput VCLInput, Database VCLDatabase)
	{

	}
	void LineSolver::WriteResults()
	{
		int i, j;
		ofstream Out("Phase fractions.txt");
		Out << "P = " << LCond.P
			<< ", N = " << LCond.N;;
		//cout << "   Compositions:";
		for (i = 2; i < LCond.nele; i++)
			Out << ", X(" << LCond.ele[i] << ") = " << LCond.x[i];
		Out << endl;
		Out << Phn << endl;
		Out << Lineid - 1 << endl;
		Out << 0 << "  " << "Temperature" << endl;
		for (i = 0;i < Phn; i++)
		{
			Out << i + 1 << "  " << LPhases[i].name << endl;;
		}
		for (i = 0;i < Lineid;i++)
		{
			Out << LV[i];
			for (j = 0;j < Phn;j++)
			{
				 Out<<"  "<< LPhfrac[i][j] ;
			}
			Out << endl;
		}
		Out.close();
	}


	int LineSolver::FindPhase(string phasename, vector<Phase> Phases)
	{
		int i = 0;
		int phaseid = -1;
		for (auto x : Phases)   // find in DBPhases
		{
			if (x.name == phasename) phaseid = i;
			i++;
		}
		// if phase not be found, return -1;
		return phaseid;
	}
}