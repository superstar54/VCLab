#include "VCLInput.h"
#include "CalcEquilibrium.h"
namespace VCLab
{
	class LineSolver
	{
	private:
	protected:
	public:
		CalcEquilibrium				VCLCalEq;
		Condition					LCond;

		int Phn, Lineid;
		vector<Phase> oldSPhases;
		vector<Phase> LPhases;
		vector<Phase> LPhasesChage;

		double LV[1000], CLV[100];
		double LPhfrac[1000][10], CLPhfrac[100][10];
		double LPhx[1000][10][5], CLPhx[100][10][5];
		int LPhn[1000];
		int LPhid[1000][10];
		int LPhb[1000][10];

		void Solver(ProjectInput VCLInput, vector<Phase> SysPhases);
		void SaveEq(vector<Phase> SPhases);
		int PhasesChange(vector<Phase> Phases);
		void WriteResults();
		void WriteConditions();
		void SortResults(ProjectInput VCLInput, Database VCLDatabase);
		int FindPhase(string phasename, vector<Phase> Phases);
	};
}