#include "VCLInput.h"
#include "CalcEquilibrium.h"
namespace VCLab
{
	class SinglePoint
	{
	private:
	protected:
	public:
		CalcEquilibrium				SPCalEq;
		Condition					SPCond;
		GlobalGrid					SPGGP;

		int Phn;
		int nele;
		vector<Phase> SPPhases;

		void Solver(ProjectInput VCLInput, vector<Phase> SysPhases);
		void SaveEq(vector<Phase> SPhases);
		int PhasesChange(vector<Phase> Phases);
		void WriteResults();
		void SortResults(ProjectInput VCLInput, Database VCLDatabase);
		int FindPhase(string phasename, vector<Phase> Phases);
	};
}