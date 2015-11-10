/*
*   This file is a part of the VCLab software project.
*   For more details please contact xingwang@csu.edu.cn
*
*   Authors:    Xing Wang
*
*   Copyright (c) 2015-2015 Phase Diagram Center. Central South University. China
*
*   This program is free software: you can redistribute it and/or modify
*   it under the terms of the GNU General Public License as published by
*   the Free Software Foundation, either version 3 of the License, or
*   (at your option) any later version.
*
*   This program is distributed in the hope that it will be useful,
*   but WITHOUT ANY WARRANTY; without even the implied warranty of
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*   GNU General Public License for more details.
*
*/
//

#ifndef MinEnergy_H_
#define MinEnergy_H_

#include "VCLInput.h"
#include "Database.h"
#include "Tool.h"

namespace VCLab
{
	class MinEnergy
	{
	private:
		string classname;
		int linenumber;
		//
		int vn = 0, mvn = 0;  // totla varibles
		double var[MDim3], mvar[MDim];
		double eF[MDim3], meF[MDim3];
		double eFJ[MDim3][MDim3], meFJ[MDim3][MDim3];
		double deF[MDim3], dmeF[MDim3];

		int nph = 2;
		// chemical potential
		double Chp[MDim];
		int Chpidce[MDim];
		int Chpb[MDim];
		
		int nele = 5;
		// T
		int Tidce;

	public:
		//
		Condition MEConditions;

		//
		~MinEnergy();

		//
		void Solver(vector<Phase> & MEPhases, Condition Conditions, int nele, double(&Chp)[MDim], int(&chb)[MDim]);
		// 
		void SetupConditions(Condition Conditions);
		// 
		void SetupVariables(vector<Phase> & MEPhases);
		// 
		void AssignVariables(vector<Phase> & MEPhases);
		// Initialization phase, y, m, lamda
		void Initialization(vector<Phase> & MEPhases);
		// 
		void SetupHillert(vector<Phase> & MEPhases);
		void MendHillert(vector<Phase> & MEPhases);
		// mend Hillert according conditions
		void reMendHillert(vector<Phase> & MEPhases);
		//
		void SetupDrivingForce();
		//
		int Iteration(vector<Phase> & MEPhases);
		double CalcToteF(vector<Phase> & MEPhases);
		//
		void CalcCompositions(vector<Phase> & MEPhases);

		// Calc energy and derivate
		void CalcGE(Phase & Phases, double T);
		double CalcGF(Phase & phase, double T); 
		void CalcGPGF(Phase & phase, double T); 
		double CalcGFdx(Phase & Phases, double T, int idy);		
		double CalcGFdx2(Phase & Phases, double T, int idy1, int idy2);	
		// parameter value
		double CalcPara(int Tseg, Parameter para, double T);
		void CalcPara(double pvalue[3], int Tnseg, Parameter para, double T);
		double CalcGtao(double tao, double p, int order); 

		void ShowEquilibrium();
		void ShowCondition();
	};

} // end of VCLab

#endif