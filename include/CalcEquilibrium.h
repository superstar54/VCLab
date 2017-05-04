<<<<<<< HEAD
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

#ifndef CalcEquilibrium_H_
#define CalcEquilibrium_H_

#include "VCLInput.h"
#include "Database.h"
#include "MinEnergy.h"
#include "Tool.h"
#include "FindConvexHull.h"

namespace VCLab
{
	using namespace std;

	class CalcEquilibrium
	{
	private:
	public:
		string classname;

		// class
		MinEnergy MinE;
		FindCH    CHull;
		GlobalGrid GGP;

		Condition CEConditions;
		vector<Phase> CEPhases;
		vector<Phase> SPhases; // stable phases
		vector<Phase> MPhases; // metastable phases

		double GPdx = 0.10; // Global grid interval

		// Phase sets
		int Phn = 3;	// Phase number
		int Phns = 1;	// Stable Phases number
		int Phnm = 2;	// MetaStable Phases number
		int Phids[10];	// Stable phase id
		int Phidm[100]; //MetaStable phase id
		// chemical potential
		int nele = 2;		// elements number
		double Chp[MDim];	// chemical potential
		int Chpb[MDim] = { 1 };
		// Driving Force
		double DF[100]; 

		int GPmode;  // 1, use Grid Point, otherwise 0.

		clock_t start, finish;
		double totaltime;

		//
		~CalcEquilibrium();

		// setup global grid point
		void SetupGP();
		// setup conditions
		void SetupConditions(Condition VCLCondition);
		// setup phase sets
		void SetupPhaseSets(vector<Phase> SysPhases);
		// initialization for y, chp
		void Initialization(vector<Phase> SysPhases);
		// calc driving force
		void CalcDrivingForce();
		// iteration
		int Iteration();
		// solver phase set
		void Solver(int GPmod, vector<Phase> SysPhases, Condition CLCondition);
		//
		void Initialization();

		// Find convex face
		void FindConvexFace();
		// 
		void GetConvexFace(int Dim, Face &f, Point PA);
		// merge
		void MergeConvexFace(int Dim, Face &f);

		// print on screen
		void ShowEquilibrium();
		void ShowCondition();
		//
		int FindPhase(string phasename, vector<Phase> SPhases);
	};

} // end of VCLab

=======
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

#ifndef CalcEquilibrium_H_
#define CalcEquilibrium_H_

#include "VCLInput.h"
#include "Database.h"
#include "MinEnergy.h"
#include "Tool.h"
#include "FindConvexHull.h"

namespace VCLab
{
	using namespace std;

	class CalcEquilibrium
	{
	private:
		//class
		combination comb[10]; //最多10个亚点阵

		string classname;
		int linenumber;
		// Grid Point
		double GPdx = 0.10; // Global grid interval
		int GPn; 
		Point *GP;
		//Cartesian production
		int **input;
		int **output;
	public:
		// class
		MinEnergy MinE;
		FindCH    CHull;

		Condition CEConditions;
		vector<Species> CESpecies;
		vector<Species> CEElements;
		vector<Function> CEFunctions;
		vector<Phase> CEPhases;
		vector<Phase> SPhases; // stable phases
		vector<Phase> MPhases; // metastable phases

		// Phase sets
		int Phn = 3;	// Phase number
		int Phns = 1;	// Stable Phases number
		int Phnm = 2;	// MetaStable Phases number
		int Phids[10];	// Stable phase id
		int Phidm[100]; //MetaStable phase id
		// chemical potential
		int nele = 2;		// elements number
		double Chp[MDim];	// chemical potential
		int Chpb[MDim] = { 1 };
		// Driving Force
		double DF[100]; 

		int GPmode;  // 1, use Grid Point, otherwise 0.

		clock_t start, finish;
		double totaltime;

		//
		~CalcEquilibrium();

		// setup global grid point
		void SetupGP();
		// setup conditions
		void SetupConditions(Condition VCLCondition);
		// setup phase sets
		void SetupPhaseSets(vector<Phase> SysPhases);
		// initialization for y, chp
		void Initialization(vector<Phase> SysPhases);
		// calc driving force
		void CalcDrivingForce();
		// iteration
		int Iteration();
		// solver phase set
		void Solver(int GPmod, vector<Phase> SysPhases, Condition CLCondition);
		//
		void Initialization();

		// Global grid point
		// total GP number
		void CalcGPn();
		//Grid Point, include x, y, energy
		void GeneGP();
		// Find convex face
		void FindConvexFace();
		// 
		void GetConvexFace(int Dim, Face &f, Point PA);
		// merge
		void MergeConvexFace(int Dim, Face &f);

		// print on screen
		void ShowEquilibrium();
		void ShowCondition();
		//
		int FindPhase(string phasename, vector<Phase> SPhases);
	};

} // end of VCLab

>>>>>>> b34a88e731ed6049f6c364b37b469ae02b913009
#endif