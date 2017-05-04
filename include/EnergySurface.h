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

#ifndef EnergySurface_H_
#define EnergySurface_H_

#include "VCLInput.h"
#include "Database.h"
#include "MinEnergy.h"
#include "Tool.h"
#include "FindConvexHull.h"

namespace VCLab
{
	using namespace std;

	class EnergySurface
	{
	private:
		
	public:
		string classname;

		// class
		MinEnergy MinE;
		FindCH    CHull;
		GlobalGrid GGP;

		Condition CEConditions;
		vector<Phase> ESPhases;
		vector<Phase> MinEPhases;

		int GPdxn = 0;
		int Bn = 1;
		int **BGPn;
		double GPdx = 0.10; // Global grid interval
		int *GPidB;
		Point **BGP;
		int *EqBGPPhn;
		Point **EqBGP;
		int *TEqBGPn;
		Point **TEqBGP;

		//
		int GPn;
		int fn;		//
		int GPfn;
		int TGPfn;
		Face *GPf;		//
		Face *TGPf;	//
		Face f;

		double mat[MDim3][MDim3];

		int Phn = 0;	// Phase number
		
		int nele = 0;		// elements number
		double Chp[MDim];	// chemical potential
		int Chpb[MDim] = { 1 };  // chemical potential

		int GPmode = 1;  // 1, use Grid Point, otherwise 0.


		//
		~EnergySurface();

		// setup conditions
		void SetupConditions(Condition VCLCondition);
		// solver phase set
		void Solver(int GPmod, vector<Phase> SysPhases, Condition CLCondition);

		// divide GP into block
		void BlockGP(int Phid, Phase Phases);
		
		//
		void FindNetFace(int in_GPn, Point *GP);
		void FindGFace(int GPn, Point *GP, Face &f);
		void RecurNetFace(int GPn, Point *GP, Face f);
		int FindNewP(int GPn, Point *GP, Face f);
		void GetNetFace();

		bool inside(Point p, Face f);
		void inside_edge(int &subedge, int parray[MDim], Point p, Face f);
		double distance(int mode, int n, Point p, Face f);
		double calcdet(int n, double array[MDim3][MDim3]);
		void deletep();
		// print on screen
		void ShowPoint(string str, int n, Point p);
		void ShowFace(string str, int n, int m, Face f);
		void ShowNetFace();
		void WriteResults();
		void WriteFaces();
	};

} // end of VCLab

#endif