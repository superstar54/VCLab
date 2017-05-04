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

#ifndef GlobalGrid_H_
#define GlobalGrid_H_

#include "VCLInput.h"
#include "Database.h"
#include "MinEnergy.h"
#include "Tool.h"

namespace VCLab
{
	using namespace std;

	

	class GlobalGrid
	{
	private:
		// combinations
		combination GGComb[10]; //最多10个亚点阵
		//Cartesian production
		
	public:
		string classname;

		// class
		MinEnergy MinE;

		// Grid Point
		int GPn;
		Point *GP;

		int Phn;	// Phase number
		int nele;
		double GPdx = 0.05; // Global grid interval

		//
		~GlobalGrid();
		
		//Grid Point, include x, y, 
		void GeneGPy(vector<Phase> &CEPhases, Condition CEConditions);
		// total GP number
		void CalcGPn(vector<Phase> &CEPhases);
		//energy
		void GeneGPGF(vector<Phase> &CEPhases, Condition CEConditionss);

		//
		void ShowPoint(string str, int n, Point p);
	};

} // end of VCLab

#endif