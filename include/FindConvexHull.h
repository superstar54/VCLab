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
#ifndef FINDCONVEXHULL_H_
#define FINDCONVEXHULL_H_

#include"Tool.h"
#include"MinEnergy.h"

namespace VCLab
{
	struct Point
	{
		double x[MDim + 1];	// compositions
		double y[MDim * 3];	// yfractions
		int phid;			// belong to id phase
		double phfrac;		// phase farctions
		int fid;			// id in the convex face
		int fb;				// ture or not
	};

	struct Face
	{
		int b;
		int id;
		Point p[MDim];		// point
		double cp[MDim];	// chemical potential

	};

	class FindCH
	{
	private:
		int Dim;  // dimension == nele
		int GPn;  // 
		int new_GP;
		//Point *new_GP; 
		int GPfn = 0; //
		Face *GPf;  // convex face

		int ncountp = 0;
		Face subf;
		Point p;
		double mat[MDim3][MDim3];
		double mat1[MDim3];
	protected:
	public:
		combination comb;
		MinEnergy MinE;

		Face f;    // convex face
		int TGPfn;  //
		Face *TGPf;  //
		double Chp[MDim3]; // chemical potential

		// Find Convex Face
		void FindConvexFace(int Dim, int GPn, Point *GP, Point  PA);
		// Global Face in zero iteration
		void FindCH::FindGFace(int Dim, int GPn, Point *GP, Face &f);
		// iteration
		void IterConvexFace(int GPn, Point *GP, Face &f, Point  PA);
		// 
		void GetConvexFace(Face f);
		//
		void ShowConvexFace(string OutFileName);
		//calc chemical potential
		void CalcChemicalpotential(Face f, double (&Chp) [MDim3]);

		//Find Convex Hull
		void FindConvexHull(int Dim, int GPn, Point *GP);
		//
		void RecurConvexHull(int GPn, Point *GP, Face f);
		//
		void GetConvexHull();
		//
		void ShowConvexHull(string OutFileName);
		void deletep();
		//
		void ShowFace(string str, int n, int m, Face f);
		//
		void ShowPoint(string str, int n, Point p);

		//Find the lowest energy point, based on current chmical potetnial surface
		int FindCH::FindMaxD(int GPn, Point *GP, Face f, double Chp[MDim3]);
		//jugde point inside face f
		bool inside(Point p, Face f);
		//inside the edge
		int inside_edge(int n, int &subedge, int parray[MDim], Point p, Face f);
		//distance between point and face
		double distance(int mode, int n, Point p, Face f);
		//calc deteminant
		double calcdet(int n, double array[MDim3][MDim3]);
	};
} // end of VCLab
#endif
