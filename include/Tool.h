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
#ifndef TOOL_H
#define TOOL_H

#include "VCLInput.h"

namespace VCLab
{
	using namespace std;
	
	//const int MDim = 10;
	const int MDim3 = 60;


	///read a line, skip empty line, and the comment line,change to capital
	string readline(fstream& Inp, int &linecount);
	// read the first word before endmark
	std::string ReadFirWor(std::string & tline, std::string endmark, int flag);
    // find the endmark forward
	int findendmarkforward(std::string tline, std::string endmark, std::string endmark_array[], int n);
	// find the endmark backward
	int findendmarkbackward(std::string tline, std::string endmark, std::string endmark_array[], int n);
	// mend the express, eliminate the space
	std::string MendExp(std::string firstword);  // delete spacing

	//change to capital, delete space before and behind
	std::string capital(std::string tline);   // change low to upper

	double StrToDou(string strval);
	int StrToInt(string strval);

	//
	int gauss_elimination(int n, double(&a)[MDim3][MDim3], double(&b)[MDim3]);


	// cartesian production
	// Input matrix n by m; output matrix om by on, with on = m;
	void cart_product(int **output, int om, int **input, int n[], int m);

	class combination
	{
	private:
	protected:
	public:
		int toln = 1;
		int ncount = 0;
		double **comb;
		int combn;
		int combm;
		int mode = 0; 
		double dy;

		void combination::combinations_recursive(int *pos, int depth, int margin);
		void combinations(int mod, double d, int n, int m);
		void show();
		void deletecomb();
	};

	//
	// Print somethings
	void PrintTitle();
	void PrintLine(string str);
	void PrintMiddle(string str);
	void PrintWarning(string Message);
	void PrintError(string Message);
	void PrintVstrstr(string vstr, string str);
	void PrintVstrint(string vstr, int i);
}
#endif

