<<<<<<< HEAD
/*
*   This file is a part of the VCLab (Virtual Calphad Laboratory) software project.
*   For more details please contact xingwang@csu.edu.cn
*
*   Authors:    Xing Wang
*
*   Copyright (c) 2015- Phase Diagram Center. Central South University. China
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

#ifndef VCLInput_H
#define VCLInput_H

#include <cstring>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <strstream>
#include <cmath>
#include <algorithm>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include "Tool.h"


#define eps100 1e-100
#define eps 1e-15
#define eps9 1e-9
#define eps6 1e-6
#define eps3 1e-3

const int MDim = 5;
//const int MDim3 = 60;

const int debug = 0;
using namespace std;


namespace VCLab
{
	struct Condition
	{
		double	N;					// Mole of total system
		double	T;					// Temperature
		int		Tb = 0;				// Default: not a varibles
		double	P;					// Pressure
		int		nele;				// Number of elements
		string	ele[MDim + 2];		// Elements list
		double	x[MDim];			// Compositions list (mole)
		double	m[MDim];			// Compositions list (mass)

		double	dx = 0.01;			// Global grid interval

	};

	class ProjectInput
	{
	private:
	protected:
	public:

	Condition Cond;
	/*
	Fil: File
	Ele: Element
	Com: Composition
	*/

	string			FilPro			= "VCLInput.txt";	// Input file of project
	
	string			Dimension		= "0D";				// Point, Line, Map
	string			Mode			= "Equilibrium";	// Point, Line, Map
	string			FilTDB			= "";				// TDB file
	string			FilOut			= "";				// Output file
	string			EleInp			= "";				// Elements input
	vector<string>	EleList			= { "/-", "VA" };	// Elements list
	string			ComInp			= "";				// Compositions input
	vector<double>	ComList			= { 0.0, 0.0 };		// Compositions list

	string			SelPhaInp = "";						// Fixed phase input
	vector<string>	SelPhaList = {};					// Fixed phase list
	string			RejPhaInp = "";						// Fixed phase input
	vector<string>	RejPhaList = {};					// Fixed phase list

	string			FixPhaInp		= "";				// Fixed phase input
	vector<string>	FixPhaList		= {};				// Fixed phase list
	vector<double>	FixPhaFraList	= {};				// Fixed phases fractions list


	int				nV = 0;			// Number of Varibles
	string			Var;			// Varibles
	double			Var_s;		// start value
	double			Var_e;		// end value
	double			Var_d;		// interval

	void ReadVCLInput();
	void SortEle(vector<string> &EleList, vector<double> &composilist);
	double StrToDou(string strval);
};
}
#endif
=======
/*
*   This file is a part of the VCLab (Virtual Calphad Laboratory) software project.
*   For more details please contact xingwang@csu.edu.cn
*
*   Authors:    Xing Wang
*
*   Copyright (c) 2015- Phase Diagram Center. Central South University. China
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

#ifndef VCLInput_H
#define VCLInput_H

#include <cstring>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <strstream>
#include <cmath>
#include <algorithm>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include "Tool.h"

#define eps100 1e-100
#define eps 1e-15
#define eps9 1e-9
#define eps6 1e-6
#define eps3 1e-3

const int MDim = 6;
//const int MDim3 = 60;

const int debug = 0;
using namespace std;


namespace VCLab
{


	struct Condition
	{
		double	N = 1.0;		// Mole of total system
		double	T = 298.15;	// Temperature
		int		Tb = 0;		// Default: not a varibles
		double	P = 101325.0; // Pressure
		int		nele = 2;		// Number of elements
		string	ele[MDim + 2];			// Elements list
		double	x[MDim];			// Compositions list (mole)
		double	m[MDim];			// Compositions list (mass)

		double	dx = 0.01;			// Global grid interval

	};

	class ProjectInput
	{
	private:
	protected:
	public:

	Condition Cond;
	/*
	Fil: File
	Ele: Element
	Com: Composition
	*/

	string			FilPro			= "VCLInput.txt";	// Input file of project

	string			VCLmode			= "Point";			// Point, Line, Map
	string			FilTDB			= ".TDB";			// TDB file
	string			FilOut			= "Results.txt";	// Output file
	string			EleInp			= "";				// Elements input
	vector<string>	EleList			= { "/-", "VA" };	// Elements list
	string			ComInp			= "";				// Compositions input
	vector<double>	ComList			= { 0.0, 0.0 };		// Compositions list
	string			FixPhaInp		= "";				// Fixed phase input
	vector<string>	FixPhaList		= {};				// Fixed phase list
	vector<double>	FixPhaFraList	= {};				// Fixed phases fractions list


	int				nV = 0;			// Number of Varibles
	string			Var;			// Varibles
	double			Var_s;		// start value
	double			Var_e;		// end value
	double			Var_d;		// interval

	void ReadVCLInput();
	void SortEle(vector<string> &EleList, vector<double> &composilist);
	double StrToDou(string strval);
};
}
#endif
>>>>>>> b34a88e731ed6049f6c364b37b469ae02b913009
