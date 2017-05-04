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

#ifndef DATABASE_H_
#define DATABASE_H_

#include "VCLInput.h"
#include "Tool.h"

namespace VCLab
{
	using namespace std;

	struct Element
	{
		string name;			// "AL"
		string ref_state;		// "FCC"
		double mass;			//
		double h298k;			//
		double s298k;			//
	};

	struct Species
	{
		string name;			//
		string formula;			//
		double charge;			//
	};

	struct Function
	{
		string name;				// "AL"
		vector<double> T;			//
		vector<string> express;		//

	};

	struct express_digit
	{
		vector<int> powerT;			// 
		vector<double> coffT;		//
		vector<double> coffTLNT;	//
	};

	struct Parameter
	{
		int order;
		/* para types
		1 end member;2 binary interaction parameters; 
		3 ternary interaction parameters
		4 reciprocal interaction parameter;
		*/
		int kind;   
		int nsub2 = 0;  // binary para
		int nsub3 = 0;  // ternary para
		int idsub2[2];  // interaction ele id in binary para
		int idsub3[3];  // interaction ele id in ternary para
		int vidsub2[2];
		int vidsub3[3];

		string phasename;
		string type;
		
		vector<string> con[10]; // 10 sublattices, element in each sublattice is a vector
		// ele id in phase constitution£¬eg, constituent :Al,Mg,Zn:Zn,Va:
		// para G(BCC,Al:Zn,0),become G(BCC,1:4,0)
		int yn = 0;		// constituent number
		int yidc[MDim*3];	// constituent id
		int yids[MDim*3];	// sublattice id

		int vyn = 0;	 // varible y number
		int vyidc[MDim*3]; // constituent id
		int vyids[MDim*3]; // sublattice id
		int vyidv[MDim*3]; // para's vy in Phase's vy id
		vector<double> T;		//
		vector<string> express;			//
		express_digit express_digit[10]; // 10 T segment, terms in each segment is a vector 

	};

	struct Type_definition
	{
		string label;			//
		string model;			//
		string command;			//
		string phasename;		//
		string property;		//
		string disname;
		double value1;			//
		double value2;			//
	};

	struct Point
	{
		double x[MDim + 1];	// compositions
		double y[MDim * 3];	// yfractions
		int phid;			// belong to id phase
		double phfrac;		// phase farctions
		int fid;			// id in the convex face
		int fb;				// ture or not
	};

	struct Phase
	{
		int id;					//
		string name;			//
		string label;			//
		string type_definition;	//

		//sublattices
		int subln = 0;			//number
		double sublper[10];		//
		double sublpersum;

		//Constituent
		vector<string>	con[10];		// ele£¬   eg, (Al, Zn)(Zn,Va)
		vector<int>		conid[10];		// con id, eg, (0,  1)  (2, 3)
		vector<int>		conidele[10];	// ele id, eg, (0,  2)  (2, 1)
		int				conn[10];		// sub num, eg,  (2)     (2)
		// Parameters
		vector<Parameter> Parameters;


		//Grid Point
		int GPn;			// n
		int GPsn[10];		// GP n in sub
		double GPsdy[10];	// dy in sub
		int GPsndy[10];		// int(1/dy)
		Point *GP;
		double *GPGF;		// Gibbs energy, per mole atom;
		double *GPSF;		// entropy
		double *GPTC;		// TC for magnetism
		double *GPBMAGN;	// BMAGN
		double *GPtao;		// tao = T/TC

		// Gibbs energy, first, second derivative
		double GF, SF;
		double dGFx[MDim*3], dSFx[MDim * 3];
		double dGFt, dSFt, d2GFt, d2SFt;
		double d2GFxt[MDim * 3], d2SFxt[MDim * 3];
		double d2GFx[MDim * 3][MDim * 3], d2SFx[MDim * 3][MDim * 3];

		//x
		double x[MDim];
		//y
		int yn = 0;			// n
		double y[MDim * 3];		// yfrac(0.6, 0.4) (0.5, 0.5)
		int yidc[MDim * 3];		// id in phase conds
		int yide[MDim * 3];		// id in sys eles, eg, (0,  2)  (2, 1)
		int yids[MDim * 3];		// id in sublattices
		double ysp[MDim * 3];	// sublattices per
		int yidv[MDim * 3];		// id in vy
		//
		double T = 300;
		double P = 1e5;
		int mgap = 1;
		double magp = 0;	// magp, fcc hcp -3.0, bcc -1.0
		double TC;
		double BMAGN;
		//
		// Calc Equilibrium
		int b = 1; //c-e
		double DF = 0;  //Driving Force

		int vn = 0;      // overall varibles numbers
		// phase amount
		int vmb = 1;
		double vm = 1;	// default, is a varibles
		int vmidce;		// id in c-e
		double vmbound[2] = { 0,1 }; // boundary
		//lamda
		double vla[10] = { 0 };
		int vlan = 0;     
		int vlaids[10];    // id in sublattices
		int vlaidce[10];
		//yfractions
		double vy[MDim];
		int vyn = 0;		// 
		int vyids[MDim * 3];	// id in sub
		double vysp[MDim * 3];	//
		int vyidla[MDim * 3];	// id in lamda
		int vyidc[MDim * 3];	// id in cons
		int vyidce[MDim * 3];	// id in c-e
		int vyide[MDim * 3];	// id in ele
	};

	struct Reference
	{
		string number;				//
		string source;				//
	};

	class Database
	{
	private:
		string classname;
		int linenumber;
		vector<Element> DBElements;
		vector<Species> DBSpecies;
		vector<Function> DBFunctions;
		vector<Phase> DBPhases;
		vector<Reference> DBReferences;
		vector<Type_definition> DBType_definitions;
	public:
		vector<Element> SysElements;
		vector<Species> SysSpecies;
		vector<Function> SysFunctions;
		vector<Phase> SysPhases;
		vector<Reference> SysReferences;

		//
		~Database();
		// read database
		void readDatabase(string FilTDB);
		int readElement(string mlines);
		int readSpecies(string mlines);
		int readFunction(string mlines);
		int readPhase(string mlines);
		int readConstituent(string mlines);
		int readType_definition(string mlines);
		int readDefine_system_default(string mlines);
		int readDefault_command(string mlines);
		int readReference(string mlines);

		// define system
		void defineSystem(vector<string> EleList, vector<string> SelPhaList,
			vector<string> RejPhaList, string FilTDB);
		void getElement(vector<string> EleList);
		void getSpecies(vector<string> EleList);
		void getFunction(vector<string> EleList);
		void getPhase(vector<string> EleList);
		void regetPhase(vector<string> SelPhaList,
			vector<string> RejPhaList);

		void readParameter(string mlines, vector<string> EleList);
		void getParameter(vector<string> EleList, string FilTDB);
		void getReference(vector<string> elementlis);
		express_digit getexpress(string express, double T_start);

		// search
		int findPhase(string phasename, int n);
		int findElement(string phasename, vector<Element> Elements);

		// print on screen
		void showDatabase();
		void showSystem();
		void showElements(vector<Element> Elements);
		void showPhases(vector<Phase> Phases);
		void showConstituents(vector<Phase> Phase);

		// sort function for Phase struct
		bool SortPhase(const Phase& phase1, const Phase& phase2);

		// overloading function "<<" for Phase sturct
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

#ifndef DATABASE_H_
#define DATABASE_H_

#include "VCLInput.h"
#include "Tool.h"

namespace VCLab
{
	using namespace std;

	struct Element
	{
		string name;			// "AL"
		string ref_state;		// "FCC"
		double mass;			//
		double h298k;			//
		double s298k;			//
	};

	struct Species
	{
		string name;			//
		string formula;			//
		double charge;			//
	};

	struct Function
	{
		string name;				// "AL"
		vector<double> T;			//
		vector<string> express;		//

	};

	struct express_digit
	{
		vector<int> powerT;			// 
		vector<double> coffT;		//
		vector<double> coffTLNT;	//
	};

	struct Parameter
	{
		int order;
		/* para types
		1 end member;2 binary interaction parameters; 
		3 ternary interaction parameters
		4 reciprocal interaction parameter;
		*/
		int kind;   
		int nsub2 = 0;  // binary para
		int nsub3 = 0;  // ternary para
		int idsub2[2];  // interaction ele id in binary para
		int idsub3[3];  // interaction ele id in ternary para
		int vidsub2[2];
		int vidsub3[3];

		string phasename;
		string type;
		
		vector<string> con[10]; // 10 sublattices, element in each sublattice is a vector
		// ele id in phase constitution£¬eg, constituent :Al,Mg,Zn:Zn,Va:
		// para G(BCC,Al:Zn,0),become G(BCC,1:4,0)
		int yn = 0;		// constituent number
		int yidc[MDim*3];	// constituent id
		int yids[MDim*3];	// sublattice id

		int vyn = 0;	 // varible y number
		int vyidc[MDim*3]; // constituent id
		int vyids[MDim*3]; // sublattice id
		int vyidv[MDim*3]; // para's vy in Phase's vy id
		vector<double> T;		//
		vector<string> express;			//
		express_digit express_digit[10]; // 10 T segment, terms in each segment is a vector 

	};

	struct Type_definition
	{
		string label;			//
		string model;			//
		string command;			//
		string phasename;		//
		string property;		//
		string disname;
		double value1;			//
		double value2;			//
	};

	struct Phase
	{
		int id;					//
		string name;			//
		string label;			//
		string type_definition;	//

		//sublattices
		int subln = 0;			//number
		double sublper[10];		//
		double sublpersum;

		//Constituent
		vector<string>	con[10];		// ele£¬   eg, (Al, Zn)(Zn,Va)
		vector<int>		conid[10];		// con id, eg, (0,  1)  (2, 3)
		vector<int>		conidele[10];	// ele id, eg, (0,  2)  (2, 1)
		int				conn[10];		// sub num, eg,  (2)     (2)
		// Parameters
		vector<Parameter> Parameters;


		//Grid Point
		int GPn;			// n
		int GPsn[10];		// GP n in sub
		double GPsdy[10];	// dy in sub
		int GPsndy[10];		// int(1/dy)
		double **GPy;		// yfractions
		double *GPGF;		// Gibbs energy, per mole atom;
		double *GPSF;		// entropy
		double *GPTC;		// TC for magnetism
		double *GPBMAGN;	// BMAGN
		double *GPtao;		// tao = T/TC

		// Gibbs energy, first, second derivative
		double GF, SF;
		double dGFx[MDim*3], dSFx[MDim * 3];
		double dGFt, dSFt, d2GFt, d2SFt;
		double d2GFxt[MDim * 3], d2SFxt[MDim * 3];
		double d2GFx[MDim * 3][MDim * 3], d2SFx[MDim * 3][MDim * 3];

		//x
		double x[MDim];
		//y
		int yn = 0;			// n
		double y[MDim * 3];		// yfrac(0.6, 0.4) (0.5, 0.5)
		int yidc[MDim * 3];		// id in phase conds
		int yide[MDim * 3];		// id in sys eles, eg, (0,  2)  (2, 1)
		int yids[MDim * 3];		// id in sublattices
		double ysp[MDim * 3];	// sublattices per
		int yidv[MDim * 3];		// id in vy
		//
		double T = 300;
		double P = 1e5;
		int mgap = 1;
		double magp = 0;	// magp, fcc hcp -3.0, bcc -1.0
		double TC;
		double BMAGN;
		//
		// Calc Equilibrium
		int b = 1; //c-e
		double DF = 0;  //Driving Force

		int vn = 0;      // overall varibles numbers
		// phase amount
		int vmb = 1;
		double vm = 1;	// default, is a varibles
		int vmidce;		// id in c-e
		double vmbound[2] = { 0,1 }; // boundary
		//lamda
		double vla[10] = { 0 };
		int vlan = 0;     
		int vlaids[10];    // id in sublattices
		int vlaidce[10];
		//yfractions
		double vy[MDim];
		int vyn = 0;		// 
		int vyids[MDim * 3];	// id in sub
		double vysp[MDim * 3];	//
		int vyidla[MDim * 3];	// id in lamda
		int vyidc[MDim * 3];	// id in cons
		int vyidce[MDim * 3];	// id in c-e
		int vyide[MDim * 3];	// id in ele
	};

	struct Reference
	{
		string number;				//
		string source;				//
	};

	class Database
	{
	private:
		string classname;
		int linenumber;
		vector<Element> DBElements;
		vector<Species> DBSpecies;
		vector<Function> DBFunctions;
		vector<Phase> DBPhases;
		vector<Reference> DBReferences;
		vector<Type_definition> DBType_definitions;
	public:
		vector<Element> SysElements;
		vector<Species> SysSpecies;
		vector<Function> SysFunctions;
		vector<Phase> SysPhases;
		vector<Reference> SysReferences;

		//
		~Database();
		// read database
		void readDatabase(string FilTDB);
		int readElement(string mlines);
		int readSpecies(string mlines);
		int readFunction(string mlines);
		int readPhase(string mlines);
		int readConstituent(string mlines);
		int readType_definition(string mlines);
		int readDefine_system_default(string mlines);
		int readDefault_command(string mlines);
		int readReference(string mlines);

		// define system
		void defineSystem(vector<string> EleList, string FilTDB);
		void getElement(vector<string> EleList);
		void getSpecies(vector<string> EleList);
		void getFunction(vector<string> EleList);
		void getPhase(vector<string> EleList);
		void readParameter(string mlines, vector<string> EleList);
		void getParameter(vector<string> EleList, string FilTDB);
		void getReference(vector<string> elementlis);
		express_digit getexpress(string express, double T_start);

		// search
		int findPhase(string phasename, int n);
		int findElement(string phasename, vector<Element> Elements);

		// print on screen
		void showDatabase();
		void showSystem();
		void showElements(vector<Element> Elements);
		void showPhases(vector<Phase> Phases);
		void showConstituents(vector<Phase> Phase);

		// sort function for Phase struct
		bool SortPhase(const Phase& phase1, const Phase& phase2);

		// overloading function "<<" for Phase sturct
	};

} // end of VCLab

>>>>>>> b34a88e731ed6049f6c364b37b469ae02b913009
#endif