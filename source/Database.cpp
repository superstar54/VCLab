#include "..\include\Database.h"

namespace VCLab
{
	using namespace std;

	Database::~Database()
	{
		if(debug>=2)
			cout << "Database modulue exit." << endl;
	}

	void Database::readDatabase(string FilTDB)
	{
		int		flag		= 0;
		int		checkline	= 0;
		string	line;		// read a line,
		string	mlines;		// merge lines to mlines
		string	keyword;

		linenumber = 0;
		classname = "Database";
		PrintLine(">");
		PrintMiddle("Read Database");
		PrintVstrstr("Database file", FilTDB);
		fstream Inp(FilTDB.c_str(), ios::in);
		if (!Inp)
		{
			PrintError("File " + FilTDB + " could not be opened");
		}
		//start to read TDB line by line
		while (!Inp.eof())
		{					
			mlines = readline(Inp, linenumber);  //read a line, skip empty line, and the comment line,change to capital
			keyword = ReadFirWor(mlines,{ ' ' },0);		//read key word
			if (keyword=="ELEMENT")
				checkline = readElement(mlines);
			else if (keyword == "SPECIES")
				checkline = readSpecies(mlines);
			else if(keyword=="FUNCTION")
				checkline = readFunction(mlines);
			else if (keyword=="PHASE")
				checkline = readPhase(mlines);
			else if (keyword=="CONSTITUENT")
				checkline = readConstituent(mlines);
			else if (keyword == "TYPE_DEFINITION")
				checkline = readType_definition(mlines);
			else if(keyword=="DEFINE_SYSTEM_DEFAULT")
				checkline = readDefine_system_default(mlines);
			else if (keyword == "DEFAULT_COMMAND ")
				checkline = readDefault_command(mlines);
			else if (keyword.find("LIST_OF_REFERENCES")!=string::npos)
				checkline = readReference(mlines);
			
			if (checkline != 0)  // if not sucess, print the error line
			{
				cerr << "Error in line: " << linenumber << endl;
				cerr << mlines<<endl;
				return;
			}
		}
		Inp.close();
	}

	int Database::readElement(string mlines)
	{
		int		checkline	= -1;
		string	endmark		= { " " };
		Element	temp;

		// ELEMENT AL   FCC_A1  2.6982E+01  4.5773E+03  2.8322E+01
		temp.name		= ReadFirWor(mlines, endmark, 0);			//name
		temp.ref_state	= ReadFirWor(mlines, endmark, 0);		//ref_state
		temp.mass		= StrToDou(ReadFirWor(mlines, endmark, 0));	//mass
		temp.h298k		= StrToDou(ReadFirWor(mlines, endmark, 0));	//H298K
		temp.s298k		= StrToDou(ReadFirWor(mlines, endmark, 0));	//S298K

		DBElements.push_back(temp);								//add a new element
		checkline = 0;
		return checkline;
	}

	int Database::readFunction(string mlines)
	{
		int			checkline	= -1;
		int			nseg		= 0;		// number of T segment
		string		endmark1	= { ' ' };
		string		endmark2	= { ';' };
		string		flag		= "Y";
		Function	temp;

		/* FUNCTION GHSERAL    2.98150E+02 
		-79.15+137.01*T-24.3*T*LN(T)-.001*T**2;  7.00E+02  Y
		- 118.378 + 188.6*T - 31.72*T*LN(T) - 1.24E+28*T**(-9);  2.90E+03  N 
		*/
		temp.name = ReadFirWor(mlines, endmark1, 0);						//name
		temp.T.push_back(StrToDou(ReadFirWor(mlines, endmark1, 0)));		//start T
		while (!flag.compare("Y"))
		{
			//read express, and mend it, delete all spacing, and add "+" if needed.
			temp.express.push_back(MendExp(ReadFirWor(mlines, endmark2, 0)));	
			temp.T.push_back(StrToDou(ReadFirWor(mlines, endmark1, 0))); //end T
			flag = ReadFirWor(mlines, endmark1, 0);						//end "Y" or "N"
		}
		DBFunctions.push_back(temp);									// add a new function
		checkline = 0;
		return checkline;
	}

	int Database::readPhase(string mlines)
	{
		int		checkline	= -1;
		string	endmark1	= { ' ' };
		string	endmark2	= { ':' };
		Phase	temp;
		string	firstword;
		string	templine;
		
		//       name   type_defition  subln  sublper 
		//PHASE LIQUID:L     %&			 2			 1   3 
		templine	= ReadFirWor(mlines, endmark1, 0);
		firstword	= ReadFirWor(templine, endmark2, 0);
		temp.name	= firstword;
		if (firstword == templine)
			temp.label = { ' ' };//BCC_A2
		else
			temp.label = templine; //LIQUID:L
		temp.type_definition	= ReadFirWor(mlines, endmark1, 0);  //type_definition
		temp.subln = StrToInt(ReadFirWor(mlines, endmark1, 0));  //number of sublattices
		for (int i = 0; i < temp.subln;i++)
			temp.sublper[i] = StrToDou(ReadFirWor(mlines, endmark1, 0));  //percentage of sublattices
		DBPhases.push_back(temp);  // add a new phase
		checkline = 0;
		return checkline;
	}

	int Database::readConstituent(string mlines)
	{
		int		checkline = -1;
		int		phaseid;
		string	phasename;
		string	firstword;
		string	templine;
		string	endmark1 = { ' ' };
		string	endmark2 = { ':' };
		string	endmark3 = { ',' };
		string	endmark4 = { "%" };
		//            phasename   constituent
		//CONSTITUENT BCC_A2  :AL,Mg,VA : VA :  !
		templine	= ReadFirWor(mlines, endmark1, 0);
		phasename	= ReadFirWor(templine, endmark2, 0);	//phasename
		phaseid		= findPhase(phasename, 0);				//find phase in DBPhases, and return id 
		// :Al%,Cu:Al,Cu%:Cu:
		firstword	= ReadFirWor(mlines, endmark2, 0);		// delelte the first ":"
		for (int i = 0;i < DBPhases[phaseid].subln;i++)
		{
			templine = ReadFirWor(mlines, endmark2, 0);
			// Al%,Cu:
			while (templine.length()!=0)					//add a new element in sublattice i
			{
				if (templine.find(endmark4) != string::npos)
				{
					firstword = ReadFirWor(templine, endmark3, 0);							// Al%,
					DBPhases[phaseid].con[i].push_back(ReadFirWor(firstword, endmark4, 0)); //Al
				}
				else //Al,
				{
					DBPhases[phaseid].con[i].push_back(ReadFirWor(templine, endmark3, 0));
				}
			}
		}
		checkline = 0;
		return checkline;
	}

	int Database::readSpecies(string mlines)
	{
		Species temp;
		int		checkline	= -1;
		string	endmark		= { " " };
		// SPECIES  CU2                 CU2!
		// SPECIES  MG2SN               MG2SN1 !
		temp.name		= ReadFirWor(mlines, endmark, 0);	// name
		temp.formula	= ReadFirWor(mlines, endmark, 0);	// compositions
		DBSpecies.push_back(temp);							//add a new element
		checkline = 0;
		return checkline;
	}

	int Database::readType_definition(string mlines)
	{
		Type_definition temp;
		int		checkline	= -1;
		string	endmark		= { " " };
		// ELEMENT AL   FCC_A1  2.6982E+01  4.5773E+03  2.8322E+01
		temp.label		= ReadFirWor(mlines, endmark, 0);		//name
		temp.model		= ReadFirWor(mlines, endmark, 0);		//ref_state
		temp.command	= ReadFirWor(mlines, endmark, 0);		//mass
		if (temp.command == "*")
		{
			DBType_definitions.push_back(temp);					//add a new element
			checkline = 0;
			return checkline;
		}
		else if (temp.command == "A_P_D")
		{
			//temp.phasename = ReadFirWor(mlines, endmark, 0);	//H298K
			//temp.property = ReadFirWor(mlines, endmark, 0);	//S298K
			//temp.value1 = StrToDou(ReadFirWor(mlines, endmark, 0));	//S298K
			//temp.value2 = StrToDou(ReadFirWor(mlines, endmark, 0));	//S298K
		}
		else if (temp.command == "AMEND_PHASE_DESCRIPTION")
		{
			//temp.phasename = ReadFirWor(mlines, endmark, 0);	//H298K
			//temp.property = ReadFirWor(mlines, endmark, 0);	//S298K
			//temp.disname = ReadFirWor(mlines, endmark, 0);	//S298K
		}
		DBType_definitions.push_back(temp);					//add a new element
		checkline = 0;
		return checkline;
	}

	int Database::readReference(string mlines)
	{
		Reference temp;
		int checkline	= -1;
		string endmark	= { " " };
		string endmark2 = { "'" };
		string endmark3 = { "SOURCE" };
		/*
		LIST_OF_REFERENCES
		NUMBER  SOURCE
		REF1   'PURE4 - SGTE Pure Elements' ! '*/
		temp.number = ReadFirWor(mlines, endmark3, 0);			// "SOURCE"
		while (mlines.length() != 0)
		{
			temp.number = ReadFirWor(mlines, endmark, 0);		// "REF1"
			temp.source = ReadFirWor(mlines, endmark2, 0);
			temp.source = ReadFirWor(mlines, endmark2, 0);		// "PURE4 ..."
			DBReferences.push_back(temp);						//add a new element

		}
		checkline = 0;
		return checkline;
	}

	int Database::readDefine_system_default(string mlines)
	{
		int checkline = -1;
		checkline = 0;
		return checkline;
	}

	int Database::readDefault_command(string mlines)
	{
		int checkline = -1;
		checkline = 0;
		return checkline;
	}

	void Database::showDatabase()
	{
		cout << endl;
		PrintVstrint("Elements: ", DBElements.size());
		for (auto x : DBElements) cout << x.name << "  ";
		cout << "\n\n";
		cout << "Phases: "<<DBPhases.size()<<endl;
		for (auto x : DBPhases) cout << x.name << "  ";
		cout << endl;
	}

	void Database::showElements(vector<Element> Elements)
	{
		cout << endl;
		cout << "Element" << "  Mass g/mol" << "  Reference state"
			<< "  H298K" << "  S298K" << endl;
		for (auto x : DBElements) 
		{
			cout <<"  "<< x.name << "  " << x.mass << "  " << x.ref_state
				<< "  "<< x.h298k << "  " << x.s298k << endl;
		}
	}
	//
	void Database::showPhases(vector<Phase> Phases)
	{
		
	}
	//
	void Database::showConstituents(vector<Phase> Phases)
	{

	}
	// find phase in DBPhases base on phase name, and return the id;
	int Database::findPhase(string phasename, int n)
	{
		int i		= 0;
		int phaseid = -1;
		if (n == 0)
		{
			for (auto x : DBPhases)   // find in DBPhases
			{
				if (x.name == phasename) phaseid = i;
				i++;
			}
			// if phase not be found, exit, and print error
			if (phaseid == -1)
			{
				cerr << "line " << linenumber << ": Phase: " << phasename << "doesn't exit\n";
				return -1;
			}
		}
		else   // find in SysPhases
		{
			for (auto x : SysPhases)
			{
				if (x.name == phasename) phaseid = i;
				i++;
			}
		}
		return phaseid;
	}
	
} //end of VCLab

/* On going work
0. Species
1. exp term in express
2. magnetism
3. order and disorder
4. print results
*/

