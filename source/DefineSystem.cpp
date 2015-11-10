
#include "..\include\Database.h"
#include "..\include\Tool.h"

namespace VCLab
{
	using namespace std;
	//
	void Database::defineSystem(vector<string> EleList, string FilTDB)
	{
		cout << endl;
		PrintLine(">");
		PrintMiddle("Define system");

		//get data
		getElement(EleList);
		getPhase(EleList);
		getParameter(EleList, FilTDB);
		
		//sortPhase(EleList, SysPhases, FilTDB);
		//sort(SysPhases.begin(), SysPhases.end(), SortPhase);  // SortPhase: compare function for Phase struct
		//for (vector<Phase>::iterator iter = SysPhases.begin(); iter != SysPhases.end(); ++iter)
		//{
		//	cout << iter << endl;  // over load the "<<" for Phase struct
		//}
	}

	void Database::getElement(vector<string> EleList)
	{
		int elementid = 0;
		for (auto ele : EleList)
		{
			elementid = findElement(ele, DBElements);
			if (elementid == -1)
			{
				cerr << "Element: " << ele << " doesn't exit\n";
				//exit(1);
			}
			else
			SysElements.push_back(DBElements[elementid]);
		}
	}

	void Database::getSpecies(vector<string> EleList)
	{

	}

	void Database::getPhase(vector<string> EleList)
	{
		int i, j;
		int elementid = 0;
		int checkphase = 0;
		int checksublattice = 0;
		Phase temp;
		for (auto phase : DBPhases)
		{
			temp = phase;
			temp.sublpersum = 0;
			checkphase = 1;  // phase exist
			int idcon = 0;
			for (i = 0;i < phase.subln;i++)
			{
				j = 0;
				temp.con[i].clear();
				checksublattice = 0;
				for (auto ele : phase.con[i])
				{
					elementid = findElement(ele, SysElements);
					if (elementid != -1)
					{
						checksublattice = 1;
						temp.con[i].push_back(ele);
						temp.conid[i].push_back(idcon);
						temp.conidele[i].push_back(elementid);
						//
						temp.yids[temp.yn] = i;				// sublattice
						temp.ysp[temp.yn] = phase.sublper[i]; // sublattice
						temp.yidc[temp.yn] = idcon;			// constitution
						temp.yide[temp.yn] = elementid;		// element
						temp.yn++;
						//
						temp.vyids[temp.vyn] = i;			//vy in sublattice
						temp.vysp[temp.vyn] = phase.sublper[i]; //sublattice
						temp.vyidla[temp.vyn] = temp.vlan;	// varibles lamda for sub
						temp.vyidc[temp.vyn] = idcon;		// constitution
						temp.vyide[temp.vyn] = elementid;	//element
						temp.vyn++;
						idcon++;
					}
					j++;
				}
				temp.conn[i] = temp.con[i].size();
				if (checksublattice == 0)
				{
					checkphase = 0; // this sublattice has no system element
					break;			// this phase is not exist.
				}
				if (temp.conn[i] > 1)
				{
					temp.vlaids[temp.vlan] = i;   // add a lamda
					temp.vlan++;
					temp.sublpersum += phase.sublper[i];
				}
				else if (temp.conn[i] == 1)
				{
					temp.y[temp.yn - 1] = 1;	//only one y in this sublattices，therefore y==1;
					temp.yidv[temp.yn - 1] = -1; //this y is not a varibles
					temp.vyn--; //not a varible , remove this one.
					if (elementid != 1)  //skip va
						temp.sublpersum += phase.sublper[i];
				}
				//sort(temp.con[i].begin(), temp.con[i].end());
			}
			if (checkphase == 0)
				continue;  //this phase is not exist
			else
			{
				for (i = 0;i < temp.vyn;i++)
				{
					temp.yidv[temp.vyidc[i]] = i;
				}
				SysPhases.push_back(temp);
			}
		}
	}

	void Database::getFunction(vector<string> EleList)
	{

	}

	void Database::getParameter(vector<string> EleList,  string FilTDB)
	{
		
		string line;    
		string mlines;  
		string keyword; 

		linenumber = 0;
		//start read file line by line
		fstream Inp(FilTDB.c_str(), ios::in);
		while (!Inp.eof())
		{
			if (debug >= 2)
				cout << "line: " << linenumber << endl;
			mlines = readline(Inp, linenumber);  //read a line, skip empty line, and the comment line,change to capital
			keyword = ReadFirWor(mlines, { ' ' }, 0);					//定位文件末尾

			if (keyword == "PARAMETER")
			{
				readParameter(mlines, EleList);
			}
		}
		Inp.close();
	}

	void Database::readParameter(string mlines, vector<string> EleList)
	{
		Parameter temp;
		int i;
		int phaseid;
		int elementid;
		int tempid;
		int nele;   // ele number in one sublattice
		int nseg = 0;
		string templine;
		string constituenmlines;
		string flag("Y");
		string Referencename;
		//                 0   1    2    3    4   5 
		string endmark[] = { " ", ":", ",","(",")",";" };
		templine		= ReadFirWor(mlines, endmark[4], 0);
		temp.type		= ReadFirWor(templine, endmark[3], 0);
		temp.phasename	= ReadFirWor(templine, endmark[2], 0);
		phaseid			= findPhase(temp.phasename, 1);
		if (phaseid != -1)
		{
			constituenmlines = ReadFirWor(templine, endmark[5], 0);
			temp.order = StrToInt(templine);
			// ele in para belong to phase cons
			int checkparameter = 1;
			for (i = 0;i < SysPhases[phaseid].subln;i++) // loop sub
			{
				nele = 0;
				templine = ReadFirWor(constituenmlines, endmark[1], 0);
				while (templine.length() != 0)  // loop ele in sub
				{
					string ele = ReadFirWor(templine, endmark[2], 0);
					tempid = 0;
					elementid = -1;
					for (auto constele : SysPhases[phaseid].con[i])   // find ele in sys ele
					{
						if (constele == ele) elementid = tempid;
						tempid++;
					}
					if (elementid==-1)
					{
						checkparameter = 0;//there ele in this para doesn't belong to phase
						break;  //not necessary to compare other ele
					}
					temp.con[i].push_back(ele);   // add new ele
					temp.yidc[temp.yn] = SysPhases[phaseid].conid[i][elementid];  // id in cons

					temp.vyidc[temp.vyn] = temp.yidc[temp.yn]; // vy in cons
					if (SysPhases[phaseid].yidv[temp.yidc[temp.yn]] >= 0)
					{
						temp.vyidv[temp.vyn] = SysPhases[phaseid].yidv[temp.yidc[temp.yn]];
						temp.vyn++;
					}
					temp.yn++;
				}
				switch (temp.con[i].size())
				{
				case 1:
					
					break;
				case 2:
					temp.nsub2++;  // binary interaction parameter
					temp.idsub2[0]	= temp.yidc[temp.yn - 2];
					temp.vidsub2[0] = SysPhases[phaseid].yidv[temp.yidc[temp.yn - 2]];
					temp.idsub2[1]	= temp.yidc[temp.yn - 1];
					temp.vidsub2[1] = SysPhases[phaseid].yidv[temp.yidc[temp.yn - 1]];
					break;
				case 3:
					temp.kind = 3;  // ternary interaction parameter
					temp.idsub3[0] = temp.yidc[temp.yn - 3];
					temp.idsub3[1] = temp.yidc[temp.yn - 2];
					temp.idsub3[2] = temp.yidc[temp.yn - 1];
					temp.vidsub3[0] = SysPhases[phaseid].yidv[temp.yidc[temp.yn - 3]];
					temp.vidsub3[0] = SysPhases[phaseid].yidv[temp.yidc[temp.yn - 2]];
					temp.vidsub3[0] = SysPhases[phaseid].yidv[temp.yidc[temp.yn - 1]];
					break;
				}
				if (checkparameter == 0)
					break;  
			}

			if (checkparameter == 1) //read express
			{
				// parameter kind
				if (temp.nsub2 == 1)
					temp.kind = 2;	// binary interaction parameter
				else if (temp.nsub2 == 2)
					temp.kind = 4;	// sublattice interaction parameter
				else if (temp.nsub3 == 1)
					temp.kind = 3;	// ternary interaction parameter
				else
					temp.kind = 1;	// endmember parameter

				// the same as readFunction, which can be relize using a function.
				temp.T.push_back(StrToDou(ReadFirWor(mlines, endmark[0], 0)));
				int nseg = 0;
				while (!flag.compare("Y"))
				{
					temp.express.push_back(MendExp(ReadFirWor(mlines, endmark[5], 0)));
					temp.T.push_back(StrToDou(ReadFirWor(mlines, endmark[0], 0)));
					// read express into digit form
					temp.express_digit[nseg] = getexpress(temp.express[nseg], temp.T[nseg]);
					// careful: add " ", when one line merge another line, done before
					flag = ReadFirWor(mlines, endmark[0], 0);
					nseg++;
				}
				// reference
				if (mlines.length() != 0)
				{
					Referencename = mlines;
					//find reference id

					//add ref
				}
				SysPhases[phaseid].Parameters.push_back(temp);
			}
		}
	}

	void Database::getReference(vector<string> EleList)
	{

	}

	express_digit Database::getexpress(string express, double T_start)
	{
		//+3028.879+125.251171*T  -24.3671976*T*LN(T) - .001884662*T**2 - 8.77664E-07*T**3 + 74092 * T**(-1)+ 7.9337E-20*T**7;
		//                 0   1    2    3    4   5 
		string endmark1[]	= { "+", "-"};
		string endmark2		= { "#" };
		string endmark3[]	= { "*T**(","*T**"};
		string endmark4		= { "*T*LN(T)"};
		string endmark5		= { "*T" };
		string endmark6[]	= { "(",")" };
		string endmark7		= {"*"};
		string endmark8		= { "#*T" };
		
		int coff_sign = 1;
		int functionid = 0;
		double coff = 1;
		int power = 0;
		express_digit express_digit1;
		express_digit express_digit2;
		string firstword;
		string tempexpress;
		string insertplus = {"+"};
		string functionname;

		//sign: "+","-"
		coff_sign = 1;
		if(express[0]=='-')
			coff_sign = -1;
		// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		while (express.size()!=0)
		{
			if (express.find(endmark8) != string::npos)   // # function
			{
				//find "+","-" before function
				int i = findendmarkforward(express, endmark8, endmark1, 2);
				if (i == -1)
				{
					cerr << "line " << linenumber << ": Cannot find " << endmark2 << " in" << express << endl;
					exit(1);
				}
				firstword = express.substr(i, express.find(endmark2) - i); //already delete "#"
				express.erase(i, firstword.length() + 1); //  delete the part already found
														  //firstword.erase(firstword.length());  
				power = 0;
				if (firstword.find(endmark5) != string::npos)  // GHSERAL*T
				{
					power = 1;
					firstword.erase(firstword.length()-3,3);  //mode 0: delete "*"
				}
				if (firstword.find(endmark7) != string::npos)  // 2*GHSERAL
				{
					coff = StrToDou(ReadFirWor(firstword, endmark7, 0)); //mode 0: delete "*"
					functionname = firstword;
				}
				else
				{
					if (firstword[0] == '+')
						coff = 1;
					else
						coff = -1;
					functionname = firstword.substr(1, firstword.length());
				}
				//find function id in DBfunctions
				functionid = -1;
				i = 0;
				for (auto func : DBFunctions)
				{
					if (func.name == functionname)
					{
						functionid = i;
						break;
					}
					i++;
				}
				if (functionid == -1)
				{
					cerr << "line " << linenumber << ": Function: " << functionname << "not found!" << endl;
					exit(1);
				}
				// temperature segmentation
				int T_id = -1;
				int Tsize;
				Tsize = DBFunctions[functionid].T.size();
				for (int i = 0; i < Tsize; i++)
				{
					if (T_start >= (DBFunctions[functionid].T[i] - 0.2) &&  // T - 0.02, avoid case: 298.14 < 298.15
						T_start < (DBFunctions[functionid].T[i + 1] + 0.2))
					{
						T_id = i;
						break; // find the T interval
					}
				}
				if (T_id == -1)
				{
					cerr << "line " << linenumber << ": T interval of " << functionname << " not found!" << endl;
					exit(1);
				}
				else
				{
					string expressfunc = DBFunctions[functionid].express[T_id];
					express_digit2 = getexpress(expressfunc, T_start);
				}
				// add digit express, and multiply the coff
				for (auto coffT : express_digit2.coffT) express_digit1.coffT.push_back(coffT*coff);
				for (auto powerT : express_digit2.powerT) express_digit1.powerT.push_back(powerT + power);
				for (auto coffTLNT : express_digit2.coffTLNT) express_digit1.coffTLNT.push_back(coffTLNT*coff);
			}
			else if (express.find(endmark2) != string::npos)   // # function
			{
				//find "+","-" before function
				int i = findendmarkforward(express, endmark2, endmark1, 2);
				if (i == -1)
				{
					cerr << "line "<<linenumber<<": Cannot find " << endmark2 << " in" << express << endl;
					exit(1);
				}
				firstword = express.substr(i, express.find(endmark2)-i); //already delete "#"
				express.erase(i, firstword.length()+1); //  delete the part already found
				//firstword.erase(firstword.length());  
				if (firstword.find(endmark7) != string::npos)  // 2*GHSERAL
				{
					coff = StrToDou(ReadFirWor(firstword, endmark7, 0)); //mode 0: delete "*"
					functionname = firstword;
				}
				else
				{
					if (firstword[0] == '+')
						coff = 1;
					else
						coff = -1;
					functionname = firstword.substr(1, firstword.length());
				}
				//find function id in DBfunctions
				functionid = -1;
				i = 0;
				for (auto func : DBFunctions)
				{
					if (func.name == functionname)
					{
						functionid = i;
						break;
					}
					i++;
				}
				if (functionid == -1)
				{
					cerr << "line " << linenumber << ": Function: " << functionname << "not found!" << endl;
					exit(1);
				}
				// temperature segmentation
				int T_id = -1;
				int Tsize;
				Tsize = DBFunctions[functionid].T.size();
				for (int i = 0; i < Tsize; i++)
				{
					if (T_start >= (DBFunctions[functionid].T[i] - 0.2) &&  // T - 0.02, avoid case: 298.14 < 298.15
						T_start < (DBFunctions[functionid].T[i+1] + 0.2))
					{
						T_id = i;
						break; // find the T interval
					}
				}
				if (T_id == -1)
				{
					cerr << "line " << linenumber << ": T interval of " << functionname << " not found!" << endl;
					exit(1);
				}
				else
				{
					string expressfunc = DBFunctions[functionid].express[T_id];
					express_digit2 = getexpress(expressfunc, T_start);
				}
				// add digit express, and multiply the coff
				for (auto coffT : express_digit2.coffT) express_digit1.coffT.push_back(coffT*coff);
				for (auto powerT : express_digit2.powerT) express_digit1.powerT.push_back(powerT);
				for (auto coffTLNT : express_digit2.coffTLNT) express_digit1.coffTLNT.push_back(coffTLNT*coff);
			}
			else if (express.find(endmark3[0]) != string::npos) //*T**(
			{
				//2.1E+01*T**(-1)
				int k = express.find(endmark3[0]);
				int i = findendmarkforward(express, endmark3[0], endmark1, 2);
				if (i == -1)
				{
					cerr << "line " << linenumber << ": Cannot find " << endmark3[0] << " in" << express << endl;
					exit(1);
				}
				firstword = express.substr(i, k-i);
				express_digit1.coffT.push_back(StrToDou(firstword));
				string tempendmark[] = { ")" };
				int j = findendmarkbackward(express,endmark3[0], tempendmark, 1);
				if (j == -1)
				{
					cerr << "line " << linenumber << ": Cannot find "<<endmark6[1]<<" in" << express << endl;
					exit(1);
				}
				firstword = express.substr(k+5, j-k-5);
				express_digit1.powerT.push_back(StrToInt(firstword));
				express.erase(i, j-i+1);
			}
			else if (express.find(endmark3[1]) != string::npos) //*T**
			{
				//2.1E+01*T**2
				int k = express.find(endmark3[1]);
				int i = findendmarkforward(express, endmark3[1],endmark1, 2);
				if (i == -1)
				{
					cerr << "line " << linenumber << ": Cannot find " << endmark3[1] << " in" << express << endl;
					exit(1);
				}
				firstword = express.substr(i, k-i);
				express_digit1.coffT.push_back(StrToDou(firstword));
				int j = findendmarkbackward(express, endmark3[1], endmark1, 2);
				if (j == -1)
				{
					cerr << "line " << linenumber << ": Cannot find " << endmark6[1] << " in" << express << endl;
					exit(1);
				}
				firstword = express.substr(k + 4, j - k - 4);
				express_digit1.powerT.push_back(StrToInt(firstword));
				express.erase(i, j - i);
			}
			else if (express.find(endmark4) != string::npos) //*T*LN(T)
			{
				//2.1*T*LN(T)
				int k = express.find(endmark4);
				int i = findendmarkforward(express, endmark4, endmark1, 2);
				if (i == -1)
				{
					cerr << "line " << linenumber << ": Cannot find " << endmark4 << " in" << express << endl;
					exit(1);
				}
				firstword = express.substr(i, k - i);
				express_digit1.coffTLNT.push_back(StrToDou(firstword));
				express.erase(i, k-i + 8);
			}
			else if (express.find(endmark5) != string::npos) //*,   *T
			{
				//2.1E+01*T
				int k = express.find(endmark5);
				int i = findendmarkforward(express, endmark5, endmark1, 2);
				if (i == -1)
				{
					cerr << "line " << linenumber << ": Cannot find " << endmark5 << " in" << express << endl;
					exit(1);
				}
				firstword = express.substr(i, k - i);
				express_digit1.coffT.push_back(StrToDou(firstword));
				express_digit1.powerT.push_back(1);
				express.erase(i, k - i + 2);
			}
			else  // constant
			{
				express_digit1.coffT.push_back(StrToDou(express));
				express_digit1.powerT.push_back(0);
				express = { "" };
			}
		}
		return express_digit1;
	}
	// find Element in DBElements base on element name, and return the id;
	int Database::findElement(string elementname, vector<Element> Elements)
	{
		int i = 0;
		int elementid = -1;
		for (auto x : Elements)
		{
			if (x.name == elementname) elementid = i;
			i++;
		}
		// if phase not be found, exit, and print error
		//if (elementid == -1)
		//{
		//	cerr << "Element: " << elementname << "doesn't exit\n";
		//	return elementid;
			//exit(1);
		//}
		return elementid;
	}

	//
	void Database::showSystem()
	{
		cout << endl;
		cout << "Elements: " << SysElements.size() << endl;
		for (auto x : SysElements) cout << x.name << "  ";
		cout << "\n\n";
		cout << "Phases: " << SysPhases.size() << endl;
		for (auto x : SysPhases) cout << x.name << "  ";
		cout << endl;
	}

	// sort function for phase
	bool Database::SortPhase(const Phase& phase1, const Phase& phase2)  
	{
		return phase1.name < phase2.name;// ascending order
	}

} //end of VCLab