#include "..\include\VCLInput.h"

namespace VCLab
{

void ProjectInput::ReadVCLInput()
{
	int		i;
	int		linecount	= 0;
	double	sumx		= 0;
	string	mlines;				// merge lines until "!"
	string	keyword;			//  
	string	strvalue;			// value of the key word
	string	endmark		= { "=" };

	fstream Inp(FilPro.c_str(), ios::in);
	if (!Inp)
		PrintError("File " + FilPro + "could not be opened.");
	
	// read file line by line
	while (!Inp.eof())
	{
		mlines = readline(Inp, linecount);					// read and merge lines until "!"
		keyword = ReadFirWor(mlines, endmark, 0);		// read first word of mlines
		if (keyword == "VCLMODE")
		{
			VCLmode = ReadFirWor(mlines, endmark, 0);
			continue;
		}
		else if (keyword == "DATABASE_FILENAME")
		{
			
			FilTDB = ReadFirWor(mlines, endmark, 0); 	
			continue;
		}
		else if (keyword == "ELEMENTS")
		{
			EleInp = ReadFirWor(mlines, endmark, 0); 
			continue;
		}
		else if (keyword == "COMPOSITIONS")
		{
			ComInp = ReadFirWor(mlines, endmark, 0); 
			continue;
		}
		else if (keyword == "TEMPERATURE")
		{
			Cond.T = StrToDou(ReadFirWor(mlines, endmark, 0)); 
			continue;
		}
		else if (keyword == "PRESSURE")
		{
			Cond.P = StrToDou(ReadFirWor(mlines, endmark, 0));
			continue;
		}
		else if (keyword == "RESULT_OUTPUT_FILE")
		{
			FilOut = ReadFirWor(mlines, endmark, 0);
			continue;
		}
		else if (keyword == "GLOBAL_GRID_INTERVAL")
		{
			Cond.dx = StrToDou(ReadFirWor(mlines, endmark, 0)); 
			continue;
		}
		if (keyword == "VARIBLES")
		{
			Var = ReadFirWor(mlines, endmark, 0); 
			continue;
		}
		else if (keyword == "V_START")
		{
			Var_s = StrToDou(ReadFirWor(mlines, endmark, 0));
			continue;
		}
		else if (keyword == "V_END")
		{
			Var_e = StrToDou(ReadFirWor(mlines, endmark, 0));
			continue;
		}
		else if (keyword == "V_INTERVAL")
		{
			Var_d = StrToDou(ReadFirWor(mlines, endmark, 0));
			continue;
		}

		else if (keyword == "PHASES_FIXED")
		{
			FixPhaInp = ReadFirWor(mlines, endmark, 0);
			continue;
		}
	}
	while (EleInp.length() != 0)	// read elements one by one
		EleList.push_back(ReadFirWor(EleInp, { "," }, 0));
	while (ComInp.length() != 0)	// read composition one by one
		ComList.push_back(StrToDou(ReadFirWor(ComInp, { "," }, 0)));
	SortEle(EleList, ComList);		// Sort elements, composition in alphabetical order
	while (FixPhaInp.length() != 0) // read fixed phases one by one
		FixPhaList.push_back(ReadFirWor(FixPhaInp, { "," }, 0));

	// nele == ncom, this is not ture, when compositions are varibles.
	// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>on going work
	if (EleList.size() != ComList.size())
	{
		PrintError("number of elements is not equal to compositions!");
	}
	//
	Cond.nele = EleList.size();
	sumx = 0;
	for (i = 0; i < Cond.nele; i++)
	{
		Cond.ele[i] = EleList[i];
		Cond.x[i] = ComList[i];
		sumx += Cond.x[i];
	}
	//
	if (fabs(sumx-1.0)>eps3)
	{
		cout << sumx;
		PrintError("overall composition is not equal to 1!");
	}
	Inp.close();
}

void ProjectInput::SortEle(vector<string> &EleList, vector<double> &composilist)
{
	int i;
	double tempc;
	string tempe;
	for (i = 3; i < EleList.size(); i++)
	{
		if (EleList[i] < EleList[i-1])
		{
			tempe = EleList[i];
			EleList[i] = EleList[i - 1];
			EleList[i - 1] = tempe;
			tempc = composilist[i];
			composilist[i] = composilist[i - 1];
			composilist[i - 1] = tempc;
		}
	}
}

double ProjectInput::StrToDou(string strval)
{
	double dval;
	stringstream strout;  //建立输出字符串流,与数组c建立关联,缓冲区长
	strout << strval;
	strout >> dval;
	return dval;
}

} // end of VCLab