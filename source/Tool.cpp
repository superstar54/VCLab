<<<<<<< HEAD

#include "..\include\Tool.h"

namespace VCLab
{
	//*****************************   string tool    ****************************//
	//
	string readline(fstream &Inp, int & linenumber)
	{

		string line, mlines = "";
		while (!Inp.eof())
		{
			getline(Inp, line);	linenumber++;	// read a line
												//if (debug >= 1) cout << "line: " << linenumber << endl;
			if (line.length() == 0) continue;		// skip the empty line
			line = capital(line);				// change to capital£¬change "Tab" to "space",remove the space before and behind
			if (line[0] == '$') continue;		// skip the comment line
			mlines = line;
			while (line.find('!') == string::npos&&!Inp.eof())		//if end line mark "!" not found
			{
				getline(Inp, line);	linenumber++; // read another line
												  //if (debug >= 1)	cout << "line: " << linenumber << endl;
				if (line.length() == 0) continue;	// skip the empty line
				line = capital(line);				// change to capital£¬change "Tab" to "space",remove the space before and behind
				if (line[0] == '$') continue;		// skip the comment line
				mlines += " ";
				mlines += line;  // merge lines
			}
			mlines = mlines.substr(0, mlines.length() - 1);	//delete "!"
			return mlines;
		}
		return mlines;
	}


	// this function is very important!!!
	// extract the first part base on given endmark;
	string ReadFirWor(string & tline, string endmark, int flag)
	{
		string firstword;
		string templine;
		// delete spacing before the string
		while (tline[0] == ' ')
		{
			tline = tline.substr(1, tline.length());
		}
		//find the first endmark
		firstword = tline.substr(0, tline.find(endmark));
		//delete the firstword in the tline, and creat a new tline.
		if (flag != 2)
			if (tline.find(endmark) != string::npos)
			{
				if (flag == 0)   // delete firstword and endmark
					tline = tline.substr(firstword.length() + endmark.length(), tline.length());
				else if (flag == 1)  // only delete firstword
					tline = tline.substr(firstword.length(), tline.length());
			}
			else
			{
				tline = { "" };   //if endmark are not found in the tline, mean reach the end of the tline, after that, tline become null
			}
		// eliminate the space
		int nspace = 0;
		for (int i = firstword.length() - 1;i >= 0;i--)
		{
			if (firstword[i] == ' ')
				nspace++;
			else
				break;
		}
		firstword = firstword.substr(0, firstword.length() - nspace);
		return firstword;
	}
	//mend the express for cases: .018531982*T**2, -.001884662*T**2
	string MendExp(string firstword)
	{
		string express_new;
		if (firstword[0] != '-' && firstword[0] != '+')
		{
			string insertstring = { "+" };
			firstword.insert(0, insertstring);  // add "0", -.3 became -0.3
		}
		express_new = firstword;
		int j = 0;
		int flength;
		flength = firstword.length();
		// eliminate the space
		for (int i = 0; i < flength; i++)
		{
			if (firstword[i] == ' ')
			{
				express_new.erase(i - j, 1);
				j++;
			}
		}
		return express_new;
	}

	//change lower to capital, change "\t" to " ", that is tab to spacing.
	string capital(string tline)
	{
		int tlength;
		tlength = tline.length();
		for (int i = 0;i < tlength;i++)
		{
			if (tline[i] == '\t') tline[i] = ' ';
			tline[i] = toupper(tline[i]);
		}
		// delete spacing before the string
		while (tline[0] == ' ')
		{
			tline.erase(0,1);
		}
		if (tline.length() != 0)
		{
			while (tline[tline.length() - 1] == ' ')
			{
				tline.erase(tline.length() - 1, 1);
			}
		}
		
		return tline;
	}
	//
	int findendmarkforward(string express, string endmark, string endmark_array[], int n)
	{
		int i;
		string firstwrod;
		for (i = express.find(endmark);i >= 0;i--)
		{
			for (int j = 0;j < n;j++)
				if (express[i] == endmark_array[j][0])
				{
					if (i > 0)
					{
						if (express[i - 1] != 'E' && express[i - 1] != '(')
						{
							return i;
						}
					}
					else
					{
						return i;
					}

				}
		}
		return -1;
	}
	//
	int findendmarkbackward(string express, string endmark, string endmark_array[], int n)
	{
		int i;
		int endi = express.length();
		for (i = express.find(endmark);i<endi;i++)
		{
			if (i != endi - 1)
				for (int j = 0;j < n;j++)
				{
					if (express[i] == endmark_array[j][0])
					{
						return i;
					}
				}

		}
		if (i == express.length()) return i; // reach the end
		return -1;
	}

	double StrToDou(string strval)
	{
		double dval;
		stringstream strout;
		strout << strval;
		strout >> dval;
		return dval;
	}

	int StrToInt(string strval)
	{
		int dval;
		stringstream strout; 
		strout << strval;
		strout >> dval;
		return dval;
	}

	//*****************************   matrix tool    ****************************//
	int gauss_elimination(int n, double (&a)[MDim3][MDim3], double (&b)[MDim3]) 
	{
		int i, j, k, row;
		double maxp, t;
		for (k = 0; k<n; k++) 
		{
			for (maxp = 0, i = k; i<n; i++)
				if (fabs(a[i][k])>fabs(maxp)) // find the maxium element
					maxp = a[row = i][k];
			if (fabs(maxp) < eps100)
			{
				//maxp = eps;
				//a[i][k] = eps;
				return 0;
			}
			if (row != k) 
			{
				for (j = k; j<n; j++)
					t = a[k][j], a[k][j] = a[row][j], a[row][j] = t;
				t = b[k], b[k] = b[row], b[row] = t;
			}
			for (j = k + 1; j<n; j++)
				for (a[k][j] /= maxp, i = k + 1; i<n; i++)
					a[i][j] -= a[i][k] * a[k][j];
			for (b[k] /= maxp, i = k + 1; i<n; i++)
				b[i] -= b[k] * a[i][k];
		}
		for (i = n - 1; i >= 0; i--)
			for (j = i + 1; j<n; j++)
				b[i] -= a[i][j] * b[j];
		return 1;
	}

	//
	void combination::combinations_recursive(int *pos, int depth, int margin)
	{
		// choice m element
		if (depth >= combm)
		{
			for (int ii = 0; ii < combm; ii++)
			{
				comb[ncount][ii] = pos[ii];
			}
			ncount++;
			//cout << endl;
		}
		else
		{
			// add a new element in right side
			for (int ii = margin; ii < combn; ii++) 
			{
				pos[depth] = ii;
				combinations_recursive(pos, depth + 1, ii + 1);
			}
		}
	}

	void combination::combinations(int mod, double d, int n, int m)
	{
		
		int i;

		combn = n;
		combm = m;
		mode = mod;
		dy = d;

		if (m<0 || m>n)
		{
			cout << "Error in combinations: n = " << n << " m = " << m << endl;
			exit(1);
		}
		else if(m == 0)
		{
			toln = 1;
			comb = new double*[toln];
			comb[0] = new double[m + 1];
			comb[0][0] = 1;
			return;
		}

		int *pos;
		pos = new int[m];
		
		// for multiplication, should initialize to 1
		toln = 1;
		ncount = 0;
		for (i = 0; i < m; i++)
		{
			toln = toln*(n - i) / (i + 1);
		}
		comb = new double*[toln];
		for (i = 0; i < toln; i++)
		{
			if (mode == 1)
			{
				comb[i] = new double[m + 1]; 
			}
			else
			{
				comb[i] = new double[m];
			}
		}
		combinations_recursive(pos, 0, 0);
		if (mode == 1)
		{
			for (i = 0; i < toln; i++)
			{
				for (int ii = 0; ii < combm; ii++)
				{
					comb[i][ii] -= ii;
				}
				comb[i][combm] = (combn - combm - comb[i][combm - 1])*dy + eps;
				for (int ii = (combm - 1); ii > 0; ii--)
				{
					comb[i][ii] = (comb[i][ii] - comb[i][ii - 1])*dy + eps;
				}
				comb[i][0] = comb[i][0] * dy + eps;
			}
		}
		delete[] pos;
	}
	void combination::show()
	{
		cout << "Total combinations: " << toln << endl;
		for (int i = 0; i < toln; i++)
		{
			if (mode == 1)
			{
				for (int j = 0; j < combm + 1; j++)
					cout << comb[i][j] << "  ";
			}
			else
			{
				for (int j = 0; j < combm; j++)
					cout << comb[i][j] << "  ";
			}

			cout << endl;
		}
	}
	void combination::deletecomb()
	{
		for (int i = 0; i < toln; i++)
			delete[] comb[i];
		delete[] comb;
	}

	void cart_product(int **output, int om, int **input, int n[], int m)
	{
		int i, j, k, ki;
		int nc;
		int on = m;
		int om1, on1;
		int **temp;
		nc = 0;
		temp = new int*[om];
		for (i = 0; i < om; i++)
			temp[i] = new int[on];
		for (i = 0; i < m; i++)
		{
			om1 = 1;
			for (j = 0; j < i; j++)
				om1 = om1*n[j];
			on1 = i;
			//
			for (j = 0; j < om1; j++)
			{
				for (k = 0; k < on1; k++)
				{
					temp[j][k] = output[j][k];
				}
			}
			//
			nc = 0;
			for (j = 0; j < om1; j++)
			{
				for (k = 0; k < n[i]; k++)
				{
					for (ki = 0; ki < on1; ki++)
					{
						output[nc][ki] = temp[j][ki];
						//cout << output[nc][ki] << "  ";
					}
					output[nc][i] = input[i][k];
					//cout << output[nc][i] << endl;
					nc++;
				}
			}
		}
		for (i = 0; i < om; i++)
			delete[] temp[i];
		delete[] temp;

	}

	void readline(fstream& Inp, string &tline, int &linecount)
	{
		string line;
		while (!Inp.eof())
		{
			getline(Inp, line);
			linecount++;
			line = capital(line);	//
			if (line.length() == 0) continue;
			if (line[0] == '$') continue;
			tline = line;
			while (line.find('!') == string::npos&&!Inp.eof())		//merge "!"
			{
				getline(Inp, line);
				linecount++;
				line = capital(line);								//
				if (line.length() == 0) continue;					//
				if (line[0] == '$') continue;						//
				tline += line;
			}
			tline = tline.substr(0, tline.length() - 1);			//delete "!"
			break;
		}
	}
	//======================================================================================
	void PrintTitle()
	{
		//software information
		PrintLine("=");
		PrintMiddle("VCLab (Virtual CALPHAD Laboratory)");
		cout << "\n  An open-source software for CALPHAD Calculations\n\n";
		cout << "  Contact: Xing Wang,  xingwang@csu.edu.cn \n";
		cout << "  Central South University, China \n\n";
		//cout<< "   Empa - Swiss Federal Laboratories for Materials Science and Technology" << endl
		cout << "  Update in https://github.com/superstar54/VCLab \n\n";
		cout << "  2015 - , Licensed under GNU \n";
		PrintLine("=");
		cout << endl << endl;
	}
	// Print a line with str
	void PrintLine(string str)
	{
		cout << setfill(str.c_str()[0]) << setw(60) << "" << endl;
	}
	void PrintLineFile(string str, fstream &Out)
	{
		Out << setfill(str.c_str()[0]) << setw(60) << "" << endl;
	}
	// Print str in the middle position of the line
	void PrintMiddle(string str)
	{
		int ns = (60 - str.length()) / 2;
		cout << setfill(' ') << setw(ns) << "" <<str << endl;
	}
	// Print Vstr : str
	void PrintVstrstr(string vstr, string str)
	{
		int ns = 20 - vstr.length();
		cout <<vstr << setfill(' ') << setw(ns) << "" <<": " <<str<< endl;
	}
	// Print Vstr : int
	void PrintVstrint(string vstr, int i)
	{
		int ns = 20 - vstr.length();
		cout << vstr << setfill(' ') << setw(ns) << "" << ": " << i << endl;
	}
	// Print warning information
	void PrintWarning(string str)
	{
		cout << endl;
		PrintLine("-");
		cout << "    W       W       W       A        RRRRRR     NN     N  !!!\n";
		cout << "     W      W      W       A A       R     R    N N    N  !!!\n";
		cout << "      W    WWW    W       A   A      R     R    N  N   N  !!!\n";
		cout << "       W  W   W  W       AAAAAAA     RRRRRR     N   N  N  !!!\n";
		cout << "        WW     WW       A      A     R    R     N    N N  !!!\n";
		cout << "        WW     WW      A        A    R     R    N     NN   ! \n";
		cout << "        WW     WW     A          A   R      R   N      N  !!!\n\n";

		cout << "\n**: " << str << endl;
	}
	// Print error information
	void PrintError(string str)
	{

		cout << endl;
		PrintLine("-");
		cout << "   EEEEEEE RRRRRR   RRRRRR     OOOOO   RRRRRR    !!!\n";
		cout << "   E       R     R  R     R   O     O  R     R   !!!\n";
		cout << "   E       R     R  R     R  O       O R     R   !!!\n";
		cout << "   EEEEEE  RRRRRR   RRRRRR   O       O RRRRRR    !!!\n";
		cout << "   E       R    R   R    R   O       O R    R    !!!\n";
		cout << "   E       R     R  R     R   O     O  R     R    !\n";
		cout << "   EEEEEEE R      R R      R   OOOOO   R      R  !!!\n\n";
		cout << "   " << str << endl;
		cout << "   Calculation terminated! \n";
		PrintLine("-");
		exit(1);
	}
}
=======

#include "..\include\Tool.h"

namespace VCLab
{
	//*****************************   string tool    ****************************//
	//
	string readline(fstream &Inp, int & linenumber)
	{

		string line, mlines = "";
		while (!Inp.eof())
		{
			getline(Inp, line);	linenumber++;	// read a line
												//if (debug >= 1) cout << "line: " << linenumber << endl;
			if (line.length() == 0) continue;		// skip the empty line
			line = capital(line);				// change to capital£¬change "Tab" to "space",remove the space before and behind
			if (line[0] == '$') continue;		// skip the comment line
			mlines = line;
			while (line.find('!') == string::npos&&!Inp.eof())		//if end line mark "!" not found
			{
				getline(Inp, line);	linenumber++; // read another line
												  //if (debug >= 1)	cout << "line: " << linenumber << endl;
				if (line.length() == 0) continue;	// skip the empty line
				line = capital(line);				// change to capital£¬change "Tab" to "space",remove the space before and behind
				if (line[0] == '$') continue;		// skip the comment line
				mlines += " ";
				mlines += line;  // merge lines
			}
			mlines = mlines.substr(0, mlines.length() - 1);	//delete "!"
			return mlines;
		}
		return mlines;
	}


	// this function is very important!!!
	// extract the first part base on given endmark;
	string ReadFirWor(string & tline, string endmark, int flag)
	{
		string firstword;
		string templine;
		// delete spacing before the string
		while (tline[0] == ' ')
		{
			tline = tline.substr(1, tline.length());
		}
		//find the first endmark
		firstword = tline.substr(0, tline.find(endmark));
		//delete the firstword in the tline, and creat a new tline.
		if (flag != 2)
			if (tline.find(endmark) != string::npos)
			{
				if (flag == 0)   // delete firstword and endmark
					tline = tline.substr(firstword.length() + endmark.length(), tline.length());
				else if (flag == 1)  // only delete firstword
					tline = tline.substr(firstword.length(), tline.length());
			}
			else
			{
				tline = { "" };   //if endmark are not found in the tline, mean reach the end of the tline, after that, tline become null
			}
		// eliminate the space
		int nspace = 0;
		for (int i = firstword.length() - 1;i >= 0;i--)
		{
			if (firstword[i] == ' ')
				nspace++;
			else
				break;
		}
		firstword = firstword.substr(0, firstword.length() - nspace);
		return firstword;
	}
	//mend the express for cases: .018531982*T**2, -.001884662*T**2
	string MendExp(string firstword)
	{
		string express_new;
		if (firstword[0] != '-' && firstword[0] != '+')
		{
			string insertstring = { "+" };
			firstword.insert(0, insertstring);  // add "0", -.3 became -0.3
		}
		express_new = firstword;
		int j = 0;
		int flength;
		flength = firstword.length();
		// eliminate the space
		for (int i = 0; i < flength; i++)
		{
			if (firstword[i] == ' ')
			{
				express_new.erase(i - j, 1);
				j++;
			}
		}
		return express_new;
	}

	//change lower to capital, change "\t" to " ", that is tab to spacing.
	string capital(string tline)
	{
		int tlength;
		tlength = tline.length();
		for (int i = 0;i < tlength;i++)
		{
			if (tline[i] == '\t') tline[i] = ' ';
			tline[i] = toupper(tline[i]);
		}
		// delete spacing before the string
		while (tline[0] == ' ')
		{
			tline.erase(0,1);
		}
		if (tline.length() != 0)
		{
			while (tline[tline.length() - 1] == ' ')
			{
				tline.erase(tline.length() - 1, 1);
			}
		}
		
		return tline;
	}
	//
	int findendmarkforward(string express, string endmark, string endmark_array[], int n)
	{
		int i;
		string firstwrod;
		for (i = express.find(endmark);i >= 0;i--)
		{
			for (int j = 0;j < n;j++)
				if (express[i] == endmark_array[j][0])
				{
					if (i > 0)
					{
						if (express[i - 1] != 'E' && express[i - 1] != '(')
						{
							return i;
						}
					}
					else
					{
						return i;
					}

				}
		}
		return -1;
	}
	//
	int findendmarkbackward(string express, string endmark, string endmark_array[], int n)
	{
		int i;
		int endi = express.length();
		for (i = express.find(endmark);i<endi;i++)
		{
			if (i != endi - 1)
				for (int j = 0;j < n;j++)
				{
					if (express[i] == endmark_array[j][0])
					{
						return i;
					}
				}

		}
		if (i == express.length()) return i; // reach the end
		return -1;
	}

	double StrToDou(string strval)
	{
		double dval;
		stringstream strout;
		strout << strval;
		strout >> dval;
		return dval;
	}

	int StrToInt(string strval)
	{
		int dval;
		stringstream strout; 
		strout << strval;
		strout >> dval;
		return dval;
	}

	//*****************************   matrix tool    ****************************//
	int gauss_elimination(int n, double (&a)[MDim3][MDim3], double (&b)[MDim3]) 
	{
		int i, j, k, row;
		double maxp, t;
		for (k = 0; k<n; k++) 
		{
			for (maxp = 0, i = k; i<n; i++)
				if (fabs(a[i][k])>fabs(maxp)) // find the maxium element
					maxp = a[row = i][k];
			if (fabs(maxp) < eps100)
			{
				//maxp = eps;
				//a[i][k] = eps;
				return 0;
			}
			if (row != k) 
			{
				for (j = k; j<n; j++)
					t = a[k][j], a[k][j] = a[row][j], a[row][j] = t;
				t = b[k], b[k] = b[row], b[row] = t;
			}
			for (j = k + 1; j<n; j++)
				for (a[k][j] /= maxp, i = k + 1; i<n; i++)
					a[i][j] -= a[i][k] * a[k][j];
			for (b[k] /= maxp, i = k + 1; i<n; i++)
				b[i] -= b[k] * a[i][k];
		}
		for (i = n - 1; i >= 0; i--)
			for (j = i + 1; j<n; j++)
				b[i] -= a[i][j] * b[j];
		return 1;
	}

	//
	void combination::combinations_recursive(int *pos, int depth, int margin)
	{
		// choice m element
		if (depth >= combm)
		{
			for (int ii = 0; ii < combm; ii++)
			{
				comb[ncount][ii] = pos[ii];
			}
			ncount++;
			//cout << endl;
		}
		else
		{
			// add a new element in right side
			for (int ii = margin; ii < combn; ii++) 
			{
				pos[depth] = ii;
				combinations_recursive(pos, depth + 1, ii + 1);
			}
		}
	}

	void combination::combinations(int mod, double d, int n, int m)
	{
		
		int i;

		combn = n;
		combm = m;
		mode = mod;
		dy = d;

		if (m<0 || m>n)
		{
			cout << "Error in combinations: n = " << n << " m = " << m << endl;
			exit(1);
		}
		else if(m == 0)
		{
			toln = 1;
			comb = new double*[toln];
			comb[0] = new double[m + 1];
			comb[0][0] = 1;
			return;
		}

		int *pos;
		pos = new int[m];
		
		// for multiplication, should initialize to 1
		toln = 1;
		ncount = 0;
		for (i = 0; i < m; i++)
		{
			toln = toln*(n - i) / (i + 1);
		}
		comb = new double*[toln];
		for (i = 0; i < toln; i++)
		{
			if (mode == 1)
			{
				comb[i] = new double[m + 1]; 
			}
			else
			{
				comb[i] = new double[m];
			}
		}
		combinations_recursive(pos, 0, 0);
		if (mode == 1)
		{
			for (i = 0; i < toln; i++)
			{
				for (int ii = 0; ii < combm; ii++)
				{
					comb[i][ii] -= ii;
				}
				comb[i][combm] = (combn - combm - comb[i][combm - 1])*dy + eps;
				for (int ii = (combm - 1); ii > 0; ii--)
				{
					comb[i][ii] = (comb[i][ii] - comb[i][ii - 1])*dy + eps;
				}
				comb[i][0] = comb[i][0] * dy + eps;
			}
		}
		delete[] pos;
	}
	void combination::show()
	{
		cout << "Total combinations: " << toln << endl;
		for (int i = 0; i < toln; i++)
		{
			if (mode == 1)
			{
				for (int j = 0; j < combm + 1; j++)
					cout << comb[i][j] << "  ";
			}
			else
			{
				for (int j = 0; j < combm; j++)
					cout << comb[i][j] << "  ";
			}

			cout << endl;
		}
	}
	void combination::deletecomb()
	{
		for (int i = 0; i < toln; i++)
			delete[] comb[i];
		delete[] comb;
	}

	void cart_product(int **output, int om, int **input, int n[], int m)
	{
		int i, j, k, ki;
		int nc;
		int on = m;
		int om1, on1;
		int **temp;
		nc = 0;
		temp = new int*[om];
		for (i = 0; i < om; i++)
			temp[i] = new int[on];
		for (i = 0; i < m; i++)
		{
			om1 = 1;
			for (j = 0; j < i; j++)
				om1 = om1*n[j];
			on1 = i;
			//
			for (j = 0; j < om1; j++)
			{
				for (k = 0; k < on1; k++)
				{
					temp[j][k] = output[j][k];
				}
			}
			//
			nc = 0;
			for (j = 0; j < om1; j++)
			{
				for (k = 0; k < n[i]; k++)
				{
					for (ki = 0; ki < on1; ki++)
					{
						output[nc][ki] = temp[j][ki];
						//cout << output[nc][ki] << "  ";
					}
					output[nc][i] = input[i][k];
					//cout << output[nc][i] << endl;
					nc++;
				}
			}
		}
		for (i = 0; i < om; i++)
			delete[] temp[i];
		delete[] temp;

	}

	void readline(fstream& Inp, string &tline, int &linecount)
	{
		string line;
		while (!Inp.eof())
		{
			getline(Inp, line);
			linecount++;
			line = capital(line);	//
			if (line.length() == 0) continue;
			if (line[0] == '$') continue;
			tline = line;
			while (line.find('!') == string::npos&&!Inp.eof())		//merge "!"
			{
				getline(Inp, line);
				linecount++;
				line = capital(line);								//
				if (line.length() == 0) continue;					//
				if (line[0] == '$') continue;						//
				tline += line;
			}
			tline = tline.substr(0, tline.length() - 1);			//delete "!"
			break;
		}
	}
	//======================================================================================
	void PrintTitle()
	{
		//software information
		PrintLine("=");
		PrintMiddle("VCLab (Virtual CALPHAD Laboratory)");
		cout << "\n  An open-source software for CALPHAD Calculations\n\n";
		cout << "  Contact: Xing Wang,  xingwang@csu.edu.cn \n";
		cout << "  Central South University, China \n\n";
		//cout<< "   Empa - Swiss Federal Laboratories for Materials Science and Technology" << endl
		cout << "  Update in https://github.com/superstar54/VCLab \n\n";
		cout << "  2015 - , Licensed under GNU \n";
		PrintLine("=");
		cout << endl << endl;
	}
	// Print a line with str
	void PrintLine(string str)
	{
		cout << setfill(str.c_str()[0]) << setw(60) << "" << endl;
	}
	// Print str in the middle position of the line
	void PrintMiddle(string str)
	{
		int ns = (60 - str.length()) / 2;
		cout << setfill(' ') << setw(ns) << "" <<str << endl;
	}
	// Print Vstr : str
	void PrintVstrstr(string vstr, string str)
	{
		int ns = 20 - vstr.length();
		cout <<vstr << setfill(' ') << setw(ns) << "" <<": " <<str<< endl;
	}
	// Print Vstr : int
	void PrintVstrint(string vstr, int i)
	{
		int ns = 20 - vstr.length();
		cout << vstr << setfill(' ') << setw(ns) << "" << ": " << i << endl;
	}
	// Print warning information
	void PrintWarning(string str)
	{

		cout << "\n***Warning: " << str << endl;
	}
	// Print error information
	void PrintError(string str)
	{

		cout << endl;
		PrintLine("-");
		cout << "   EEEEEEE RRRRRR   RRRRRR     OOOOO   RRRRRR    !!!\n";
		cout << "   E       R     R  R     R   O     O  R     R   !!!\n";
		cout << "   E       R     R  R     R  O       O R     R   !!!\n";
		cout << "   EEEEEE  RRRRRR   RRRRRR   O       O RRRRRR    !!!\n";
		cout << "   E       R    R   R    R   O       O R    R    !!!\n";
		cout << "   E       R     R  R     R   O     O  R     R    !\n";
		cout << "   EEEEEEE R      R R      R   OOOOO   R      R  !!!\n\n";
		cout << "   " << str << endl;
		cout << "   Calculation terminated! \n";
		PrintLine("-");
		exit(1);
	}
}
>>>>>>> b34a88e731ed6049f6c364b37b469ae02b913009
