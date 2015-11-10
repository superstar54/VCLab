#include "..\include\MinEnergy.h"

namespace VCLab
{
	using namespace std;

	MinEnergy::~MinEnergy()
	{
		//cout << "Minimize energy modulue exit." << endl;
	}

	void MinEnergy::Solver(vector<Phase> & MEPhases, Condition Conditions, int ns, double (&ch)[MDim], int(&chb)[MDim])
	{
		int i;
		nph = MEPhases.size();
		nele = ns;
		for (int i = 0; i < nele; i++)
		{
			Chp[i] = ch[i];
			Chpb[i] = chb[i];
		}
		//
		SetupConditions(Conditions);
		//
		SetupVariables(MEPhases);
		//
		Initialization(MEPhases);
		if (debug >= 1)
		{
			for (i = 0; i < vn; i++)
				cout << var[i] << "  " << endl;
		}
		//
		Iteration(MEPhases);
		//
		CalcCompositions(MEPhases);
		// return chemical potential
		for (int i = 0; i < nele; i++)
			ch[i] = Chp[i];
	}

	void MinEnergy::SetupConditions(Condition Conditions)
	{
		MEConditions = Conditions;
	}

	void MinEnergy::SetupVariables(vector<Phase> & MEPhases)
	{
		int i, j;
  
		vn = 0; 
		for (i = 0; i < nph; i++)
		{
			//Phases amount
			MEPhases[i].vmidce = vn;
			vn++;
			//亚点阵拉格朗日乘子
			for (j = 0; j < MEPhases[i].vlan; j++)
			{
				MEPhases[i].vlaidce[j] = vn;
				vn++;
			}
			for (j = 0; j < MEPhases[i].vyn; j++)
			{
				MEPhases[i].vyidce[j] = vn;
				vn++;
			}
		}
		//Chemical potential
		for (i = 2; i < nele; i++)
		{
			Chpidce[i] = vn;
			vn++;
		}
		// temperature
		Tidce = vn;
		vn++;
	}

	void MinEnergy::AssignVariables(vector<Phase> & MEPhases)
	{
		int i, j;
		vector<Phase>::iterator iter;
		

		for (iter = MEPhases.begin(); iter != MEPhases.end();)
		{

			if (var[iter->vmidce] < eps && iter->vmb!=0)
			{
				iter->vmb = 0;
				iter->vm = 0;
				var[iter->vmidce] = 0;
				// remove this phase
				iter = MEPhases.erase(iter);
				nph--;
				// update the varibles
				SetupVariables(MEPhases);
				//
				Initialization(MEPhases);
				break;

			}
			else
			{
				if (var[iter->vmidce] > 1)
				{
					iter->vm = 1 - eps;
					var[iter->vmidce] = 1 - eps;
				}
				else
					iter->vm = var[iter->vmidce];

				for (j = 0; j < iter->vlan; j++)
				{
					iter->vla[j] = var[iter->vlaidce[j]];
				}

				for (j = 0; j < iter->vyn; j++)
				{
					//boundary
					if (var[iter->vyidce[j]] < eps)
					{
						iter->vy[j] = eps;
						var[iter->vyidce[j]] = eps;
						iter->y[iter->vyidc[j]] = eps;
					}
					else if (var[iter->vyidce[j]]>1)
					{
						iter->vy[j] = 1;
						var[iter->vyidce[j]] = 1;
						iter->y[iter->vyidc[j]] = 1;
					}
					else
					{
						iter->vy[j] = var[iter->vyidce[j]];
						iter->y[iter->vyidc[j]] = iter->vy[j];
					}
				}
				iter++;
			}
		}
		// chemical potential
		for (i = 2; i < nele; i++)
			Chp[i] = var[Chpidce[i]];
		// temperature
		MEConditions.T = var[Tidce];
	}

	void MinEnergy::Initialization(vector<Phase> & MEPhases)
	{
		int i, j;

		for (i = 2; i < nele; i++)
		{
			var[Chpidce[i]] = Chp[i];
		}
		// temperature
		var[Tidce] = MEConditions.T;

		for (i = 0; i < nph; i++)
		{
			//CalcGE(MEPhases[i], MEConditions.T);
			var[MEPhases[i].vmidce] = MEPhases[i].vm;
			for (j = 0; j < MEPhases[i].vyn; j++)
			{
				MEPhases[i].vy[j] = MEPhases[i].y[MEPhases[i].vyidc[j]];
				var[MEPhases[i].vyidce[j]] = MEPhases[i].vy[j];
				//MEPhases[i].vla[MEPhases[i].vyidla[j]] = (MEPhases[i].vla[MEPhases[i].vyidla[j]] 
				//	-MEPhases[i].dGFx[j]
				//	- Chp[MEPhases[i].vyide[j]] * MEPhases[i].vysp[j])/2;
			}

			for (j = 0; j < MEPhases[i].vlan; j++)
			{
				var[MEPhases[i].vlaidce[j]] = MEPhases[i].vla[j];
			}
		}
	}

	double MinEnergy::CalcToteF(vector<Phase> & MEPhases)
	{
		int i, j;
		double toteF = 0;

		for (i = 0; i < nph; i++)
		{
			//CalcGE(MEPhases[i], MEConditions.T);
			toteF += eF[MEPhases[i].vmidce];

			for (j = 0; j < MEPhases[i].vlan; j++)
			{
				toteF += eF[MEPhases[i].vlaidce[j]];
			}
		}
		// chemical potential
		for (i = 2; i < nele; i++)
		{
			if(Chpb[i]==1)
				toteF += eF[Chpidce[i]];
		}
		// temperature
		if(MEConditions.Tb==1)
			toteF += eF[Tidce];
		return toteF;
	}

	int MinEnergy::Iteration(vector<Phase> & MEPhases)
	{
		int i, j;
		int niter = 0;
		double totvar = 1e5;
		double toteF = 1e5, otoleF;
		double var0[MDim3];
		//
		while (totvar > eps9 || toteF > eps3)
		{
			niter++;
			totvar = 0;
			toteF = 0;
			//cout << "iter: " << niter << endl;
			if (niter > 5000)
			{
				cout << " Too much iteration" << endl; 
				return 0;
			}
			//
			for (i = 0; i < vn; i++)
				var0[i] = var[i];
			if (debug >= 1)
			{
				if (MEConditions.Tb == 1)
				{
					cout << "   0MinE: " << MEConditions.T << endl;
					cout << "    phases: " << MEPhases.size() << endl;;
				}
			}

			SetupHillert(MEPhases);
			toteF = CalcToteF(MEPhases);
			MendHillert(MEPhases);
			if (debug >= 2)
			{
				for (i = 0; i < mvn; i++)
					cout << setprecision(1) << scientific << meF[i] << "   ";
				cout << endl;
				for (i = 0; i < mvn; i++)
				{
					for (j = 0; j < mvn; j++)
						cout << setprecision(1) << scientific << meFJ[i][j] << "   ";
					cout << endl;
				}
			}
			

			//求解线性方程组得到dvar
			/*ofstream Out("A.txt");
			ofstream Outb("b.txt");

			for (i = 0;i < mvn;i++)
			{
				Outb << meF[i]<<"  ";
				for (j = 0;j < mvn;j++)
				{
					Out << meFJ[i][j]<<"  ";
				}
				Out << endl;
			}
			Out.close();
			Outb.close();
			*/
			if (gauss_elimination(mvn, meFJ, meF))
			{
				for (i = 0; i < mvn; i++)
				{
					dmeF[i] = meF[i];
					totvar += dmeF[i] * dmeF[i];
				}
				//
				reMendHillert(MEPhases);
				for (i = 0; i < vn; i++)
				{
					deF[i] = eF[i];
					var[i] = var0[i] - deF[i];
				}
			}
			else
			{
				PrintError("\n matrix singular, calculation failed...");
				return 0;
			}
			if (debug >= 1)
			{
				cout << "tolvar: " << totvar << endl;
				cout << "toleF: " << toteF << endl;
			}
			//
			AssignVariables(MEPhases);
		}
		/*
		for (i = 0; i < vn; i++)
			cout << var[i] << endl;
			*/
		return 1;
	}

	void MinEnergy::SetupHillert(vector<Phase> & MEPhases)
	{
		int i, j, k;
		int nF = 0;
		int tempi;
		double tempd, tempvm, sumx, tempdGFdx;

		for (i = 0; i < vn; i++)
		{
			eF[i] = 0;				
			for (j = 0; j < vn; j++)
				eFJ[i][j] = 0;		// jacobian matrix
		}
		// Phase matrix	
		for (i = 0; i < nph; i++)
		{
			CalcGE(MEPhases[i], MEConditions.T);
			tempvm = MEPhases[i].vm;
			//=======================================================================================
			//GE = sum(x*chp)
			nF = MEPhases[i].vmidce;
			eF[nF] += MEPhases[i].GF;  // G[i];
			for (j = 0; j < MEPhases[i].yn; j++)
			{
				if (MEPhases[i].yide[j] >= 2) // skip vacancy
				{
					tempd = MEPhases[i].ysp[j] * MEPhases[i].y[j];
					eF[nF] += Chp[MEPhases[i].yide[j]] * tempd;   // -Chp*sp*y
					eFJ[nF][Chpidce[MEPhases[i].yide[j]]] += tempd; //sp*y
				}
			}
			// to vy
			for (j = 0; j < MEPhases[i].vyn; j++)
			{
				tempi = MEPhases[i].vyidce[j];   //vyj对应c-e中的id
				eFJ[nF][tempi] += MEPhases[i].dGFx[j]
					+ Chp[MEPhases[i].vyide[j]] * MEPhases[i].vysp[j];   // Chp*sp
			}
			//to T
			eFJ[nF][Tidce] += MEPhases[i].dGFt;
			//=======================================================================================
			// sum(y) == 1
			for (j = 0; j < MEPhases[i].vyn; j++)
			{
				nF = MEPhases[i].vlaidce[MEPhases[i].vyidla[j]];
				eF[nF] += MEPhases[i].vy[j];
				eFJ[nF][MEPhases[i].vyidce[j]] += 1;
			}
			for (j = 0; j < MEPhases[i].vlan; j++)
			{
				/*% y1 + ... + yj + ... + y5 + ... = 0;
				% 0   ...    1   ...   0  */
				eF[MEPhases[i].vlaidce[j]] -= 1;  
			}
			//=======================================================================================
			//dGm/dy
			for (j = 0; j < MEPhases[i].vyn; j++)
			{
				nF = MEPhases[i].vyidce[j];   
				//m[i]*pG/py - m[i]*Chp[]*subper + lam, eliminate m
				tempdGFdx = MEPhases[i].dGFx[j];
				eF[nF] += tempdGFdx
					+ Chp[MEPhases[i].vyide[j]] * MEPhases[i].vysp[j] + MEPhases[i].vla[MEPhases[i].vyidla[j]];
				for (k = 0; k < MEPhases[i].vyn; k++)   //to vy_j
				{
					eFJ[nF][MEPhases[i].vyidce[k]] += MEPhases[i].d2GFx[j][k];
				}
				if (MEPhases[i].vyide[j] > 1 ) // skip vacancy
					eFJ[nF][Chpidce[MEPhases[i].vyide[j]]] += MEPhases[i].vysp[j]; //sp

				eFJ[nF][MEPhases[i].vmidce] += MEPhases[i].dGFx[j]
					+ Chp[MEPhases[i].vyide[j]] * MEPhases[i].vysp[j];
				eFJ[nF][MEPhases[i].vlaidce[MEPhases[i].vyidla[j]]] += 1;
				// to T
				if (MEConditions.Tb == 1)
					eFJ[nF][Tidce] += MEPhases[i].d2GFxt[j];
			}
			//=======================================================================================
			// overall compositions
			for (j = 0; j < MEPhases[i].yn; j++)
			{
				if (MEPhases[i].yide[j] <= 1) continue;
				nF = Chpidce[MEPhases[i].yide[j]];
				eF[nF] += tempvm*MEPhases[i].ysp[j] * MEPhases[i].y[j];   // m*sp*y
				if (MEPhases[i].vmb == 1)
				{
					eFJ[nF][MEPhases[i].vmidce] += MEPhases[i].ysp[j] * MEPhases[i].y[j];  // to m
				}
			}
			// to vy
			for (j = 0; j < MEPhases[i].vyn; j++)
			{
				if (MEPhases[i].vyide[j] <= 1) continue;
				nF = Chpidce[MEPhases[i].vyide[j]];
				eFJ[nF][MEPhases[i].vyidce[j]] += tempvm*MEPhases[i].vysp[j]; //
			}
			//=======================================================================================
			// others
		}
		// sum(y_i) - x_i = 0
		for (i = 2; i < nele; i++)
		{
			nF = Chpidce[i];
			eF[nF] = eF[nF];
			eF[nF] -= MEConditions.x[i];
		}
	}
	// Mend Hillert according conditions
	void MinEnergy::MendHillert(vector<Phase> & MEPhases)
	{
		int i, j, k;
		int nF = 0;
		int tempi;
		double tempd, tempvm;
		mvn = 0;
		//
		for (i = 0;i < MDim3;i++)
			meF[i] = eF[i];
		for (i = 0; i < nph; i++)
		{
			//Phases amount
			if (MEPhases[i].vmb == 1)
			{
				for (j = 0;j < vn;j++)
					meFJ[j][mvn] = eFJ[j][MEPhases[i].vmidce];
				mvn++;
			}

			for (j = 0; j < MEPhases[i].vlan; j++)
			{
				for (k = 0;k < vn;k++)
					meFJ[k][mvn] = eFJ[k][MEPhases[i].vlaidce[j]];
				mvn++;
			}
			//
			for (j = 0; j < MEPhases[i].vyn; j++)
			{
				for (k = 0; k < vn; k++)
					meFJ[k][mvn] = eFJ[k][MEPhases[i].vyidce[j]];
				mvn++;
			}
		}
		for (i = 2; i < nele; i++)
		{
			if (Chpb[i] == 1)
			{
				for (j = 0;j < vn;j++)
					meFJ[j][mvn] = eFJ[j][Chpidce[i]];
				mvn++;
			}
			
		}
		// temperature
		if (MEConditions.Tb == 1)
		{
			for (j = 0;j < vn;j++)
				meFJ[j][mvn] = eFJ[j][Tidce];
			mvn++;
		}
	}
	//
	//根据条件修改Hillert矩阵
	void MinEnergy::reMendHillert(vector<Phase> & MEPhases)
	{
		int i, j, k;
		int nF = 0;
		int tempi;
		double tempd, tempvm;
		mvn = 0;
		for (i = 0; i < nph; i++)
		{
			//Phases amount
			if (MEPhases[i].vmb == 1)
			{
				eF[MEPhases[i].vmidce] = meF[mvn];
				mvn++;
			}
			else
				eF[MEPhases[i].vmidce] = 0;
			for (j = 0; j < MEPhases[i].vlan; j++)
			{
				eF[MEPhases[i].vlaidce[j]] = meF[mvn];
				mvn++;
			}
			//
			for (j = 0; j < MEPhases[i].vyn; j++)
			{
				eF[MEPhases[i].vyidce[j]] = meF[mvn];
				mvn++;
			}
		}
		for (i = 2; i < nele; i++)
		{
			if (Chpb[i] == 1)
			{
				eF[Chpidce[i]] = meF[mvn];
				mvn++;
			}
			else
				eF[Chpidce[i]] = 0;

		}
		// temperature
		if (MEConditions.Tb == 1)
		{
			eF[Tidce] = meF[mvn];
			mvn++;
		}
		else
			eF[Tidce] = 0;
	}
	//
	void MinEnergy::ShowEquilibrium()
	{

	}
	void MinEnergy::ShowCondition()
	{

	}
	//
	void MinEnergy::CalcCompositions(vector<Phase> & MEPhases)
	{
		int i, j;
		double tempx;

		for (i = 0; i < nph; i++)
		{
			for (j = 0; j < nele; j++)
				MEPhases[i].x[j] = 0;
			for (j = 0; j < MEPhases[i].yn; j++)
			{
				MEPhases[i].x[MEPhases[i].yide[j] - 2] += MEPhases[i].ysp[j] * MEPhases[i].y[j];
			}

			tempx = 0;
			for (j = 0; j < nele; j++)
			{
				tempx += MEPhases[i].x[j];
			}
			for (j = 0; j < nele; j++)
			{
				MEPhases[i].x[j] = MEPhases[i].x[j] / tempx;
			}
		}
	}
}// end of namespace