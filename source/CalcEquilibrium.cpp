#include "..\include\CalcEquilibrium.h"

namespace VCLab
{
	using namespace std;

	CalcEquilibrium::~CalcEquilibrium()
	{
		//cout << "Calculate Equilibrium modulue exit." << endl;
	}
	// 求解给定的 phase set
	void CalcEquilibrium::Solver(int GPmod, vector<Phase> SysPhases, Condition CLCondition)
	{
		int i;
		//清除数据
		CEPhases.clear();
		SPhases.clear();
		MPhases.clear();
		for (i = 0;i < MDim;i++)
			Chpb[i] = 1;
		//读取数据
		CEConditions = CLCondition;
		CEPhases = SysPhases;
		Phn = CEPhases.size();
		GPdx = CLCondition.dx;
		nele = CLCondition.nele - 2;
		GPmode = GPmod;
		if (GPmode == 1)
		{
			//建立Grid point
			SetupGP();
		}
		// 确定Phase sets
		SetupPhaseSets(SysPhases);
		// 迭代Phase set
		//Iteration();
		//求解平衡相集的平衡
		start = clock();
		MinE.Solver(SPhases, CEConditions, CLCondition.nele, Chp, Chpb);
		//update Phns
		Phns = SPhases.size();
		CEConditions.T = MinE.MEConditions.T;
		finish = clock();
		totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
		if (debug >= 1)
		{
			cout << "Newton solution in " << totaltime << "s." << endl;
		}
		// nomalization phase amount
		Initialization();
		//ShowEquilibrium();
	}
	//
	void CalcEquilibrium::SetupConditions(Condition CLCondition)
	{
		CEConditions = CLCondition;
	}
	// 初始化，相选择，点阵分数，化学势
	void CalcEquilibrium::Initialization(vector<Phase> SysPhases)
	{
		//从数据库的SysPhases中提取相
		
	}
	// 确定Phase sets
	void CalcEquilibrium::SetupPhaseSets(vector<Phase> SysPhases)
	{
		int i, j;
		int Dim, phaseid;
		double natom;
		if (GPmode == 1)
		{
			// use Grid Hull
			Dim = CEConditions.nele - 2;
			Phns = 0;
			for (i = 0; i < Dim; i++)
			{
				if (CHull.f.p[i].fb != 0 && CHull.f.p[i].phfrac>eps9)
				{
					phaseid = FindPhase(CEPhases[CHull.f.p[i].phid].name, SPhases);
					if (phaseid != -1)
					{
						if (CEPhases[CHull.f.p[i].phid].x[0] > SPhases[phaseid].x[0])
						{
							CEPhases[CHull.f.p[i].phid].name += "#2";
							//SPhases[phaseid].name += "#1";
						}
						else
						{
							SPhases[phaseid].name += "#2";
							//CEPhases[CHull.f.p[i].phid].name += "#1";
						}
					}
					SPhases.push_back(CEPhases[CHull.f.p[i].phid]);
					SPhases[Phns].vm = CHull.f.p[i].phfrac;
					for (j = 0; j < SPhases[Phns].yn;j++)
					{
						if (CHull.f.p[i].y[j]<eps)
							SPhases[Phns].y[j] = eps;
						else
							SPhases[Phns].y[j] = CHull.f.p[i].y[j];
					}
					Phns++;
				}
			}
			//Phns--;
		}
		else
		{
			//
			Phns = CEPhases.size();
			
			SPhases = CEPhases;

			vector<Phase>::iterator iter;

			for (iter = SPhases.begin(); iter != SPhases.end();)
			{
				//相含量m
				//相含量m的上下界
				if (iter->vmb == 0)
				{
					//
					break;
				}
				else
				{
					iter++;
				}
			}

		}
		// 
		for (i = 0;i < Phns;i++)
		{
			natom = 0;
			for (j = 0; j < SPhases[i].yn;j++)
			{
				if (SPhases[i].yide[j] > 1)  // skip Va
					natom += SPhases[i].ysp[j] * SPhases[i].y[j];
			}
			SPhases[i].vm = SPhases[i].vm / natom;
		}
		//
	}
	
	// 求解亚稳相的驱动力
	void CalcEquilibrium::CalcDrivingForce()
	{

	}
	// 迭代求解平衡相
	int CalcEquilibrium::Iteration()
	{
		//
		int i;
		int maxDFid;
		double maxDF;
		//求解平衡相集的平衡
		MinE.Solver(SPhases, CEConditions, nele, Chp, Chpb);
		/*
		maxDF = -1;
		while (maxDF < 0)
		{
			//change pahse sets
			//求解平衡相集的平衡
			MinE.Solver(1, SPhases, CEConditions, nele, Chp);
			//求解亚稳相的驱动力
			//MinE.Solver(0, MPhases, CEConditions, nele, Chp);
			//判断是否需要更改Phases Sets
			//驱动力判断
			maxDF = MPhases[0].DF;
			maxDFid = 0;
			for (i = 1; i < Phnm; i++)
			{
				if (MPhases[i].DF < maxDF)
				{
					maxDFid = i;
					maxDF = MPhases[i].DF;
				}
			}
			// 相含量判断
			for (i = 1; i < Phns; i++)
			{
				if (MPhases[i].vm < 0)
				{
					maxDF = MPhases[i].DF;
				}
		}
		*/
		return 0;
	}
	
	//Grid point
	void CalcEquilibrium::SetupGP()
	{
		string OutFileName = "output.txt";
		start = clock();
		CalcGPn();  //计算grid points数目
		GeneGP(); //生成Grid Point，成分，并计算能量
		FindConvexFace();  //Find convex face
		finish = clock();
		totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
		if (debug >= 1)
		{
			PrintLine("-");
			cout << "\n\nGlobal minimization " << GPn << " grid points in " << totaltime << "s." << endl;
		}
		//merge grid point

	}
	//
	void CalcEquilibrium::Initialization()
	{
		int i, j;
		double natom;
		for (i = 0;i < Phns;i++ )
		{
			natom = 0;
			for (j = 0; j < SPhases[i].yn;j++)
			{
				if (SPhases[i].yide[j]>1)  // skip vacancy
					natom += SPhases[i].ysp[j] * SPhases[i].y[j];
			}
			SPhases[i].vm = SPhases[i].vm * natom;
		}
	}
	//
	void CalcEquilibrium::ShowEquilibrium()
	{
		int i, j, kk;
		cout << endl << endl;
		PrintLine("=");
		cout << "Conditions:";
		cout << "   T = " << CEConditions.T
			<< "   P = " << CEConditions.P
			<< "   N = " << CEConditions.N << endl;;
		//cout << "   Compositions:";
		for (i = 2; i < CEConditions.nele; i++)
			cout << "   X(" << CEConditions.ele[i] << ") = " << CEConditions.x[i] << endl;
		cout << "\nEquilibrium: " << endl;
		cout << "\nChemical potential: " << endl;
		for (i = 2; i < CEConditions.nele; i++)
			cout << CEConditions.ele[i]<<"   "<< Chp[i] << endl;
		cout << "\nPhase amount: " << endl;
		for (i = 0; i < SPhases.size(); i++)
		{
			cout << SPhases[i].name << endl;
			cout << "Moles: " << SPhases[i].vm << endl;
			for (j = 0; j < nele; j++)
			{
				if(SPhases[i].x[j]>0)
					cout << "X(" << CEConditions.ele[j+2] << ") = " << SPhases[i].x[j] << endl;
			}
			if (debug >= 1)
			{
				cout << endl;
				kk = SPhases[i].yids[0];
				for (j = 0; j < SPhases[i].yn; j++)
				{
					if (SPhases[i].yids[j] != kk) cout << "   :   ";
					cout << "  " << SPhases[i].y[j];
					kk = SPhases[i].yids[j];
				}
			}
			cout << endl;
		}
	}
	void CalcEquilibrium::ShowCondition()
	{

	}
	int CalcEquilibrium::FindPhase(string phasename, vector<Phase> SPhases)
	{
		int i = 0;
		int phaseid = -1;
		for (auto x : SPhases)   // find in DBPhases
		{
			if (x.name == phasename) phaseid = i;
			i++;
		}
		// if phase not be found, return -1;
		return phaseid;
	}
} // end of VCLab