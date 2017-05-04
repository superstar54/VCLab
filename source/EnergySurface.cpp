#include "..\include\EnergySurface.h"

namespace VCLab
{
	using namespace std;

	EnergySurface::~EnergySurface()
	{
		//cout << "Calculate Equilibrium modulue exit." << endl;
	}
	// 求解给定的 phase set
	void EnergySurface::Solver(int GPmod, vector<Phase> SysPhases, Condition CLCondition)
	{
		int i, j, k;
		//清除数据
		ESPhases.clear();
		for (i = 0;i < MDim;i++)
			Chpb[i] = 1;
		//读取数据
		CEConditions = CLCondition;
		ESPhases = SysPhases;
		Phn = ESPhases.size();
		

		GPdx = CLCondition.dx;
		CEConditions.dx = CEConditions.dx / 2.000001;
		nele = CLCondition.nele - 2;
		//
		Bn = 1;
		GPdxn = floor(1 / GPdx) + 1;
		for (i = 0;i < nele - 1;i++)
			Bn = Bn*GPdxn;
		//
		EqBGP		= new Point*[Phn];
		TEqBGP = new Point*[Phn];
		TEqBGPn = new int [Phn];

		BGPn		= new int*[Phn];
		EqBGPPhn	= new int[Phn];
		for (i = 0;i < Phn;i++)
		{
			EqBGP[i]	= new Point[Bn];
			TEqBGP[i] = new Point[Bn];

			BGPn[i]		= new int[Bn];
			for (j = 0;j < Bn;j++)
				BGPn[i][j] = 0;

			EqBGPPhn[i] = 0;
		}

		//
		GPmode = GPmod;
		if (GPmode == 1)
		{
			//Grid Point
			GGP.GeneGPy(ESPhases, CEConditions);
			//建立Grid point
			GGP.GeneGPGF(ESPhases, CEConditions);
		}
		// loop phases
		for (i = 0;i < Phn;i++)
		{
			// block GP
			BlockGP(i, ESPhases[i]);

			for (j = 0;j < Bn;j++)
			{
				if (BGPn[i][j] > 0)
				{
					//
					if (debug >= 1)
					{
						ShowPoint(ESPhases[i].name, nele, EqBGP[i][j]);
						cout << endl;
					}
					cout << " j = " << j << "  ";
					
					for (k = 0;k < ESPhases[i].yn; k++)
						ESPhases[i].y[k] = EqBGP[i][j].y[k];
					for (k = 0;k < nele;k++)
					{
						CEConditions.x[k + 2] = EqBGP[i][j].x[k];
						cout << CEConditions.x[k + 2] << "   ";
					}
					cout << endl;
					ESPhases[i].vm = 1;
					MinEPhases.clear();
					MinEPhases.push_back(ESPhases[i]);
					
					MinE.Solver(MinEPhases, CEConditions, CEConditions.nele, Chp, Chpb);
					//
					for (k = 0;k < ESPhases[i].yn; k++)
						EqBGP[i][j].y[k] = MinEPhases[0].y[k];
					EqBGP[i][j].x[nele] = MinEPhases[0].GF;
					//
					//delete[] BGP[j];
				}
			}
			// delete
			//delete[] BGP;
			//delete[] GPidB;
			//delete[] GGP.GP;
		}
		
		//ShowEquilibrium();
		WriteResults();
		//WriteFaces();
	}
	void EnergySurface::BlockGP(int Phid, Phase Phases)
	{
		int i, j, k;
		int blockid;
		int nc;
		
		//
		BGP = new Point*[Bn];
		GPidB = new int[Phases.GPn];
		
		for (i = 0;i < Phases.GPn;i++)
		{
			blockid = 0;
			for (j = 0;j < nele - 1;j++)
			{
				blockid += floor(Phases.GP[i].x[j]/GPdx)*pow(GPdxn, nele - j - 2);
			}

			GPidB[i] = blockid;
			BGPn[Phid][blockid]++;
			if(debug==1)
				cout << "blockid  " << blockid << "   BGPn  " << BGPn[Phid][blockid] << endl;
		}

		for (i = 0;i < Bn;i++)
		{
			if (BGPn[Phid][i] > 0)
			{
				BGP[i] = new Point[BGPn[Phid][i]];
				EqBGP[Phid][i].x[nele] = 1e15;
				BGPn[Phid][i] = 0;
				EqBGPPhn[Phid]++;
			}
		}

		for (i = 0;i < Phases.GPn;i++)
		{
			//cout << "i = " << i << endl;
			BGP[GPidB[i]][BGPn[Phid][GPidB[i]]] = Phases.GP[i];
			if (debug == 1)
			{
				ShowPoint(Phases.name, nele, BGP[GPidB[i]][BGPn[Phid][GPidB[i]]]);
				cout << endl;
			}
			
			if (Phases.GP[i].x[nele] < EqBGP[Phid][GPidB[i]].x[nele])
				EqBGP[Phid][GPidB[i]] = Phases.GP[i];
			BGPn[Phid][GPidB[i]]++;
			
		}
	}

	//
	void EnergySurface::SetupConditions(Condition CLCondition)
	{
		CEConditions = CLCondition;
	}
	
	void EnergySurface::WriteResults()
	{
		int i, j, k;

		ofstream ESOut("Energy Surface.txt");
		ESOut << "T = " << CEConditions.T << " K"
			<< ", P = " << CEConditions.P << " Pa"
			<< ", N = " << CEConditions.N << endl;
		ESOut << nele << endl;
		ESOut << Phn << endl;

		for (i = 0;i < Phn;i++)
		{
			TEqBGPn[i] = 0;

			ESOut << ESPhases[i].name << endl;
			ESOut << EqBGPPhn[i] << endl;
			for (j = 0;j < Bn;j++)
			{
				if (BGPn[i][j] > 0)
				{
					TEqBGP[i][TEqBGPn[i]] = EqBGP[i][j];
					TEqBGPn[i]++;
					for (k = 0;k < nele + 1;k++)
					{
						ESOut << EqBGP[i][j].x[k] << "  ";
					}
					ESOut << endl;
				}
			}
		}

		ESOut.close();
	}

	void EnergySurface::WriteFaces()
	{
		int i, j, k, ki;
		for (i = 0;i < Phn;i++)
		{
			FindNetFace(TEqBGPn[i], TEqBGP[i]);

			ofstream ESOut(ESPhases[i].name + ".txt");
			ESOut << "T = " << CEConditions.T << " K"
				<< ", P = " << CEConditions.P << " Pa"
				<< ", N = " << CEConditions.N << endl;
			ESOut << nele << endl;
			ESOut << TGPfn << endl;
			for (j = 0;j < TGPfn;j++)
			{
				for (k = 0;k < nele + 1;k++)
				{
					for (ki = 0;ki < nele;ki++)
						ESOut << TGPf[j].p[k].x[ki] << "  ";
					ESOut << endl;
				}
			}
			ESOut.close();
		}
	}


	//Find Net Face ============================================================
	void EnergySurface::FindNetFace(int in_GPn, Point *GP)
	{
		GPn = in_GPn;			//
		int fn = GPn*nele;		//
		GPf = new Face[fn];		//
		TGPf = new Face[fn];	//

		FindGFace(GPn, GP, f);
		GPfn = 0;  //
		GPf[GPfn] = f;
		GPf[GPfn].b = 1;
		GPf[GPfn].id = GPfn;
		if (debug == 2)
		{
			ShowFace("Initial face", nele, nele + 1, GPf[0]);
		}
		GPfn++;

		RecurNetFace(GPn, GP, GPf[0]);
		GetNetFace();
	}

	void EnergySurface::FindGFace(int GPn, Point *GP, Face &f)
	{
		int i, j;
		int id[MDim];

		//initial face is unit matrix with large energy
		for (i = 0; i < nele; i++)
		{
			id[i] = -1;
			for (j = 0; j < nele; j++)
				f.p[i].x[j] = eps;
			f.p[i].x[i] = 1 - 10 * eps;
			f.p[i].x[nele] = 1e15;
			f.p[i].fid = i;
			f.p[i].fb = 0;
		}
		// the i point in face have the largest i compositions, same com, lower energy
		for (i = 0; i < GPn; i++)
		{
			for (j = 0; j < nele; j++)
			{
				if (GP[i].x[j] - f.p[j].x[j] > -2.0 * eps)
				{
					id[j] = i;
					//f.p[j] = GP[i];
				}
			}
		}

		for (i = 0; i < nele; i++)
		{
			if (id[i] > -1)
			{
				f.p[i] = GP[id[i]];
				GP[id[i]].fb = 0;
			}
		}
	}

	void EnergySurface::RecurNetFace(int GPn, Point *GP, Face f)
	{
		int i, j;
		int pid;
		int subedge;
		int parray[MDim];
		//divide GP,get maximum nele new region
		int newfn;					// number new face
		int newfid[MDim];			// id of new face
		int newGPn[MDim];			// number of GP that belong to this new face
		Point **newGP;
		newGP = new Point*[nele];	//新凸面分得的点
		for (i = 0; i < nele; i++)
		{
			newGPn[i] = 0;
			newGP[i] = new Point[GPn];
		}
		pid = FindNewP(GPn, GP, f); //遍历所有的点，寻找下方，距离最远的点。
										  //找到新的凸点，更新全局的凸点和凸面，即删除当前凸面，新增加nele个凸面
		if (pid > -1)
		{
			if (debug == 2)
			{
				ShowPoint("new point: ", nele, GP[pid]);
			}
			GPf[f.id].b = 0;  //original vace is false，0
			GP[pid].fb = 0;
			//判断该凸点的投影是否属于该面的某边个投影
			subedge = 0;
			//递归寻找属于哪个维度的边
			inside_edge(subedge, parray, GP[pid], f);
			//新凸面是原图的nele-subedge-1个点，加上原图的subedge个点，加上新的凸点构成。
			int pture;
			newfn = 0;
			for (i = 0; i < nele; i++)
			{
				pture = 1;
				for (j = 0; j < subedge; j++)
				{
					if (i == parray[j]) //新的凸面必须包含iedge点
					{
						pture = 0;
						break;
					}
				}
				if (pture == 0)
					continue;
				GPf[GPfn] = f;
				//将新凸点加入新凸面，即代替原来的第i个点
				GPf[GPfn].p[i] = GP[pid];
				newfid[newfn] = GPfn;  //新凸面的id号
				GPf[GPfn].id = GPfn;
				//对面中的点进行编号
				for (j = 0; j < nele; j++)
					GPf[GPfn].p[j].fid = j;
				if (debug == 2)
				{
					ShowFace("增加一个新凸面", nele, nele + 1, GPf[GPfn]);
				}
				newfn++;  //新增加的面的个数
				GPfn++;  //增加一个凸面
			}

			//将点和新凸面投影到成分面上，判断该点是否属于该凸面，进行划分
			for (i = 0; i < GPn; i++)
			{
				if (GP[i].fb == 0)
					continue;
				//循环新凸面进行判断，要循环所有的新面，因为一个点可以属于多个凸面
				if (debug == 2)
				{
					ShowPoint("点", nele, GP[i]);
				}
				for (j = 0; j < newfn; j++)
				{
					//如果在这个区域内部，这个区域的GP增加一个，跳出循环
					if (debug == 2)
					{
						ShowFace("面的投影", nele, nele, GPf[newfid[j]]);
					}
					if (inside(GP[i], GPf[newfid[j]])) //该点在本次凸面内部
					{
						newGP[j][newGPn[j]] = GP[i];
						newGPn[j]++;
					}
				}
			}
			//新的递归查找
			for (i = 0; i < newfn; i++)
				RecurNetFace(newGPn[i], newGP[i], GPf[newfid[i]]);
		}
	}

	int EnergySurface::FindNewP(int GPn, Point *GP, Face f)
	{
		int i, j, id;
		double tempg, tempd, maxd;
		id = -1;
		for (i = 0; i < GPn; i++)
		{
			if (GP[i].fb == 0)
				continue;
			//showPoint("Point: ", Dim, GP[i]);
			id = i;
			break;
		}
		return id;
	}

	//提取真凸面函数
	void EnergySurface::GetNetFace()
	{
		int i, j, k;
		int tb;

		TGPfn = 0;
		for (i = 0; i < GPfn; i++)
		{
			tb = 1;
			if (debug > 1)
			{
				ShowFace("Face: ", nele, nele + 1, GPf[i]);
			}
			if (GPf[i].b == 1)
			{
				for (j = 0;j < nele;j++)
				{
					if (GPf[i].p[j].x[nele]>1e10)
					{
						tb = 0; 
						break;
					}
				}
				if (tb == 0)
					continue;
				TGPf[TGPfn] = GPf[i];
				TGPfn++;
			}
		}
			
	}


	

	
	void EnergySurface::ShowPoint(string str, int n, Point p)
	{
		int i;
		cout << str.c_str() << endl;
		for (i = 0; i < n; i++)
			cout << p.x[i] << "  ";
	}

	void EnergySurface::ShowFace(string str, int n, int m, Face f)
	{
		int i, j;
		cout << str.c_str() << endl;
		for (i = 0; i < n; i++)
		{
			cout << "the " << i << " Point： ";
			for (j = 0; j < m; j++)
				cout << f.p[i].x[j] << "  ";
			cout << endl;
		}
	}

	void EnergySurface::deletep()
	{
		//释放空间
		//delete[] GPf;   //递归过程中记录
		//delete[] TGPf;  //根据GPfb值，提取最终的凸面
	}
	//N维空间中，判断点k是否在N维面的下方 (在小于N维面的情况，排列组合，在inside_ndmd函数中有体系)
	//mode==0，判断方向，mode=1，在下方时，计算距离
	double EnergySurface::distance(int mode, int n, Point p, Face f)
	{
		int i, j, k;
		double det;
		//构建矩阵，并计算矩阵行列式
		/// 三维空间示例
		//| x11 x12 x13 1 |
		//| x21 x22 x23 1 |
		//| x31 x32 x33 1 |
		//| xk1 xk2 xk3 1 |
		//
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
			{
				mat[i][j] = f.p[i].x[j];
			}
			mat[n][i] = p.x[i];
			mat[i][n] = 1;
		}
		mat[n][n] = 1;

		//for (i = 0; i < n + 1; i++)
		//{
		//for (j = 0; j < n + 1; j++)
		//cout << mat[i][j] << "  ";
		//cout << endl;
		//}
		//
		//
		det = calcdet(n + 1, mat);
		if (det < -eps)  //在平面上方
		{
			return -1;
		}
		else //在平面或者平面上下方
		{

			if (mode == 0)
			{
				return det;  //只需要计算方向
			}
			//计算距离
			double det2, det3;
			det2 = 0;
			// 构建矩阵，并计算行列式，进行累加
			for (i = 0; i < n; i++)
			{
				for (j = 0; j < n; j++)
				{
					for (k = 0; k < i; k++)
						mat[j][k] = f.p[j].x[k];
					for (k = i + 1; k < n; k++)
						mat[j][k - 1] = f.p[j].x[k];
					mat[j][n - 1] = 1;
				}
				det3 = calcdet(n, mat);
				det2 += det3*det3;
			}
			double dis = fabs(det) / sqrt(det2);
			return dis;
		}
	}

	bool EnergySurface::inside(Point p, Face f)
	{
		int i, j, k, in;
		double dx, dfx;
		double eF[MDim3];
		double eFJ[MDim3][MDim3];
		// building matrix
		for (i = 0; i < nele; i++)
		{
			for (j = 0; j < nele; j++)
				eFJ[i][j] = f.p[j].x[i];
			eF[i] = p.x[i];
		}
		if (gauss_elimination(nele, eFJ, eF))
		{
			// all phase farctions should in [0,1]
			for (i = 0; i < nele; i++)
			{
				if (eF[i] < -eps || eF[i] > (1.0 + eps))
					return false;
			}
		}
		else
		{
			//cout << "\n The matrix is singular, calculations failed!!!" << endl;
			return false;
		}
		return true;
	}

	//判断点是否在面投影的边上
	void EnergySurface::inside_edge(int &subedge, int parray[MDim], Point p, Face f)
	{
		int i, j, k, in;
		double dx, dfx;
		double eF[MDim3];
		double eFJ[MDim3][MDim3];
		// building matrix
		for (i = 0; i < nele; i++)
		{
			for (j = 0; j < nele; j++)
				eFJ[i][j] = f.p[j].x[i];
			eF[i] = p.x[i];
		}
		//
		subedge = 0;
		if (gauss_elimination(nele, eFJ, eF))
		{
			// all phase farctions should in [0,1]
			for (i = 0; i < nele; i++)
			{
				if (fabs(eF[i]) < eps)
				{
					parray[subedge] = i;
					subedge++;
				}
			}
		}
	}

	//求矩阵行列式,将一个矩阵经过初等行变换化为上三角矩阵。
	double EnergySurface::calcdet(int n, double array[MDim3][MDim3])
	{
		int i, j, k, u;
		int iter = 0;  //记录行变换的次数（交换）
		double det1 = 1;
		double temp, yin;

		for (i = 0; i<n; i++)
		{
			if (fabs(array[i][i]) <2 * eps)

				for (j = i; j<n; j++)
				{
					if (fabs(array[j][i]) >eps)
					{
						//交换两行
						for (k = 0; k<n; k++)
						{
							temp = array[i][k];
							array[i][k] = array[j][k];
							array[j][k] = temp;
						}
						iter++;
					}
				}

			for (k = i + 1; k<n; k++)
			{
				if (fabs(array[i][i]) < 2 * eps)  // 如果对角元素为0，表示行列式的该列都为0，故det为0
					return 0;
				yin = -1 * array[k][i] / array[i][i];

				for (u = 0; u<n; u++)
				{
					array[k][u] = array[k][u] + array[i][u] * yin;
				}
			}
		}
		for (i = 0; i<n; i++)  //求对角线的积 即 行列式的值
			det1 = det1 * array[i][i];
		//行变换偶数次符号不变
		if (iter % 2 == 1)
			det1 = -det1;
		return (det1);
	}


} // end of VCLab