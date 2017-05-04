<<<<<<< HEAD

#include "..\include\FindConvexHull.h"
namespace VCLab
{
	//================   Find Convex Face    =====================================
	// GP, 0:dim-2 is x, dim -1 is energy, dim is x
	void FindCH::FindConvexFace(int in_Dim, int in_GPn, Point *GP, Point  PA)
	{
		int i;
		Dim = in_Dim;			//
		GPn = in_GPn;			//
		int fn = GPn*Dim;		//
		
		if (debug >= 2)
		{
			cout << "Input data: \n";
			for (i = 0; i < GPn; i++)
				ShowPoint(" ", Dim + 1, GP[i]);
		}
		

		//Found the global face
		FindGFace(Dim, GPn, GP, f);
		if (debug >= 1)
		{
			ShowPoint("alloy composition：", Dim, PA);
			ShowFace("initial global face：", Dim, Dim + 1, f);
		}
		//interation
		IterConvexFace(GPn, GP, f, PA); 

	}
	void FindCH::FindGFace(int Dim, int GPn, Point *GP, Face &f)
	{
		int i, j;
		//initial face is unit matrix with large energy
		for (i = 0; i < Dim; i++)
		{
			for (j = 0; j < Dim; j++)
				f.p[i].x[j] = eps;
			f.p[i].x[i] = 1 - 10*eps;  
			f.p[i].x[Dim] = 1e8;
			f.p[i].fid = i; 
			f.p[i].fb = 1;
		}
		// the i point in face have the largest i compositions, same com, lower energy
		for (i = 0; i < GPn; i++)
		{
			for (j = 0; j < Dim; j++)
			{
				if (GP[i].x[j] - f.p[j].x[j] > -2.0 * eps && GP[i].x[Dim] < f.p[j].x[Dim]) 
					f.p[j] = GP[i]; 
			}
		}
	}

	void FindCH::IterConvexFace(int GPn, Point *GP, Face &f, Point PA)
	{
		int i, j, maxi, newmaxi;
		double tempg, tempgmin;
		double tempgmin_old;
		int newfb = 0, pmin;
		Point tempP;

		CalcChemicalpotential(f, Chp);
		//Find the lowest energy point, based on current energy surface
		ncountp = 0;
		maxi = FindMaxD(GPn, GP, f, Chp);
		tempgmin_old = 1e30;
		while (maxi > -1)
		{
			if (ncountp > GPn)
			{
				cout << "Too much iteration in find convex point!!!";
				exit(1);
			}
			if (debug >= 1)
				ShowPoint("New point: ", Dim, GP[maxi]);
			newmaxi = -1;
			newfb = 0;
			tempgmin = 1e30;
			//replace one point in old face using new point, jugde alloy in new face or not
			for (i = 0; i < Dim; i++)
			{
				tempP = f.p[i];
				f.p[i] = GP[maxi];
				// inside new face or not
				if (inside(PA, f))
				{
					newfb = 1;
					// new energy surface
					CalcChemicalpotential(f, Chp);
					//
					tempg = 0;
					for (j = 0; j < Dim; j++)
					{
						tempg += PA.x[j] * Chp[j + 2];
					}
					if (tempg < tempgmin)
					{
						tempgmin = tempg;
						pmin = i;
					}
				}
				f.p[i] = tempP;
			}
			if (newfb == 1)
			{
				
				f.p[pmin] = GP[maxi];
				CalcChemicalpotential(f, Chp);
				if (debug >= 1)
					ShowFace("Find a new face: ", Dim, Dim + 1, f);
				if (debug >= 1)
					cout << "Current gibbs energy: " << tempgmin << endl;
				newmaxi = FindMaxD(GPn, GP, f, Chp);
				if ((tempgmin - tempgmin_old)>-1e-9)
				{
					newmaxi = -1;
					break;
				}
			}
			tempgmin_old = tempgmin;
			maxi = newmaxi;
		}
	}

	void FindCH::CalcChemicalpotential(Face f, double(&Chp)[MDim3])
	{
		int i, j;
		for (i = 0; i < Dim; i++)
		{
			for (j = 0; j < Dim; j++)
				eF[i][j] = f.p[i].x[j];
			eFJ[i] = f.p[i].x[Dim];
		}
		
		if (gauss_elimination(Dim, eF, eFJ))
		{
			for (i = 0; i < Dim; i++)
				Chp[i + 2] = eFJ[i];
		}
		else
		{
			cout << "\n The matrix is singular, calculations failed!!!" << endl;
			exit(1);
		}
		

	}

	int FindCH::FindMaxD(int GPn, Point *GP, Face f, double Chp[MDim3])
	{
		int i, j, maxi;
		double tempg, tempd, maxd;
		maxi = -1;
		maxd = eps6;
		if (ncountp < 1)
		{
			for (i = 0; i < GPn; i++)
			{
				if (GP[i].fb == 0)
					continue;
				//showPoint("Point: ", Dim, GP[i]);
				tempg = 0;
				for (j = 0; j < Dim; j++)
				{
					tempg += GP[i].x[j] * Chp[j + 2];
				}
				tempd = tempg - GP[i].x[Dim];
				if (tempd>maxd) 
				{
					maxd = tempd;
					maxi = i;
				}
				else if(tempd<-eps3)
				{
					GP[i].fb = 0;
				}
			}
		}
		else
		{
			for (i = 0; i < GPn; i++)
			{
				if (GP[i].fb == 0)
					continue;
				//showPoint("Point: ", Dim, GP[i]);
				tempg = 0;
				for (j = 0; j < Dim; j++)
				{
					tempg += GP[i].x[j] * Chp[j + 2];
				}
				tempd = tempg - GP[i].x[Dim];
				if (tempd>maxd)
				{
					maxd = tempd;
					maxi = i;
				}
			}
		}
		ncountp++;
		return maxi;
	}

	void FindCH::ShowConvexFace(string OutFileName)
	{
		int i, j;
		cout << "Face found: " << endl;
		ofstream Out(OutFileName.c_str());
		Out << Dim << endl;
		for (i = 0; i < Dim; i++)
		{
			for (j = 0; j < Dim + 1; j++)
			{
				cout << f.p[i].x[j] << "  ";
				Out << f.p[i].x[j] << "  ";
			}
			cout << endl;
		}
		Out.close();
	}
	void FindCH::ShowFace(string str, int n, int m, Face f)
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

	void FindCH::ShowPoint(string str, int n, Point p)
	{
		int i;
		cout << str.c_str() << endl;
		for (i = 0; i < n; i++)
			cout << p.x[i] << "  ";
	}

	

	bool FindCH::inside(Point p, Face f)
	{
		int i, j, k, in;
		double dx, dfx;
		double eF[MDim3];
		double eFJ[MDim3][MDim3];
		// building matrix
		for (i = 0; i < Dim; i++)
		{
			for (j = 0; j < Dim; j++)
				eFJ[i][j] = f.p[j].x[i];
			eF[i] = p.x[i];
		}
		if (gauss_elimination(Dim, eFJ, eF))
		{
			// all phase farctions should in [0,1]
			for (i = 0; i < Dim; i++)
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

	
} //end of VCLab
=======

#include "..\include\FindConvexHull.h"
namespace VCLab
{
	//================   Find Convex Face    =====================================
	// GP, 0:dim-2 is x, dim -1 is energy, dim is x
	void FindCH::FindConvexFace(int in_Dim, int in_GPn, Point *GP, Point  PA)
	{
		int i;
		Dim = in_Dim;			//
		GPn = in_GPn;			//
		int fn = GPn*Dim;		//
		/*
		cout << "Input data: \n";
		for (i = 0; i < GPn; i++)
			ShowPoint(" ", Dim, GP[i]);
		*/

		//Found the global face
		FindGFace(Dim, GPn, GP, f);
		if (debug >= 1)
		{
			ShowPoint("alloy composition：", Dim, PA);
			ShowFace("initial global face：", Dim, Dim + 1, f);
		}
		//interation
		IterConvexFace(GPn, GP, f, PA); 

		delete GP;
	}
	void FindCH::FindGFace(int Dim, int GPn, Point *GP, Face &f)
	{
		int i, j;
		//initial face is unit matrix with large energy
		for (i = 0; i < Dim; i++)
		{
			for (j = 0; j < Dim + 1; j++)
				f.p[i].x[j] = eps;
			f.p[i].x[i] = 1 - 10*eps;  
			f.p[i].x[Dim - 1] = 1e8;
			f.p[i].fid = i; 
			f.p[i].fb = 1;
		}
		f.p[Dim - 1].x[Dim] = 1 - eps;
		// the i point in face have the largest i compositions, same com, lower energy
		for (i = 0; i < GPn; i++)
		{
			for (j = 0; j < Dim - 1; j++)
			{
				if (GP[i].x[j] - f.p[j].x[j] > -2.0 * eps && GP[i].x[Dim - 1] < f.p[j].x[Dim - 1]) 
					f.p[j] = GP[i]; 
			}
			// Dim is x
			if (GP[i].x[Dim] - f.p[Dim - 1].x[Dim] > -2.0 * eps && GP[i].x[Dim - 1] < f.p[Dim - 1].x[Dim - 1])
				f.p[Dim - 1] = GP[i];
		}
	}

	void FindCH::IterConvexFace(int GPn, Point *GP, Face &f, Point PA)
	{
		int i, j, maxi, newmaxi;
		double tempg, tempgmin;
		double tempgmin_old;
		int newfb = 0, pmin;
		Point tempP;

		CalcChemicalpotential(f, Chp);
		//Find the lowest energy point, based on current energy surface
		ncountp = 0;
		maxi = FindMaxD(GPn, GP, f, Chp);
		tempgmin_old = 1e30;
		while (maxi > -1)
		{
			if (ncountp > GPn)
			{
				cout << "Too much iteration in find convex point!!!";
				exit(1);
			}
			if (debug >= 1)
				ShowPoint("New point: ", Dim, GP[maxi]);
			newmaxi = -1;
			newfb = 0;
			tempgmin = 1e30;
			//replace one point in old face using new point, jugde alloy in new face or not
			for (i = 0; i < Dim; i++)
			{
				tempP = f.p[i];
				f.p[i] = GP[maxi];
				// inside new face or not
				if (inside(PA, f))
				{
					newfb = 1;
					// new energy surface
					CalcChemicalpotential(f, Chp);
					//
					tempg = 0;
					for (j = 0; j < Dim; j++)
					{
						tempg += PA.x[j] * Chp[j + 2];
					}
					if (tempg < tempgmin)
					{
						tempgmin = tempg;
						pmin = i;
					}
				}
				f.p[i] = tempP;
			}
			if (newfb == 1)
			{
				
				f.p[pmin] = GP[maxi];
				CalcChemicalpotential(f, Chp);
				if (debug >= 1)
					ShowFace("Find a new face: ", Dim, Dim + 1, f);
				if (debug >= 1)
					cout << "Current gibbs energy: " << tempgmin << endl;
				newmaxi = FindMaxD(GPn, GP, f, Chp);
				if ((tempgmin - tempgmin_old)>-1e-9)
				{
					newmaxi = -1;
					break;
				}
			}
			tempgmin_old = tempgmin;
			maxi = newmaxi;
		}
	}

	void FindCH::CalcChemicalpotential(Face f, double(&Chp)[MDim3])
	{
		int i, j;
		double temp[MDim3];
		double omat[MDim3][MDim3];
		double omat1[MDim3];
		for (i = 0; i < Dim; i++)
		{
			for (j = 0; j < Dim - 1; j++)
				mat[i][j] = f.p[i].x[j];
			mat[i][Dim - 1] = f.p[i].x[Dim];
			mat1[i] = f.p[i].x[Dim - 1];
		}
		
		if (gauss_elimination(Dim, mat, mat1))
		{
			for (i = 0; i < Dim; i++)
				Chp[i + 2] = mat1[i];
		}
		else
		{
			cout << "\n The matrix is singular, calculations failed!!!" << endl;
			exit(1);
		}
		

	}

	int FindCH::FindMaxD(int GPn, Point *GP, Face f, double Chp[MDim3])
	{
		int i, j, maxi;
		double tempg, tempd, maxd;
		maxi = -1;
		maxd = eps6;
		if (ncountp < 1)
		{
			for (i = 0; i < GPn; i++)
			{
				if (GP[i].fb == 0)
					continue;
				//showPoint("Point: ", Dim, GP[i]);
				tempg = 0;
				for (j = 0; j < Dim - 1; j++)
				{
					tempg += GP[i].x[j] * Chp[j + 2];
				}
				tempg += GP[i].x[Dim] * Chp[Dim + 1];
				tempd = tempg - GP[i].x[Dim - 1];
				if (tempd>maxd) 
				{
					maxd = tempd;
					maxi = i;
				}
				else if(tempd<-eps3)
				{
					GP[i].fb = 0;
				}
			}
		}
		else
		{
			for (i = 0; i < GPn; i++)
			{
				if (GP[i].fb == 0)
					continue;
				//showPoint("Point: ", Dim, GP[i]);
				tempg = 0;
				for (j = 0; j < Dim - 1; j++)
				{
					tempg += GP[i].x[j] * Chp[j + 2];
				}
				tempg += GP[i].x[Dim] * Chp[Dim + 1];
				tempd = tempg - GP[i].x[Dim - 1];
				if (tempd>maxd)
				{
					maxd = tempd;
					maxi = i;
				}
			}
		}
		ncountp++;
		return maxi;
	}

	void FindCH::ShowConvexFace(string OutFileName)
	{
		int i, j;
		cout << "Face found: " << endl;
		ofstream Out(OutFileName.c_str());
		Out << Dim << endl;
		for (i = 0; i < Dim; i++)
		{
			for (j = 0; j < Dim + 1; j++)
			{
				cout << f.p[i].x[j] << "  ";
				Out << f.p[i].x[j] << "  ";
			}
			cout << endl;
		}
		Out.close();
	}
	void FindCH::ShowFace(string str, int n, int m, Face f)
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

	void FindCH::ShowPoint(string str, int n, Point p)
	{
		int i;
		cout << str.c_str() << endl;
		for (i = 0; i < n; i++)
			cout << p.x[i] << "  ";
	}

	

	bool FindCH::inside(Point p, Face f)
	{
		int i, j, k, in;
		double dx, dfx;
		double eF[MDim3];
		double eFJ[MDim3][MDim3];
		// building matrix
		for (i = 0; i < Dim - 1; i++)
		{
			for (j = 0; j < Dim; j++)
				eFJ[i][j] = f.p[j].x[i];
			eFJ[Dim - 1][i] = f.p[i].x[Dim];
			eF[i] = p.x[i];
		}
		eFJ[Dim - 1][Dim - 1] = f.p[Dim - 1].x[Dim];
		eF[Dim - 1] = p.x[Dim - 1];
		if (gauss_elimination(Dim, eFJ, eF))
		{
			// all phase farctions should in [0,1]
			for (i = 0; i < Dim; i++)
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

	//Find Convex Hull============================================================
	//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> not implement right now
	/*void FindCH::FindConvexHull(int in_Dim, int in_GPn, Point *GP)
	{
		Dim = in_Dim;			//
		GPn = in_GPn;			//
		int fn = GPn*Dim;		//
		GPf = new Face[fn];		//
		TGPf = new Face[fn];	//

		FindGFace(Dim, GPn, GP, f);
		GPfn = 0;  //
		GPf[GPfn] = f; 
		GPf[GPfn].b = 1;
		GPf[GPfn].id = GPfn;
		ShowFace("Initial face", Dim, Dim + 1, GPf[0]);
		GPfn++;

		RecurConvexHull(GPn, GP, GPf[0]); 
		GetConvexHull(); 
	}
	
	

	void FindCH::deletep()
	{
		//释放空间
		delete[] GPf;   //递归过程中记录
		delete[] TGPf;  //根据GPfb值，提取最终的凸面
	}
	//N维空间中，判断点k是否在N维面的下方 (在小于N维面的情况，排列组合，在inside_ndmd函数中有体系)
	//mode==0，判断方向，mode=1，在下方时，计算距离
	double FindCH::distance(int mode, int n, Point p, Face f)
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

	//判断点是否在面投影的边上
	int FindCH::inside_edge(int n, int &subedge, int parray[MDim], Point p, Face f)
	{
	int i, j, k, in;
	double dx;
	in = -1; //不在
	//选择fx上的一个点i
	for (i = 0; i < Dim; i++)
	{
	//构建新的面，不包含点i
	for (j = 0; j < i; j++)
	for (k = 0; k < Dim - 1; k++)  //因为是投影到成分面上，所以Dim-1能量值不要，另外一个Dim的成分值也不要
	subf.p[j].x[k] = f.p[j].x[k];
	for (j = i + 1; j < Dim; j++)
	for (k = 0; k < Dim - 1; k++)
	subf.p[j - 1].x[k] = f.p[j].x[k];
	//判断点x和fx的一个点与新的面的方向性
	dx = distance(0, Dim - 1, p, subf);
	if (fabs(dx) < eps)  //所有的都同向则在内部。
	return in = i; //在，返回边的序号
	}
	return in;
	}

	//求矩阵行列式,将一个矩阵经过初等行变换化为上三角矩阵。
	double FindCH::calcdet(int n, double array[MDim3][MDim3])
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
	*/
	
	
} //end of VCLab
>>>>>>> b34a88e731ed6049f6c364b37b469ae02b913009
