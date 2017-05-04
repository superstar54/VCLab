<<<<<<< HEAD
#include "..\include\MinEnergy.h"

namespace VCLab
{
	using namespace std;
	/*---------------------------------------------------------------------*/
	/*                       计算Gibbs自由能                               */
	/*---------------------------------------------------------------------*/
	/*
	void MinEnergy::CalcGE(Phase & Phases, double T)
	double MinEnergy::CalcGF(Phase & Phases, double T)
	double MinEnergy::CalcGFdx(Phase & Phases, double T, int idy)
	double MinEnergy::CalcGFdx2(Phase & Phases, double T, int idy1, int idy2)
	double MinEnergy::CalcPara(int Tnseg, Parameter para, double T)
	double MinEnergy::CalcGtao(double tao, double p, int order)
	*/

	// Calc the Gibbs energy, the first , second derivative of site fraction and temperature
	// mode: 0, 1, 2; 10, 11, 12;
	void MinEnergy::CalcGE(int mode, Phase & Phases, double T)
	{
		//
		int i, j;
		double R = 8.31451;
		double xpro_all, xpro, xdiff;
		double dxpro[10], dxpro2[10];
		double d2xpro[10][10], d2xpro2[10][10], d2xpro3[10][10], d2xpro4[10][10];

		double TC = 0;
		double BMAGN = 0;
		double tao;
		//
		double yfracs[MDim];
		double paravalue[3];

		for (i = 0; i < Phases.yn; i++)
		{
			if(Phases.y[i]<eps)
				yfracs[i] = eps;
			else
				yfracs[i] = Phases.y[i];
		}
		//
		Phases.GF = 0; Phases.SF = 0;
		Phases.dGFt = 0; Phases.d2GFt = 0;
		for (i = 0;i < Phases.vyn;i++)
		{
			Phases.dGFx[i] = 0;
			Phases.d2GFxt[i] = 0;
			for (j = 0;j < Phases.vyn;j++)
			{
				Phases.d2GFx[i][j] = 0;
			}
		}

		// 循环相参数
		for (auto para : Phases.Parameters)
		{
			// 判断温度区间
			int Tnseg = 0;
			for (i = 0; i < para.T.size() - 1;i++)
			{
				if ((para.T[i] - 0.02)< T && T <(para.T[i + 1] + 0.02))   // T - 0.02，防止298.14和298.15这种误差
					break;
				Tnseg++;
			}
			//
			for (i = 0;i < Phases.vyn;i++)
			{
				dxpro[i] = 0;
				d2xpro[i][i] = 0;
				d2xpro2[i][i] = 0;
				d2xpro3[i][i] = 0;
				d2xpro4[i][i] = 0;

				for (j = i + 1;j < Phases.vyn;j++)
				{
					d2xpro[i][j] = 0; d2xpro[j][i] = 0;
					d2xpro2[i][j] = 0; d2xpro2[j][i] = 0;
					d2xpro3[i][j] = 0; d2xpro3[j][i] = 0;
					d2xpro4[i][j] = 0; d2xpro4[j][i] = 0;
				}
			}
			// 计算相参数前面的系数
			xpro_all = 1;
			for (i = 0; i < para.yn; i++)
			{
				xpro_all = xpro_all*yfracs[para.yidc[i]];
			}
			// 根据参数类型，进行修正
			switch (para.kind)
			{
			case 1:
				//g
				xpro = xpro_all;
				// first derivative
				for (i = 0;i < para.vyn;i++)
				{
					dxpro[para.vyidv[i]] = xpro_all / yfracs[para.vyidc[i]];
					// second derivate
					// for i==j, in case 1, this part is zero
					for (j = i + 1;j < para.vyn;j++)
					{
						d2xpro[para.vyidv[i]][para.vyidv[j]] = dxpro[para.vyidv[i]] / yfracs[para.vyidc[j]];
						d2xpro[para.vyidv[j]][para.vyidv[i]] = d2xpro[para.vyidv[i]][para.vyidv[j]];
					}
				}
				break;
			case 2:  // 二元相互作用参数
					 //亚点阵内部二元相互作用，求(x1-x2)
				xpro = xpro_all;
				xdiff = yfracs[para.idsub2[0]] - yfracs[para.idsub2[1]];
				xpro = xpro*pow(xdiff, para.order);
				//
				//first derivative
				//亚点阵内部二元相互作用，求(x1-x2)
				for (i = 0;i < para.vyn;i++)
				{
					dxpro[para.vyidv[i]] = xpro_all/yfracs[para.vyidc[i]]*pow(xdiff, para.order);
				}
				//第二部分, idy在相互作用的亚点阵内部，Xi*Xj*Xk(Xi-Xj)^n，idy是Xi或者Xj
				dxpro2[para.vidsub2[0]] = 0;
				if (para.order > 0)  // n=0,求导之后为0，不考虑。
				{
					dxpro2[para.vidsub2[0]] = xpro_all*para.order*pow(xdiff, para.order - 1);
				}
				dxpro[para.vidsub2[0]] += dxpro2[para.vidsub2[0]];
				dxpro[para.vidsub2[1]] += -dxpro2[para.vidsub2[0]];
				//second derivative
				// i == j, inside the interation sublattice, ortherwise is zero
				if (para.order > 0) //0阶参数, 为0，跳过这个参数
				{
					d2xpro[para.vidsub2[0]][para.vidsub2[0]] = 2 * xpro_all / yfracs[para.idsub2[0]]
						* para.order*pow(xdiff, para.order - 1);
					d2xpro[para.vidsub2[1]][para.vidsub2[1]] = -2 * xpro_all / yfracs[para.idsub2[1]]
						* para.order*pow(xdiff, para.order - 1);
					if (para.order > 1)
					{
						d2xpro2[para.vidsub2[0]][para.vidsub2[0]] = xpro_all*para.order
							*(para.order - 1)*pow(xdiff, para.order - 2);
					}
					d2xpro[para.vidsub2[0]][para.vidsub2[0]] += d2xpro2[para.vidsub2[0]][para.vidsub2[0]];
					d2xpro[para.vidsub2[1]][para.vidsub2[1]] += d2xpro2[para.vidsub2[0]][para.vidsub2[0]];
				}
				// i!=j, 
				// both i and j inside the interaction sublattice
				d2xpro[para.vidsub2[0]][para.vidsub2[1]] = xpro_all
					/ yfracs[para.idsub2[0]] / yfracs[para.idsub2[1]]
					* pow(xdiff, para.order);
				d2xpro2[para.vidsub2[0]][para.vidsub2[1]] = 0;
				d2xpro3[para.vidsub2[0]][para.vidsub2[1]] = 0;
				d2xpro4[para.vidsub2[0]][para.vidsub2[1]] = 0;
				if (para.order > 0)
				{
					d2xpro2[para.vidsub2[0]][para.vidsub2[1]] = xpro_all
						/ yfracs[para.idsub2[1]] * para.order*pow(xdiff, para.order - 1);
					d2xpro3[para.vidsub2[0]][para.vidsub2[1]] = -xpro_all
						/ yfracs[para.idsub2[0]] * para.order*pow(xdiff, para.order - 1);
					if (para.order > 1)
					{
						d2xpro4[para.vidsub2[0]][para.vidsub2[1]] = -xpro_all*para.order
							*(para.order - 1)*pow(xdiff, para.order - 2);
					}
				}
				d2xpro[para.vidsub2[0]][para.vidsub2[1]] += d2xpro2[para.vidsub2[0]][para.vidsub2[1]]
					+ d2xpro3[para.vidsub2[0]][para.vidsub2[1]]
					+ d2xpro4[para.vidsub2[0]][para.vidsub2[1]];
				d2xpro[para.vidsub2[1]][para.vidsub2[0]] = d2xpro[para.vidsub2[0]][para.vidsub2[1]];
				// only i inside the interaction sublattice
				//
				for (i = 0; i < para.vyn; i++)
				{
					for (j = i+1; j < para.vyn; j++)
					{
						d2xpro[para.vyidv[i]][para.vyidv[j]] = xpro_all / yfracs[para.vyidc[i]]
							/ yfracs[para.vyidc[j]] * pow(xdiff, para.order);
						if (para.vyidc[i] == para.vidsub2[0]&& para.vyidc[j] != para.vidsub2[1])
						{
							d2xpro2[para.vyidv[i]][para.vyidv[j]] = 0;

							if (para.order > 0)
							{
								d2xpro2[para.vyidv[i]][para.vyidv[j]] = xpro_all
									/ yfracs[para.vyidc[j]] * para.order*pow(xdiff, para.order - 1);
							}
							d2xpro[para.vyidv[i]][para.vyidv[j]] +=
								d2xpro2[para.vyidv[i]][para.vyidv[j]];
							d2xpro[para.vyidv[j]][para.vyidv[i]] = d2xpro[para.vyidv[i]][para.vyidv[j]];
						}
						else if(para.vyidc[i] == para.vidsub2[1] && para.vyidc[j] != para.vidsub2[0])
						{
							d2xpro2[para.vyidv[i]][para.vyidv[j]] = 0;

							if (para.order > 0)
							{
								d2xpro2[para.vyidv[i]][para.vyidv[j]] = -xpro_all
									/ yfracs[para.vyidc[j]] * para.order*pow(xdiff, para.order - 1);
							}
							d2xpro[para.vyidv[i]][para.vyidv[j]] +=
								d2xpro2[para.vyidv[i]][para.vyidv[j]];
							d2xpro[para.vyidv[j]][para.vyidv[i]] = d2xpro[para.vyidv[i]][para.vyidv[j]];
						}
						else if ((para.vyidc[i] != para.vidsub2[1] && para.vyidc[j] != para.vidsub2[0])
							&& (para.vyidc[i] != para.vidsub2[0] && para.vyidc[j] != para.vidsub2[1]))
						{
							d2xpro2[para.vyidv[i]][para.vyidv[j]] = 0;

							//both i and j are not inside the inter sub, part2 == 0;
							d2xpro[para.vyidv[i]][para.vyidv[j]] +=
								d2xpro2[para.vyidv[i]][para.vyidv[j]];
							d2xpro[para.vyidv[j]][para.vyidv[i]] = d2xpro[para.vyidv[i]][para.vyidv[j]];
						}
						
					}
				}
				break;
			case 3:  // 三元相互作用参数
				xpro = 0;
				for (i = 0; i < para.vyn; i++)
				{
					dxpro[para.vyidv[i]] = 0;
					for (j = i + 1; j < para.vyn; j++)
					{
						d2xpro[para.vyidv[i]][para.vyidv[j]] = 0;
						d2xpro[para.vyidv[j]][para.vyidv[i]] = d2xpro[para.vyidv[i]][para.vyidv[j]];
					}
				}
				//亚点阵内部三元相互作用，求
				break;
			case 4:  // 三元相互作用参数
				xpro = 0;
				for (i = 0; i < para.vyn; i++)
				{
					dxpro[i] = 0;
					for (j = i + 1; j < para.vyn; j++)
					{
						d2xpro[i][j] = 0;
						d2xpro[j][i] = d2xpro[i][j];
					}
				}
				//亚点阵内部三元相互作用，求
				break;
			default:
				break;
			}

			// 累加到Gibbs自由能中
			CalcPara(paravalue, Tnseg, para, T);
			if (para.type == "G" || para.type == "L")
			{
				Phases.GF += xpro*paravalue[0];
				Phases.dGFt += xpro*paravalue[1];
				Phases.d2GFt += xpro*paravalue[2];
				for (i = 0;i < para.vyn;i++)
				{
					Phases.dGFx[para.vyidv[i]] += dxpro[para.vyidv[i]] * paravalue[0];
					Phases.d2GFxt[para.vyidv[i]] += dxpro[para.vyidv[i]] * paravalue[1];
					for (j = i; j < para.vyn;j++)
					{
						Phases.d2GFx[para.vyidv[i]][para.vyidv[j]] += d2xpro[para.vyidv[i]][para.vyidv[j]] * paravalue[0];
						Phases.d2GFx[para.vyidv[j]][para.vyidv[i]] = Phases.d2GFx[para.vyidv[i]][para.vyidv[j]];
					}
				}
			}
			else if (para.type == "TC")
				TC += xpro*CalcPara(Tnseg, para, T);
			else if (para.type == "BMAGN")
				BMAGN += xpro*CalcPara(Tnseg, para, T);
		}
		//存储磁性参数，用于其他性质的计算，比如一阶、二阶导数
		Phases.TC = TC;
		Phases.BMAGN = BMAGN;
		//磁性
		tao = T / TC;
		//GF += R*T*log(BMAGN + 1)*CalcGtao(tao, Phases.magp, 1);
		
		// 熵
		// 循环所有的y
		double natom;
		natom = 0;
		Phases.SF = 0;
		for (i = 0; i < Phases.yn; i++)
		{
			Phases.SF += Phases.ysp[i] * Phases.y[i] * log(Phases.y[i]);
			if (Phases.yide[i]>1)  //Calc the mole atom in this formular, skip vacancy
				natom += Phases.ysp[i] * Phases.y[i];
		}
		Phases.SF = R*T*Phases.SF;
		Phases.GF += Phases.SF;

		Phases.dSFt = Phases.SF / T;
		Phases.d2SFt = 0;
		Phases.dGFt += Phases.dSFt;
		Phases.d2GFt += Phases.d2SFt;
		for (i = 0; i < Phases.vyn; i++)
		{
			Phases.dSFx[i] = Phases.ysp[Phases.vyidc[i]] * (log(yfracs[Phases.vyidc[i]]) + 1);
			Phases.dSFx[i] = R*T*Phases.dSFx[i];
			Phases.dGFx[i] += Phases.dSFx[i];
			Phases.d2SFxt[i] = Phases.dSFx[i] / T;
			Phases.d2GFxt[i] += Phases.d2SFxt[i];
			//
			Phases.d2SFx[i][i] = Phases.ysp[Phases.vyidc[i]] /yfracs[Phases.vyidc[i]];
			Phases.d2SFx[i][i] = R*T*Phases.d2SFx[i][i];
			Phases.d2GFx[i][i] += Phases.d2SFx[i][i];
		}
		if(mode==11)
			Phases.GF = Phases.GF / natom;
	}

	//计算Gibbs自由能
	double MinEnergy::CalcGF(Phase & Phases, double T)
	{
		//
		int i;
		double R = 8.31451;
		double xpro, xdiff;
		double GF = 0;
		double SF = 0;
		double TC = 0;
		double BMAGN = 0;
		double tao;
		//
		double yfracs[MDim];
		for (i = 0; i < Phases.yn; i++)
			yfracs[i] = Phases.y[i];

		// 循环相参数
		for (auto para : Phases.Parameters)
		{
			xpro = 1;
			xdiff = 1;
			// 判断温度区间
			int Tnseg = 0;
			for (i = 0; i < para.T.size() - 1;i++)
			{
				if ((para.T[i] - 0.02)< T && T <(para.T[i+1] + 0.02))   // T - 0.02，防止298.14和298.15这种误差
					break;
				Tnseg++;
			}
			// 计算相参数前面的系数
			for (i = 0; i < para.yn; i++)
			{
				xpro = xpro*yfracs[para.yidc[i]];
			}
			// 根据参数类型，进行修正
			switch (para.kind)
			{
			case 2:  // 二元相互作用参数
					 //亚点阵内部二元相互作用，求(x1-x2)
				xdiff = yfracs[para.idsub2[0]] - yfracs[para.idsub2[1]];
				xpro = xpro*pow(xdiff, para.order);
				break;
			case 3:  // 三元相互作用参数
				xpro = 0;
					 //亚点阵内部三元相互作用，求
				break;
			case 4:  // 三元相互作用参数
				xpro = 0;
					 //亚点阵内部三元相互作用，求
				break;
			}

			// 累加到Gibbs自由能中
			if (para.type == "G" || para.type == "L")
				GF += xpro*CalcPara(Tnseg, para, T);
			else if (para.type == "TC")
				TC += xpro*CalcPara(Tnseg, para, T);
			else if (para.type == "BMAGN")
				BMAGN += xpro*CalcPara(Tnseg, para, T);
		}
		//存储磁性参数，用于其他性质的计算，比如一阶、二阶导数
		Phases.TC = TC;
		Phases.BMAGN = BMAGN;
		//磁性
		tao = T / TC;
		GF += R*T*log(BMAGN + 1)*CalcGtao(tao, Phases.magp, 1);

		// 熵
		// 循环所有的y
		double natom;
		natom = 0;
		for (i = 0; i < Phases.yn; i++)
		{
			if (Phases.y[i]<eps)
				SF += Phases.ysp[i] *eps * log(eps);
			else
				SF += Phases.ysp[i] * Phases.y[i] * log(Phases.y[i]);

			if (Phases.yide[i]>1)  //不包括空位
				natom += Phases.ysp[i] * Phases.y[i];
		}
		SF = R*T*SF;
		GF += SF;
		//GF = GF / natom;
		return GF;
	}
	//
	void MinEnergy::CalcGPGF(Phase & Phases, double T)
	{
		//
		int i, j, k;
		int GPn;
		double R = 8.31451;
		double *xpro, *xdiff;
		double paravalue;

		GPn = Phases.GPn;
		//申请空间
		
		xpro = new double[GPn];
		xdiff = new double[GPn];
		//初始化
		for (i = 0; i < GPn; i++)
		{
			Phases.GPGF[i] = 0;
			Phases.GPSF[i] = 0;
			Phases.GPTC[i] = 0;
			Phases.GPBMAGN[i] = 0;
		}
		if (Phases.magp != 0)
		{
			for (i = 0; i < GPn; i++)
			{
				Phases.GPTC[i] = 0;
				Phases.GPBMAGN[i] = 0;
			}
		}
		

		// 循环相参数
		for (auto para : Phases.Parameters)
		{
			for (i = 0; i < GPn; i++)
				xpro[i] = 1;
			// 判断温度区间
			int Tnseg = 0;
			for (i = 0; i < para.T.size() - 1; i++)
			{
				if ((para.T[i] - 0.02)< T && T <(para.T[i + 1] + 0.02))   // T - 0.02，防止298.14和298.15这种误差
					break;
				Tnseg++;
			}
			// 计算相参数前面的系数
			for (i = 0; i < para.yn; i++)
				for (j = 0; j < Phases.GPn; j++)
					xpro[j] = xpro[j] * Phases.GP[j].y[para.yidc[i]];
			// 根据参数类型，进行修正
			switch (para.kind)
			{
			case 2:  // 二元相互作用参数
					 //亚点阵内部二元相互作用，求(x1-x2)
				for (j = 0; j < Phases.GPn; j++)
				{
					xdiff[j] = Phases.GP[j].y[para.idsub2[0]] - Phases.GP[j].y[para.idsub2[1]];
					xpro[j] = xpro[j] * pow(xdiff[j], para.order);
				}
				break;
			case 3:  // 三元相互作用参数
				xpro = 0;

					 //亚点阵内部三元相互作用，求
				break;
			case 4:  // 三元相互作用参数
				xpro = 0;

					 //亚点阵内部三元相互作用，求
				break;
			}

			// 累加到Gibbs自由能中
			paravalue = CalcPara(Tnseg, para, T);
			if (para.type == "G" || para.type == "L")
			{
				for (j = 0; j < GPn; j++)
					Phases.GPGF[j] += xpro[j] * paravalue;
			}
			else if (para.type == "TC")
			{
				for (j = 0; j < GPn; j++)
					Phases.GPTC[j] += xpro[j] * paravalue;
			}
			else if (para.type == "BMAGN")
			{
				for (j = 0; j < GPn; j++)
					Phases.GPBMAGN[j] += xpro[j] * paravalue;
			}
		}
		//存储磁性参数，用于其他性质的计算，比如一阶、二阶导数
		//磁性
		if (Phases.magp != 0)
		{
			for (j = 0; j < GPn; j++)
			{
				Phases.GPtao[j] = T / Phases.GPTC[j];
				Phases.GPGF[j] = R*T*log(Phases.GPBMAGN[j] + 1)*CalcGtao(Phases.GPtao[j], Phases.magp, 1);
			}
		}

		// 熵
		//循环亚点阵
		double natom;
		for (j = 0; j < GPn; j++)
		{
			natom = 0;
			for (k = 0; k < Phases.yn; k++)
			{
				if(Phases.GP[j].y[k]<eps)
					Phases.GPSF[j] += Phases.ysp[k] * eps * log(eps);
				else
					Phases.GPSF[j] += Phases.ysp[k] * Phases.GP[j].y[k] * log(Phases.GP[j].y[k]);

				if (Phases.yide[k]>1)  //不包括空位
					natom += Phases.ysp[k] * Phases.GP[j].y[k];
			}
			Phases.GPSF[j] = R*T*Phases.GPSF[j];
			Phases.GPGF[j] += Phases.GPSF[j];
			//计算1mol化学式中包含多少摩尔原子数，不包括空位
			Phases.GPGF[j] = Phases.GPGF[j]/natom;
		}

		//delete[] xpro;
		//delete[] xdiff;
	}

	// 计算Gibbs自由能一阶导数
	double MinEnergy::CalcGFdx(Phase & Phases, double T, int idy)
	{
		//
		int i;
		double R = 8.31451;
		double xpro, xpro2, xdiff;
		double dGF = 0;
		double dSF = 0;
		double dTC = 0;
		double dBMAGN = 0;
		double tao;

		//
		double yfracs[MDim];
		for (i = 0; i < Phases.yn; i++)
		{
			
			yfracs[i] = Phases.y[i];

		}

		//G和L
		for (auto para : Phases.Parameters)
		{
			//判断idy是否在参数中，不在，导数为0
			int checkcon = 0;
			int naddcon = 0;
			for (i = 0; i < para.yn; i++)
			{
				if (para.yidc[i] == idy)
					checkcon = 1;
			}
			if (checkcon == 0) //idy不在参数列表中，该参数对的一阶导数为0，跳过。
				continue;

			// 判断温度区间
			int Tnseg = 0;
			for (i = 0; i < para.T.size() - 1; i++)
			{
				if ((para.T[i] - 0.02)< T && T <(para.T[i + 1] + 0.02))   // T - 0.02，防止298.14和298.15这种误差
					break;
				Tnseg++;
			}
			xpro = 1;
			xpro2 = 1;
			xdiff = 1;
			// 计算相参数前面的系数
			for (i = 0; i < para.yn; i++)
			{
				xpro = xpro*yfracs[para.yidc[i]];
			}
			xpro = xpro / yfracs[idy]; // %求导数之后为1,（yfracs[idy]/yfracs[idy] =1）
									   // 根据参数类型，进行修正
			switch (para.kind)
			{

			case 2:  // 二元相互作用参数
					 //第一部分
					 //亚点阵内部二元相互作用，求(x1-x2)
				xdiff = yfracs[para.idsub2[0]] - yfracs[para.idsub2[1]];
				xpro = xpro*pow(xdiff, para.order);
				//第二部分, idy在相互作用的亚点阵内部，Xi*Xj*Xk(Xi-Xj)^n，idy是Xi或者Xj
				if (para.order > 0)  // n=0,求导之后为0，不考虑。
				{
					if (idy == para.idsub2[0])  // (Xi-Xj)^n, idy ==Xi
					{
						for (i = 0; i < Phases.subln + 1; i++)
						{
							xpro2 = xpro2*yfracs[para.yidc[i]];
						}
						xpro2 = xpro2*para.order*pow(xdiff, para.order - 1);
					}
					else if ((idy == para.idsub2[1]))  // (Xi-Xj)^n, idy ==Xj
					{
						for (i = 0; i < Phases.subln + 1; i++)
						{
							xpro2 = xpro2*yfracs[para.yidc[i]];
						}
						xpro2 = -xpro2*para.order*pow(xdiff, para.order - 1);
					}
					xpro += xpro2;
				}
				break;
			case 3:  // 三元相互作用参数
				xpro = 0;

					 //亚点阵内部三元相互作用，求
				break;
			case 4:
				xpro = 0;
				break;
			}

			// 累加到Gibbs自由能中

			if (para.type == "G" || para.type == "L")
				dGF += xpro*CalcPara(Tnseg, para, T);
			else if (para.type == "TC")
				dTC += xpro*CalcPara(Tnseg, para, T);
			else if (para.type == "BMAGN")
				dBMAGN += xpro*CalcPara(Tnseg, para, T);
		}
		//磁性
		if (Phases.magp != 0)
		{
			tao = T / Phases.TC;
			dGF += R*T*(CalcGtao(tao, Phases.magp, 0) / (Phases.BMAGN + 1)*dBMAGN -
				T / Phases.TC / Phases.TC*log(Phases.BMAGN + 1)*CalcGtao(tao, Phases.magp, 1)*dTC);
		}

		// 熵
		
			dSF += Phases.ysp[idy] * (log(yfracs[idy]) + 1);
		dSF = R*T*dSF;
		dGF += dSF;
		//计算摩尔原子的G
		double natom;
		natom = 0;
		for (i = 0; i < Phases.yn; i++)
		{
			if (Phases.yide[i]>1)  //不包括空位
				natom += Phases.ysp[i] * Phases.y[i];
		}
		//dGF = dGF / natom;

		return dGF;
	}
	// 计算Gibbs自由能二阶导数
	double MinEnergy::CalcGFdx2(Phase & Phases, double T, int idy1, int idy2)
	{
		//
		int i;
		double R = 8.31451;
		double xpro, xpro2, xpro3, xpro4, xdiff;
		int signx1, signx2;
		double dGF2 = 0;
		double dSF2 = 0;
		int checkcon1, checkcon2;

		//
		double yfracs[MDim];
		for (i = 0; i < Phases.yn; i++)
			yfracs[i] = Phases.y[i];

		//G和L
		for (auto para : Phases.Parameters)
		{
			//判断idy是否在参数中，不在，导数为0
			checkcon1 = -1;
			checkcon2 = -1;
			for (i = 0; i < para.yn; i++)
			{
				if (para.yidc[i] == idy1)
					checkcon1 = 1;
				if (para.yidc[i] == idy2)
					checkcon2 = 1;
			}
			if (checkcon1 == -1 || checkcon2 == -1) //idy1和idy2中只要有一个不在列表中，该参数对的二阶导数为0，跳过。
				continue;
			// 判断温度区间
			int Tnseg = 0;
			for (i = 0; i < para.T.size() - 1; i++)
			{
				if ((para.T[i] - 0.02)< T && T <(para.T[i + 1] + 0.02))   // T - 0.02，防止298.14和298.15这种误差
					break;
				Tnseg++;
			}
			// 成分系数
			xpro = 1;
			xpro2 = 0;
			xpro3 = 0;
			xpro4 = 0;
			// 判断二阶导数类型
			if (idy1 == idy2)
			{
				switch (para.kind)
				{
				case 1:
					continue;  // end member， 二阶导数为0，跳过这个参数
					break;
				case 2:
					if (para.order == 0)
						continue;  //0阶参数, 为0，跳过这个参数
								   //求导后的符号
					signx1 = 0;
					if (idy1 == para.idsub2[0])
						signx1 = 1;
					else
						signx1 = -1;
					if (signx1 == 0)
						continue; //不在相互作用的亚点阵中, 为0，跳过这个参数
								  //第一部分
					for (i = 0; i < para.yn; i++)
					{
						xpro = xpro*yfracs[para.yidc[i]];
					}
					xdiff = yfracs[para.idsub2[0]] - yfracs[para.idsub2[1]];  //亚点阵内部二元相互作用，求(x1-x2)
					if (para.order > 1)
						xpro2 = xpro* para.order*(para.order - 1)*pow(xdiff, para.order - 2);  // yi*yj*yk*n*(n-1)*(yi-yk)^(n-2)
					xpro = signx1 * 2 * xpro / yfracs[idy1] * para.order*pow(xdiff, para.order - 1); // 2*yj*yk*n*(yi-yk)^(n-1)
					xpro = xpro + xpro2;
					break;
				case 3:
					xpro = 0;
					break;
				case 4:
					xpro = 0;
					break;
				}
			}
			else    // idy1 != idy2
			{
				switch (para.kind)
				{
				case 1:   // end member  yk*G(i:j:k)
					for (i = 0; i < para.yn; i++)
					{
						xpro = xpro*yfracs[para.yidc[i]];
					}
					xpro = xpro / yfracs[idy1] / yfracs[idy2];; // %求导数之后为1
					break;
				case 2:  // 二元相互作用参数
						 //第一部分
					for (i = 0; i < para.yn; i++)
					{
						xpro = xpro*yfracs[para.yidc[i]];
					}
					xdiff = yfracs[para.idsub2[0]] - yfracs[para.idsub2[1]];
					//第二部分
					//求导后的符号
					signx1 = 0;
					signx2 = 0;
					if (idy1 == para.idsub2[0])
						signx1 = 1;
					else if (idy1 == para.idsub2[1])
						signx1 = -1;
					if (idy2 == para.idsub2[0])
						signx2 = 1;
					else if (idy2 == para.idsub2[1])
						signx2 = -1;
					if (signx1 != 0 && signx2 != 0) // 第一种情况：都在相互作用亚点阵内部
					{
						if (para.order > 0)
						{
							xpro2 = signx1*xpro / yfracs[idy2] * para.order*pow(xdiff, para.order - 1); // yi*yk*n*(yi-yj)^(n-1)
							xpro3 = signx2*xpro / yfracs[idy1] * para.order*pow(xdiff, para.order - 1); // -yj*yk*n*(yi-yj)^(n-1)
							if (para.order > 1)
								xpro4 = -xpro*para.order*(para.order - 1)*pow(xdiff, para.order - 2); // -yi*yj*yk*n*(n-1)*(yi-yj)^(n-2)
							xpro = xpro / yfracs[idy1] / yfracs[idy2] * pow(xdiff, para.order);  // yk*n*(yi-yj)^n
						}
						else
							xpro = xpro / yfracs[idy1] / yfracs[idy2] * pow(xdiff, para.order);  // yk*n*(yi-yj)^n
						xpro = xpro + xpro2 + xpro3 + xpro4;
					}
					else if (signx1 != 0 && signx2 == 0) // 第二种情况：只有一个在相互作用的亚点阵内部
					{
						if (para.order > 0)
							xpro2 = signx1*xpro / yfracs[idy2] * para.order*pow(xdiff, para.order - 1);
						xpro = xpro / yfracs[idy1] / yfracs[idy2] * pow(xdiff, para.order);
						xpro = xpro + xpro2;
					}
					else if (signx1 == 0 && signx2 != 0) // 第二种情况：只有一个在相互作用的亚点阵内部
					{
						if (para.order > 0)
							xpro2 = signx2*xpro*yfracs[idy1] * para.order*pow(xdiff, para.order - 1);
						xpro = xpro / yfracs[idy1] / yfracs[idy2] * pow(xdiff, para.order);
						xpro = xpro + xpro2;
					}
					else if (signx1 != 0 && signx2 != 0) // 第四种情况：都不在相互作用的亚点阵内部
					{
						xpro = xpro / yfracs[idy1] / yfracs[idy2] * pow(xdiff, para.order);
					}
					break;
				case 3:
					for (i = 0; i < para.yn; i++)
					{
						xpro = xpro*yfracs[para.yidc[i]];
					}
					xpro = xpro / yfracs[idy1] / yfracs[idy2];; // %求导数之后为1
					break;
				case 4:
					for (i = 0; i < para.yn; i++)
					{
						xpro = xpro*yfracs[para.yidc[i]];
					}
					xpro = xpro / yfracs[idy1] / yfracs[idy2];; // %求导数之后为1
					break;
				}

			}
			// 累加到Gibbs自由能中
			dGF2 += xpro*CalcPara(Tnseg, para, T);
		}
		// 熵

		if (idy1 == idy2)
		{

			dSF2 += Phases.ysp[idy1] / yfracs[idy1];
		}
		dSF2 = R*T*dSF2;
		dGF2 += dSF2;
		//计算摩尔原子的G
		double natom;
		natom = 0;
		for (i = 0; i < Phases.yn; i++)
		{
			if (Phases.yide[i]>1)  //不包括空位
				natom += Phases.ysp[i] * Phases.y[i];
		}
		//dGF2 = dGF2 / natom;
		//if (dGF2 > 1e9)
		//dGF2 = 1e9;
		return dGF2;
	}
	//
	//
	double MinEnergy::CalcPara(int Tnseg, Parameter para, double T)
	{
		int k = 0;
		double paravalue = 0;
		for (auto coffT : para.express_digit[Tnseg].coffT)
		{
			paravalue += coffT*pow(T, para.express_digit[Tnseg].powerT[k]);
			k++;
		}
		// 计算T*LN(T)
		for (auto coffTLNT : para.express_digit[Tnseg].coffTLNT)
		{
			paravalue += coffTLNT*T*log(T);
		}
		return paravalue;
	}
	//
	void MinEnergy::CalcPara(double pvalue[3], int Tnseg, Parameter para, double T)
	{
		int k = 0;
		double temp0, temp1;
		pvalue[0] = 0;
		pvalue[1] = 0;
		pvalue[2] = 0;
		for (auto coffT : para.express_digit[Tnseg].coffT)
		{
			temp0 = coffT*pow(T, para.express_digit[Tnseg].powerT[k]);
			pvalue[0] += temp0;
			if (para.express_digit[Tnseg].powerT[k] != 0)
			{
				temp1 = para.express_digit[Tnseg].powerT[k] * temp0 / T;
				pvalue[1] += temp1;
				if (para.express_digit[Tnseg].powerT[k] != 1)
				{
					pvalue[2] += (para.express_digit[Tnseg].powerT[k] - 1) *temp1 / T;
				}
			}
			k++;
		}
		// 计算T*LN(T)
		for (auto coffTLNT : para.express_digit[Tnseg].coffTLNT)
		{
			pvalue[0] += coffTLNT*T*log(T);
			pvalue[1] += coffTLNT*(log(T) + 1);
			pvalue[2] += coffTLNT / T;
		}
	}
	double MinEnergy::CalcGtao(double tao, double p, int order)
	{
		double  Gtao;
		double A = 0.460444444444444 + 0.731893583724570*(1 / p - 1);
		if (order == 0)
		{
			if (tao < 1)
				Gtao = 1 - 1 / A*(0.564285714285714 / tao / p + 0.953722334004024*
					(1 / p - 1)*(pow(tao, 3) / 6 + pow(tao, 9) / 135 + pow(tao, 15) / 600));
			else
				Gtao = -1 / A*(pow(tao, -5) / 10 + pow(tao, -15) / 315 + pow(tao, -25) / 1500);
		}
		else if (order == 1)
		{
			if (tao < 1)
				Gtao = 1 - 1 / A*(0.564285714285714 / tao / p + 0.953722334004024*
					(1 / p - 1)*(pow(tao, 3) / 6 + pow(tao, 9) / 135 + pow(tao, 15) / 600));
			else
				Gtao = -1 / A*(pow(tao, -5) / 10 + pow(tao, -15) / 315 + pow(tao, -25) / 1500);
		}
		else if (order == 2)
		{
			if (tao < 1)
				Gtao = -1 / A*(0.564285714285714 / tao / tao / p + 0.953722334004024*
					(1 / p - 1)*(pow(tao, 2) / 3 + pow(tao, 8) / 15 + pow(tao, 14) / 40));
			else
				Gtao = 1 / A*(pow(tao, -6) / 2 + pow(tao, -16) / 21 + pow(tao, -26) / 60);
		}
		return Gtao;
	}

=======
#include "..\include\MinEnergy.h"

namespace VCLab
{
	using namespace std;
	/*---------------------------------------------------------------------*/
	/*                       计算Gibbs自由能                               */
	/*---------------------------------------------------------------------*/
	/*
	void MinEnergy::CalcGE(Phase & Phases, double T)
	double MinEnergy::CalcGF(Phase & Phases, double T)
	double MinEnergy::CalcGFdx(Phase & Phases, double T, int idy)
	double MinEnergy::CalcGFdx2(Phase & Phases, double T, int idy1, int idy2)
	double MinEnergy::CalcPara(int Tnseg, Parameter para, double T)
	double MinEnergy::CalcGtao(double tao, double p, int order)
	*/

	// Calc the Gibbs energy, the first , second derivative of site fraction and temperature
	void MinEnergy::CalcGE(Phase & Phases, double T)
	{
		//
		int i, j;
		double R = 8.31451;
		double xpro_all, xpro, xdiff;
		double dxpro[10], dxpro2[10];
		double d2xpro[10][10], d2xpro2[10][10], d2xpro3[10][10], d2xpro4[10][10];

		double TC = 0;
		double BMAGN = 0;
		double tao;
		//
		double yfracs[MDim];
		double paravalue[3];

		for (i = 0; i < Phases.yn; i++)
		{
			if(Phases.y[i]<eps)
				yfracs[i] = eps;
			else
				yfracs[i] = Phases.y[i];
		}
		//
		Phases.GF = 0; Phases.SF = 0;
		Phases.dGFt = 0; Phases.d2GFt = 0;
		for (i = 0;i < Phases.vyn;i++)
		{
			Phases.dGFx[i] = 0;
			Phases.d2GFxt[i] = 0;
			for (j = 0;j < Phases.vyn;j++)
			{
				Phases.d2GFx[i][j] = 0;
			}
		}

		// 循环相参数
		for (auto para : Phases.Parameters)
		{
			// 判断温度区间
			int Tnseg = 0;
			for (i = 0; i < para.T.size() - 1;i++)
			{
				if ((para.T[i] - 0.02)< T && T <(para.T[i + 1] + 0.02))   // T - 0.02，防止298.14和298.15这种误差
					break;
				Tnseg++;
			}
			//
			for (i = 0;i < Phases.vyn;i++)
			{
				dxpro[i] = 0;
				d2xpro[i][i] = 0;
				d2xpro2[i][i] = 0;
				d2xpro3[i][i] = 0;
				d2xpro4[i][i] = 0;

				for (j = i + 1;j < Phases.vyn;j++)
				{
					d2xpro[i][j] = 0; d2xpro[j][i] = 0;
					d2xpro2[i][j] = 0; d2xpro2[j][i] = 0;
					d2xpro3[i][j] = 0; d2xpro3[j][i] = 0;
					d2xpro4[i][j] = 0; d2xpro4[j][i] = 0;
				}
			}
			// 计算相参数前面的系数
			xpro_all = 1;
			for (i = 0; i < para.yn; i++)
			{
				xpro_all = xpro_all*yfracs[para.yidc[i]];
			}
			// 根据参数类型，进行修正
			switch (para.kind)
			{
			case 1:
				//g
				xpro = xpro_all;
				// first derivative
				for (i = 0;i < para.vyn;i++)
				{
					dxpro[para.vyidv[i]] = xpro_all / yfracs[para.vyidc[i]];
					// second derivate
					// for i==j, in case 1, this part is zero
					for (j = i + 1;j < para.vyn;j++)
					{
						d2xpro[para.vyidv[i]][para.vyidv[j]] = dxpro[para.vyidv[i]] / yfracs[para.vyidc[j]];
						d2xpro[para.vyidv[j]][para.vyidv[i]] = d2xpro[para.vyidv[i]][para.vyidv[j]];
					}
				}
				break;
			case 2:  // 二元相互作用参数
					 //亚点阵内部二元相互作用，求(x1-x2)
				xpro = xpro_all;
				xdiff = yfracs[para.idsub2[0]] - yfracs[para.idsub2[1]];
				xpro = xpro*pow(xdiff, para.order);
				//
				//first derivative
				//亚点阵内部二元相互作用，求(x1-x2)
				for (i = 0;i < para.vyn;i++)
				{
					dxpro[para.vyidv[i]] = xpro_all/yfracs[para.vyidc[i]]*pow(xdiff, para.order);
				}
				//第二部分, idy在相互作用的亚点阵内部，Xi*Xj*Xk(Xi-Xj)^n，idy是Xi或者Xj
				dxpro2[para.vidsub2[0]] = 0;
				if (para.order > 0)  // n=0,求导之后为0，不考虑。
				{
					dxpro2[para.vidsub2[0]] = xpro_all*para.order*pow(xdiff, para.order - 1);
				}
				dxpro[para.vidsub2[0]] += dxpro2[para.vidsub2[0]];
				dxpro[para.vidsub2[1]] += -dxpro2[para.vidsub2[0]];
				//second derivative
				// i == j, inside the interation sublattice, ortherwise is zero
				if (para.order > 0) //0阶参数, 为0，跳过这个参数
				{
					d2xpro[para.vidsub2[0]][para.vidsub2[0]] = 2 * xpro_all / yfracs[para.idsub2[0]]
						* para.order*pow(xdiff, para.order - 1);
					d2xpro[para.vidsub2[1]][para.vidsub2[1]] = -2 * xpro_all / yfracs[para.idsub2[1]]
						* para.order*pow(xdiff, para.order - 1);
					if (para.order > 1)
					{
						d2xpro2[para.vidsub2[0]][para.vidsub2[0]] = xpro_all*para.order
							*(para.order - 1)*pow(xdiff, para.order - 2);
					}
					d2xpro[para.vidsub2[0]][para.vidsub2[0]] += d2xpro2[para.vidsub2[0]][para.vidsub2[0]];
					d2xpro[para.vidsub2[1]][para.vidsub2[1]] += d2xpro2[para.vidsub2[0]][para.vidsub2[0]];
				}
				// i!=j, 
				// both i and j inside the interaction sublattice
				d2xpro[para.vidsub2[0]][para.vidsub2[1]] = xpro_all
					/ yfracs[para.idsub2[0]] / yfracs[para.idsub2[1]]
					* pow(xdiff, para.order);
				d2xpro2[para.vidsub2[0]][para.vidsub2[1]] = 0;
				d2xpro3[para.vidsub2[0]][para.vidsub2[1]] = 0;
				d2xpro4[para.vidsub2[0]][para.vidsub2[1]] = 0;
				if (para.order > 0)
				{
					d2xpro2[para.vidsub2[0]][para.vidsub2[1]] = xpro_all
						/ yfracs[para.idsub2[1]] * para.order*pow(xdiff, para.order - 1);
					d2xpro3[para.vidsub2[0]][para.vidsub2[1]] = -xpro_all
						/ yfracs[para.idsub2[0]] * para.order*pow(xdiff, para.order - 1);
					if (para.order > 1)
					{
						d2xpro4[para.vidsub2[0]][para.vidsub2[1]] = -xpro_all*para.order
							*(para.order - 1)*pow(xdiff, para.order - 2);
					}
				}
				d2xpro[para.vidsub2[0]][para.vidsub2[1]] += d2xpro2[para.vidsub2[0]][para.vidsub2[1]]
					+ d2xpro3[para.vidsub2[0]][para.vidsub2[1]]
					+ d2xpro4[para.vidsub2[0]][para.vidsub2[1]];
				d2xpro[para.vidsub2[1]][para.vidsub2[0]] = d2xpro[para.vidsub2[0]][para.vidsub2[1]];
				// only i inside the interaction sublattice
				//
				for (i = 0; i < para.vyn; i++)
				{
					for (j = i+1; j < para.vyn; j++)
					{
						d2xpro[para.vyidv[i]][para.vyidv[j]] = xpro_all / yfracs[para.vyidc[i]]
							/ yfracs[para.vyidc[j]] * pow(xdiff, para.order);
						if (para.vyidc[i] == para.vidsub2[0]&& para.vyidc[j] != para.vidsub2[1])
						{
							d2xpro2[para.vyidv[i]][para.vyidv[j]] = 0;

							if (para.order > 0)
							{
								d2xpro2[para.vyidv[i]][para.vyidv[j]] = xpro_all
									/ yfracs[para.vyidc[j]] * para.order*pow(xdiff, para.order - 1);
							}
							d2xpro[para.vyidv[i]][para.vyidv[j]] +=
								d2xpro2[para.vyidv[i]][para.vyidv[j]];
							d2xpro[para.vyidv[j]][para.vyidv[i]] = d2xpro[para.vyidv[i]][para.vyidv[j]];
						}
						else if(para.vyidc[i] == para.vidsub2[1] && para.vyidc[j] != para.vidsub2[0])
						{
							d2xpro2[para.vyidv[i]][para.vyidv[j]] = 0;

							if (para.order > 0)
							{
								d2xpro2[para.vyidv[i]][para.vyidv[j]] = -xpro_all
									/ yfracs[para.vyidc[j]] * para.order*pow(xdiff, para.order - 1);
							}
							d2xpro[para.vyidv[i]][para.vyidv[j]] +=
								d2xpro2[para.vyidv[i]][para.vyidv[j]];
							d2xpro[para.vyidv[j]][para.vyidv[i]] = d2xpro[para.vyidv[i]][para.vyidv[j]];
						}
						else if ((para.vyidc[i] != para.vidsub2[1] && para.vyidc[j] != para.vidsub2[0])
							&& (para.vyidc[i] != para.vidsub2[0] && para.vyidc[j] != para.vidsub2[1]))
						{
							d2xpro2[para.vyidv[i]][para.vyidv[j]] = 0;

							//both i and j are not inside the inter sub, part2 == 0;
							d2xpro[para.vyidv[i]][para.vyidv[j]] +=
								d2xpro2[para.vyidv[i]][para.vyidv[j]];
							d2xpro[para.vyidv[j]][para.vyidv[i]] = d2xpro[para.vyidv[i]][para.vyidv[j]];
						}
						
					}
				}
				break;
			case 3:  // 三元相互作用参数
				xpro = 0;
				for (i = 0; i < para.vyn; i++)
				{
					dxpro[para.vyidv[i]] = 0;
					for (j = i + 1; j < para.vyn; j++)
					{
						d2xpro[para.vyidv[i]][para.vyidv[j]] = 0;
						d2xpro[para.vyidv[j]][para.vyidv[i]] = d2xpro[para.vyidv[i]][para.vyidv[j]];
					}
				}
				//亚点阵内部三元相互作用，求
				break;
			case 4:  // 三元相互作用参数
				xpro = 0;
				for (i = 0; i < para.vyn; i++)
				{
					dxpro[i] = 0;
					for (j = i + 1; j < para.vyn; j++)
					{
						d2xpro[i][j] = 0;
						d2xpro[j][i] = d2xpro[i][j];
					}
				}
				//亚点阵内部三元相互作用，求
				break;
			default:
				break;
			}

			// 累加到Gibbs自由能中
			CalcPara(paravalue, Tnseg, para, T);
			if (para.type == "G" || para.type == "L")
			{
				Phases.GF += xpro*paravalue[0];
				Phases.dGFt += xpro*paravalue[1];
				Phases.d2GFt += xpro*paravalue[2];
				for (i = 0;i < para.vyn;i++)
				{
					Phases.dGFx[para.vyidv[i]] += dxpro[para.vyidv[i]] * paravalue[0];
					Phases.d2GFxt[para.vyidv[i]] += dxpro[para.vyidv[i]] * paravalue[1];
					for (j = i; j < para.vyn;j++)
					{
						Phases.d2GFx[para.vyidv[i]][para.vyidv[j]] += d2xpro[para.vyidv[i]][para.vyidv[j]] * paravalue[0];
						Phases.d2GFx[para.vyidv[j]][para.vyidv[i]] = Phases.d2GFx[para.vyidv[i]][para.vyidv[j]];
					}
				}
			}
			else if (para.type == "TC")
				TC += xpro*CalcPara(Tnseg, para, T);
			else if (para.type == "BMAGN")
				BMAGN += xpro*CalcPara(Tnseg, para, T);
		}
		//存储磁性参数，用于其他性质的计算，比如一阶、二阶导数
		Phases.TC = TC;
		Phases.BMAGN = BMAGN;
		//磁性
		tao = T / TC;
		//GF += R*T*log(BMAGN + 1)*CalcGtao(tao, Phases.magp, 1);
		
		// 熵
		// 循环所有的y
		double natom;
		natom = 0;
		Phases.SF = 0;
		for (i = 0; i < Phases.yn; i++)
		{
			Phases.SF += Phases.ysp[i] * Phases.y[i] * log(Phases.y[i]);
			if (Phases.yide[i]>1)  //Calc the mole atom in this formular, skip vacancy
				natom += Phases.ysp[i] * Phases.y[i];
		}
		Phases.SF = R*T*Phases.SF;
		Phases.GF += Phases.SF;

		Phases.dSFt = Phases.SF / T;
		Phases.d2SFt = 0;
		Phases.dGFt += Phases.dSFt;
		Phases.d2GFt += Phases.d2SFt;
		for (i = 0; i < Phases.vyn; i++)
		{
			Phases.dSFx[i] = Phases.ysp[Phases.vyidc[i]] * (log(yfracs[Phases.vyidc[i]]) + 1);
			Phases.dSFx[i] = R*T*Phases.dSFx[i];
			Phases.dGFx[i] += Phases.dSFx[i];
			Phases.d2SFxt[i] = Phases.dSFx[i] / T;
			Phases.d2GFxt[i] += Phases.d2SFxt[i];
			//
			Phases.d2SFx[i][i] = Phases.ysp[Phases.vyidc[i]] /yfracs[Phases.vyidc[i]];
			Phases.d2SFx[i][i] = R*T*Phases.d2SFx[i][i];
			Phases.d2GFx[i][i] += Phases.d2SFx[i][i];
		}
		//GF = GF / natom;
	}

	//计算Gibbs自由能
	double MinEnergy::CalcGF(Phase & Phases, double T)
	{
		//
		int i;
		double R = 8.31451;
		double xpro, xdiff;
		double GF = 0;
		double SF = 0;
		double TC = 0;
		double BMAGN = 0;
		double tao;
		//
		double yfracs[MDim];
		for (i = 0; i < Phases.yn; i++)
			yfracs[i] = Phases.y[i];

		// 循环相参数
		for (auto para : Phases.Parameters)
		{
			xpro = 1;
			xdiff = 1;
			// 判断温度区间
			int Tnseg = 0;
			for (i = 0; i < para.T.size() - 1;i++)
			{
				if ((para.T[i] - 0.02)< T && T <(para.T[i+1] + 0.02))   // T - 0.02，防止298.14和298.15这种误差
					break;
				Tnseg++;
			}
			// 计算相参数前面的系数
			for (i = 0; i < para.yn; i++)
			{
				xpro = xpro*yfracs[para.yidc[i]];
			}
			// 根据参数类型，进行修正
			switch (para.kind)
			{
			case 2:  // 二元相互作用参数
					 //亚点阵内部二元相互作用，求(x1-x2)
				xdiff = yfracs[para.idsub2[0]] - yfracs[para.idsub2[1]];
				xpro = xpro*pow(xdiff, para.order);
				break;
			case 3:  // 三元相互作用参数
				xpro = 0;
					 //亚点阵内部三元相互作用，求
				break;
			case 4:  // 三元相互作用参数
				xpro = 0;
					 //亚点阵内部三元相互作用，求
				break;
			}

			// 累加到Gibbs自由能中
			if (para.type == "G" || para.type == "L")
				GF += xpro*CalcPara(Tnseg, para, T);
			else if (para.type == "TC")
				TC += xpro*CalcPara(Tnseg, para, T);
			else if (para.type == "BMAGN")
				BMAGN += xpro*CalcPara(Tnseg, para, T);
		}
		//存储磁性参数，用于其他性质的计算，比如一阶、二阶导数
		Phases.TC = TC;
		Phases.BMAGN = BMAGN;
		//磁性
		tao = T / TC;
		GF += R*T*log(BMAGN + 1)*CalcGtao(tao, Phases.magp, 1);

		// 熵
		// 循环所有的y
		double natom;
		natom = 0;
		for (i = 0; i < Phases.yn; i++)
		{
			if (Phases.y[i]<eps)
				SF += Phases.ysp[i] *eps * log(eps);
			else
				SF += Phases.ysp[i] * Phases.y[i] * log(Phases.y[i]);

			if (Phases.yide[i]>1)  //不包括空位
				natom += Phases.ysp[i] * Phases.y[i];
		}
		SF = R*T*SF;
		GF += SF;
		//GF = GF / natom;
		return GF;
	}
	//
	void MinEnergy::CalcGPGF(Phase & Phases, double T)
	{
		//
		int i, j, k;
		int GPn;
		double R = 8.31451;
		double *xpro, *xdiff;
		double paravalue;
		GPn = Phases.GPn;
		//申请空间
		Phases.GPGF = new double[GPn];
		Phases.GPSF = new double[GPn];
		Phases.GPTC = new double[GPn];
		Phases.GPBMAGN = new double[GPn];
		Phases.GPtao = new double[GPn];
		xpro = new double[GPn];
		xdiff = new double[GPn];
		//初始化
		for (i = 0; i < GPn; i++)
		{
			Phases.GPGF[i] = 0;
			Phases.GPSF[i] = 0;
			Phases.GPTC[i] = 0;
			Phases.GPBMAGN[i] = 0;
		}
		if (Phases.magp != 0)
		{
			for (i = 0; i < GPn; i++)
			{
				Phases.GPTC[i] = 0;
				Phases.GPBMAGN[i] = 0;
			}
		}
		

		// 循环相参数
		for (auto para : Phases.Parameters)
		{
			for (i = 0; i < GPn; i++)
				xpro[i] = 1;
			// 判断温度区间
			int Tnseg = 0;
			for (i = 0; i < para.T.size() - 1; i++)
			{
				if ((para.T[i] - 0.02)< T && T <(para.T[i + 1] + 0.02))   // T - 0.02，防止298.14和298.15这种误差
					break;
				Tnseg++;
			}
			// 计算相参数前面的系数
			for (i = 0; i < para.yn; i++)
				for (j = 0; j < Phases.GPn; j++)
					xpro[j] = xpro[j] * Phases.GPy[j][para.yidc[i]];
			// 根据参数类型，进行修正
			switch (para.kind)
			{
			case 2:  // 二元相互作用参数
					 //亚点阵内部二元相互作用，求(x1-x2)
				for (j = 0; j < Phases.GPn; j++)
				{
					xdiff[j] = Phases.GPy[j][para.idsub2[0]] - Phases.GPy[j][para.idsub2[1]];
					xpro[j] = xpro[j] * pow(xdiff[j], para.order);
				}
				break;
			case 3:  // 三元相互作用参数
				xpro = 0;

					 //亚点阵内部三元相互作用，求
				break;
			case 4:  // 三元相互作用参数
				xpro = 0;

					 //亚点阵内部三元相互作用，求
				break;
			}

			// 累加到Gibbs自由能中
			paravalue = CalcPara(Tnseg, para, T);
			if (para.type == "G" || para.type == "L")
			{
				for (j = 0; j < GPn; j++)
					Phases.GPGF[j] += xpro[j] * paravalue;
			}
			else if (para.type == "TC")
			{
				for (j = 0; j < GPn; j++)
					Phases.GPTC[j] += xpro[j] * paravalue;
			}
			else if (para.type == "BMAGN")
			{
				for (j = 0; j < GPn; j++)
					Phases.GPBMAGN[j] += xpro[j] * paravalue;
			}
		}
		//存储磁性参数，用于其他性质的计算，比如一阶、二阶导数
		//磁性
		if (Phases.magp != 0)
		{
			for (j = 0; j < GPn; j++)
			{
				Phases.GPtao[j] = T / Phases.GPTC[j];
				Phases.GPGF[j] = R*T*log(Phases.GPBMAGN[j] + 1)*CalcGtao(Phases.GPtao[j], Phases.magp, 1);
			}
		}

		// 熵
		//循环亚点阵
		double natom;
		for (j = 0; j < GPn; j++)
		{
			natom = 0;
			for (k = 0; k < Phases.yn; k++)
			{
				if(Phases.GPy[j][k]<eps)
					Phases.GPSF[j] += Phases.ysp[k] * eps * log(eps);
				else
					Phases.GPSF[j] += Phases.ysp[k] * Phases.GPy[j][k] * log(Phases.GPy[j][k]);

				if (Phases.yide[k]>1)  //不包括空位
					natom += Phases.ysp[k] * Phases.GPy[j][k];
			}
			Phases.GPSF[j] = R*T*Phases.GPSF[j];
			Phases.GPGF[j] += Phases.GPSF[j];
			//计算1mol化学式中包含多少摩尔原子数，不包括空位
			Phases.GPGF[j] = Phases.GPGF[j]/natom;
		}
		delete[] Phases.GPSF;
		delete[] Phases.GPTC;
		delete[] Phases.GPBMAGN;
		delete[] Phases.GPtao;
		delete[] xpro;
		delete[] xdiff;
	}

	// 计算Gibbs自由能一阶导数
	double MinEnergy::CalcGFdx(Phase & Phases, double T, int idy)
	{
		//
		int i;
		double R = 8.31451;
		double xpro, xpro2, xdiff;
		double dGF = 0;
		double dSF = 0;
		double dTC = 0;
		double dBMAGN = 0;
		double tao;

		//
		double yfracs[MDim];
		for (i = 0; i < Phases.yn; i++)
		{
			
			yfracs[i] = Phases.y[i];

		}

		//G和L
		for (auto para : Phases.Parameters)
		{
			//判断idy是否在参数中，不在，导数为0
			int checkcon = 0;
			int naddcon = 0;
			for (i = 0; i < para.yn; i++)
			{
				if (para.yidc[i] == idy)
					checkcon = 1;
			}
			if (checkcon == 0) //idy不在参数列表中，该参数对的一阶导数为0，跳过。
				continue;

			// 判断温度区间
			int Tnseg = 0;
			for (i = 0; i < para.T.size() - 1; i++)
			{
				if ((para.T[i] - 0.02)< T && T <(para.T[i + 1] + 0.02))   // T - 0.02，防止298.14和298.15这种误差
					break;
				Tnseg++;
			}
			xpro = 1;
			xpro2 = 1;
			xdiff = 1;
			// 计算相参数前面的系数
			for (i = 0; i < para.yn; i++)
			{
				xpro = xpro*yfracs[para.yidc[i]];
			}
			xpro = xpro / yfracs[idy]; // %求导数之后为1,（yfracs[idy]/yfracs[idy] =1）
									   // 根据参数类型，进行修正
			switch (para.kind)
			{

			case 2:  // 二元相互作用参数
					 //第一部分
					 //亚点阵内部二元相互作用，求(x1-x2)
				xdiff = yfracs[para.idsub2[0]] - yfracs[para.idsub2[1]];
				xpro = xpro*pow(xdiff, para.order);
				//第二部分, idy在相互作用的亚点阵内部，Xi*Xj*Xk(Xi-Xj)^n，idy是Xi或者Xj
				if (para.order > 0)  // n=0,求导之后为0，不考虑。
				{
					if (idy == para.idsub2[0])  // (Xi-Xj)^n, idy ==Xi
					{
						for (i = 0; i < Phases.subln + 1; i++)
						{
							xpro2 = xpro2*yfracs[para.yidc[i]];
						}
						xpro2 = xpro2*para.order*pow(xdiff, para.order - 1);
					}
					else if ((idy == para.idsub2[1]))  // (Xi-Xj)^n, idy ==Xj
					{
						for (i = 0; i < Phases.subln + 1; i++)
						{
							xpro2 = xpro2*yfracs[para.yidc[i]];
						}
						xpro2 = -xpro2*para.order*pow(xdiff, para.order - 1);
					}
					xpro += xpro2;
				}
				break;
			case 3:  // 三元相互作用参数
				xpro = 0;

					 //亚点阵内部三元相互作用，求
				break;
			case 4:
				xpro = 0;
				break;
			}

			// 累加到Gibbs自由能中

			if (para.type == "G" || para.type == "L")
				dGF += xpro*CalcPara(Tnseg, para, T);
			else if (para.type == "TC")
				dTC += xpro*CalcPara(Tnseg, para, T);
			else if (para.type == "BMAGN")
				dBMAGN += xpro*CalcPara(Tnseg, para, T);
		}
		//磁性
		if (Phases.magp != 0)
		{
			tao = T / Phases.TC;
			dGF += R*T*(CalcGtao(tao, Phases.magp, 0) / (Phases.BMAGN + 1)*dBMAGN -
				T / Phases.TC / Phases.TC*log(Phases.BMAGN + 1)*CalcGtao(tao, Phases.magp, 1)*dTC);
		}

		// 熵
		
			dSF += Phases.ysp[idy] * (log(yfracs[idy]) + 1);
		dSF = R*T*dSF;
		dGF += dSF;
		//计算摩尔原子的G
		double natom;
		natom = 0;
		for (i = 0; i < Phases.yn; i++)
		{
			if (Phases.yide[i]>1)  //不包括空位
				natom += Phases.ysp[i] * Phases.y[i];
		}
		//dGF = dGF / natom;

		return dGF;
	}
	// 计算Gibbs自由能二阶导数
	double MinEnergy::CalcGFdx2(Phase & Phases, double T, int idy1, int idy2)
	{
		//
		int i;
		double R = 8.31451;
		double xpro, xpro2, xpro3, xpro4, xdiff;
		int signx1, signx2;
		double dGF2 = 0;
		double dSF2 = 0;
		int checkcon1, checkcon2;

		//
		double yfracs[MDim];
		for (i = 0; i < Phases.yn; i++)
			yfracs[i] = Phases.y[i];

		//G和L
		for (auto para : Phases.Parameters)
		{
			//判断idy是否在参数中，不在，导数为0
			checkcon1 = -1;
			checkcon2 = -1;
			for (i = 0; i < para.yn; i++)
			{
				if (para.yidc[i] == idy1)
					checkcon1 = 1;
				if (para.yidc[i] == idy2)
					checkcon2 = 1;
			}
			if (checkcon1 == -1 || checkcon2 == -1) //idy1和idy2中只要有一个不在列表中，该参数对的二阶导数为0，跳过。
				continue;
			// 判断温度区间
			int Tnseg = 0;
			for (i = 0; i < para.T.size() - 1; i++)
			{
				if ((para.T[i] - 0.02)< T && T <(para.T[i + 1] + 0.02))   // T - 0.02，防止298.14和298.15这种误差
					break;
				Tnseg++;
			}
			// 成分系数
			xpro = 1;
			xpro2 = 0;
			xpro3 = 0;
			xpro4 = 0;
			// 判断二阶导数类型
			if (idy1 == idy2)
			{
				switch (para.kind)
				{
				case 1:
					continue;  // end member， 二阶导数为0，跳过这个参数
					break;
				case 2:
					if (para.order == 0)
						continue;  //0阶参数, 为0，跳过这个参数
								   //求导后的符号
					signx1 = 0;
					if (idy1 == para.idsub2[0])
						signx1 = 1;
					else
						signx1 = -1;
					if (signx1 == 0)
						continue; //不在相互作用的亚点阵中, 为0，跳过这个参数
								  //第一部分
					for (i = 0; i < para.yn; i++)
					{
						xpro = xpro*yfracs[para.yidc[i]];
					}
					xdiff = yfracs[para.idsub2[0]] - yfracs[para.idsub2[1]];  //亚点阵内部二元相互作用，求(x1-x2)
					if (para.order > 1)
						xpro2 = xpro* para.order*(para.order - 1)*pow(xdiff, para.order - 2);  // yi*yj*yk*n*(n-1)*(yi-yk)^(n-2)
					xpro = signx1 * 2 * xpro / yfracs[idy1] * para.order*pow(xdiff, para.order - 1); // 2*yj*yk*n*(yi-yk)^(n-1)
					xpro = xpro + xpro2;
					break;
				case 3:
					xpro = 0;
					break;
				case 4:
					xpro = 0;
					break;
				}
			}
			else    // idy1 != idy2
			{
				switch (para.kind)
				{
				case 1:   // end member  yk*G(i:j:k)
					for (i = 0; i < para.yn; i++)
					{
						xpro = xpro*yfracs[para.yidc[i]];
					}
					xpro = xpro / yfracs[idy1] / yfracs[idy2];; // %求导数之后为1
					break;
				case 2:  // 二元相互作用参数
						 //第一部分
					for (i = 0; i < para.yn; i++)
					{
						xpro = xpro*yfracs[para.yidc[i]];
					}
					xdiff = yfracs[para.idsub2[0]] - yfracs[para.idsub2[1]];
					//第二部分
					//求导后的符号
					signx1 = 0;
					signx2 = 0;
					if (idy1 == para.idsub2[0])
						signx1 = 1;
					else if (idy1 == para.idsub2[1])
						signx1 = -1;
					if (idy2 == para.idsub2[0])
						signx2 = 1;
					else if (idy2 == para.idsub2[1])
						signx2 = -1;
					if (signx1 != 0 && signx2 != 0) // 第一种情况：都在相互作用亚点阵内部
					{
						if (para.order > 0)
						{
							xpro2 = signx1*xpro / yfracs[idy2] * para.order*pow(xdiff, para.order - 1); // yi*yk*n*(yi-yj)^(n-1)
							xpro3 = signx2*xpro / yfracs[idy1] * para.order*pow(xdiff, para.order - 1); // -yj*yk*n*(yi-yj)^(n-1)
							if (para.order > 1)
								xpro4 = -xpro*para.order*(para.order - 1)*pow(xdiff, para.order - 2); // -yi*yj*yk*n*(n-1)*(yi-yj)^(n-2)
							xpro = xpro / yfracs[idy1] / yfracs[idy2] * pow(xdiff, para.order);  // yk*n*(yi-yj)^n
						}
						else
							xpro = xpro / yfracs[idy1] / yfracs[idy2] * pow(xdiff, para.order);  // yk*n*(yi-yj)^n
						xpro = xpro + xpro2 + xpro3 + xpro4;
					}
					else if (signx1 != 0 && signx2 == 0) // 第二种情况：只有一个在相互作用的亚点阵内部
					{
						if (para.order > 0)
							xpro2 = signx1*xpro / yfracs[idy2] * para.order*pow(xdiff, para.order - 1);
						xpro = xpro / yfracs[idy1] / yfracs[idy2] * pow(xdiff, para.order);
						xpro = xpro + xpro2;
					}
					else if (signx1 == 0 && signx2 != 0) // 第二种情况：只有一个在相互作用的亚点阵内部
					{
						if (para.order > 0)
							xpro2 = signx2*xpro*yfracs[idy1] * para.order*pow(xdiff, para.order - 1);
						xpro = xpro / yfracs[idy1] / yfracs[idy2] * pow(xdiff, para.order);
						xpro = xpro + xpro2;
					}
					else if (signx1 != 0 && signx2 != 0) // 第四种情况：都不在相互作用的亚点阵内部
					{
						xpro = xpro / yfracs[idy1] / yfracs[idy2] * pow(xdiff, para.order);
					}
					break;
				case 3:
					for (i = 0; i < para.yn; i++)
					{
						xpro = xpro*yfracs[para.yidc[i]];
					}
					xpro = xpro / yfracs[idy1] / yfracs[idy2];; // %求导数之后为1
					break;
				case 4:
					for (i = 0; i < para.yn; i++)
					{
						xpro = xpro*yfracs[para.yidc[i]];
					}
					xpro = xpro / yfracs[idy1] / yfracs[idy2];; // %求导数之后为1
					break;
				}

			}
			// 累加到Gibbs自由能中
			dGF2 += xpro*CalcPara(Tnseg, para, T);
		}
		// 熵

		if (idy1 == idy2)
		{

			dSF2 += Phases.ysp[idy1] / yfracs[idy1];
		}
		dSF2 = R*T*dSF2;
		dGF2 += dSF2;
		//计算摩尔原子的G
		double natom;
		natom = 0;
		for (i = 0; i < Phases.yn; i++)
		{
			if (Phases.yide[i]>1)  //不包括空位
				natom += Phases.ysp[i] * Phases.y[i];
		}
		//dGF2 = dGF2 / natom;
		//if (dGF2 > 1e9)
		//dGF2 = 1e9;
		return dGF2;
	}
	//
	//
	double MinEnergy::CalcPara(int Tnseg, Parameter para, double T)
	{
		int k = 0;
		double paravalue = 0;
		for (auto coffT : para.express_digit[Tnseg].coffT)
		{
			paravalue += coffT*pow(T, para.express_digit[Tnseg].powerT[k]);
			k++;
		}
		// 计算T*LN(T)
		for (auto coffTLNT : para.express_digit[Tnseg].coffTLNT)
		{
			paravalue += coffTLNT*T*log(T);
		}
		return paravalue;
	}
	//
	void MinEnergy::CalcPara(double pvalue[3], int Tnseg, Parameter para, double T)
	{
		int k = 0;
		double temp0, temp1;
		pvalue[0] = 0;
		pvalue[1] = 0;
		pvalue[2] = 0;
		for (auto coffT : para.express_digit[Tnseg].coffT)
		{
			temp0 = coffT*pow(T, para.express_digit[Tnseg].powerT[k]);
			pvalue[0] += temp0;
			if (para.express_digit[Tnseg].powerT[k] != 0)
			{
				temp1 = para.express_digit[Tnseg].powerT[k] * temp0 / T;
				pvalue[1] += temp1;
				if (para.express_digit[Tnseg].powerT[k] != 1)
				{
					pvalue[2] += (para.express_digit[Tnseg].powerT[k] - 1) *temp1 / T;
				}
			}
			k++;
		}
		// 计算T*LN(T)
		for (auto coffTLNT : para.express_digit[Tnseg].coffTLNT)
		{
			pvalue[0] += coffTLNT*T*log(T);
			pvalue[1] += coffTLNT*(log(T) + 1);
			pvalue[2] += coffTLNT / T;
		}
	}
	double MinEnergy::CalcGtao(double tao, double p, int order)
	{
		double  Gtao;
		double A = 0.460444444444444 + 0.731893583724570*(1 / p - 1);
		if (order == 0)
		{
			if (tao < 1)
				Gtao = 1 - 1 / A*(0.564285714285714 / tao / p + 0.953722334004024*
					(1 / p - 1)*(pow(tao, 3) / 6 + pow(tao, 9) / 135 + pow(tao, 15) / 600));
			else
				Gtao = -1 / A*(pow(tao, -5) / 10 + pow(tao, -15) / 315 + pow(tao, -25) / 1500);
		}
		else if (order == 1)
		{
			if (tao < 1)
				Gtao = 1 - 1 / A*(0.564285714285714 / tao / p + 0.953722334004024*
					(1 / p - 1)*(pow(tao, 3) / 6 + pow(tao, 9) / 135 + pow(tao, 15) / 600));
			else
				Gtao = -1 / A*(pow(tao, -5) / 10 + pow(tao, -15) / 315 + pow(tao, -25) / 1500);
		}
		else if (order == 2)
		{
			if (tao < 1)
				Gtao = -1 / A*(0.564285714285714 / tao / tao / p + 0.953722334004024*
					(1 / p - 1)*(pow(tao, 2) / 3 + pow(tao, 8) / 15 + pow(tao, 14) / 40));
			else
				Gtao = 1 / A*(pow(tao, -6) / 2 + pow(tao, -16) / 21 + pow(tao, -26) / 60);
		}
		return Gtao;
	}

>>>>>>> b34a88e731ed6049f6c364b37b469ae02b913009
}// end of namespace