#include "..\include\CalcEquilibrium.h"

namespace VCLab
{
	using namespace std;

	void CalcEquilibrium::CalcGPn()
	{
		int i, j, k;
		int nperm, yn;
		double sublsum;
		double dy;

		GPn = 0;
		for (i = 0; i < Phn; i++)
		{
			CEPhases[i].GPn = 1;
			sublsum = CEPhases[i].sublpersum;
			for (j = 0; j < CEPhases[i].subln; j++)
			{
				if (CEPhases[i].conn[j]>1)
				{
					//dy = dx/(a/sum(a)) in sublattice 
					dy = GPdx / CEPhases[i].sublper[j]*sublsum;
					yn = int(1 / dy);
					CEPhases[i].GPsdy[j] = 1.0 / yn;
					CEPhases[i].GPsndy[j] = yn;
					//combination, put n ball into m box£¬box allow to be empty£¬total: C[(n+m-1),(m-1)]
					yn += CEPhases[i].conn[j] - 1;
					nperm = 1;
					for (k = 0; k < (CEPhases[i].conn[j] - 1); k++)
						nperm = nperm*(yn - k) / (k + 1);
					CEPhases[i].GPsn[j] = nperm;
					//remove identical GP, due to symmetry
					//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>on going work
					CEPhases[i].GPn = CEPhases[i].GPn*nperm;
				}
				else
				{
					CEPhases[i].GPsndy[j] = 1;
					CEPhases[i].GPsn[j] = 1;
				}
			}
			GPn += CEPhases[i].GPn;
		}
	}
	// generate grid point, permutation of x, yfractions and energy
	void CalcEquilibrium::GeneGP()
	{
		int i, j, k, ki, nc;
		int ystart;
		int combn, combm;
		double tempx;
		
		GP = new Point[GPn];
		nc = 0;
		for (i = 0; i < Phn; i++)
		{
			for (j = 0; j < CEPhases[i].subln; j++) //loop sublattices
			{
				//generate combinations in this sublattice=========================
				combm = CEPhases[i].conn[j] - 1;
				combn = CEPhases[i].GPsndy[j] + combm;
				//start = clock();
				comb[j].combinations(1, CEPhases[i].GPsdy[j], combn, combm);
				//finish = clock();
				//totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
				//cout << "\n combinations running time " << totaltime << " s!" << endl;
				//cout << "combn, combm" << combn << "  " << combm << endl;
			}
			//cartesian production of sublattices to get all permuntations==========
			input = new int*[CEPhases[i].subln];
			for (j = 0; j < CEPhases[i].subln; j++)
			{
				input[j] = new int[CEPhases[i].GPsn[j]];
				for (k = 0; k < CEPhases[i].GPsn[j]; k++)
					input[j][k] = k;
			}
			output = new int*[CEPhases[i].GPn];
			for (j = 0; j < CEPhases[i].GPn; j++)
				output[j] = new int[CEPhases[i].subln];
			//cartesian production
			//start = clock();
			cart_product(output, CEPhases[i].GPn, input, CEPhases[i].GPsn, CEPhases[i].subln);
			//finish = clock();
			//totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
			//cout << "\n cart_product running times " << totaltime << "s!" << endl;
			//use the permutations to get GPy
			//start = clock();
			CEPhases[i].GPy = new double*[CEPhases[i].GPn];
			for (j = 0; j < CEPhases[i].GPn; j++)  
			{
				CEPhases[i].GPy[j] = new double[CEPhases[i].yn];
				ystart = 0;
				for (k = 0; k < CEPhases[i].subln; k++)
				{
					for (ki = 0; ki < CEPhases[i].conn[k]; ki++)
					{
						CEPhases[i].GPy[j][ki + ystart] = comb[k].comb[output[j][k]][ki];
					}
					ystart += CEPhases[i].conn[k];
				}
			}
			//finish = clock();
			//totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
			//cout << "\n running time " << totaltime << "s!" << endl;
			//delete
			for (j = 0; j < CEPhases[i].subln; j++) 
			{
				comb[j].deletecomb();
			}
			// Calculate Gibbs energy, per mole atom======================================
			//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> on going work: parallel 
			MinE.CalcGPGF(CEPhases[i], CEConditions.T);
			//finish = clock();
			//totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
			//cout << "\nCalcGPGF running time " << totaltime << "s!" << endl;

			for (j = 0; j < CEPhases[i].GPn; j++)
			{
				for (k = 0; k < nele; k++)
					GP[nc].x[k] = 0;  
				// Calculate the compositions===========================================
				for (k = 0; k < CEPhases[i].yn; k++)
				{
					GP[nc].y[k] = CEPhases[i].GPy[j][k]; 
					GP[nc].x[CEPhases[i].yide[k] - 2] += CEPhases[i].ysp[k] * CEPhases[i].GPy[j][k]; 
				}
				// calc the overall compositions
				tempx = 0;
				for (k = 0; k < nele; k++)
				{
					tempx += GP[nc].x[k];
				}
				for (k = 0; k < nele; k++)
				{
					GP[nc].x[k] = GP[nc].x[k]/tempx;
				}
				tempx = GP[nc].x[nele - 1];
				GP[nc].x[nele- 1] = CEPhases[i].GPGF[j];
				GP[nc].x[nele] = tempx;
				GP[nc].phid = i;
				GP[nc].fb= 1;
				nc++;
			}

			// Print GP in screen
			if (debug >= 2)
			{
				cout << "Grid points in yfracs for  " << CEPhases[i].name << ",    total: " << CEPhases[i].GPn << endl;
				for (j = 0; j < CEPhases[i].GPn; j++)
				{
					for (k = 0; k < CEPhases[i].yn; k++)
					{
						cout << CEPhases[i].GPy[j][k] << "  ";
					}
					cout << endl;
				}
			}
			for (j = 0; j < CEPhases[i].GPn; j++)
			{
				delete[] CEPhases[i].GPy[j];
				delete[] output[j];
			}
			delete[] CEPhases[i].GPy;
			delete[] CEPhases[i].GPGF;
			delete[] output;
			for (j = 0; j < CEPhases[i].subln; j++)
			{
				delete[] input[j];
			}
			delete[] input;
		}
		//
		if (debug >= 2)
		{
			cout << "Grid points, total: " <<GPn << endl;
			for (j = 0; j < GPn; j++)
			{
				cout <<setfill(' ') << "Phases: " << CEPhases[GP[j].phid].name <<setw(10);
				for (k = 0; k < nele + 1; k++)
				{
					cout  << GP[j].x[k] << "  ";
				}
				cout << endl;
			}
		}
	}
	//fine the convex hull
	void CalcEquilibrium::FindConvexFace()
	{
		string OutFileName = "output.txt";
		Point PA;
		int i, Dim;
		Dim = CEConditions.nele - 2;
		for (i = 0; i < Dim; i++)
			PA.x[i] = CEConditions.x[i + 2];
		CHull.FindConvexFace(Dim, GPn, GP, PA);  // 
		if(debug >= 1)
			CHull.ShowConvexFace(OutFileName);  //
		GetConvexFace(Dim, CHull.f, PA); // Calc the chemical potetnail and phase fractions
		
		MergeConvexFace(Dim, CHull.f);	// merge GP
	}
	
	// Calc the chemical potetnail and phase fractions
	void CalcEquilibrium::GetConvexFace(int Dim, Face &f, Point PA)
	{
		int i, j, k;
		double eF[MDim3], eFJ[MDim3][MDim3];
		double temp;
		Point aveGP;
		// building the martic
		for (i = 0; i < Dim; i++)
		{
			for (j = 0; j < Dim - 1; j++)
				eFJ[i][j] = f.p[i].x[j];
			eFJ[i][Dim - 1] = f.p[i].x[Dim];
			eF[i] = f.p[i].x[Dim - 1];
		}
		//gauss
		gauss_elimination(Dim, eFJ, eF);
		for (i = 0; i < Dim; i++)
			Chp[i + 2] = -eF[i];
		//Phases fractions
		for (i = 0; i < Dim - 1; i++)
		{
			for (j = 0; j < Dim; j++)
				eFJ[i][j] = f.p[j].x[i];
			eFJ[Dim - 1][i] = f.p[i].x[Dim];
			eF[i] = PA.x[i];
		}
		eFJ[Dim - 1][Dim - 1] = f.p[Dim - 1].x[Dim];
		eF[Dim - 1] = PA.x[Dim - 1];
		gauss_elimination(Dim, eFJ, eF);
		for (i = 0; i < Dim; i++)
			f.p[i].phfrac = eF[i];
		if (debug >= 1)
		{
			cout << "Chemical potential: " << endl;
			for (i = 0; i < Dim; i++)
				cout << Chp[i + 2] << "  ";
			cout << endl;
			cout << "Phase fraction: " << endl;
			for (i = 0; i < Dim; i++)
			{
				cout << CEPhases[f.p[i].phid].name << "  " << f.p[i].phfrac << endl;
				for (j = 0;j < CEPhases[f.p[i].phid].vyn;j++)
				{
					cout << f.p[i].y[j] << "   ";
				}
				cout << endl;
			}
			cout << endl;
		}
		
	}

	void CalcEquilibrium::MergeConvexFace(int Dim, Face &f)
	{
		int i, j, k;
		int fphn[MDim];
		int fphid[MDim];
		Point aveGP;
		
		for (i = 0; i < Dim; i++)
		{
			for (j = i + 1; j < Dim; j++)
			{
				if (f.p[i].phid == f.p[j].phid)
				{
					aveGP = f.p[j];
					CEPhases[f.p[j].phid].GPn = 1;
					CEPhases[f.p[j].phid].GPy = new double*[1];
					CEPhases[f.p[j].phid].GPy[0] = new double[CEPhases[f.p[j].phid].yn];
					//calc the average compositions
					for (k = 0; k < CEPhases[f.p[j].phid].yn; k++)
					{
						aveGP.y[k] = (f.p[i].phfrac*f.p[i].y[k] + f.p[j].phfrac*f.p[j].y[k]) /
							(f.p[i].phfrac + f.p[j].phfrac);
						CEPhases[f.p[j].phid].GPy[0][k] = aveGP.y[k];
					}
					MinE.CalcGPGF(CEPhases[f.p[j].phid], CEConditions.T);  
					//if ture, or miscibility gap
					if (CEPhases[f.p[j].phid].GPGF[0]<=(f.p[i].phfrac*f.p[i].x[Dim-1]+ 
						f.p[j].phfrac*f.p[j].x[Dim - 1])/(f.p[i].phfrac + f.p[j].phfrac))
					{
						for (k = 0; k < Dim + 1; k++)
							aveGP.x[k] = (f.p[i].phfrac*f.p[i].x[k] + f.p[j].phfrac*f.p[j].x[k]) 
							/ (f.p[i].phfrac + f.p[j].phfrac);
						f.p[i].fb = 0;
						f.p[j] = aveGP;
						f.p[j].phfrac += f.p[i].phfrac;
						break;
					}
				}
			}
		}
	}

} // end of VCLab