<<<<<<< HEAD
/*
*   This file is a part of the VCLab software project.
*   For more details please contact xingwang@csu.edu.cn
*
*   Authors:    Xing Wang
*
*   Copyright (c) 2015-2015 Phase Diagram Center. Central South University. China
*
*   This program is free software: you can redistribute it and/or modify
*   it under the terms of the GNU General Public License as published by
*   the Free Software Foundation, either version 3 of the License, or
*   (at your option) any later version.
*
*   This program is distributed in the hope that it will be useful,
*   but WITHOUT ANY WARRANTY; without even the implied warranty of
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*   GNU General Public License for more details.
*
*/
//

#ifndef DRIVINGFORCE_H_
#define DRIVINGFORCE_H_

#include "VCLInput.h"
#include "Database.h"
#include "Tool.h"

namespace VCLab
{
	using namespace std;
	//

	class DrivingForce
	{
	private:
		string classname;
		int linenumber;
	public:
		//
		Condition DFConditions;
		vector<Species> DFSpecies;
		vector<Phase> DFPhases;

		//
		int vn = 0;  //总变量数
		double DF[MDim3];
		int npha = 2;  //相数目
		double Chp[MDim]; //化学势
		//int nsp = 4;  //组元数

					  // mode = -1 确定计算亚稳相Driving Force 的Hillert's 平衡方程组及牛顿法需要的值及导数
		int mode = 1;

		//
		~DrivingForce();

		// 计算平衡
		// 建立计算条件
		void SetupConditions(Condition CLCondition);
		// 确定Phase sets
		void SetupPhaseSets(Database CLDatabase);
		// 确定系统变量，并编号
		void SetupVariables();
		// 根据系统变量编号，赋值给每个相及化学势
		void AssignVariables();
		// 初始化，相选择，点阵分数，化学势
		void Initialization();
		// 确定Hillert's 平衡方程组,牛顿法，值及导数
		void SetupHillert();
		// 确定计算亚稳相Driving Force 的Hillert's 平衡方程组及牛顿法需要的值及导数
		void SetupDrivingForce();
		// 迭代计算
		int Iteration();
		// 求解给定的 phase set
		void Solver(Database CLDatabase, Condition CLCondition);

		// print on screen
		void ShowEquilibrium();
		void ShowCondition();

	};

} // end of VCLab

=======
/*
*   This file is a part of the VCLab software project.
*   For more details please contact xingwang@csu.edu.cn
*
*   Authors:    Xing Wang
*
*   Copyright (c) 2015-2015 Phase Diagram Center. Central South University. China
*
*   This program is free software: you can redistribute it and/or modify
*   it under the terms of the GNU General Public License as published by
*   the Free Software Foundation, either version 3 of the License, or
*   (at your option) any later version.
*
*   This program is distributed in the hope that it will be useful,
*   but WITHOUT ANY WARRANTY; without even the implied warranty of
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*   GNU General Public License for more details.
*
*/
//

#ifndef DRIVINGFORCE_H_
#define DRIVINGFORCE_H_

#include "VCLInput.h"
#include "Database.h"
#include "Tool.h"

namespace VCLab
{
	using namespace std;
	//

	class DrivingForce
	{
	private:
		string classname;
		int linenumber;
	public:
		//
		Condition DFConditions;
		vector<Species> DFSpecies;
		vector<Phase> DFPhases;

		//
		int vn = 0;  //总变量数
		double DF[MDim3];
		int npha = 2;  //相数目
		double Chp[MDim]; //化学势
		//int nsp = 4;  //组元数

					  // mode = -1 确定计算亚稳相Driving Force 的Hillert's 平衡方程组及牛顿法需要的值及导数
		int mode = 1;

		//
		~DrivingForce();

		// 计算平衡
		// 建立计算条件
		void SetupConditions(Condition CLCondition);
		// 确定Phase sets
		void SetupPhaseSets(Database CLDatabase);
		// 确定系统变量，并编号
		void SetupVariables();
		// 根据系统变量编号，赋值给每个相及化学势
		void AssignVariables();
		// 初始化，相选择，点阵分数，化学势
		void Initialization();
		// 确定Hillert's 平衡方程组,牛顿法，值及导数
		void SetupHillert();
		// 确定计算亚稳相Driving Force 的Hillert's 平衡方程组及牛顿法需要的值及导数
		void SetupDrivingForce();
		// 迭代计算
		int Iteration();
		// 求解给定的 phase set
		void Solver(Database CLDatabase, Condition CLCondition);

		// print on screen
		void ShowEquilibrium();
		void ShowCondition();

	};

} // end of VCLab

>>>>>>> b34a88e731ed6049f6c364b37b469ae02b913009
#endif