#include "Cell.h"

int main()
{
	//并行线程数
	omp_set_num_threads(12);
	//动态线程开关
	omp_set_dynamic(0);
	//定义时间变量 初始化
	time_t start, stop;
	start = time(NULL);

	//生成random_device对象sd做种子
	random_device sd;
	//使用种子初始化linear_congruential_engine对象
	mt19937_64 linearRan(sd());
	uniform_real_distribution<double>dis1(0.0, 1.0);
	uniform_int_distribution<int>dis2(1, 89);

	//初始化元胞的二维结构体并分配动态内存
	//初始化储存上一迭代元胞状态的二维结构体并分配动态内存
	Cell** cell, ** precell;
	cell = new Cell * [ROWS + 4];
	precell = new Cell * [ROWS + 4];
	for (int i = 0; i < ROWS + 4; i++)
	{
		cell[i] = new Cell[COLS + 4];
		precell[i] = new Cell[COLS + 4];
	}

	//初始化元胞状态
	initCell(cell, precell);

	//初始化元胞种子
	initseed(cell, precell);

	//预览图框生成
	initgraph(SPACE * ROWS, SPACE * COLS, EW_SHOWCONSOLE);
	setbkcolor(0xAA0000);
	cleardevice();

	//初始化迭代时间与时间间隔
	//初始化迭代次数变量
	double t = 0, delta_t = dt;
	int ite = 0;
	string filename = "CA";			//存储文件前缀

	double vel_max = DBL_MIN;		//获取一个浮点数的最小值，防止分母为零
	double Ds_max = DBL_MIN, Dl_max = DBL_MIN;

	//迭代开始
	while (ite < num_ite)
	{
		//确定界面元胞对应偏心正方形顶点坐标
		ESVCCell(cell);

		//保存上一次迭代后的元胞状态 1
		copyCell_1(cell, precell);

		//确定界面元胞判断是否捕捉周围液相胞	液相胞→界面胞
		LtoICell(cell, precell);

		//更新元胞状态	液相胞→界面胞
		LtoICell_1(cell);

		//更新元胞状态	界面胞→固面胞
		ItoSCell(cell);

		//更新总浓度
		updateCe(cell);

		//计算曲率与各项异性（双线性插值法）
		Curvature(cell);

		//界面胞界面两侧溶质浓度
		Cl8Cs8(cell);

		//计算界面胞生长速度
		vel(cell);
	
		//计算几何修正因子
		gf(cell);

		//计算偏心正方形半对角线长度
		Ldia(cell, delta_t);

		//计算固相增长
		deltafs(cell, delta_t);

		//保存上一次迭代后的元胞状态2
		copyCell_2(cell, precell);

		//更新元胞固相分数
		updatefs(cell, precell);

		//保存上一次迭代后的元胞状态3
		copyCell_3(cell, precell);

		//-----------------------------------------
		//温度场迭代
		
		//-----------------------------------------
		//溶质场迭代
		//计算液相溶质场求解系数
		Cla(cell, delta_t);

		//计算固相溶质场求解系数
		Csa(cell, delta_t);

		//计算域边界条件
		Cl0(cell);

		//计算溶质场求解系数
		aCl0(cell);

		//液相溶质显式迭代
		Cl(cell, precell);

		//固相溶质显式迭代
		Cs(cell, precell);

		//-----------------------------------------
		//数据导出
		//if (t >= 0.212)
		//if (ite % 1000 == 0.0)
		//if (ite > num_ite - 2)
		//if (true)
		if (false)
		{
			filename = "./data/" + to_string(ite) + "s CA.csv";
			rData(cell, filename);
		}

		//绘制晶粒	
		//if (false)
		//if (true)
		if (ite % 100 == 0.0)
		//if (ite > num_ite - 2)
		{
			drawCell(cell);
		}

		//if (ite % 100 == 0.0)
		if (false)
		{
			char buf[100] = { 0 };
			sprintf_s(buf, "./data/%d.png", ite);//字符串里换上你的文件路径
			saveimage(buf);
		}

		//求迭代时间间隔
		vel_max = DBL_MIN;
		Ds_max = DBL_MIN;
		Dl_max = DBL_MIN;
		delta_t = fun_delta_t(cell, vel_max, Ds_max, Dl_max);

		//计算总时间
		t = t + delta_t;
		ite = ite + 1;
		//if (true)
		if (ite % 100 == 0.0)

		{
			cout << "迭代时间：" << t << endl;
			cout << "时间间隔：" << delta_t << endl;
			cout << "第" << ite << "次迭代" << endl;
		}

		//cout << "第" << ite << "次迭代" << endl;
	}
	//outFile.close();

	cout << "计算时间" << t << endl;
	cout << "计算完成" << endl;
	stop = time(NULL);
	stop = stop - start;
	cout << "运行时间" << stop << endl;
	system("pause");

	//释放内存
	for (int i = 0; i < ROWS + 4; i++)
	{
		delete[] cell[i];
		delete[] precell[i];
	}
	delete[] cell;
	delete[] precell;

	return 0;
}