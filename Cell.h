#pragma once

#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <graphics.h>
#include <omp.h>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <random>

using namespace std;

constexpr auto SPACE = 2;						//元胞像素

constexpr auto ROWS = 201;						//元胞行列数
constexpr auto COLS = 201;
constexpr auto num_ite = 15001;					//CA迭代总次数
constexpr auto maxTstep = 12000;				//温度场最大迭代时间步
constexpr auto maxClstep = 10000;				//液相溶质场最大迭代时间步
constexpr auto Neb = 4;							//邻居数目

constexpr auto len = 1e-6;						//单元格长度 1μm
constexpr auto dy = 1e-6;						//单元格y方向长度
constexpr auto dx = 1e-6;						//单元格x方向长度

constexpr auto iang= M_PI / 6;					//初始取向角

constexpr auto T0 = 1809.15;					//液相线温度
constexpr auto C0 = 0.0082;						//初始溶质分数
constexpr auto ml = -7800;
//液相线斜率

constexpr auto k = 0.34;						//溶质分配系数

constexpr auto dT = 15;							//直接给定的常数过冷度
constexpr auto dTn = 15;						//形核过冷度的平均值
constexpr auto dTσ = 1.5;						//形核过冷度的标准偏差
constexpr auto nmax = 2.224e6;					//最大形核密度
constexpr auto dTdt = 0;						//降温速率

constexpr auto dt = 0.00001;					//时间步长

constexpr auto GT = 1.9e-7;						//Gibbs–Thomson系数

constexpr auto angd = 0.04;						//角度的各向异性系数

constexpr auto Conductivity = 33;				//材料热导率
constexpr auto Density_s = 7400;				//固相密度		
constexpr auto Density_l = 7020;				//液相密度
constexpr auto Shc_s = 648;						//固相比热
constexpr auto Shc_l = 824;						//液相比热
constexpr auto Shc_m = 770;
constexpr auto Latent = 2.72e5;					//熔化潜热

//申明结构体 元胞
struct Cell
{
	//申明元胞状态变量
	int sta = 0, TL = 0;
	double T = 0.0, Tb = 0.0,
		x = 0.0, y = 0.0,													//元胞温度 过冷度 元胞坐标
		MDCS_x = 0.0, MDCS_y = 0.0, L_dia = 0.0,							//偏心正方形坐标 半对角线长度
		P1_x = 0.0, P1_y = 0.0, P2_x = 0.0, P2_y = 0.0,
		P3_x = 0.0, P3_y = 0.0, P4_x = 0.0, P4_y = 0.0,						//偏心正方形顶点坐标（按象限分）
		vel = 0.0, vel_x = 0.0, vel_y = 0.0,								//元胞法向速度，x方向速度，y方向速度
		Cl = 0.0, Cs = 0.0, Cl8 = 0.0, Cs8 = 0.0, Ce = 0.0,					//元胞液相溶质浓度，固相溶质浓度，固液界面前沿液相浓度，固液界面前沿固相浓度
		fss = 0.0, fl = 0.0, fs = 0.0, delta_fs = 0.0,						//界面元胞液相固相分数
		cur = 0.0, ani = 0.0, ang = 0.0, angn = 0.0, test = 0.0,			//元胞曲率，各向异性系数，取向角
		nx = 0.0, ny = 0.0,													//界面法向分量
		gf = 0.0,															//形状因子
		nuc_p = 0.0, P = 0.0,												//非均质连续形核概率
		//温度场浓度场求解要素
		Tae1 = 0.0, Taw1 = 0.0, Tan1 = 0.0, Tas1 = 0.0, Tap1 = 0.0, Tap0 = 0.0, Tbp = 0.0,
		Tn0 = 0.0, Tresi = 0.0,
		Clae0 = 0.0, Claw0 = 0.0, Clan0 = 0.0, Clas0 = 0.0,
		Clae1 = 0.0, Claw1 = 0.0, Clan1 = 0.0, Clas1 = 0.0, Clap1 = 0.0, Clap0 = 0.0, Clbp = 0.0,
		Cln0 = 0.0, Clresi = 0.0,
		Csae0 = 0.0, Csaw0 = 0.0, Csan0 = 0.0, Csas0 = 0.0, Csap1 = 0.0, Csap0 = 0.0, Csbp = 0.0,
		Csn0 = 0.0;

	//对赋值运算符进行重载，以快速复制元胞数组，保存前一迭代元胞状态
	Cell operator=(Cell& cal)
	{
		sta = cal.sta; TL = cal.TL;
		T = cal.T; Tb = cal.Tb;
		vel = cal.vel; vel_x = cal.vel_x; vel_y = cal.vel_y;
		Cl = cal.Cl; Cs = cal.Cs; Cl8 = cal.Cl8; Cs8 = cal.Cs8;
		fss = fs; fl = cal.fl; fs = cal.fs;
		cur = cal.cur; ani = cal.ani; ang = cal.ang;
		delta_fs = cal.delta_fs;
		return*this;
	}
};

//固相扩散系数
double Ds(double TT);

//液相扩散系数
double Dl(double TT);
//晶粒绘图
void drawRect(int i, int j, int cell_color);

//初始化元胞参数
//对元胞坐标（元胞中心）进行赋值
void initCell(Cell** cell, Cell** precell);

//初始化元胞种子
void initseed(Cell** cell, Cell** precell);

//确定界面元胞对应偏心正方形顶点坐标
void ESVCCell(Cell** cell);

//保存上一次迭代后的元胞状态 1
void copyCell_1(Cell** cell, Cell** precell);

//确定界面元胞判断是否捕捉周围液相胞	液相胞→界面胞
void LtoICell(Cell** cell, Cell** precell);

//更新元胞状态	液相胞→界面胞
void LtoICell_1(Cell** cell);

//更新元胞状态	界面胞→固面胞
void ItoSCell(Cell** cell);

//更新总浓度
void updateCe(Cell** cell);

//计算曲率与各项异性（双线性插值法）
void Curvature(Cell** cell);

//界面胞界面两侧溶质浓度
void Cl8Cs8(Cell** cell);

//计算界面胞生长速度
void vel(Cell** cell);

//计算几何修正因子
void gf(Cell** cell);

//计算偏心正方形半对角线长度
void Ldia(Cell** cell, double delta_t);

//计算固相增长
void deltafs(Cell**cell, double delta_t);

//保存上一次迭代后的元胞状态2
void copyCell_2(Cell** cell, Cell** precell);

//更新元胞固相分数
void updatefs(Cell** cell, Cell** precell);

//保存上一次迭代后的元胞状态3
void copyCell_3(Cell** cell, Cell** precell);

//计算液相溶质场求解系数
void Cla(Cell** cell, double delta_t);

//计算固相溶质场求解系数
void Csa(Cell** cell, double delta_t);

//计算域边界条件
void Cl0(Cell** cell);

//计算溶质场求解系数
void aCl0(Cell** cell);

//液相溶质显式迭代
void Cl(Cell** cell, Cell** precell);

//固相溶质显式迭代
void Cs(Cell** cell, Cell** precell);

//数据写入
void rData(Cell** cell, string filename);

//晶粒绘制
void drawCell(Cell** cell);

//求迭代时间间隔
double fun_delta_t(Cell** cell, double vell, double Dss, double Dll);