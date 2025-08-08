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

constexpr auto SPACE = 2;						//Ԫ������

constexpr auto ROWS = 201;						//Ԫ��������
constexpr auto COLS = 201;
constexpr auto num_ite = 15001;					//CA�����ܴ���
constexpr auto maxTstep = 12000;				//�¶ȳ�������ʱ�䲽
constexpr auto maxClstep = 10000;				//Һ�����ʳ�������ʱ�䲽
constexpr auto Neb = 4;							//�ھ���Ŀ

constexpr auto len = 1e-6;						//��Ԫ�񳤶� 1��m
constexpr auto dy = 1e-6;						//��Ԫ��y���򳤶�
constexpr auto dx = 1e-6;						//��Ԫ��x���򳤶�

constexpr auto iang= M_PI / 6;					//��ʼȡ���

constexpr auto T0 = 1809.15;					//Һ�����¶�
constexpr auto C0 = 0.0082;						//��ʼ���ʷ���
constexpr auto ml = -7800;
//Һ����б��

constexpr auto k = 0.34;						//���ʷ���ϵ��

constexpr auto dT = 15;							//ֱ�Ӹ����ĳ��������
constexpr auto dTn = 15;						//�κ˹���ȵ�ƽ��ֵ
constexpr auto dT�� = 1.5;						//�κ˹���ȵı�׼ƫ��
constexpr auto nmax = 2.224e6;					//����κ��ܶ�
constexpr auto dTdt = 0;						//��������

constexpr auto dt = 0.00001;					//ʱ�䲽��

constexpr auto GT = 1.9e-7;						//Gibbs�CThomsonϵ��

constexpr auto angd = 0.04;						//�Ƕȵĸ�������ϵ��

constexpr auto Conductivity = 33;				//�����ȵ���
constexpr auto Density_s = 7400;				//�����ܶ�		
constexpr auto Density_l = 7020;				//Һ���ܶ�
constexpr auto Shc_s = 648;						//�������
constexpr auto Shc_l = 824;						//Һ�����
constexpr auto Shc_m = 770;
constexpr auto Latent = 2.72e5;					//�ۻ�Ǳ��

//�����ṹ�� Ԫ��
struct Cell
{
	//����Ԫ��״̬����
	int sta = 0, TL = 0;
	double T = 0.0, Tb = 0.0,
		x = 0.0, y = 0.0,													//Ԫ���¶� ����� Ԫ������
		MDCS_x = 0.0, MDCS_y = 0.0, L_dia = 0.0,							//ƫ������������ ��Խ��߳���
		P1_x = 0.0, P1_y = 0.0, P2_x = 0.0, P2_y = 0.0,
		P3_x = 0.0, P3_y = 0.0, P4_x = 0.0, P4_y = 0.0,						//ƫ�������ζ������꣨�����޷֣�
		vel = 0.0, vel_x = 0.0, vel_y = 0.0,								//Ԫ�������ٶȣ�x�����ٶȣ�y�����ٶ�
		Cl = 0.0, Cs = 0.0, Cl8 = 0.0, Cs8 = 0.0, Ce = 0.0,					//Ԫ��Һ������Ũ�ȣ���������Ũ�ȣ���Һ����ǰ��Һ��Ũ�ȣ���Һ����ǰ�ع���Ũ��
		fss = 0.0, fl = 0.0, fs = 0.0, delta_fs = 0.0,						//����Ԫ��Һ��������
		cur = 0.0, ani = 0.0, ang = 0.0, angn = 0.0, test = 0.0,			//Ԫ�����ʣ���������ϵ����ȡ���
		nx = 0.0, ny = 0.0,													//���淨�����
		gf = 0.0,															//��״����
		nuc_p = 0.0, P = 0.0,												//�Ǿ��������κ˸���
		//�¶ȳ�Ũ�ȳ����Ҫ��
		Tae1 = 0.0, Taw1 = 0.0, Tan1 = 0.0, Tas1 = 0.0, Tap1 = 0.0, Tap0 = 0.0, Tbp = 0.0,
		Tn0 = 0.0, Tresi = 0.0,
		Clae0 = 0.0, Claw0 = 0.0, Clan0 = 0.0, Clas0 = 0.0,
		Clae1 = 0.0, Claw1 = 0.0, Clan1 = 0.0, Clas1 = 0.0, Clap1 = 0.0, Clap0 = 0.0, Clbp = 0.0,
		Cln0 = 0.0, Clresi = 0.0,
		Csae0 = 0.0, Csaw0 = 0.0, Csan0 = 0.0, Csas0 = 0.0, Csap1 = 0.0, Csap0 = 0.0, Csbp = 0.0,
		Csn0 = 0.0;

	//�Ը�ֵ������������أ��Կ��ٸ���Ԫ�����飬����ǰһ����Ԫ��״̬
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

//������ɢϵ��
double Ds(double TT);

//Һ����ɢϵ��
double Dl(double TT);
//������ͼ
void drawRect(int i, int j, int cell_color);

//��ʼ��Ԫ������
//��Ԫ�����꣨Ԫ�����ģ����и�ֵ
void initCell(Cell** cell, Cell** precell);

//��ʼ��Ԫ������
void initseed(Cell** cell, Cell** precell);

//ȷ������Ԫ����Ӧƫ�������ζ�������
void ESVCCell(Cell** cell);

//������һ�ε������Ԫ��״̬ 1
void copyCell_1(Cell** cell, Cell** precell);

//ȷ������Ԫ���ж��Ƿ�׽��ΧҺ���	Һ����������
void LtoICell(Cell** cell, Cell** precell);

//����Ԫ��״̬	Һ����������
void LtoICell_1(Cell** cell);

//����Ԫ��״̬	������������
void ItoSCell(Cell** cell);

//������Ũ��
void updateCe(Cell** cell);

//����������������ԣ�˫���Բ�ֵ����
void Curvature(Cell** cell);

//�����������������Ũ��
void Cl8Cs8(Cell** cell);

//�������������ٶ�
void vel(Cell** cell);

//���㼸����������
void gf(Cell** cell);

//����ƫ�������ΰ�Խ��߳���
void Ldia(Cell** cell, double delta_t);

//�����������
void deltafs(Cell**cell, double delta_t);

//������һ�ε������Ԫ��״̬2
void copyCell_2(Cell** cell, Cell** precell);

//����Ԫ���������
void updatefs(Cell** cell, Cell** precell);

//������һ�ε������Ԫ��״̬3
void copyCell_3(Cell** cell, Cell** precell);

//����Һ�����ʳ����ϵ��
void Cla(Cell** cell, double delta_t);

//����������ʳ����ϵ��
void Csa(Cell** cell, double delta_t);

//������߽�����
void Cl0(Cell** cell);

//�������ʳ����ϵ��
void aCl0(Cell** cell);

//Һ��������ʽ����
void Cl(Cell** cell, Cell** precell);

//����������ʽ����
void Cs(Cell** cell, Cell** precell);

//����д��
void rData(Cell** cell, string filename);

//��������
void drawCell(Cell** cell);

//�����ʱ����
double fun_delta_t(Cell** cell, double vell, double Dss, double Dll);