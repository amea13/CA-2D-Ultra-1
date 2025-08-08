#include "Cell.h"

//固相扩散系数
double Ds(double TT)
{
	double ds;
	ds = 7.61e-6 * exp(-16185.23 / TT);
	return ds;
}

//液相扩散系数
double Dl(double TT)
{
	double dl;
	dl = 7.67e-6 * exp(-12749.58 / TT);
	return dl;
}

//晶粒绘图
void drawRect(int i, int j, int cell_color)
{
	setfillcolor(cell_color);
	solidrectangle(j * SPACE, i * SPACE, j * SPACE + SPACE, i * SPACE + SPACE);
}

//初始化元胞参数
//对元胞坐标（元胞中心）进行赋值
void initCell(Cell** cell, Cell** precell)
{
#pragma omp parallel for
	for (int i = 0; i <= ROWS + 3; i++)
	{
		for (int j = 0; j <= COLS + 3; j++)
		{
			//坐标值初始化
			cell[i][j].x = j * len + len / 2.0;
			cell[i][j].y = i * len + len / 2.0;
			//液相初始化 -1液相胞 0界面胞 1固相胞
			cell[i][j].sta = -1;
			cell[i][j].Cl = C0;
			cell[i][j].Cs = 0.0;
			cell[i][j].Cl8 = 0.0;
			cell[i][j].Cs8 = 0.0;
			cell[i][j].fl = 1.0;
			//温度场初始化（待补充函数）
			cell[i][j].T = T0 - dT;
			precell[i][j].T = T0 - dT;
		}
	}
}

//初始化元胞种子
void initseed(Cell** cell, Cell** precell)
{
	precell[(ROWS + 1) / 2 + 1][(COLS + 1) / 2 + 1].fs = 1.0;
	cell[(ROWS + 1) / 2 + 1][(COLS + 1) / 2 + 1].fs = 1.0;
	cell[(ROWS + 1) / 2 + 1][(COLS + 1) / 2 + 1].Cs = C0 * k;
	cell[(ROWS + 1) / 2 + 1][(COLS + 1) / 2 + 1].Cl = 0;
	cell[(ROWS + 1) / 2 + 1][(COLS + 1) / 2 + 1].MDCS_x = cell[(ROWS + 1) / 2 + 1][(COLS + 1) / 2 + 1].x;
	cell[(ROWS + 1) / 2 + 1][(COLS + 1) / 2 + 1].MDCS_y = cell[(ROWS + 1) / 2 + 1][(COLS + 1) / 2 + 1].y;
	cell[(ROWS + 1) / 2 + 1][(COLS + 1) / 2 + 1].ang = iang;
	//cell[(ROWS + 1) / 2+1][(COLS + 1) / 2+1].ang = 0.0;
	cell[(ROWS + 1) / 2 + 1][(COLS + 1) / 2 + 1].L_dia = len / 2 / cos(cell[(ROWS + 1) / 2 + 1][(COLS + 1) / 2 + 1].ang) + 1e-10;
	//cell[(ROWS + 1) / 2 + 1][(COLS + 1) / 2 + 1].L_dia = len * sqrt(2) / 2 + 2e-10;
	//cell[(ROWS + 1) / 2 + 1][(COLS + 1) / 2 + 1].L_dia = len;
	cell[(ROWS + 1) / 2 + 1][(COLS + 1) / 2 + 1].sta = 0;
}

//确定界面元胞对应偏心正方形顶点坐标
void ESVCCell(Cell** cell)
{
#pragma omp parallel for
	for (int i = 2; i <= ROWS + 1; i++)
	{
		for (int j = 2; j <= COLS + 1; j++)
		{
			if (cell[i][j].sta == 0 && cell[i][j].TL == 0)
			{
				//jp方向
				cell[i][j].P1_x = cell[i][j].MDCS_x + cell[i][j].L_dia * cos(cell[i][j].ang);
				cell[i][j].P1_y = cell[i][j].MDCS_y - cell[i][j].L_dia * sin(cell[i][j].ang);
				//im方向
				cell[i][j].P2_x = cell[i][j].MDCS_x - cell[i][j].L_dia * sin(cell[i][j].ang);
				cell[i][j].P2_y = cell[i][j].MDCS_y - cell[i][j].L_dia * cos(cell[i][j].ang);
				//jm
				cell[i][j].P3_x = cell[i][j].MDCS_x - cell[i][j].L_dia * cos(cell[i][j].ang);
				cell[i][j].P3_y = cell[i][j].MDCS_y + cell[i][j].L_dia * sin(cell[i][j].ang);
				//ip
				cell[i][j].P4_x = cell[i][j].MDCS_x + cell[i][j].L_dia * sin(cell[i][j].ang);
				cell[i][j].P4_y = cell[i][j].MDCS_y + cell[i][j].L_dia * cos(cell[i][j].ang);

			}
		}
	}
}

//保存上一次迭代后的元胞状态 1
void copyCell_1(Cell** cell, Cell** precell)
{
#pragma omp parallel for
	for (int i = 0; i <= ROWS + 3; i++)
	{
		for (int j = 0; j <= COLS + 3; j++)
		{
			precell[i][j].sta = cell[i][j].sta;
			precell[i][j].TL = cell[i][j].TL;
		}
	}
}

//确定界面元胞判断是否捕捉周围液相胞	液相胞→界面胞
void LtoICell(Cell** cell, Cell** precell)
{
#pragma omp parallel for
	for (int i = 3; i <= ROWS; i++)
	{
		for (int j = 3; j <= COLS; j++)
		{
			if (precell[i][j].sta == 0 && precell[i][j].TL == 0)											//TL为1是带锁的 界面胞转化为固相胞后 周围若还有液相胞 则该液相胞转变为带锁界面胞 带锁界面胞不能正常参与感染 需要被感染后才能参与
			{
				int im = i - 1, ip = i + 1, jm = j - 1, jp = j + 1;
				if ((abs(cell[i][j].P1_x - cell[i][j].x) >= (len / 2)										//判断MDCS顶点的落点是否在MDCS形心之外的元胞上
					|| abs(cell[i][j].P1_y - cell[i][j].y) >= (len / 2)))
				{
					if (precell[i][jp].sta == -1 || precell[i][jp].TL == 1)									//液相胞与带锁界面胞可以被感染		
					{
						if (abs(cell[i][j].P1_x - cell[i][jp].x) <= (len / 2)								//确定该MDCS顶点感染的位置
							&& abs(cell[i][j].P1_y - cell[i][jp].y) <= (len / 2))
						{
							if (cell[i + 1][jp].fs == 1.0 || cell[i - 1][jp].fs == 1.0 ||
								cell[i][jp + 1].fs == 1.0 || cell[i][jp - 1].fs == 1.0 ||
								cell[i + 1][jp + 1].fs == 1.0 || cell[i + 1][jp - 1].fs == 1.0 ||
								cell[i - 1][jp + 1].fs == 1.0 || cell[i - 1][jp - 1].fs == 1.0)				//周围无固相胞者不可被感染 防止出现多层界面
							{
								if (cell[i][jp].sta == 0 && cell[i][jp].TL == 0)								//若此轮已被感染 则平均所得形心坐标
								{
									cell[i][jp].MDCS_x = (cell[i][j].P1_x + cell[i][jp].MDCS_x) / 2;
									cell[i][jp].MDCS_y = (cell[i][j].P1_y + cell[i][jp].MDCS_y) / 2;
								}
								else
								{
									cell[i][jp].MDCS_x = cell[i][j].P1_x;
									cell[i][jp].MDCS_y = cell[i][j].P1_y;
								}
								cell[i][jp].sta = 0;
								cell[i][jp].fs = 0.00001;
								cell[i][jp].fl = 1 - cell[i][jp].fs;
								//cell[i][jp].MDCS_x = cell[i][j].P1_x;
								//cell[i][jp].MDCS_y = cell[i][j].P1_y;
								cell[i][jp].ang = cell[i][j].ang;
								cell[i][jp].Cs = cell[i][j].Cs;
								cell[i][jp].TL = 0;
							}
						}
					}
					if (precell[im][j].sta == -1 || precell[im][j].TL == 1)
					{
						if (abs(cell[i][j].P1_y - cell[im][j].y) <= (len / 2)
							&& abs(cell[i][j].P1_x - cell[im][j].x) <= (len / 2))
						{
							if (cell[im + 1][j].fs == 1.0 || cell[im - 1][j].fs == 1.0 ||
								cell[im][j + 1].fs == 1.0 || cell[im][j - 1].fs == 1.0 ||
								cell[im + 1][j + 1].fs == 1.0 || cell[im + 1][j - 1].fs == 1.0 ||
								cell[im - 1][j + 1].fs == 1.0 || cell[im - 1][j - 1].fs == 1.0)
							{
								if (cell[im][j].sta == 0 && cell[im][j].TL == 0)
								{
									cell[im][j].MDCS_x = (cell[i][j].P1_x + cell[im][j].MDCS_x) / 2;
									cell[im][j].MDCS_y = (cell[i][j].P1_y + cell[im][j].MDCS_y) / 2;
								}
								else
								{
									cell[im][j].MDCS_x = cell[i][j].P1_x;
									cell[im][j].MDCS_y = cell[i][j].P1_y;
								}
								cell[im][j].sta = 0;
								cell[im][j].fs = 0.00001;
								cell[im][j].fl = 1 - cell[im][j].fs;
								//cell[im][j].MDCS_x = cell[i][j].P1_x;
								//cell[im][j].MDCS_y = cell[i][j].P1_y;
								cell[im][j].ang = cell[i][j].ang;
								cell[im][j].Cs = cell[i][j].Cs;
								cell[im][j].TL = 0;
							}
						}
					}
					if (precell[i][jm].sta == -1 || precell[i][jm].TL == 1)
					{
						if (abs(cell[i][j].P1_x - cell[i][jm].x) <= (len / 2)
							&& abs(cell[i][j].P1_y - cell[i][jm].y) <= (len / 2))
						{
							if (cell[i + 1][jm].fs == 1.0 || cell[i - 1][jm].fs == 1.0 ||
								cell[i][jm + 1].fs == 1.0 || cell[i][jm - 1].fs == 1.0 ||
								cell[i + 1][jm + 1].fs == 1.0 || cell[i + 1][jm - 1].fs == 1.0 ||
								cell[i - 1][jm + 1].fs == 1.0 || cell[i - 1][jm - 1].fs == 1.0)
							{
								if (cell[i][jm].sta == 0 && cell[i][jm].TL == 0)
								{
									cell[i][jm].MDCS_x = (cell[i][j].P1_x + cell[i][jm].MDCS_x) / 2;
									cell[i][jm].MDCS_y = (cell[i][j].P1_y + cell[i][jm].MDCS_y) / 2;
								}
								else
								{
									cell[i][jm].MDCS_x = cell[i][j].P1_x;
									cell[i][jm].MDCS_y = cell[i][j].P1_y;
								}
								cell[i][jm].sta = 0;
								cell[i][jm].fs = 0.00001;
								cell[i][jm].fl = 1 - cell[i][jm].fs;
								//cell[i][jm].MDCS_x = cell[i][j].P1_x;
								//cell[i][jm].MDCS_y = cell[i][j].P1_y;
								cell[i][jm].ang = cell[i][j].ang;
								cell[i][jm].Cs = cell[i][j].Cs;
								cell[i][jm].TL = 0;
							}
						}
					}
					if (precell[ip][j].sta == -1 || precell[ip][j].TL == 1)
					{
						if (abs(cell[i][j].P1_y - cell[ip][j].y) <= (len / 2)
							&& abs(cell[i][j].P1_x - cell[ip][j].x) <= (len / 2))
						{
							if (cell[ip + 1][j].fs == 1.0 || cell[ip - 1][j].fs == 1.0 ||
								cell[ip][j + 1].fs == 1.0 || cell[ip][j - 1].fs == 1.0 ||
								cell[ip + 1][j + 1].fs == 1.0 || cell[ip + 1][j - 1].fs == 1.0 ||
								cell[ip - 1][j + 1].fs == 1.0 || cell[ip - 1][j - 1].fs == 1.0)
							{
								if (cell[ip][j].sta == 0 && cell[ip][j].TL == 0)
								{
									cell[ip][j].MDCS_x = (cell[i][j].P1_x + cell[ip][j].MDCS_x) / 2;
									cell[ip][j].MDCS_y = (cell[i][j].P1_y + cell[ip][j].MDCS_y) / 2;
								}
								else
								{
									cell[ip][j].MDCS_x = cell[i][j].P1_x;
									cell[ip][j].MDCS_y = cell[i][j].P1_y;
								}
								cell[ip][j].sta = 0;
								cell[ip][j].fs = 0.00001;
								cell[ip][j].fl = 1 - cell[ip][j].fs;
								//cell[ip][j].MDCS_x = cell[i][j].P1_x;
								//cell[ip][j].MDCS_y = cell[i][j].P1_y;
								cell[ip][j].ang = cell[i][j].ang;
								cell[ip][j].Cs = cell[i][j].Cs;
								cell[ip][j].TL = 0;
							}
						}
					}

				}
				if ((abs(cell[i][j].P3_x - cell[i][j].x) >= (len / 2)
					|| abs(cell[i][j].P3_y - cell[i][j].y) >= (len / 2)))
				{
					if (precell[i][jp].sta == -1 || precell[i][jp].TL == 1)
					{
						if (abs(cell[i][j].P3_x - cell[i][jp].x) <= (len / 2)
							&& abs(cell[i][j].P3_y - cell[i][jp].y) <= (len / 2))
						{
							if (cell[i + 1][jp].fs == 1.0 || cell[i - 1][jp].fs == 1.0 ||
								cell[i][jp + 1].fs == 1.0 || cell[i][jp - 1].fs == 1.0 ||
								cell[i + 1][jp + 1].fs == 1.0 || cell[i + 1][jp - 1].fs == 1.0 ||
								cell[i - 1][jp + 1].fs == 1.0 || cell[i - 1][jp - 1].fs == 1.0)
							{
								if (cell[i][jp].sta == 0 && cell[i][jp].TL == 0)
								{
									cell[i][jp].MDCS_x = (cell[i][j].P3_x + cell[i][jp].MDCS_x) / 2;
									cell[i][jp].MDCS_y = (cell[i][j].P3_y + cell[i][jp].MDCS_y) / 2;
								}
								else
								{
									cell[i][jp].MDCS_x = cell[i][j].P3_x;
									cell[i][jp].MDCS_y = cell[i][j].P3_y;
								}
								cell[i][jp].sta = 0;
								cell[i][jp].fs = 0.00001;
								cell[i][jp].fl = 1 - cell[i][jp].fs;
								//cell[i][jp].MDCS_x = cell[i][j].P3_x;
								//cell[i][jp].MDCS_y = cell[i][j].P3_y;
								cell[i][jp].ang = cell[i][j].ang;
								cell[i][jp].Cs = cell[i][j].Cs;
								cell[i][jp].TL = 0;
							}
						}
					}
					if (precell[im][j].sta == -1 || precell[im][j].TL == 1)
					{
						if (abs(cell[i][j].P3_y - cell[im][j].y) <= (len / 2)
							&& abs(cell[i][j].P3_x - cell[im][j].x) <= (len / 2))
						{
							if (cell[im + 1][j].fs == 1.0 || cell[im - 1][j].fs == 1.0 ||
								cell[im][j + 1].fs == 1.0 || cell[im][j - 1].fs == 1.0 ||
								cell[im + 1][j + 1].fs == 1.0 || cell[im + 1][j - 1].fs == 1.0 ||
								cell[im - 1][j + 1].fs == 1.0 || cell[im - 1][j - 1].fs == 1.0)
							{
								if (cell[im][j].sta == 0 && cell[im][j].TL == 0)
								{
									cell[im][j].MDCS_x = (cell[i][j].P3_x + cell[im][j].MDCS_x) / 2;
									cell[im][j].MDCS_y = (cell[i][j].P3_y + cell[im][j].MDCS_y) / 2;
								}
								else
								{
									cell[im][j].MDCS_x = cell[i][j].P3_x;
									cell[im][j].MDCS_y = cell[i][j].P3_y;
								}
								cell[im][j].sta = 0;
								cell[im][j].fs = 0.00001;
								cell[im][j].fl = 1 - cell[im][j].fs;
								//cell[im][j].MDCS_x = cell[i][j].P3_x;
								//cell[im][j].MDCS_y = cell[i][j].P3_y;
								cell[im][j].ang = cell[i][j].ang;
								cell[im][j].Cs = cell[i][j].Cs;
								cell[im][j].TL = 0;
							}
						}
					}
					if (precell[i][jm].sta == -1 || precell[i][jm].TL == 1)
					{
						if (abs(cell[i][j].P3_x - cell[i][jm].x) <= (len / 2)
							&& abs(cell[i][j].P3_y - cell[i][jm].y) <= (len / 2))
						{
							if (cell[i + 1][jm].fs == 1.0 || cell[i - 1][jm].fs == 1.0 ||
								cell[i][jm + 1].fs == 1.0 || cell[i][jm - 1].fs == 1.0 ||
								cell[i + 1][jm + 1].fs == 1.0 || cell[i + 1][jm - 1].fs == 1.0 ||
								cell[i - 1][jm + 1].fs == 1.0 || cell[i - 1][jm - 1].fs == 1.0)
							{
								if (cell[i][jm].sta == 0 && cell[i][jm].TL == 0)
								{
									cell[i][jm].MDCS_x = (cell[i][j].P3_x + cell[i][jm].MDCS_x) / 2;
									cell[i][jm].MDCS_y = (cell[i][j].P3_y + cell[i][jm].MDCS_y) / 2;
								}
								else
								{
									cell[i][jm].MDCS_x = cell[i][j].P3_x;
									cell[i][jm].MDCS_y = cell[i][j].P3_y;
								}
								cell[i][jm].sta = 0;
								cell[i][jm].fs = 0.00001;
								cell[i][jm].fl = 1 - cell[i][jm].fs;
								//cell[i][jm].MDCS_x = cell[i][j].P3_x;
								//cell[i][jm].MDCS_y = cell[i][j].P3_y;
								cell[i][jm].ang = cell[i][j].ang;
								cell[i][jm].Cs = cell[i][j].Cs;
								cell[i][jm].TL = 0;
							}
						}
					}
					if (precell[ip][j].sta == -1 || precell[ip][j].TL == 1)
					{
						if (abs(cell[i][j].P3_y - cell[ip][j].y) <= (len / 2)
							&& abs(cell[i][j].P3_x - cell[ip][j].x) <= (len / 2))
						{
							if (cell[ip + 1][j].fs == 1.0 || cell[ip - 1][j].fs == 1.0 ||
								cell[ip][j + 1].fs == 1.0 || cell[ip][j - 1].fs == 1.0 ||
								cell[ip + 1][j + 1].fs == 1.0 || cell[ip + 1][j - 1].fs == 1.0 ||
								cell[ip - 1][j + 1].fs == 1.0 || cell[ip - 1][j - 1].fs == 1.0)
							{
								if (cell[ip][j].sta == 0 && cell[ip][j].TL == 0)
								{
									cell[ip][j].MDCS_x = (cell[i][j].P3_x + cell[ip][j].MDCS_x) / 2;
									cell[ip][j].MDCS_y = (cell[i][j].P3_y + cell[ip][j].MDCS_y) / 2;
								}
								else
								{
									cell[ip][j].MDCS_x = cell[i][j].P3_x;
									cell[ip][j].MDCS_y = cell[i][j].P3_y;
								}
								cell[ip][j].sta = 0;
								cell[ip][j].fs = 0.00001;
								cell[ip][j].fl = 1 - cell[ip][j].fs;
								//cell[ip][j].MDCS_x = cell[i][j].P3_x;
								//cell[ip][j].MDCS_y = cell[i][j].P3_y;
								cell[ip][j].ang = cell[i][j].ang;
								cell[ip][j].Cs = cell[i][j].Cs;
								cell[ip][j].TL = 0;
							}
						}
					}

				}
				if ((abs(cell[i][j].P2_x - cell[i][j].x) >= (len / 2)
					|| abs(cell[i][j].P2_y - cell[i][j].y) >= (len / 2)))
				{
					if (precell[i][jp].sta == -1 || precell[i][jp].TL == 1)
					{
						if (abs(cell[i][j].P2_x - cell[i][jp].x) <= (len / 2)
							&& abs(cell[i][j].P2_y - cell[i][jp].y) <= (len / 2))
						{
							if (cell[i + 1][jp].fs == 1.0 || cell[i - 1][jp].fs == 1.0 ||
								cell[i][jp + 1].fs == 1.0 || cell[i][jp - 1].fs == 1.0 ||
								cell[i + 1][jp + 1].fs == 1.0 || cell[i + 1][jp - 1].fs == 1.0 ||
								cell[i - 1][jp + 1].fs == 1.0 || cell[i - 1][jp - 1].fs == 1.0)
							{
								if (cell[i][jp].sta == 0 && cell[i][jp].TL == 0)
								{
									cell[i][jp].MDCS_x = (cell[i][j].P2_x + cell[i][jp].MDCS_x) / 2;
									cell[i][jp].MDCS_y = (cell[i][j].P2_y + cell[i][jp].MDCS_y) / 2;
								}
								else
								{
									cell[i][jp].MDCS_x = cell[i][j].P2_x;
									cell[i][jp].MDCS_y = cell[i][j].P2_y;
								}
								cell[i][jp].sta = 0;
								cell[i][jp].fs = 0.00001;
								cell[i][jp].fl = 1 - cell[i][jp].fs;
								//cell[i][jp].MDCS_x = cell[i][j].P2_x;
								//cell[i][jp].MDCS_y = cell[i][j].P2_y;
								cell[i][jp].ang = cell[i][j].ang;
								cell[i][jp].Cs = cell[i][j].Cs;
								cell[i][jp].TL = 0;
							}
						}
					}
					if (precell[im][j].sta == -1 || precell[im][j].TL == 1)
					{
						if (cell[im + 1][j].fs == 1.0 || cell[im - 1][j].fs == 1.0 ||
							cell[im][j + 1].fs == 1.0 || cell[im][j - 1].fs == 1.0 ||
							cell[im + 1][j + 1].fs == 1.0 || cell[im + 1][j - 1].fs == 1.0 ||
							cell[im - 1][j + 1].fs == 1.0 || cell[im - 1][j - 1].fs == 1.0)
						{
							if (abs(cell[i][j].P2_y - cell[im][j].y) <= (len / 2)
								&& abs(cell[i][j].P2_x - cell[im][j].x) <= (len / 2))
							{
								if (cell[im][j].sta == 0 && cell[im][j].TL == 0)
								{
									cell[im][j].MDCS_x = (cell[i][j].P2_x + cell[im][j].MDCS_x) / 2;
									cell[im][j].MDCS_y = (cell[i][j].P2_y + cell[im][j].MDCS_y) / 2;
								}
								else
								{
									cell[im][j].MDCS_x = cell[i][j].P2_x;
									cell[im][j].MDCS_y = cell[i][j].P2_y;
								}
								cell[im][j].sta = 0;
								cell[im][j].fs = 0.00001;
								cell[im][j].fl = 1 - cell[im][j].fs;
								//cell[im][j].MDCS_x = cell[i][j].P2_x;
								//cell[im][j].MDCS_y = cell[i][j].P2_y;
								cell[im][j].ang = cell[i][j].ang;
								cell[im][j].Cs = cell[i][j].Cs;
								cell[im][j].TL = 0;
							}
						}
					}
					if (precell[i][jm].sta == -1 || precell[i][jm].TL == 1)
					{
						if (abs(cell[i][j].P2_x - cell[i][jm].x) <= (len / 2)
							&& abs(cell[i][j].P2_y - cell[i][jm].y) <= (len / 2))
						{
							if (cell[i + 1][jm].fs == 1.0 || cell[i - 1][jm].fs == 1.0 ||
								cell[i][jm + 1].fs == 1.0 || cell[i][jm - 1].fs == 1.0 ||
								cell[i + 1][jm + 1].fs == 1.0 || cell[i + 1][jm - 1].fs == 1.0 ||
								cell[i - 1][jm + 1].fs == 1.0 || cell[i - 1][jm - 1].fs == 1.0)
							{
								if (cell[i][jm].sta == 0 && cell[i][jm].TL == 0)
								{
									cell[i][jm].MDCS_x = (cell[i][j].P2_x + cell[i][jm].MDCS_x) / 2;
									cell[i][jm].MDCS_y = (cell[i][j].P2_y + cell[i][jm].MDCS_y) / 2;
								}
								else
								{
									cell[i][jm].MDCS_x = cell[i][j].P2_x;
									cell[i][jm].MDCS_y = cell[i][j].P2_y;
								}
								cell[i][jm].sta = 0;
								cell[i][jm].fs = 0.00001;
								cell[i][jm].fl = 1 - cell[i][jm].fs;
								//cell[i][jm].MDCS_x = cell[i][j].P2_x;
								//cell[i][jm].MDCS_y = cell[i][j].P2_y;
								cell[i][jm].ang = cell[i][j].ang;
								cell[i][jm].Cs = cell[i][j].Cs;
								cell[i][jm].TL = 0;
							}
						}
					}
					if (precell[ip][j].sta == -1 || precell[ip][j].TL == 1)
					{
						if (abs(cell[i][j].P2_y - cell[ip][j].y) <= (len / 2)
							&& abs(cell[i][j].P2_x - cell[ip][j].x) <= (len / 2))
						{
							if (cell[ip + 1][j].fs == 1.0 || cell[ip - 1][j].fs == 1.0 ||
								cell[ip][j + 1].fs == 1.0 || cell[ip][j - 1].fs == 1.0 ||
								cell[ip + 1][j + 1].fs == 1.0 || cell[ip + 1][j - 1].fs == 1.0 ||
								cell[ip - 1][j + 1].fs == 1.0 || cell[ip - 1][j - 1].fs == 1.0)
							{
								if (cell[ip][j].sta == 0 && cell[ip][j].TL == 0)
								{
									cell[ip][j].MDCS_x = (cell[i][j].P2_x + cell[ip][j].MDCS_x) / 2;
									cell[ip][j].MDCS_y = (cell[i][j].P2_y + cell[ip][j].MDCS_y) / 2;
								}
								else
								{
									cell[ip][j].MDCS_x = cell[i][j].P2_x;
									cell[ip][j].MDCS_y = cell[i][j].P2_y;
								}
								cell[ip][j].sta = 0;
								cell[ip][j].fs = 0.00001;
								cell[ip][j].fl = 1 - cell[ip][j].fs;
								//cell[ip][j].MDCS_x = cell[i][j].P2_x;
								//cell[ip][j].MDCS_y = cell[i][j].P2_y;
								cell[ip][j].ang = cell[i][j].ang;
								cell[ip][j].Cs = cell[i][j].Cs;
								cell[ip][j].TL = 0;
							}
						}
					}
				}
				if ((abs(cell[i][j].P4_x - cell[i][j].x) >= (len / 2)
					|| abs(cell[i][j].P4_y - cell[i][j].y) >= (len / 2)))
				{
					if (precell[i][jp].sta == -1 || precell[i][jp].TL == 1)
					{
						if (abs(cell[i][j].P4_x - cell[i][jp].x) <= (len / 2)
							&& abs(cell[i][j].P4_y - cell[i][jp].y) <= (len / 2))
						{
							if (cell[i + 1][jp].fs == 1.0 || cell[i - 1][jp].fs == 1.0 ||
								cell[i][jp + 1].fs == 1.0 || cell[i][jp - 1].fs == 1.0 ||
								cell[i + 1][jp + 1].fs == 1.0 || cell[i + 1][jp - 1].fs == 1.0 ||
								cell[i - 1][jp + 1].fs == 1.0 || cell[i - 1][jp - 1].fs == 1.0)
							{
								if (cell[i][jp].sta == 0 && cell[i][jp].TL == 0)
								{
									cell[i][jp].MDCS_x = (cell[i][j].P4_x + cell[i][jp].MDCS_x) / 2;
									cell[i][jp].MDCS_y = (cell[i][j].P4_y + cell[i][jp].MDCS_y) / 2;
								}
								else
								{
									cell[i][jp].MDCS_x = cell[i][j].P4_x;
									cell[i][jp].MDCS_y = cell[i][j].P4_y;
								}
								cell[i][jp].sta = 0;
								cell[i][jp].fs = 0.00001;
								cell[i][jp].fl = 1 - cell[i][jp].fs;
								//cell[i][jp].MDCS_x = cell[i][j].P4_x;
								//cell[i][jp].MDCS_y = cell[i][j].P4_y;
								cell[i][jp].ang = cell[i][j].ang;
								cell[i][jp].Cs = cell[i][j].Cs;
								cell[i][jp].TL = 0;
							}
						}
					}
					if (precell[im][j].sta == -1 || precell[im][j].TL == 1)
					{
						if (abs(cell[i][j].P4_y - cell[im][j].y) <= (len / 2)
							&& abs(cell[i][j].P4_x - cell[im][j].x) <= (len / 2))
						{
							if (cell[im + 1][j].fs == 1.0 || cell[im - 1][j].fs == 1.0 ||
								cell[im][j + 1].fs == 1.0 || cell[im][j - 1].fs == 1.0 ||
								cell[im + 1][j + 1].fs == 1.0 || cell[im + 1][j - 1].fs == 1.0 ||
								cell[im - 1][j + 1].fs == 1.0 || cell[im - 1][j - 1].fs == 1.0)
							{
								if (cell[im][j].sta == 0 && cell[im][j].TL == 0)
								{
									cell[im][j].MDCS_x = (cell[i][j].P4_x + cell[im][j].MDCS_x) / 2;
									cell[im][j].MDCS_y = (cell[i][j].P4_y + cell[im][j].MDCS_y) / 2;
								}
								else
								{
									cell[im][j].MDCS_x = cell[i][j].P4_x;
									cell[im][j].MDCS_y = cell[i][j].P4_y;
								}
								cell[im][j].sta = 0;
								cell[im][j].fs = 0.00001;
								cell[im][j].fl = 1 - cell[im][j].fs;
								//cell[im][j].MDCS_x = cell[i][j].P4_x;
								//cell[im][j].MDCS_y = cell[i][j].P4_y;
								cell[im][j].ang = cell[i][j].ang;
								cell[im][j].Cs = cell[i][j].Cs;
								cell[im][j].TL = 0;
							}
						}
					}
					if (precell[i][jm].sta == -1 || precell[i][jm].TL == 1)
					{
						if (abs(cell[i][j].P4_x - cell[i][jm].x) <= (len / 2)
							&& abs(cell[i][j].P4_y - cell[i][jm].y) <= (len / 2))
						{
							if (cell[i + 1][jm].fs == 1.0 || cell[i - 1][jm].fs == 1.0 ||
								cell[i][jm + 1].fs == 1.0 || cell[i][jm - 1].fs == 1.0 ||
								cell[i + 1][jm + 1].fs == 1.0 || cell[i + 1][jm - 1].fs == 1.0 ||
								cell[i - 1][jm + 1].fs == 1.0 || cell[i - 1][jm - 1].fs == 1.0)
							{
								if (cell[i][jm].sta == 0 && cell[i][jm].TL == 0)
								{
									cell[i][jm].MDCS_x = (cell[i][j].P4_x + cell[i][jm].MDCS_x) / 2;
									cell[i][jm].MDCS_y = (cell[i][j].P4_y + cell[i][jm].MDCS_y) / 2;
								}
								else
								{
									cell[i][jm].MDCS_x = cell[i][j].P4_x;
									cell[i][jm].MDCS_y = cell[i][j].P4_y;
								}
								cell[i][jm].sta = 0;
								cell[i][jm].fs = 0.00001;
								cell[i][jm].fl = 1 - cell[i][jm].fs;
								//cell[i][jm].MDCS_x = cell[i][j].P4_x;
								//cell[i][jm].MDCS_y = cell[i][j].P4_y;
								cell[i][jm].ang = cell[i][j].ang;
								cell[i][jm].Cs = cell[i][j].Cs;
								cell[i][jm].TL = 0;
							}
						}
					}
					if (precell[ip][j].sta == -1 || precell[ip][j].TL == 1)
					{
						if (abs(cell[i][j].P4_y - cell[ip][j].y) <= (len / 2)
							&& abs(cell[i][j].P4_x - cell[ip][j].x) <= (len / 2))
						{
							if (cell[ip + 1][j].fs == 1.0 || cell[ip - 1][j].fs == 1.0 ||
								cell[ip][j + 1].fs == 1.0 || cell[ip][j - 1].fs == 1.0 ||
								cell[ip + 1][j + 1].fs == 1.0 || cell[ip + 1][j - 1].fs == 1.0 ||
								cell[ip - 1][j + 1].fs == 1.0 || cell[ip - 1][j - 1].fs == 1.0)
							{
								if (cell[ip][j].sta == 0 && cell[ip][j].TL == 0)
								{
									cell[ip][j].MDCS_x = (cell[i][j].P4_x + cell[ip][j].MDCS_x) / 2;
									cell[ip][j].MDCS_y = (cell[i][j].P4_y + cell[ip][j].MDCS_y) / 2;
								}
								else
								{
									cell[ip][j].MDCS_x = cell[i][j].P4_x;
									cell[ip][j].MDCS_y = cell[i][j].P4_y;
								}
								cell[ip][j].sta = 0;
								cell[ip][j].fs = 0.00001;
								cell[ip][j].fl = 1 - cell[ip][j].fs;
								//cell[ip][j].MDCS_x = cell[i][j].P4_x;
								//cell[ip][j].MDCS_y = cell[i][j].P4_y;
								cell[ip][j].ang = cell[i][j].ang;
								cell[ip][j].Cs = cell[i][j].Cs;
								cell[ip][j].TL = 0;
							}
						}
					}
				}
			}
		}
	}
}

//更新元胞状态	液相胞→界面胞
void LtoICell_1(Cell** cell)
{
#pragma omp parallel for
	for (int i = 2; i <= ROWS + 1; i++)
	{
		for (int j = 2; j <= COLS + 1; j++)
		{
			if (cell[i][j].fs == 1.0)
			{
				cell[i][j].vel_x = 0.0;
				cell[i][j].vel_y = 0.0;
				cell[i][j].vel = 0.0;
				cell[i][j].sta = 1;
				cell[i][j].fl = 0.0;
				cell[i][j].Cl = 0.0;
				cell[i][j].Cl8 = 0.0;
				cell[i][j].Cs8 = 0.0;
				cell[i][j].cur = 0.0;
				cell[i][j].ani = 0.0;
				cell[i][j].delta_fs = 0.0;

				if (cell[i + 1][j].sta == -1)
				{
					cell[i + 1][j].TL = 1;
					cell[i + 1][j].sta = 0;
					cell[i + 1][j].fs = 0.00001;
					cell[i + 1][j].fl = 1 - cell[i + 1][j].fs;
					cell[i + 1][j].MDCS_x = 0;
					cell[i + 1][j].MDCS_y = 0;
					cell[i + 1][j].ang = cell[i][j].ang;
					cell[i + 1][j].Cs = cell[i][j].Cs;
				}
				if (cell[i - 1][j].sta == -1)
				{
					cell[i - 1][j].TL = 1;
					cell[i - 1][j].sta = 0;
					cell[i - 1][j].fs = 0.00001;
					cell[i - 1][j].fl = 1 - cell[i - 1][j].fs;
					cell[i - 1][j].MDCS_x = 0;
					cell[i - 1][j].MDCS_y = 0;
					cell[i - 1][j].ang = cell[i][j].ang;
					cell[i - 1][j].Cs = cell[i][j].Cs;
				}
				if (cell[i][j - 1].sta == -1)
				{
					cell[i][j - 1].TL = 1;
					cell[i][j - 1].sta = 0;
					cell[i][j - 1].fs = 0.00001;
					cell[i][j - 1].fl = 1 - cell[i][j - 1].fs;
					cell[i][j - 1].MDCS_x = 0;
					cell[i][j - 1].MDCS_y = 0;
					cell[i][j - 1].ang = cell[i][j].ang;
					cell[i][j - 1].Cs = cell[i][j].Cs;
				}
				if (cell[i][j + 1].sta == -1)
				{
					cell[i][j + 1].TL = 1;
					cell[i][j + 1].sta = 0;
					cell[i][j + 1].fs = 0.00001;
					cell[i][j + 1].fl = 1 - cell[i][j + 1].fs;
					cell[i][j + 1].MDCS_x = 0;
					cell[i][j + 1].MDCS_y = 0;
					cell[i][j + 1].ang = cell[i][j].ang;
					cell[i][j + 1].Cs = cell[i][j].Cs;
				}
			}
		}
	}
}

//更新元胞状态	界面胞→固面胞
void ItoSCell(Cell** cell)
{
#pragma omp parallel for
	for (int i = 2; i <= ROWS + 1; i++)
	{
		for (int j = 2; j <= COLS + 1; j++)
		{
			if (cell[i][j].sta == 0)
			{
				if (cell[i + 1][j].sta != -1 && cell[i - 1][j].sta != -1 &&
					cell[i][j + 1].sta != -1 && cell[i][j - 1].sta != -1 &&
					cell[i + 1][j + 1].sta != -1 && cell[i + 1][j - 1].sta != -1 &&
					cell[i - 1][j + 1].sta != -1 && cell[i - 1][j - 1].sta != -1)
				{
					cell[i][j].fs = 1.0;
					cell[i][j].vel_x = 0.0;
					cell[i][j].vel_y = 0.0;
					cell[i][j].vel = 0.0;
					cell[i][j].sta = 1;
					cell[i][j].fl = 0.0;
					cell[i][j].Cl = 0.0;
					cell[i][j].Cl8 = 0.0;
					cell[i][j].Cs8 = 0.0;
					cell[i][j].cur = 0.0;
					cell[i][j].ani = 0.0;
					cell[i][j].delta_fs = 0.0;
				}
			}
		}
	}
}

//更新总浓度
void updateCe(Cell** cell)
{
#pragma omp parallel for
	for (int i = 2; i <= ROWS + 1; i++)
	{
		for (int j = 2; j <= COLS + 1; j++)
		{
			cell[i][j].Ce = cell[i][j].fs * cell[i][j].Cs + (1 - cell[i][j].fs) * cell[i][j].Cl;
		}
	}
}

//计算曲率与各项异性（双线性插值法）
void Curvature(Cell** cell)
{
#pragma omp parallel for
	for (int i = 2; i <= ROWS + 1; i++)
	{
		for (int j = 2; j <= COLS + 1; j++)
		{
			if (cell[i][j].sta == 0 && cell[i][j].TL == 0)
			{
				double dfsdx, dfsdy, dfs2dx, dfs2dy, dfsdxdy;
				double F1, F2, F3, F4, F5, F6, F7, F8, F9, F10, F11, F12;
				double fN, fS, fW, fE, fP;

				F1 = (cell[i - 2][j - 1].fs + cell[i - 2][j].fs + cell[i - 1][j - 1].fs + cell[i - 1][j].fs) / 4;
				F2 = (cell[i - 2][j + 1].fs + cell[i - 2][j].fs + cell[i - 1][j + 1].fs + cell[i - 1][j].fs) / 4;
				F3 = (cell[i - 1][j - 2].fs + cell[i - 1][j - 1].fs + cell[i][j - 2].fs + cell[i][j - 1].fs) / 4;
				F4 = (cell[i - 1][j - 1].fs + cell[i - 1][j].fs + cell[i][j - 1].fs + cell[i][j].fs) / 4;
				F5 = (cell[i - 1][j + 1].fs + cell[i - 1][j].fs + cell[i][j + 1].fs + cell[i][j].fs) / 4;
				F6 = (cell[i - 1][j + 1].fs + cell[i - 1][j + 2].fs + cell[i][j + 1].fs + cell[i][j + 2].fs) / 4;
				F7 = (cell[i + 1][j - 2].fs + cell[i + 1][j - 1].fs + cell[i][j - 2].fs + cell[i][j - 1].fs) / 4;
				F8 = (cell[i + 1][j - 1].fs + cell[i + 1][j].fs + cell[i][j - 1].fs + cell[i][j].fs) / 4;
				F9 = (cell[i + 1][j + 1].fs + cell[i + 1][j].fs + cell[i][j + 1].fs + cell[i][j].fs) / 4;
				F10 = (cell[i + 1][j + 1].fs + cell[i + 1][j + 2].fs + cell[i][j + 1].fs + cell[i][j + 2].fs) / 4;
				F11 = (cell[i + 1][j - 1].fs + cell[i + 1][j].fs + cell[i + 2][j - 1].fs + cell[i + 2][j].fs) / 4;
				F12 = (cell[i + 1][j + 1].fs + cell[i + 1][j].fs + cell[i + 2][j + 1].fs + cell[i + 2][j].fs) / 4;

				fN = (F1 + F2 + F4 + F5) / 4;
				fS = (F8 + F9 + F11 + F12) / 4;
				fW = (F3 + F4 + F7 + F8) / 4;
				fE = (F5 + F6 + F9 + F10) / 4;
				fP = (F4 + F5 + F8 + F9) / 4;

				dfsdx = (fE - fW) / 2 / len;
				dfsdy = (fS - fN) / 2 / len;
				dfs2dx = (fE + fW - 2 * fP) / len / len;
				dfs2dy = (fS + fN - 2 * fP) / len / len;
				dfsdxdy = (F4 + F9 - F5 - F8) / len / len;

				cell[i][j].nx = dfsdx / sqrt(dfsdx * dfsdx + dfsdy * dfsdy);
				cell[i][j].ny = dfsdy / sqrt(dfsdx * dfsdx + dfsdy * dfsdy);

				//if (abs(cell[i][j].nx) < 0.0001)cell[i][j].nx = 0.0;
				//if (abs(cell[i][j].ny) < 0.0001)cell[i][j].ny = 0.0;

				cell[i][j].cur = (2 * dfsdx * dfsdy * dfsdxdy - dfs2dx * dfsdy * dfsdy - dfs2dy * dfsdx * dfsdx)
					/ sqrt((dfsdx * dfsdx + dfsdy * dfsdy) * (dfsdx * dfsdx + dfsdy * dfsdy) * (dfsdx * dfsdx + dfsdy * dfsdy));

				if (cell[i][j].ny >= 0)
				{
					cell[i][j].angn = acos(-dfsdx / sqrt(dfsdx * dfsdx + dfsdy * dfsdy));
					cell[i][j].ani = 1 - 15 * angd * cos(4 * (acos(-dfsdx / sqrt(dfsdx * dfsdx + dfsdy * dfsdy)) - cell[i][j].ang));
				}
				else if (cell[i][j].ny < 0)
				{
					cell[i][j].angn = 2 * M_PI - acos(-dfsdx / sqrt(dfsdx * dfsdx + dfsdy * dfsdy));
					cell[i][j].ani = 1 - 15 * angd * cos(4 * (2 * M_PI - acos(-dfsdx / sqrt(dfsdx * dfsdx + dfsdy * dfsdy)) - cell[i][j].ang));
				}
			}
		}
	}
}

//界面胞界面两侧溶质浓度
void Cl8Cs8(Cell** cell)
{
#pragma omp parallel for
	for (int i = 2; i <= ROWS + 1; i++)
	{
		for (int j = 2; j <= COLS + 1; j++)
		{
			if (cell[i][j].sta == 0 && cell[i][j].TL == 0)
			{
				cell[i][j].Cl8 = C0 + (cell[i][j].T - T0 + GT * cell[i][j].cur * cell[i][j].ani) / ml;
				cell[i][j].Cs8 = k * cell[i][j].Cl8;
			}
		}
	}
}

//计算界面胞生长速度
void vel(Cell** cell)
{
#pragma omp parallel for
	for (int i = 3; i <= ROWS; i++)
	{
		for (int j = 3; j <= COLS; j++)
		{
			if (cell[i][j].sta == 0 && cell[i][j].TL == 0)
			{
				int im = i - 1, ip = i + 1, jm = j - 1, jp = j + 1;

				cell[i][j].vel_x = Dl(cell[i][j].T) / len / (1 - k)
					* ((1 - cell[i][jm].Cl / cell[i][j].Cl8) * cell[i][jm].fl
						+ (1 - cell[i][jp].Cl / cell[i][j].Cl8) * cell[i][jp].fl)
					+ k * Ds(cell[i][j].T) / len / (1 - k)
					* ((1 - cell[i][jm].Cs / cell[i][j].Cs8) * cell[i][jm].fs
						+ (1 - cell[i][jp].Cs / cell[i][j].Cs8) * cell[i][jp].fs);

				cell[i][j].vel_y = Dl(cell[i][j].T) / len / (1 - k)
					* ((1 - cell[im][j].Cl / cell[i][j].Cl8) * cell[im][j].fl
						+ (1 - cell[ip][j].Cl / cell[i][j].Cl8) * cell[ip][j].fl)
					+ k * Ds(cell[i][j].T) / len / (1 - k)
					* ((1 - cell[im][j].Cs / cell[i][j].Cs8) * cell[im][j].fs
						+ (1 - cell[ip][j].Cs / cell[i][j].Cs8) * cell[ip][j].fs);

				if (cell[i][j].vel_x >= 0 && cell[i][j].vel_y >= 0)
				{
					cell[i][j].vel = cell[i][j].vel_x * abs(cell[i][j].nx) + cell[i][j].vel_y * abs(cell[i][j].ny);
				}
				else if (cell[i][j].vel_x >= 0 && cell[i][j].vel_y < 0)
				{
					cell[i][j].vel = cell[i][j].vel_x * abs(cell[i][j].nx);
				}
				else if (cell[i][j].vel_x < 0 && cell[i][j].vel_y >= 0)
				{
					cell[i][j].vel = cell[i][j].vel_y * abs(cell[i][j].ny);
				}
				else
				{
					cell[i][j].vel = 0.0;
				}
			}
		}
	}
}

//计算几何修正因子
void gf(Cell** cell)
{
#pragma omp parallel for
	for (int i = 2; i <= ROWS + 1; i++)
	{
		for (int j = 2; j <= COLS + 1; j++)
		{
			cell[i][j].gf = 0.0;
			if (cell[i][j].sta == 0)
			{
				int im = i - 1, ip = i + 1, jm = j - 1, jp = j + 1;
				double S = 0.0, S1 = 0.0, S2 = 0.0;
				if (cell[i][j].Cl8 >= cell[i][j].Cl)	//凝固
				{
					if (cell[im][j].fs == 1.0)S1 = S1 + 1;
					if (cell[ip][j].fs == 1.0)S1 = S1 + 1;
					if (cell[i][jm].fs == 1.0)S1 = S1 + 1;
					if (cell[i][jp].fs == 1.0)S1 = S1 + 1;

					if (cell[im][jm].fs == 1.0)S2 = S2 + 1;
					if (cell[im][jp].fs == 1.0)S2 = S2 + 1;
					if (cell[ip][jm].fs == 1.0)S2 = S2 + 1;
					if (cell[ip][jp].fs == 1.0)S2 = S2 + 1;
				}
				S = (S1 + S2 / sqrt(2.0)) / 2.0;
				cell[i][j].gf = double(min(S, 1.0));
			}
		}
	}
}

//计算偏心正方形半对角线长度
void Ldia(Cell** cell, double delta_t)
{
#pragma omp parallel for
	for (int i = 2; i <= ROWS + 1; i++)
	{
		for (int j = 2; j <= COLS + 1; j++)
		{
			if (cell[i][j].sta == 0 && cell[i][j].fs != 1.0)
			{
				cell[i][j].L_dia = cell[i][j].L_dia + cell[i][j].gf * cell[i][j].vel * delta_t;
			}
		}
	}
}

//计算固相增长
void deltafs(Cell** cell, double delta_t)
{
#pragma omp parallel for
	for (int i = 2; i <= ROWS + 1; i++)
	{
		for (int j = 2; j <= COLS + 1; j++)
		{
			cell[i][j].delta_fs = 0.0;
			if (cell[i][j].sta == 0)
			{
				cell[i][j].delta_fs = (cell[i][j].gf * cell[i][j].vel * delta_t)
					/ (len * (abs(sin(cell[i][j].ang)) + cos(cell[i][j].ang)));
			}
		}
	}
}

//保存上一次迭代后的元胞状态2
void copyCell_2(Cell** cell, Cell** precell)
{
#pragma omp parallel for
	for (int i = 0; i <= ROWS + 3; i++)
	{
		for (int j = 0; j <= COLS + 3; j++)
		{
			precell[i][j].Cs = cell[i][j].Cs;
			precell[i][j].Cl = cell[i][j].Cl;
			precell[i][j].fs = cell[i][j].fs;
		}
	}
}

//更新元胞固相分数
void updatefs(Cell** cell, Cell** precell)
{
#pragma omp parallel for
	for (int i = 2; i <= ROWS + 1; i++)
	{
		for (int j = 2; j <= COLS + 1; j++)
		{
			if (cell[i][j].sta == 0 && cell[i][j].TL == 0)
			{
				if (cell[i][j].delta_fs > 1 - cell[i][j].fs)
				{
					cell[i][j].delta_fs = 1 - cell[i][j].fs;
				}
				cell[i][j].fs = cell[i][j].fs + cell[i][j].delta_fs;
				cell[i][j].fl = 1 - cell[i][j].fs;

				cell[i][j].Cs = (precell[i][j].Cs * precell[i][j].fs + k * precell[i][j].Cl * cell[i][j].delta_fs)
					/ (precell[i][j].fs + cell[i][j].delta_fs);
			}
		}
	}
}

//保存上一次迭代后的元胞状态3
void copyCell_3(Cell** cell, Cell** precell)
{
#pragma omp parallel for
	for (int i = 0; i <= ROWS + 3; i++)
	{
		for (int j = 0; j <= COLS + 3; j++)
		{
			precell[i][j].Cs = cell[i][j].Cs;
			precell[i][j].Cl = cell[i][j].Cl;
		}
	}
}

//计算液相溶质场求解系数
void Cla(Cell** cell, double delta_t)
{
#pragma omp parallel for
	for (int i = 2; i <= ROWS + 1; i++)
	{
		for (int j = 2; j <= COLS + 1; j++)
		{
			if (cell[i][j].sta == -1 || cell[i][j].sta == 0)
			{
				cell[i][j].Clae0 = Dl(cell[i][j].T) * dy / dx;
				cell[i][j].Claw0 = Dl(cell[i][j].T) * dy / dx;
				cell[i][j].Clan0 = Dl(cell[i][j].T) * dx / dy;
				cell[i][j].Clas0 = Dl(cell[i][j].T) * dx / dy;
				if (cell[i + 1][j].sta == 1)cell[i][j].Clas0 = 0.0;
				if (cell[i - 1][j].sta == 1)cell[i][j].Clan0 = 0.0;
				if (cell[i][j + 1].sta == 1)cell[i][j].Clae0 = 0.0;
				if (cell[i][j - 1].sta == 1)cell[i][j].Claw0 = 0.0;
				cell[i][j].Clap1 = dx * dy / delta_t;
				cell[i][j].Clap0 = cell[i][j].Clap1 - cell[i][j].Clae0 - cell[i][j].Claw0 - cell[i][j].Clan0 - cell[i][j].Clas0;
				cell[i][j].Clbp = cell[i][j].Cl * (1 - k) * cell[i][j].delta_fs * dx * dy / delta_t;
			}
		}
	}
}

//计算固相溶质场求解系数
void Csa(Cell** cell, double delta_t)
{
#pragma omp parallel for
	for (int i = 2; i <= ROWS + 1; i++)
	{
		for (int j = 2; j <= COLS + 1; j++)
		{
			if (cell[i][j].sta == 1 || cell[i][j].sta == 0)
			{
				cell[i][j].Csae0 = Ds(cell[i][j].T) * dy / dx;
				cell[i][j].Csaw0 = Ds(cell[i][j].T) * dy / dx;
				cell[i][j].Csan0 = Ds(cell[i][j].T) * dx / dy;
				cell[i][j].Csas0 = Ds(cell[i][j].T) * dx / dy;
				if (cell[i + 1][j].sta == -1)cell[i][j].Csas0 = 0.0;
				if (cell[i - 1][j].sta == -1)cell[i][j].Csan0 = 0.0;
				if (cell[i][j + 1].sta == -1)cell[i][j].Csae0 = 0.0;
				if (cell[i][j - 1].sta == -1)cell[i][j].Csaw0 = 0.0;
				cell[i][j].Csap1 = dx * dy / delta_t;
				cell[i][j].Csap0 = cell[i][j].Csap1 - cell[i][j].Csae0 - cell[i][j].Csaw0 - cell[i][j].Csan0 - cell[i][j].Csas0;
			}
		}
	}
}

//计算域边界条件
void Cl0(Cell** cell)
{
#pragma omp parallel for 
	for (int i = 2; i <= ROWS + 1; i++)
	{
		cell[i][1].Cln0 = C0;
		cell[i][2].Claw0 = 0;		//无扩散条件
	}
#pragma omp parallel for 
	for (int i = 2; i <= ROWS + 1; i++)
	{
		cell[i][COLS + 2].Cln0 = C0;
		cell[i][COLS + 1].Clae0 = 0;
	}
#pragma omp parallel for 
	for (int j = 2; j <= COLS + 1; j++)
	{
		cell[1][j].Cln0 = C0;
		cell[2][j].Clan0 = 0;
	}
#pragma omp parallel for 
	for (int j = 2; j <= COLS + 1; j++)
	{
		cell[ROWS + 2][j].Cln0 = C0;
		cell[ROWS + 1][j].Clas0 = 0;
	}
}

//计算溶质场求解系数
void aCl0(Cell** cell)
{
#pragma omp parallel for 
	for (int i = 2; i <= ROWS + 1; i++)
	{
		for (int j = 2; j <= COLS + 1; j++)
		{
			cell[i][j].Clap0 = cell[i][j].Clap1 - cell[i][j].Clae0 - cell[i][j].Claw0 - cell[i][j].Clan0 - cell[i][j].Clas0;
		}
	}
}

//液相溶质显式迭代
void Cl(Cell** cell, Cell** precell)
{
#pragma omp parallel for
	for (int i = 2; i <= ROWS + 1; i++)
	{
		for (int j = 2; j <= COLS + 1; j++)
		{
			if (cell[i][j].sta == -1 || cell[i][j].sta == 0)
			{
				cell[i][j].Cl = (cell[i][j].Clae0 * precell[i][j + 1].Cl + cell[i][j].Claw0 * precell[i][j - 1].Cl
					+ cell[i][j].Clan0 * precell[i - 1][j].Cl + cell[i][j].Clas0 * precell[i + 1][j].Cl
					+ cell[i][j].Clap0 * precell[i][j].Cl + cell[i][j].Clbp) / cell[i][j].Clap1;
			}
		}
	}
}

//固相溶质显式迭代
void Cs(Cell** cell, Cell** precell)
{
#pragma omp parallel for
	for (int i = 2; i <= ROWS + 1; i++)
	{
		for (int j = 2; j <= COLS + 1; j++)
		{
			if (cell[i][j].sta == 1 || cell[i][j].sta == 0)
			{
				cell[i][j].Cs = (cell[i][j].Csae0 * precell[i][j + 1].Cs + cell[i][j].Csaw0 * precell[i][j - 1].Cs
					+ cell[i][j].Csan0 * precell[i - 1][j].Cs + cell[i][j].Csas0 * precell[i + 1][j].Cs
					+ cell[i][j].Csap0 * precell[i][j].Cs) / cell[i][j].Csap1;
			}
		}
	}
}

//数据写入
void rData(Cell** cell, string filename) {
	//定义输出文件流
	ofstream outFile;
	//打开文件用于写入
	outFile.open(filename, ios::out | ios::trunc); // 打开模式可省略
	//outFile << t << endl;
	outFile << "x" << ',';
	outFile << "y" << ',';
	outFile << "z" << ',';
	outFile << "sta" << ',';
	outFile << "Ce" << ',';
	outFile << "Cl8" << ',';
	outFile << "Cs8" << ',';
	outFile << "Cs" << ',';
	outFile << "Cl" << ',';
	outFile << "vel_x" << ',';
	outFile << "vel_y" << ',';
	outFile << "cur" << ',';
	outFile << "fs" << ',';
	outFile << "fl" << ',';
	outFile << "ani" << ',';
	outFile << "T" << ',';
	outFile << "nuc_p" << ',';
	outFile << "delta_fs" << ',';
	outFile << "L_dia" << ',';
	outFile << "gf" << ',';
	outFile << "vel" << ',';
	outFile << "nx" << ',';
	outFile << "ny" << ',';
	outFile << "Px" << ',';
	outFile << "Py" << ',';
	outFile << "P1x" << ',';
	outFile << "P1y" << ',';
	outFile << "P2x" << ',';
	outFile << "P2y" << ',';
	outFile << "P3x" << ',';
	outFile << "P3y" << ',';
	outFile << "P4x" << ',';
	outFile << "P4y" << ',';
	outFile << "x" << ',';
	outFile << "y" << ',';
	outFile << "test" << ',';
	outFile << endl;

	for (int i = 2; i <= ROWS + 1; i++)
	{
		for (int j = 2; j <= COLS + 1; j++)
		{
			outFile << cell[i][j].x << ',';
			outFile << cell[i][j].y << ',';
			outFile << 0 << ',';
			outFile << cell[i][j].sta << ',';
			outFile << cell[i][j].Ce << ',';
			outFile << cell[i][j].Cl8 << ',';
			outFile << cell[i][j].Cs8 << ',';
			outFile << cell[i][j].Cs << ',';
			outFile << cell[i][j].Cl << ',';
			outFile << cell[i][j].vel_x << ',';
			outFile << cell[i][j].vel_y << ',';
			outFile << cell[i][j].cur << ',';
			outFile << cell[i][j].fs << ',';
			outFile << cell[i][j].fl << ',';
			outFile << cell[i][j].ani << ',';
			outFile << setprecision(11) << cell[i][j].T << ',';
			outFile << cell[i][j].nuc_p << ',';
			outFile << cell[i][j].delta_fs << ',';
			outFile << cell[i][j].L_dia << ',';
			outFile << cell[i][j].gf << ',';
			outFile << cell[i][j].vel << ',';
			outFile << cell[i][j].nx << ',';
			outFile << cell[i][j].ny << ',';
			outFile << cell[i][j].MDCS_x << ',';
			outFile << cell[i][j].MDCS_y << ',';
			outFile << cell[i][j].P1_x << ',';
			outFile << cell[i][j].P1_y << ',';
			outFile << cell[i][j].P2_x << ',';
			outFile << cell[i][j].P2_y << ',';
			outFile << cell[i][j].P3_x << ',';
			outFile << cell[i][j].P3_y << ',';
			outFile << cell[i][j].P4_x << ',';
			outFile << cell[i][j].P4_y << ',';
			outFile << cell[i][j].x << ',';
			outFile << cell[i][j].y << ',';
			outFile << cell[i][j].test << ',';
			outFile << endl;
		}
	}
	//关闭文件流
	outFile.close();
}

//晶粒绘制
void drawCell(Cell** cell) {
	for (int i = 2; i <= ROWS + 1; i++)
	{
		for (int j = 2; j <= COLS + 1; j++)
		{
			double r = 0, g = 0, b = 0;
			if (cell[i][j].sta == 1)
			{
				if ((cell[i][j].ang + M_PI_4) / M_PI_2 <= (4.0 / 7.0))
				{
					r = 255; g = (cell[i][j].ang + M_PI_4) / M_PI_2 * 7.0 * (255.0 / 4.0); b = 0;
				}
				else if ((cell[i][j].ang + M_PI_4) / M_PI_2 <= (6.0 / 7.0))
				{
					//r = 255; g = 0; b = 0;
					r = 255 - ((cell[i][j].ang + M_PI_4) / M_PI_2 - (4.0 / 7.0)) * 7.0 * (255.0 / 2.0); g = 255; b = 0;
				}
				else
				{
					r = 0; g = 255; b = ((cell[i][j].ang + M_PI_4) / M_PI_2 - (6.0 / 7.0)) * 7.0 * (255.0 / 2.0);
				}
				r = 255; g = 0; b = 0;
			}
			else if (cell[i][j].sta == 0)
			{
				r = 255; g = 255; b = 255;
			}
			else if (cell[i][j].sta == -1)
			{
				r = 0; g = 0; b = 255;
			}

			drawRect(i - 2, j - 2, RGB(r, g, b));
		}
	}
}

//求迭代时间间隔
double fun_delta_t(Cell** cell, double vell, double Dss, double Dll)
{
	double delta_t;
	for (int i = 2; i <= ROWS + 1; i++)
	{
		for (int j = 2; j <= COLS + 1; j++)
		{
			vell = max(vell, cell[i][j].vel);
			Dss = max(Dss, Ds(cell[i][j].T));
			Dll = max(Dll, Dl(cell[i][j].T));
		}
	}
	delta_t = 0.25 * min(min((len / vell), (len * len / Dll)), (len * len / Dss));
	//if (delta_t > dt)delta_t = dt;
	return delta_t;
}