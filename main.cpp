#include "Cell.h"

int main()
{
	//�����߳���
	omp_set_num_threads(12);
	//��̬�߳̿���
	omp_set_dynamic(0);
	//����ʱ����� ��ʼ��
	time_t start, stop;
	start = time(NULL);

	//����random_device����sd������
	random_device sd;
	//ʹ�����ӳ�ʼ��linear_congruential_engine����
	mt19937_64 linearRan(sd());
	uniform_real_distribution<double>dis1(0.0, 1.0);
	uniform_int_distribution<int>dis2(1, 89);

	//��ʼ��Ԫ���Ķ�ά�ṹ�岢���䶯̬�ڴ�
	//��ʼ��������һ����Ԫ��״̬�Ķ�ά�ṹ�岢���䶯̬�ڴ�
	Cell** cell, ** precell;
	cell = new Cell * [ROWS + 4];
	precell = new Cell * [ROWS + 4];
	for (int i = 0; i < ROWS + 4; i++)
	{
		cell[i] = new Cell[COLS + 4];
		precell[i] = new Cell[COLS + 4];
	}

	//��ʼ��Ԫ��״̬
	initCell(cell, precell);

	//��ʼ��Ԫ������
	initseed(cell, precell);

	//Ԥ��ͼ������
	initgraph(SPACE * ROWS, SPACE * COLS, EW_SHOWCONSOLE);
	setbkcolor(0xAA0000);
	cleardevice();

	//��ʼ������ʱ����ʱ����
	//��ʼ��������������
	double t = 0, delta_t = dt;
	int ite = 0;
	string filename = "CA";			//�洢�ļ�ǰ׺

	double vel_max = DBL_MIN;		//��ȡһ������������Сֵ����ֹ��ĸΪ��
	double Ds_max = DBL_MIN, Dl_max = DBL_MIN;

	//������ʼ
	while (ite < num_ite)
	{
		//ȷ������Ԫ����Ӧƫ�������ζ�������
		ESVCCell(cell);

		//������һ�ε������Ԫ��״̬ 1
		copyCell_1(cell, precell);

		//ȷ������Ԫ���ж��Ƿ�׽��ΧҺ���	Һ����������
		LtoICell(cell, precell);

		//����Ԫ��״̬	Һ����������
		LtoICell_1(cell);

		//����Ԫ��״̬	������������
		ItoSCell(cell);

		//������Ũ��
		updateCe(cell);

		//����������������ԣ�˫���Բ�ֵ����
		Curvature(cell);

		//�����������������Ũ��
		Cl8Cs8(cell);

		//�������������ٶ�
		vel(cell);
	
		//���㼸����������
		gf(cell);

		//����ƫ�������ΰ�Խ��߳���
		Ldia(cell, delta_t);

		//�����������
		deltafs(cell, delta_t);

		//������һ�ε������Ԫ��״̬2
		copyCell_2(cell, precell);

		//����Ԫ���������
		updatefs(cell, precell);

		//������һ�ε������Ԫ��״̬3
		copyCell_3(cell, precell);

		//-----------------------------------------
		//�¶ȳ�����
		
		//-----------------------------------------
		//���ʳ�����
		//����Һ�����ʳ����ϵ��
		Cla(cell, delta_t);

		//����������ʳ����ϵ��
		Csa(cell, delta_t);

		//������߽�����
		Cl0(cell);

		//�������ʳ����ϵ��
		aCl0(cell);

		//Һ��������ʽ����
		Cl(cell, precell);

		//����������ʽ����
		Cs(cell, precell);

		//-----------------------------------------
		//���ݵ���
		//if (t >= 0.212)
		//if (ite % 1000 == 0.0)
		//if (ite > num_ite - 2)
		//if (true)
		if (false)
		{
			filename = "./data/" + to_string(ite) + "s CA.csv";
			rData(cell, filename);
		}

		//���ƾ���	
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
			sprintf_s(buf, "./data/%d.png", ite);//�ַ����ﻻ������ļ�·��
			saveimage(buf);
		}

		//�����ʱ����
		vel_max = DBL_MIN;
		Ds_max = DBL_MIN;
		Dl_max = DBL_MIN;
		delta_t = fun_delta_t(cell, vel_max, Ds_max, Dl_max);

		//������ʱ��
		t = t + delta_t;
		ite = ite + 1;
		//if (true)
		if (ite % 100 == 0.0)

		{
			cout << "����ʱ�䣺" << t << endl;
			cout << "ʱ������" << delta_t << endl;
			cout << "��" << ite << "�ε���" << endl;
		}

		//cout << "��" << ite << "�ε���" << endl;
	}
	//outFile.close();

	cout << "����ʱ��" << t << endl;
	cout << "�������" << endl;
	stop = time(NULL);
	stop = stop - start;
	cout << "����ʱ��" << stop << endl;
	system("pause");

	//�ͷ��ڴ�
	for (int i = 0; i < ROWS + 4; i++)
	{
		delete[] cell[i];
		delete[] precell[i];
	}
	delete[] cell;
	delete[] precell;

	return 0;
}