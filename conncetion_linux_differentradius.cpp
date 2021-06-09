// graph_connection.cpp : ���ļ����� "main" ����������ִ�н��ڴ˴���ʼ��������
//
#include <iostream>
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sstream>
#pragma warning(disable : 4996)
using namespace std;
double three_jie_chaji(double **a);
double vector_angle(double**a);
int inone_face(double **a);                          //���淵��1�������淵��0
double det(double **a, int num);                     //����������ʽ��ֵ
void getastar(double **a, int num, double **b);      //��������������ʽ����䵽b����
int reverse_matrix(double **a, int num, double **b); //�����������󣬳ɹ����1 ���������0
const int cengshu = 5;                               //��������������������ھ�����չʱ�Ĳ���
const int yanshen = cengshu * cengshu * cengshu; 
void buble_plus(double* a, int *b, int num);
const double pi = 3.1415926;//˵�����˶��پ���
const char a[120][3] = { " ", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg" };
double dist[120][120];
int meatal_xuhao[88] = { 3, 4, 11, 12, 13, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,37, 38, 39, 40, 41, 42, 43, 44, 46, 47, 48, 49, 50,51, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111 };
void random_choose_two(int*org, int num, int**in);
double dis(double *p1, double *p2)
{
	return sqrt(pow(p1[0] - p2[0], 2) + pow(p1[1] - p2[1], 2) + pow(p1[2] - p2[2], 2));
}
void new_get_style(char *style)
{
	ifstream fin;
	fin.open(style, ios::in);
	if (!fin.is_open())
	{
		cout << "i can not find the file:" << style << endl;
		cin.get();
	}
	while (fin.good())
	{
		for (int i = 1; i < 112; i++)
		{
			for (int j = 1; j < 112; j++)
			{
				fin >> dist[i][j];
			}
		}
	}
	fin.close();
	return;
}
class cell
{
public:
	cell(char *jiegou_name);
	int num = 0;
	double **letice;
	double **p;
	double ***p_real;
	double ***real_position;
	double *ridus;
	int *type;
	int type_num;
	int* type_save;
};

cell::cell(char *name)
{
	int i, j, k;
	//cout << "expand the :" << cengshu << "layer" << endl;
	char temp[300];
	double x_pian = 0.0;
	double y_pian = 0.0;
	double z_pian = 0.0;
	//strcpy(wenjian, "atom1.config");
	FILE *in;
	in = fopen(name, "rt");
	//system("pause");
	if (in == NULL)
	{
		printf("error of rading atom.config!\n");
		printf("the filename is :%s\n", name);
		//cin.get();
		return;
	}
	fscanf(in, "%d", &num);

	type = (int *)malloc(num * sizeof(int));
	letice = (double **)malloc(3 * sizeof(double *));
	for (i = 0; i < 3; i++)
	{
		letice[i] = (double *)malloc(3 * sizeof(double));
	}
	p = (double **)malloc(num * sizeof(double *));
	for (i = 0; i < num; i++)
	{
		p[i] = (double *)malloc(3 * sizeof(double));
	}
	real_position = (double ***)malloc(yanshen * sizeof(double **));
	for (i = 0; i < yanshen; i++)
	{
		real_position[i] = (double **)malloc(num * sizeof(double *));
		for (k = 0; k < num; k++)
			real_position[i][k] = (double *)malloc(3 * sizeof(double));
	}

	p_real = new double **[yanshen];
	for (i = 0; i < yanshen; i++)
	{
		p_real[i] = new double *[num];
		for (k = 0; k < num; k++)
		{
			p_real[i][k] = new double[3];
		}
	}
	ridus = new double[num];
	while (fgets(temp, 300, in) != NULL)
	{
		if (strstr(temp, "VECTOR") != NULL)
			break;
	}
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			fscanf(in, "%lf", &letice[i][j]);
		}
	}

	fgets(temp, 300, in);
	//cout << temp << endl;
	fgets(temp, 300, in);
	//cout << temp << endl;
	for (i = 0; i < num; i++)
	{

		fscanf(in, "%d", &type[i]);
		fscanf(in, "%lf", &p[i][0]);
		fscanf(in, "%lf", &p[i][1]);
		fscanf(in, "%lf", &p[i][2]);
		fscanf(in, "%lf", &ridus[i]);
		fgets(temp, 300, in);
	}
	//int x_xishu = 0;
	//int y_xishu = 0;
	//int z_zishu = 0;

	//增加记录type的记录
	int temp_save[120] = { 0 };
	type_num = 0;
	for (i = 0; i < num; i++)
	{
		temp_save[type[i]]++;
	}
	for (i = 0; i < 120; i++)
	{
		if (temp_save[i] != 0)
			type_num++;
	}
	type_save = new int[type_num];
	j = 0;
	for (i = 0; i < 120; i++)
	{
		if (temp_save[i] != 0)
			type_save[j++] = i;
	}


	for (i = 0; i < yanshen; i++)
	{

		for (j = 0; j < num; j++)
		{
			//x_xishu = i/3;
			//y_xishu = i / 3;
			//z_zishu = (i % 9) / 3;

			real_position[i][j][0] = letice[0][0] * p[j][0] + letice[1][0] * p[j][1] + letice[2][0] * p[j][2] + (i % cengshu - ((cengshu - 1) / 2)) * letice[0][0] + (i % (cengshu * cengshu) / cengshu - ((cengshu - 1) / 2)) * letice[1][0] - (i / (cengshu * cengshu) - ((cengshu - 1) / 2)) * letice[2][0];
			real_position[i][j][1] = letice[0][1] * p[j][0] + letice[1][1] * p[j][1] + letice[2][1] * p[j][2] + (i % cengshu - ((cengshu - 1) / 2)) * letice[0][1] + (i % (cengshu * cengshu) / cengshu - ((cengshu - 1) / 2)) * letice[1][1] - (i / (cengshu * cengshu) - ((cengshu - 1) / 2)) * letice[2][1];
			real_position[i][j][2] = letice[0][2] * p[j][0] + letice[1][2] * p[j][1] + letice[2][2] * p[j][2] + (i % cengshu - ((cengshu - 1) / 2)) * letice[0][2] + (i % (cengshu * cengshu) / cengshu - ((cengshu - 1) / 2)) * letice[1][2] - (i / (cengshu * cengshu) - ((cengshu - 1) / 2)) * letice[2][2];
		}
	}
	for (i = 0; i < yanshen; i++)
	{
		for (j = 0; j < num; j++)
		{
			p_real[i][j][0] = (i % cengshu - ((cengshu - 1) / 2)) + p[j][0];
			p_real[i][j][1] = (i % (cengshu * cengshu) / cengshu - ((cengshu - 1) / 2)) + p[j][1];
			p_real[i][j][2] = -(i / (cengshu * cengshu) - ((cengshu - 1) / 2)) + p[j][2];
		}
	}

	fclose(in);
}

class save
{

public:
	int sunbets; //����Ƿֳɵ���ͨ��������
	int **save_list;
	save(int num);
	~save(void);
};
save::save(int num)
{
	int i, j, k;
	sunbets = 0;
	save_list = new int *[num];
	for (i = 0; i < num; i++)
	{
		save_list[i] = new int[num];
	}
	for (i = 0; i < num; i++)
	{
		for (j = 0; j < num; j++)
		{
			save_list[i][j] = -1;
		}
	}
}
save::~save(void) //�������ǲ��걸�ģ�˼����ô���������ص�����ֵ
{
	int i = 0, j = 0;
	delete[] save_list;
	//cout << "the object is being destroied!" << endl;
}

class element
{
public:
	int atomic_num;//atomic number
	double vdw_radius_min;
	double vdw_radius_max;
	//van der waals radius,min and max
	double cov_radius;//covalence radius
	int num_met_radius;  // number of metallic 


};

void full_expand_generate_graph(cell cell_a, int **expand_matrix, int edge_num,double***);
void ex_find_connect(int edge_num, cell cell_a, int *visited, int **matrix, save *save_a);
void ex_DFS(int a, int *ex_visited, cell cell_a, int **expand_matrix, int num_edge, int subets);
void generate_testfile(save *save_a, cell cell_a, int edge_num);
int *judge_the_2d(cell cell_a, int *save_list, int edge_num, int *judge);
void generate_outfile(cell cell_a, save *ex_pt, int edge_num, int *judge);
int generate_atom(save *ex_pt, int *judge, int edge_num, cell cell_a, string name, int *list);
string int2string(int i);
double face_point_dis(double letice[2][3], double in_face_point[3], double point[3]);
double get_vacumlayer_thick(int *list, int edge_num, int orig, cell cell_a, double **letice, int *copy);
int judge_falg_file(char* path);
static char wenjian[100]="100390.config";
static char result_path[40]="100001";
static char style[100]="rule_noorg_for_pc.txt";
static char wenjian_real[100];
const double ridus_plus_factor = 2.0;
int main(int argc, const char *argv[])
{
	//��������Ҫ�����ĳ�linux�µİ汾��
	int i = 0, j = 0,k=0;
	double temp = 0.0;		
	for (i = 0; i < argc; i++)
	{
		//if (i == 1)
		//{
		//	for (j = 0; j < strlen(argv[i]); j++)
		//	{
		//		style[j] = argv[i][j];
		//		//cout<<argv[i][j]<<endl;
		//		//cout<<wenjian[j]<<endl;
		//	}
		//	style[strlen(argv[i])] = '\0';
		//}

		if (i == 1)
		{
			for (j = 0; j < strlen(argv[i]); j++)
			{
				wenjian[j] = argv[i][j];
				//cout<<argv[i][j]<<endl;
				//cout<<wenjian[j]<<endl;
			}
			wenjian[strlen(argv[i])] = '\0';
		}
		if (i == 2)
		{
			for (j = 0; j < strlen(argv[i]); j++)
			{
				result_path[j] = argv[i][j];
			}
			result_path[strlen(argv[i])] = '\0';
		}

		/*if (i == 4)
		{
			for (j = 0; j < strlen(argv[i]); j++)
			{
				wenjian_real[j] = argv[i][j];
			}
			wenjian_real[strlen(argv[i])] = '\0';
		}*/
	}	
	//new_get_style(style);
	//测试修改：对于距离的判断我们采用半径方法
	char path_flag[200];
	strcpy(path_flag, wenjian);
	strcat(wenjian, "/atom.config");
	cell cell_a(wenjian);
	/*for (i = 0; i < cell_a.num; i++)
	{
		cout << i<<","<<cell_a.ridus[i] << endl;
	}
	cin.get();*/
	if (cell_a.num == 0 || judge_falg_file(path_flag)==0)
		return 0;

	//建立有关成键可能性的记录
	//第一个数是
	double*** possiblty = new double**[cell_a.type_num];
	for (i = 0; i < cell_a.type_num; i++)
	{
		possiblty[i] = new double*[cell_a.type_num];
		for (j = 0; j < cell_a.type_num; j++)
		{
			possiblty[i][j] = new double[3];
			possiblty[i][j][0] = -999;
			possiblty[i][j][1] = 100;
		}		
	}

	//test p_real

	/*cout << "start test the p_real!" << endl;
	for (i = 0; i < yanshen; i++)
	{
		for (j = 0; j < cell_a.num; j++)
		{
			cout << cell_a.p_real[i][j][0] << "," << cell_a.p_real[i][j][1] << "," << cell_a.p_real[i][j][2];
		}
		cout << endl;
	}	*/
	//还是先跳过有机结构
	int org_flag = 0;
	for (i = 0; i < cell_a.num; i++)
	{
		if (cell_a.type[i] == 6)
		{
			for (j = i + 1; j < cell_a.num; j++)
			{
				if (cell_a.type[j] == 1)
				{
					org_flag = 1;
					break;
				}
			}
			if (j == cell_a.num)
				break;
			if (org_flag == 1)
				break;
		}
		if (cell_a.type[i] == 1)
		{
			for (j = i + 1; j < cell_a.num; j++)
			{
				if (cell_a.type[j] == 6)
				{
					org_flag = 1; 
					break;
				}
			}
			if (j == cell_a.num)
				break;
			if (org_flag == 1)
				break;
		}
	}
	if (org_flag == 1)
		return 0;
	i = 0;
	int exit_falg = 1; //������֤�ǲ��ǺϽ�
	while (i++ < cell_a.num)
	{
		for (j = 0; j < 88; j++)
		{
			if (cell_a.type[i] == meatal_xuhao[j])
			{
				break;
			}
		}
		if (j == 88)
		{
			exit_falg = 0;
			break;
		}
	}
	if (exit_falg == 1)
	{
		//cout << "this is meatal!" << endl;
		//cin.get();
		return 0;
	}

	int num = cell_a.num * yanshen;
	save save_a(num);
	save *pt = &save_a; //ָ��ָ��save_a
	//cout << "start the expand graph processdure!" << endl;
	int *ex_visited; //�����Ƿ�������������
	ex_visited = new int[num];
	int **expand_graph = new int *[num]; //��ά���󣬴��������֮�����ϵ
	for (i = 0; i < num; i++)
	{
		expand_graph[i] = new int[num];
		ex_visited[i] = 0;
	}
	for (i = 0; i < num; i++)
	{
		for (j = 0; j < num; j++)
		{
			expand_graph[i][j] = -1;
		}
	}
	//��ʼ������ͨ�Ե��ж�
	full_expand_generate_graph(cell_a, expand_graph, num,possiblty);
	//补充完最后一项,最后记录
	int poss_save[2];
	double min_temp = 200;
	for (i = 0; i < cell_a.type_num; i++)
	{
		for (j = 0; j < cell_a.type_num; j++)
		{
			if (i != j)
			{
				possiblty[i][j][2] = abs(possiblty[i][j][0] - possiblty[i][j][1]);
				if (possiblty[i][j][2] < min_temp)
				{
					min_temp = possiblty[i][j][2];
					poss_save[0] = i;
					poss_save[1] = j;
				}
			}
		}
	}
	











	ex_find_connect(num, cell_a, ex_visited, expand_graph, pt);
	//cout << "the subnets is :" << pt->sunbets << endl;
	int success_ful = 0;
	if (save_a.sunbets == 1)
	{
		//cout << "this atom do not have the graph_connection!";
		//generate_testfile(pt, cell_a, num);
		delete[] ex_visited;
		for (i = 0; i < num; i++)
		{
			delete[] expand_graph[i];
		}
		delete[] expand_graph;
		//cout << "all work done" << endl;
		//cin.get();
		return success_ful;
	}
	//�õ�����ͨ����Ϊ1˵���϶��ǲ������������ˣ�ֱ��������
	//�������ÿ�ֵõ�����ͨ���������ж�

	else
	{
		//cout << "this atom has :" << pt->sunbets << " graph_connection!,please check!" << endl;
		int **judge = new int *[pt->sunbets];
		for (int i = 0; i < pt->sunbets; i++)
		{
			judge[i] = new int[2];
		}
		//现在建立一个标志位用来标志二维的连通分量是不是两头的，0表示起始，1表示中间，2表示末尾，-1表示不是2d
		//我们现在做这样的事情，建立二维连通分量只针对中间的来建立
		int * middle_flag = new int[pt->sunbets];
		for (i = 0; i < pt->sunbets; i++)
		{
			middle_flag[i] = 0;
		}
		//���������ͨ�����������жϽ��
		//cout << "show the connections xuhao!" << endl;
		//for (i = 0; i < pt->sunbets; i++)
		//{
		//	int num_point = 0;
		//	for (j = 0; j < num; j++)
		//	{
		//		if (pt->save_list[i][j] != -1)
		//		{
		//			num_point++;
		//			//cout << pt->save_list[i][j]<<",";
		//		}
		//	}	
		//	cout << num_point << endl;
		//	cout << endl;
		//}		
		//cin.get();
		for (j = 0; j < pt->sunbets; j++)
		{

			//cout << "now it's judge the :" << j << endl;
			judge[j][0] = 0;
			judge[j] = judge_the_2d(cell_a, pt->save_list[j], num, judge[j]);
			if (judge[j][0] == 2)
			{		
				
				success_ful = 1;
			}
			else
			{
				middle_flag[j] = -1;
			}
			//cout << "the judge result is :" << judge[j][0] << judge[j][1] << endl;
		}
		if (success_ful == 0)
			return 0;

		//这里我们重写一下关于两边的判断
		int twod_num = 0;//先看一下有多少个2d的连通分量		
		j = 0;

		for (i = 0; i < pt->sunbets; i++)
		{
			if (judge[i][0] == 2)
			{
				twod_num++;
			}
		}
		if (twod_num == pt->sunbets)//只有在全2d的情况下才决定两边不要
		{
			//先建立好该需要的空间
			double *angle_face = new double[pt->sunbets];//储存每个连通分量与其他分量角度的最小值
			int* xuhao_face = new int[pt->sunbets];//储存对应序号
			double **temp_angel_save = new double*[2];//临时储存用来计算角度的中间变量
			for (i = 0; i < 2; i++)
			{
				temp_angel_save[i] = new double[3];
			}			
			int*xuhao_a = new int[pt->sunbets];//临时储存序号
			for (i = 0; i < pt->sunbets; i++)
			{
				xuhao_a[i] = i;
			}
			int**random_choose = new int*[pt->sunbets*(pt->sunbets - 1) / 2];//储存任选两个的序号
			for (i = 0; i < pt->sunbets*(pt->sunbets - 1) / 2; i++)
			{
				random_choose[i] = new int[2];
			}
			random_choose_two(xuhao_a, pt->sunbets, random_choose);

			for (i = 0; i < pt->sunbets; i++)//每次需要找到这堆角度的最大值
			{
				//每次我需要记录三个序号
				int now_xuhao =pt->save_list[i][0];;
				int other1_xuhao;
				int oth2_xuhao;
				double temp_min_angel = -2;
				double every_angel;
				for (j = 0; j < pt->sunbets*(pt->sunbets - 1) / 2; j++)
				{
					if (random_choose[j][0] != i && random_choose[j][1] != i)//取出两个另外的来
					{
						other1_xuhao = pt->save_list[random_choose[j][0]][0];
						oth2_xuhao = pt->save_list[random_choose[j][1]][0];
						//cout << random_choose[j][0] << "," << random_choose[j][1] << endl;
					}
					else
						continue;
					//选出来之后计算角度
					for (k = 0; k < 3; k++)
					{
						temp_angel_save[0][k] = cell_a.real_position[other1_xuhao/cell_a.num][other1_xuhao%cell_a.num][k] - cell_a.real_position[now_xuhao/cell_a.num][now_xuhao%cell_a.num][k];
						temp_angel_save[1][k] = cell_a.real_position[oth2_xuhao/cell_a.num][oth2_xuhao%cell_a.num][k] - cell_a.real_position[now_xuhao / cell_a.num][now_xuhao%cell_a.num][k];
					}
					every_angel = vector_angle(temp_angel_save);
					//cout << every_angel << endl;
					if (every_angel > temp_min_angel)//需要的是最大值
					{
						temp_min_angel = every_angel;
					}
				}
				angle_face[i] = temp_min_angel;
				temp_min_angel = -2;
				//cout << angle_face[i]<< endl;
				
				
			}

			//在最大值里面找到最小值
			buble_plus(angle_face, xuhao_face, pt->sunbets);
			/*cout << angle_face[0] << angle_face[1] << endl;
			cout << xuhao_face[0] << xuhao_face[1] << endl;
			middle_flag[xuhao_face[0]] = 0;
			middle_flag[xuhao_face[1]] = 2;*/
			for (i = 2; i < pt->sunbets; i++)
			{
				middle_flag[xuhao_face[i]] = 1;
			}


			for (i = 0; i < pt->sunbets*(pt->sunbets - 1) / 2; i++)
			{
				delete[]random_choose[i];
			}
			delete[]random_choose;
			delete[]xuhao_a;			
			for (i = 0; i < 2; i++)
			{
				delete[]temp_angel_save[i];
			}
			delete[]temp_angel_save;
			delete[]angle_face;
			delete[]xuhao_face;
		}
		else
		{
			for (i = 0; i < pt->sunbets; i++)
			{
				middle_flag[i] = 1;
			}
		}

		//test middle falg
		/*for (i = 0; i < pt->sunbets; i++)
		{
			cout << middle_flag[i] << endl;
		}*/


		////然后需要输出层之间的最短距离。
		//double temp_layer_dis = 0;
		//double *layer_dis = new double[pt->sunbets];
		//string **show_ele = new string *[pt->sunbets];
		//double *duiying_style = new double[pt->sunbets];

		//for (i = 0; i < pt->sunbets; i++)
		//{
		//	show_ele[i] = new string[2];
		//}

		//for (int i = 0; i < pt->sunbets; i++)
		//{
		//	layer_dis[i] = 100.0;
		//}

		//for (j = 0; j < pt->sunbets; j++)
		//{
		//	if (judge[j][0] == 2 && middle_flag[j]==1)
		//	{
		//		for (i = 0; i < num && pt->save_list[j][i] != -1; i++)
		//		{
		//			//找到了每个2d连通分量的序号save_list[j][i]
		//			for (int m = 0; m < pt->sunbets; m++)
		//			{
		//				if (m != j && judge[m][0] == 2 && middle_flag[m] == 1)
		//				{
		//					for (int k = 0; k < num && pt->save_list[m][k] != -1; k++)
		//					{
		//						temp_layer_dis = dis(cell_a.real_position[pt->save_list[j][i] / cell_a.num][pt->save_list[j][i] % cell_a.num], cell_a.real_position[pt->save_list[m][k] / cell_a.num][pt->save_list[m][k] % cell_a.num]);
		//						if (temp_layer_dis < layer_dis[j] && temp_layer_dis > 1e-3)
		//						{
		//							layer_dis[j] = temp_layer_dis;
		//							show_ele[j][0] = a[cell_a.type[pt->save_list[j][i] % cell_a.num]];
		//							show_ele[j][1] = a[cell_a.type[pt->save_list[m][k] % cell_a.num]];
		//							duiying_style[j] = dist[cell_a.type[pt->save_list[j][i] % cell_a.num]][cell_a.type[pt->save_list[m][k] % cell_a.num]];
		//							//cout << "now the lauer:" << j << "is:" << temp_layer_dis;
		//						}
		//					}
		//				}
		//			}
		//		}
		//	}
		//}
		////到这里我们找到了最近距离
		ofstream fout;
		char dis_name[100] = ".txt";
		char dis_name_real[200];
		strcpy(dis_name_real, result_path);
		fout.open(strcat(dis_name_real, dis_name), ios::out | ios::app);
		if (!fout.is_open())
		{
			cout << "filed to open the file!" << endl;
			cout << "now the file path is :" << dis_name << endl;
			cin.get();
		}

		//然后开始根据我们的判断结果进行生成新的atom.config文件了
		//如果根据连通性获得了多个，我们不是记录了最多的，同时分别生成对应的2d文件
		

		int *generate_flag = new int[pt->sunbets]; //标记是不是需要写进文件0的话不写
		for (i = 0; i < pt->sunbets; i++)
		{
			generate_flag[i] = 1;
		}
		success_ful = 0;
		//fout << wenjian_real << "\n";
		char two_filename[100];
		int towd_num = 0;
		char temp_num[20];
		for (j = 0; j < pt->sunbets; j++)
		{
			if (judge[j][0] == 2 && generate_flag[j] == 1 && middle_flag[j] == 1)
			{
				for (int i = 0; i < pt->sunbets; i++)
				{
					if (judge[j][1] == judge[i][1])
						generate_flag[i] = 0;
				}
				//two_filename = wenjian + '-' + int2string(towd_num) + "_2d.config";
				sprintf(temp_num, "%d", towd_num);
				strcpy(two_filename, result_path);
				//strcat(two_filename, wenjian);
				strcat(two_filename, "-");
				strcat(two_filename, temp_num);
				strcat(two_filename, "_2d.config");

				if (generate_atom(pt, judge[j], num, cell_a, two_filename, pt->save_list[j]) != 0)
				{
					success_ful = 1;
					//fout << show_ele[j][0] << "\t" << show_ele[j][1] << "\t" << layer_dis[j] << "\t" << duiying_style[j];
					//fout << "\t" << double(layer_dis[j] / duiying_style[j]) << endl;
					//fout << rule_max << "\t" << no_connect << "\t" << abs(rule_max - no_connect) << endl;
					fout << possiblty[poss_save[0]][poss_save[1]][0] << "\t" << possiblty[poss_save[0]][poss_save[1]][1] << "\t" << possiblty[poss_save[0]][poss_save[1]][2] << endl;
					//break;
				}
				towd_num++;
			}
		}

		//generate_outfile(cell_a, pt, num, judge);
		//fout << endl;
		fout.close();

		for (int i = 0; i < pt->sunbets; i++)
		{
			delete[] judge[i];
			//delete[] show_ele[i];
		}
		//delete[] show_ele;
		delete[] judge;
		delete[] ex_visited;
		for (i = 0; i < num; i++)
		{
			delete[] expand_graph[i];
		}
		delete[] expand_graph;
		//delete[]subnets_num;
		delete[] generate_flag;
		//delete[] layer_dis;
		//delete[] duiying_style;
		delete[]middle_flag;
		//cin.get();
		return success_ful;
	}

	//cout << "all total work done" << endl;
	//cin.get();
	
}

void full_expand_generate_graph(cell cell_a, int **expand_matrix, int num,double*** possi) //�����ǿվ���Ͷ�Ӧ�㣬��������õľ���
{
	int i = 0, j = 0;
	int m = 0, n = 0;
	double temp = 0.0;
	for (i = 0; i < num; i++)
	{
		for (j = i; j < num; j++)
		{
			//cout << "for j,this is the :" << j / cell_a.num << "matrix ," << endl;
			//cout << "the xuhao is :" << j % cell_a.num << endl;

			temp = dis(cell_a.real_position[i / cell_a.num][i % cell_a.num], cell_a.real_position[j / cell_a.num][j % cell_a.num]);
			for (m = 0; m < cell_a.type_num; m++)
			{
				if (cell_a.type[i%cell_a.num]==cell_a.type_save[m])
				{
					break;
				}
			}
			for (n= 0; n < cell_a.type_num; n++)
			{
				if (cell_a.type[j%cell_a.num] == cell_a.type_save[n])
				{
					break;
				}
			}



			/*if (temp != 0.0 && temp < dist[cell_a.type[i % cell_a.num]][cell_a.type[j % cell_a.num]])*/
			if (temp != 0.0 && temp < (ridus_plus_factor*(cell_a.ridus[i%cell_a.num]+cell_a.ridus[j%cell_a.num])/1.2)  && cell_a.type[i%cell_a.num]!=cell_a.type[j%cell_a.num])
			{
				//if (cell_a.type[i%cell_a.num] != cell_a.type[j%cell_a.num])
				//{
				//	//cout << "for " << cell_a.type[i%cell_a.num] << "," << cell_a.type[i%cell_a.num] << endl;
				//	//cout << temp << "rule is:" << ridus_plus_factor * (cell_a.ridus[i%cell_a.num] + cell_a.ridus[j%cell_a.num]) << endl;

				//}
				if (temp > possi[m][n][0])
				{
					possi[n][m][0] =possi[m][n][0]=temp;
				}
				expand_matrix[i][j] = 1;
				expand_matrix[j][i] = 1;
			}
			else
			{
				if (temp < possi[m][n][1]&&abs(temp)>1e-3)
				{
					possi[n][m][1] =possi[m][n][1] = temp;
				}
				expand_matrix[i][j] = 0;
				expand_matrix[j][i] = 0;
			}
		}
	}

	return;
}

void ex_find_connect(int edge_num, cell cell_a, int *visited, int **matrix, save *save_a)
{
	int i = 0, ii = 0, j = 0, k = 0, jj = 0;
	int save_flag = 1;
	int subnets = 0;
	for (k = 0; k < edge_num; k++)
	{
		if (!visited[k])
		{
			subnets++;
			ex_DFS(k, visited, cell_a, matrix, edge_num, subnets);

			//cout << "the subnets is:" << subnets << endl;
			save_a->sunbets = subnets;
			for (i = 0; i < edge_num; i++)
			{
				if (visited[i] == subnets)
				{
					save_a->save_list[subnets - 1][ii] = i;
					ii++;
				}
			}
		}
		ii = 0;
	}
	//cout << "the hanshu finished before,the subnets is:" << save_a->sunbets << endl;
	//cin.get();
	return;
}

void generate_testfile(save *save_a, cell cell_a, int edge_num) //������save*ָ�룬����ǽ����txt
{
	int i = 0, j = 0;
	FILE *out;
	out = fopen("test_output.txt", "wt");
	char a[50] = "the lian tong num is:";
	char b[20];
	sprintf(b, "%d\t", save_a->sunbets);
	strcat(a, b);
	fprintf(out, "%s\t", a);
	cout << "i got the liantong geshu is:" << save_a->sunbets << endl;
	for (i = 0; i < save_a->sunbets; i++)
	{
		for (j = 0; j < edge_num; j++)
		{
			if (save_a->save_list[i][j] != -1)
			{
				fprintf(out, "%d\t", save_a->save_list[i][j]);
			}
		}
		fprintf(out, "\n");
	}
	fclose(out);
	return;
}

void generate_outfile(cell cell_a, save *ex_pt, int edge_num, int *judge) ////��������������ǻ������ͨ����֮�󣬸����жϵĽ������txt�ļ����ڹ۲�
{
	int i, j;
	FILE *out;
	out = fopen("out_judge.txt", "w");
	char a[50] = "the sunbets num is:";
	char b[20];
	sprintf(b, "%d\t", ex_pt->sunbets);
	strcat(a, b);
	fprintf(out, "%s\n", a);
	cout << "i got the liantong geshu is:" << ex_pt->sunbets << endl;
	char c[60] = "for subnets,judge out:";
	fprintf(out, "%s\n", c);
	for (i = 0; i < ex_pt->sunbets; i++)
	{
		fprintf(out, "%d\t", judge[i]);
	}
	fprintf(out, "\n");
	for (i = 0; i < ex_pt->sunbets; i++)
	{
		for (j = 0; j < edge_num; j++)
		{
			if (ex_pt->save_list[ex_pt->sunbets - 1][j] != -1)
			{
				fprintf(out, "%d\t", ex_pt->save_list[i][j]);
			}
		}
		fprintf(out, "\n");
	}

	fclose(out);
	return;
}

void ex_DFS(int a, int *ex_visited, cell cell_a, int **expand_matrix, int num_edge, int subnets) //�������������ʵ�ִӶ���a����������������δ���ʹ����ڽӽڵ�
{
	ex_visited[a] = subnets;
	//printf("now i start the dot %d\n", a);
	int i = 0, j = 0, k = 0;
	for (j = 0; j < num_edge; j++)
	{
		if (expand_matrix[a][j] == 1 && !ex_visited[j])
		{
			ex_DFS(j, ex_visited, cell_a, expand_matrix, num_edge, subnets);
			//����д���˺������,��ʵ�ܼ򵥣����Ǽ�¼���ok��
		}
	}
}

int *judge_the_2d(cell cell_a, int *save_list, int edge_num, int *judge) //���������������������õ���ÿ����ͨ�������ж����ǲ�������Һõ�2d��ͨ����
{
	int i = 0, j = 0, k = 0, m = 0;
	int *pt = judge;
	int exit_flag = 0; //�˳���־��1��ʾ1d��2��ʾ2d��3��ʾ3d,0��ʾ���ڵ�̫���ڷ�ɢ����޷���������жϵ�����
	int *copy = new int[yanshen];
	for (int i = 0; i < yanshen; i++)
	{
		copy[i] = -1;
	}
	//�������ʵ�copy
	int temp = 0;
	for (i = 0; i < edge_num; i++)
	{
		if (save_list[i] != -1)
		{
			temp = save_list[i] % cell_a.num;
			for (k = 0; k < edge_num; k++)
			{
				if (save_list[k] == -1)
					break;
				if ((save_list[k]) % cell_a.num == temp)
				{
					copy[j] = save_list[k];
					j++;
				}
			}
			if (j > 3)
				break;
			else
				j = 0;
		}
		else
			break;
	}

	if (j < 4 || j > cengshu * cengshu)
	{
		//cout << "i can not find the enough vertex to complete the judgement" << endl;
		*pt = 0;
		*(pt + 1) = -1;
		delete[] copy;
		return pt;
	}

	if (j == cengshu * cengshu)
	{
		*pt = 2;
		*(pt + 1) = temp;
		delete[] copy;
		return pt;
	}

	//cout << "i need to judge the ppints_num is:" << j << endl;
	//cin.get();
	//debug
	/*for (i = 0; i < j; i++)
	{
		cout << copy[i] << "," << endl;
	}*/
	//�ҵ��˷��������ĸ������������

	//j���Ǹ�����ĸ���

	//�������жϹ�������
	exit_flag = 1;
	double **a = new double *[3];
	for (i = 0; i < 3; i++)
	{
		a[i] = new double[3];
	}
	a[0][0] = a[0][1] = a[0][2] = 1;

	double x = cell_a.real_position[copy[1] / cell_a.num][copy[1] % cell_a.num][0] - cell_a.real_position[copy[0] / cell_a.num][copy[0] % cell_a.num][0];
	double y = cell_a.real_position[copy[1] / cell_a.num][copy[1] % cell_a.num][1] - cell_a.real_position[copy[0] / cell_a.num][copy[0] % cell_a.num][1];
	double z = cell_a.real_position[copy[1] / cell_a.num][copy[1] % cell_a.num][2] - cell_a.real_position[copy[0] / cell_a.num][copy[0] % cell_a.num][2];
	a[1][0] = x;
	a[1][1] = y;
	a[1][2] = z;

	double nx = 0;
	double ny = 0;
	double nz = 0;
	for (i = 2; i < j; i++)
	{
		if (save_list[i] != -1)
		{
			nx = cell_a.real_position[copy[i] / cell_a.num][copy[i] % cell_a.num][0] - cell_a.real_position[copy[0] / cell_a.num][copy[0] % cell_a.num][0];
			ny = cell_a.real_position[copy[i] / cell_a.num][copy[i] % cell_a.num][1] - cell_a.real_position[copy[0] / cell_a.num][copy[0] % cell_a.num][1];
			nz = cell_a.real_position[copy[i] / cell_a.num][copy[i] % cell_a.num][2] - cell_a.real_position[copy[0] / cell_a.num][copy[0] % cell_a.num][2];

			a[2][0] = nx;
			a[2][1] = ny;
			a[2][2] = nz;
			if (abs(three_jie_chaji(a)) > 1e-3) //��һ�������ߵ�����Ϳ���
			{
				exit_flag = 0;
				break;
			}
		}
	}

	for (i = 0; i < 3; i++)
	{
		delete[] a[i];
	}
	delete[] a;
	if (exit_flag == 1)
	{
		//cout << "the 1d situation" << endl;
		*pt = 1;
		*(pt + 1) = temp;
		delete[] copy;
		return pt;
	}

	//Ȼ�����Ƿ��ĵ㹲�����
	exit_flag = 2;
	double **b = new double *[4];
	for (i = 0; i < 4; i++)
	{
		b[i] = new double[3];
	}
	int number = j; //�ȰѸ���������
	int *flag = new int[number];
	for (i = 0; i < number; i++)
	{
		flag[i] = 0;
	}

	int jj = 0;
	for (i = 0; i < number; i++)
	{
		flag[i] = 1;
		for (j = i + 1; j < number; j++)
		{
			flag[j] = 1;
			for (k = j + 1; k < number; k++)
			{
				flag[k] = 1;
				for (m = k + 1; m < number; m++)
				{
					flag[m] = 1;

					for (int mm = 0; mm < number; mm++) //�����ͻ�õ�4��������
					{
						if (flag[mm] == 1)
						{
							for (int ii = 0; ii < 3; ii++)
							{
								b[jj][ii] = cell_a.real_position[copy[mm] / cell_a.num][copy[mm] % cell_a.num][ii];
							}
							jj++;
							if (jj == 4)
								break;
						}
					}
					jj = 0;
					//������д�жϹ���,���淵��1
					//ֻҪ���ַǹ������������3d���
					if (inone_face(b) == 0)
					{
						exit_flag = 3;
						for (int i = 0; i < 4; i++)
						{
							delete[] b[i];
						}
						delete[] b;
						delete[] copy;
						delete[] flag;
						*pt = 3;
						*(pt + 1) = temp;
						return pt;
					}

					flag[m] = 0;
				}
				flag[k] = 0;
			}

			flag[j] = 0;
		}
		flag[i] = 0;
	}

	*pt = 2;
	//���������Ǿ͵õ��������2d�����;
	*(pt + 1) = temp;
	delete[] copy;
	delete[] flag;
	return pt;
}


int generate_atom(save *ex_pt, int *judge, int edge_num, cell cell_a, string name, int *list)
{

	//这个的作用是根据得到的连通信息生成atom文件,输入是存储的连通分量，输出是atom文件

	int i = 0, j = 0;
	int k = 0;
	int exit_flag = 0;
	int orig = judge[1];
	//cout << "这里的对应原始点是：" << orig << endl;
	string file[2];
	int num = cell_a.num;
	file[0] = "Lattice vector";
	file[1] = "Position";
	//首先开始生成abc的三个基矢分量
	//从原型序号中找到这样的点以及这样的距离

	double **letice = new double *[3];
	for (i = 0; i < 3; i++)
	{
		letice[i] = new double[3]; //建立letice数组用来储存
	}
	double **letice_r = new double *[3];
	for (i = 0; i < 3; i++)
	{
		letice_r[i] = new double[3];
	}
	for (i = 0; i < edge_num; i++)
	{
		if (list[i] == -1)
			break;
		if (list[i] % cell_a.num == orig)
		{
			j++;
		}
	}
	int copy_num = j;
	int *copy = new int[j];
	for (i = 0; i < edge_num; i++)
	{
		if (list[i] == -1)
			break;
		if (list[i] % cell_a.num == orig)
		{
			copy[k++] = list[i];
		}
	}
	
	
	//有j个这样的复制体,并且都找到了，并且也判断了哪些原子需要放进去
	int *atom_flag_a = new int[cell_a.num];
		//这里的atom_falg_a针对的是第一种情况的二维材料
	for (i = 0; i < cell_a.num; i++)
	{
		atom_flag_a[i] = 0;
	}
	k = 0;
	for (i = 0; i < edge_num; i++)
	{
		if (list[i] == -1)
			break;
		for (j = 0; j < cell_a.num; j++)
		{
			if (list[i] % cell_a.num == j && atom_flag_a[j] == 0)
			{
				atom_flag_a[j] = 1;
				k++;
				break;
			}
		}
	}
	int atom_num_write = k;	
	double **fenshu = new double *[atom_num_write];
	for (i = 0; i < k; i++)
	{
		fenshu[i] = new double[3];
	}
	double ** fenshu_plus = new double*[edge_num];
	for (i = 0; i < edge_num; i++)
	{
		fenshu_plus[i] = new double[3];
	}
	int atom_plus = 0;
	

	// 然后根据两种情况进行分类
	if (copy_num == cengshu * cengshu)
	{	
	
		int m = 0;
		for (i = 0; i < 3; i++)
		{
			letice[0][i] = cell_a.real_position[copy[1] / cell_a.num][copy[1] % cell_a.num][i] - cell_a.real_position[copy[0] / cell_a.num][copy[0] % cell_a.num][i];
			letice[1][i] = cell_a.real_position[copy[copy_num / cengshu] / cell_a.num][copy[copy_num / cengshu] % cell_a.num][i] - cell_a.real_position[copy[0] / cell_a.num][copy[0] % cell_a.num][i];
		}
		letice[2][0] = letice[0][1] * letice[1][2] - letice[0][2] * letice[1][1];
		letice[2][1] = letice[0][2] * letice[1][0] - letice[0][0] * letice[1][2];
		letice[2][2] = letice[0][0] * letice[1][1] - letice[0][1] * letice[1][0];
		double plus = get_vacumlayer_thick(list, edge_num, orig, cell_a, letice, copy);

		//然后将c的模长度拉长，保证真空层是严格的25的距离
		double pingfanghe = pow(letice[2][0], 2) + pow(letice[2][1], 2) + pow(letice[2][2], 2);
		double a = pow((pow(25 + 2 * plus, 2) / pingfanghe), 0.5);
		for (i = 0; i < 3; i++)
		{
			letice[2][i] = a * letice[2][i];
		}
		/*double yuxian = (letice[0][0] * letice[1][0] + letice[0][1] * letice[1][1] + letice[0][2] * letice[1][2]) / (pow(letice[0][0] * letice[0][0] + letice[0][1] * letice[0][1] + letice[0][2] * letice[0][2], 0.5)*pow(letice[1][0] * letice[1][0] + letice[1][1] * letice[1][1] + letice[1][2] * letice[1][2], 0.5));
		double theta = acos(yuxian);
		double mo[3];
		for (i = 0; i < 3; i++)
		{
			mo[i] = pow(letice[i][0] * letice[i][0] + letice[i][1] * letice[i][1] + letice[i][2] * letice[i][2], 0.5);
		}

		for (i = 0; i < 3; i++)
		{
			if (i == 0)
			{
				letice[i][0] = mo[0];
				letice[i][1] = 0;
				letice[i][2] = 0;
			}
			if (i == 1)
			{
				letice[i][0] = mo[1] * cos(theta);
				letice[i][1] = mo[1] * sin(theta);
				letice[i][2] = 0;
			}
			if (i == 2)
			{
				letice[i][0] = letice[i][1] = 0;
				letice[i][2] = mo[2];
			}
		}*/
		if (reverse_matrix(letice, 3, letice_r) == 0)
			return 0; //建立了逆矩阵
		//这里我们需要做一步，就是将产生的矩阵做成是尽量0多的


		//到这里完成了基矢的获取，下面是确定需要放哪个原子，以及新原子的位置是哪里
		//这里需要注意，我们首先需要知道偏移了多少，然后才能根据偏移方向针对性的写出原子真实坐标

		//首先进行原子平移,选出c最小的原子进行平移
		j = 0;

		//改变之前的做法，这里让即使横平竖直的情况，坐标变换也是基于筛选
		for (i = 0; i < edge_num; i++)
		{
			if (list[i] == -1)
				break;
			for (m = 0; m < 3; m++)
			{
				fenshu_plus[i][m] = (cell_a.real_position[list[i] / cell_a.num][list[i] % cell_a.num][0] - cell_a.real_position[copy[0] / cell_a.num][copy[0] % cell_a.num][0])* letice_r[0][m] + (cell_a.real_position[list[i] / cell_a.num][list[i] % cell_a.num][1] - cell_a.real_position[copy[0] / cell_a.num][copy[0] % cell_a.num][1])* letice_r[1][m] + (cell_a.real_position[list[i] / cell_a.num][list[i] % cell_a.num][2] - cell_a.real_position[copy[0] / cell_a.num][copy[0] % cell_a.num][2]) * letice_r[2][m];

			}
			if (fenshu_plus[i][0] >= 1-1e-5 || fenshu_plus[i][0] < 0 || fenshu_plus[i][1] >= 1-1e-5 || fenshu_plus[i][1] < 0)
			{
				fenshu_plus[i][0] = -100;
			}
			if (fenshu_plus[i][0] != -100)
				atom_plus++;

		}

		//反变换之后开始进行c方向偏移
		double c_pingyi = 100;
		i = 0, j = 0;
		while (i < edge_num)
		{
			if (list[i] == -1)
				break;
			if (fenshu_plus[i][0] != -100)
			{
				if (fenshu_plus[i][2] < c_pingyi)
				{
					c_pingyi = fenshu_plus[i][2];
				}
			}
			i++;
		}

		for (i = 0; i < edge_num; i++)
		{
			if (list[i] == -1)
				break;
			if (fenshu_plus[i][0] != -100)
			{
				fenshu_plus[i][2] = fenshu_plus[i][2] - c_pingyi + 0.02;

			}
		}


	}

	else
	{
		//cout << "开始进入偏移的状态了哦，哈哈哈哈哈" << endl;
		//cout << "the file is:" << wenjian << endl;
		//   cin.get();
		
	
	//从这里要确认怎么生成基矢,这里生成基矢的原则是面积产生最小且最接近90°
		double** chaji_save = new double*[3];
		for (i = 0; i < 3; i++)
		{
			chaji_save[i] = new double[3];
		}
		for (i = 0; i < 3; i++)
		{
			chaji_save[1][i] = cell_a.real_position[copy[1] / cell_a.num][copy[1] % cell_a.num][i] - cell_a.real_position[copy[0] / cell_a.num][copy[0] % cell_a.num][i];
		}

		double jvli = 100;
		double temp_jvli;
		int * mianji_small = new int[cengshu*cengshu];
		for (i = 0; i < cengshu*cengshu; i++)
		{
			mianji_small[i] = -1;
		}

		for (i = 2; i < copy_num; i++)
		{
			for (j = 0; j < 3; j++)
			{
				chaji_save[2][j] = cell_a.real_position[copy[i] / cell_a.num][copy[i] % cell_a.num][j] - cell_a.real_position[copy[0] / cell_a.num][copy[0] % cell_a.num][j];
			}			
			temp_jvli = three_jie_chaji(chaji_save);
			if (abs(temp_jvli) > 1e-6 && temp_jvli < jvli)
			{
				jvli = temp_jvli;				
			}
		}

		for (i = 2,k=0; i < copy_num; i++)
		{
			for (j = 0; j < 3; j++)
			{
				chaji_save[2][j] = cell_a.real_position[copy[i] / cell_a.num][copy[i] % cell_a.num][j] - cell_a.real_position[copy[0] / cell_a.num][copy[0] % cell_a.num][j];
			}
			temp_jvli = three_jie_chaji(chaji_save);
			if (abs(temp_jvli) > 1e-6 && abs(temp_jvli-jvli)<1e-1)
			{
				mianji_small[k] = copy[i];
				k++;
			}
		}
		//先找到了最小的面积是多少
		//然后看最接近90°的情况
		double angle_cha = 100;
		double angle_cha_temp = 0;
		double **vector = new double*[2];
		for (i = 0; i < 2; i++)
		{
			vector[i] = new double[3];
		}
		for (i = 0; i < 3; i++)
		{
			vector[0][i]= cell_a.real_position[copy[1] / cell_a.num][copy[1] % cell_a.num][i] - cell_a.real_position[copy[0] / cell_a.num][copy[0] % cell_a.num][i];
		}


		int m = 0;
		for (i = 0; i < k; i++)
		{
			if (mianji_small[i] == -1)
				break;
			for (j = 0; j < 3; j++)
			{
				vector[1][j] = cell_a.real_position[mianji_small[i] / cell_a.num][mianji_small[i] % cell_a.num][j] - cell_a.real_position[copy[0] / cell_a.num][copy[0] % cell_a.num][j];
			}
			angle_cha_temp = abs(vector_angle(vector) - 90);
			if (angle_cha_temp < angle_cha)
			{
				angle_cha = angle_cha_temp;
				m = mianji_small[i];
			}

		}
		for (i = 0; i < 3; i++)
		{
			letice[0][i] = cell_a.real_position[copy[1] / cell_a.num][copy[1] % cell_a.num][i] - cell_a.real_position[copy[0] / cell_a.num][copy[0] % cell_a.num][i];
			letice[1][i] = cell_a.real_position[m / cell_a.num][m % cell_a.num][i] - cell_a.real_position[copy[0] / cell_a.num][copy[0] % cell_a.num][i];
		}
		letice[2][0] = letice[0][1] * letice[1][2] - letice[0][2] * letice[1][1];
		letice[2][1] = letice[0][2] * letice[1][0] - letice[0][0] * letice[1][2];
		letice[2][2] = letice[0][0] * letice[1][1] - letice[0][1] * letice[1][0];
		double plus = get_vacumlayer_thick(list, edge_num, orig, cell_a, letice, copy);

		//然后将c的模长度拉长，保证真空层是严格的25的距离
		double pingfanghe = pow(letice[2][0], 2) + pow(letice[2][1], 2) + pow(letice[2][2], 2);
		double a = pow((pow(25 + 2 * plus, 2) / pingfanghe), 0.5);
		for (i = 0; i < 3; i++)
		{
			letice[2][i] = a * letice[2][i];
		}
		/*for (i = 0; i < 3; i++)
		{
			cout<<letice[2][i]<<endl;
		}
*/
		/*double yuxian = (letice[0][0] * letice[1][0] + letice[0][1] * letice[1][1] + letice[0][2] * letice[1][2]) / (pow(letice[0][0] * letice[0][0] + letice[0][1] * letice[0][1] + letice[0][2] * letice[0][2], 0.5)*pow(letice[1][0] * letice[1][0] + letice[1][1] * letice[1][1] + letice[1][2] * letice[1][2], 0.5));
		double theta = acos(yuxian);
		double mo[3];
		for (i = 0; i < 3; i++)
		{
			mo[i] = pow(letice[i][0] * letice[i][0] + letice[i][1] * letice[i][1] + letice[i][2] * letice[i][2], 0.5);
		}

		for (i = 0; i < 3; i++)
		{
			if (i == 0)
			{
				letice[i][0] = mo[0];
				letice[i][1] = 0;
				letice[i][2] = 0;
			}
			if (i == 1)
			{
				letice[i][0] = mo[1] * cos(theta);
				letice[i][1] = mo[1] * sin(theta);
				letice[i][2] = 0;
			}
			if (i == 2)
			{
				letice[i][0] = letice[i][1] = 0;
				letice[i][2] = mo[2];
			}
		}*/
		if (reverse_matrix(letice, 3, letice_r) == 0)
			return 0; //建立了逆矩阵

		//现在遇到了问题，就是实际上放进去的原子个数应该是多余一个晶胞的，或者是不止一个晶胞个数，而应该要把包括的都放进去
		//这里我们这么做，先全部放进去。然后筛选
		for (i = 0; i < edge_num; i++)
		{
			if (list[i] == -1)
				break;
			for (m = 0; m < 3; m++)
			{
				fenshu_plus[i][m] = (cell_a.real_position[list[i] / cell_a.num][list[i]%cell_a.num][0]- cell_a.real_position[copy[0] / cell_a.num][copy[0] % cell_a.num][0] )* letice_r[0][m] + (cell_a.real_position[list[i] / cell_a.num][list[i] % cell_a.num][1]- cell_a.real_position[copy[0] / cell_a.num][copy[0] % cell_a.num][1])* letice_r[1][m] + (cell_a.real_position[list[i] / cell_a.num][list[i] % cell_a.num][2]- cell_a.real_position[copy[0] / cell_a.num][copy[0] % cell_a.num][2]) * letice_r[2][m];
							
			}
			if (fenshu_plus[i][0] >= (1-1e-5) || fenshu_plus[i][0] <0 || fenshu_plus[i][1] >= (1-1e-5) || fenshu_plus[i][1] < 0)
			{
				fenshu_plus[i][0] = -100;				
			}
			if (fenshu_plus[i][0] != -100)
				atom_plus++;

		}
		
		//反变换之后开始进行c方向偏移
		double c_pingyi = 100;
		i = 0, j = 0;
		while (i < edge_num)
		{
			if (list[i] == -1)
				break;
			if (fenshu_plus[i][0]!=-100)
			{
				if (fenshu_plus[i][2] < c_pingyi)
				{
					c_pingyi = fenshu_plus[i][2];
				}
			}
			i++;
		}

		for (i = 0; i < edge_num; i++)
		{
			if (list[i] == -1)
				break;
			if (fenshu_plus[i][0]!=-100)
			{
				fenshu_plus[i][2] = fenshu_plus[i][2] - c_pingyi + 0.02;
				
			}
		}

		
		for (i = 0; i < 3; i++)
		{
			delete[]chaji_save[i];
		}
		delete[]chaji_save;
		delete[]mianji_small;
		for (i = 0; i < 2; i++)
		{
			delete[]vector[i];
		}
		delete[]vector;
		//cin.get();
	}

	//然后我生成了新的坐标，开始写文件吧
	FILE *out = fopen(name.c_str(), "wt");
	char atom_head[20] = "atoms";

	if (copy_num == cengshu * cengshu)
	{
		fprintf(out, "%d\n", atom_plus);
	}
	else
	{
		fprintf(out, "%d\n", atom_plus);
	}
	




	//fprintf(out, "%s\n", atom_head);
	fprintf(out, "%s\n", file[0].c_str());
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			fprintf(out, "\t%.9lf\t", letice[i][j]);
		}
		fprintf(out, "\n");
	}
	fprintf(out, "%s\n", file[1].c_str());
	j = 0;

	

	if (copy_num == cengshu * cengshu)
	{
		for (i = 0; i < edge_num; i++)
		{
			if (list[i] == -1)
				break;
			int mm = 1;
			if (fenshu_plus[i][0] != -100)
			{
				fprintf(out, "  %d\t", cell_a.type[list[i] % cell_a.num]);
				for (int m = 0; m < 3; m++)
				{
					fprintf(out, "%.9lf\t", fenshu_plus[i][m]);
				}
				for (int m = 0; m < 3; m++)
				{
					fprintf(out, "%d\t", mm);
				}
				fprintf(out, "\n");

			}
		}
		fclose(out);
	}
	else
	{
		for (i = 0; i <edge_num; i++)
		{
			if (list[i] == -1)
				break;
			int mm = 1;
			if (fenshu_plus[i][0]!=-100)
			{
				fprintf(out, "  %d\t", cell_a.type[list[i]%cell_a.num]);
				for (int m = 0; m < 3; m++)
				{
					fprintf(out, "%.9lf\t", fenshu_plus[i][m]);
				}
				for (int m = 0; m < 3; m++)
				{
					fprintf(out, "%d\t", mm);
				}
				fprintf(out, "\n");
				
			}
		}
		fclose(out);
	}
	
	for (i = 0; i < atom_num_write; i++)
	{
		delete[] fenshu[i];
	}
	delete[] fenshu;
	for (i = 0; i < 3; i++)
	{
		delete[] letice_r[i];
	}
	delete[] letice_r;
	delete[] copy;
	for (i = 0; i < 3; i++)
	{
		delete[] letice[i];
	}
	delete[] letice;
	delete[] atom_flag_a;
	for (i = 0; i < edge_num; i++)
	{
		delete[]fenshu_plus[i];
	}
	delete[]fenshu_plus;
	return 1;
}

double get_vacumlayer_thick(int *list, int edge_num, int orig, cell cell_a, double **letice, int *copy)
{
	int i = 0, j = 0;
	int *dis_flag = new int[edge_num];
	for (int i = 0; i < edge_num; i++)
	{
		dis_flag[i] = 0;
	}
	//这个用来标记这个连通分量中哪些需要进行求距离的运算

	int *atom_flag = new int[cell_a.num];
	for (int i = 0; i < cell_a.num; i++)
	{
		atom_flag[i] = 0;
	} //标记哪些原子用到了
	for (int i = 0; i < edge_num; i++)
	{
		if (list[i] == -1)
			break;
		for (int j = 0; j < cell_a.num; j++)
		{
			if ((atom_flag[j] != 1) && (list[i] % cell_a.num == j))
			{
				atom_flag[j] = 1;
				break;
			}
		}
	}
	for (i = 0; i < edge_num; i++)
	{
		if (list[i] == -1)
			break;
		else
		{
			if (list[i] % cell_a.num != orig && atom_flag[list[i] % cell_a.num] == 1)
			{
				dis_flag[i] = 1;
				atom_flag[list[i] % cell_a.num]++;
			}
		}
	}
	//到这里完成了该计算哪些原子的点面距离了

	double plus = 0;
	double temp_dis = 0;
	//其实现在a和b就是两个向量，c就是法向量
	double xiangliang[2][3];
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			xiangliang[i][j] = letice[i][j];
		}
	}
	for (i = 0; i < edge_num; i++)
	{
		if (list[i] != -1 && dis_flag[i] == 1)
		{
			//从连通分量中找出点到面的距离
			temp_dis = face_point_dis(xiangliang, cell_a.real_position[copy[0] / cell_a.num][copy[0] % cell_a.num], cell_a.real_position[list[i] / cell_a.num][list[i] % cell_a.num]);
			if (temp_dis > plus)
			{
				plus = temp_dis;
			}
		}
	}
	delete[] dis_flag;
	delete[] atom_flag;
	return plus;
}
inline double  three_jie_chaji(double **a)
{
	int i = 0, j = 0;
	double temp = 0;
	temp = pow(pow((a[1][1] * a[2][2] - a[1][2] * a[2][1]),2) +pow( (a[1][2] * a[2][0] - a[1][0] * a[2][2]),2) + pow((a[1][0] * a[2][1] - a[1][1] * a[2][0]),2),0.5);
	return temp;
}

inline int inone_face(double **a) //���淵��1�������淵��0
{
	int i = 0, j = 0;
	double x1, x2, x3, y1, y2, y3, z1, z2, z3;
	x1 = a[0][0] - a[3][0];
	x2 = a[1][0] - a[3][0];
	x3 = a[2][0] - a[3][0];

	y1 = a[0][1] - a[3][1];
	y2 = a[1][1] - a[3][1];
	y3 = a[2][1] - a[3][1];

	z1 = a[0][2] - a[3][2];
	z2 = a[1][2] - a[3][2];
	z3 = a[2][2] - a[3][2];
	double k = 0;
	k = (x1 * y2 * z3) + (x2 * y3 * z1) + (x3 * y1 * z2) - (x3 * y2 * z1) - (y3 * z2 * x1) - (z3 * x2 * y1);

	if (-1e-1 < k && k < 1e-1)
		return 1;
	else
		return 0;
}

double det(double **a, int num)
{
	int i = 0, j = 0, k = 0;
	if (num == 1)
	{
		return a[0][0];
	}
	double ans = 0;
	double **temp_det = new double *[num];
	for (i = 0; i < num; i++)
	{
		temp_det[i] = new double[num];
	}

	for (i = 0; i < num; i++)
	{
		for (j = 0; j < num - 1; j++)
		{
			for (k = 0; k < num - 1; k++)
			{
				temp_det[j][k] = a[j + 1][(k >= i) ? k + 1 : k];
			}
		}
		double t = det(temp_det, num - 1);

		if (i % 2 == 0)
		{
			ans = ans + a[0][i] * t;
		}
		else
		{
			ans = ans - a[0][i] * t;
		}
	}
	return ans;
}

//����һ����������ĺ�����������ԭ��������Ƕ�Ӧλ�ñ�����ʽ���İ������
void getastar(double **a, int num, double **b)
{
	int i = 0, j = 0, k = 0, m = 0;
	if (num == 1)
	{
		b[0][0] = 1;
		return;
	}
	double **temp = new double *[num];
	for (i = 0; i < num; i++)
	{
		temp[i] = new double[num];
	}
	for (i = 0; i < num; i++)
	{
		for (j = 0; j < num; j++)
		{
			for (k = 0; k < num - 1; k++)
			{
				for (m = 0; m < num - 1; m++)
				{
					temp[k][m] = a[(k >= i) ? k + 1 : k][(m >= j) ? m + 1 : m];
				}
			}

			b[j][i] = det(temp, num - 1);
			if ((i + j) % 2 != 0)
			{
				b[j][i] = -b[j][i];
			}
		}
	}
	for (i = 0; i < num; i++)
	{
		delete[] temp[i];
	}
	delete[] temp;
	return;
}

int reverse_matrix(double **a, int num, double **b)
{
	int i = 0, j = 0;
	double **temp = new double *[num];
	for (i = 0; i < num; i++)
	{
		temp[i] = new double[num];
	}
	getastar(a, num, temp);
	double value = det(a, num);
	if (abs(value) < 1e-4)
	{
		cout << "can not reverse!please check!" << endl;
		return 0;
	}
	for (i = 0; i < num; i++)
	{
		for (j = 0; j < num; j++)
		{
			b[i][j] = temp[i][j] / value;
		}
	}
	return 1;
}

string int2string(int i)
{

	stringstream ss;
	ss << i;

	return ss.str();
}

inline double face_point_dis(double letice[2][3], double in_face_point[3], double point[3])
{
	int i = 0, j = 0;
	double fa[3] = { 0 };
	fa[0] = letice[0][1] * letice[1][2] - letice[0][2] * letice[1][1];
	fa[1] = letice[0][2] * letice[1][0] - letice[0][0] * letice[1][2];
	fa[2] = letice[0][0] * letice[1][1] - letice[0][1] * letice[1][0];
	double qp[3] = { 0 };
	for (i = 0; i < 3; i++)
	{
		qp[i] = in_face_point[i] - point[i];
	}
	double fenzi = qp[0] * fa[0] + qp[1] * fa[1] + qp[2] * fa[2];
	double dis = 0;
	dis = abs(fenzi / pow(pow(fa[0], 2) + pow(fa[1], 2) + pow(fa[2], 2), 0.5));
	return dis;
}


inline double vector_angle(double**a)
{
	//用来求两个向量的夹角，结果以角度制返回

	double a_mu = pow(   pow(a[0][0], 2) + pow(a[0][1], 2) + pow(a[0][2], 2) , 0.5);
	double b_mu= pow(pow(a[1][0], 2) + pow(a[1][1], 2) + pow(a[1][2], 2), 0.5);
	double diancheng = a[0][0] * a[1][0] + a[0][1] * a[1][1] + a[0][2] * a[1][2];
	double jungel_orig = acos(diancheng / (a_mu*b_mu))*180/pi;
	//cout << "the judge angle is:" << jungel_orig << endl;
	return jungel_orig;

	
}

void buble_plus(double* a, int *b, int num)
{
	int i = 0, j = 0;
	int flag = 0;
	for (i = 0; i < num; i++)
	{
		b[i] = i;
	}
	double temp;
	int temp_xuhao;
	for (i = num - 1; i >=0; i--)
	{
		flag = 0;
		for (j = 0; j < i; j++)
		{
			if (a[j] > a[j + 1])//如果是这样的话需要交换
			{
				flag = 1;
				temp = a[j];
				a[j] = a[j + 1];
				a[j + 1] = temp;

				temp_xuhao = b[j];
				b[j] = b[j + 1];
				b[j + 1] = temp_xuhao;
			}
		}
		if (flag == 0)
			break;
	}
	return;
}


//double det_a(int n, double **a)
//{
//	double **b;	/*��������b����ʼ��*/
//	int i = 0, j = 0; /*i��jΪ�����У�sumΪ����ʽ��ֵ*/
//	int x = 0, c = 0, p = 0;   /*��x�жϼ������c,pΪ�м����*/
//	double sum = 0;
//	if (n == 1)
//		return a[0][0];
//
//	for (i = 0; i < n; i++)
//	{
//		b[i] = new double[n];
//	}
//
//	for (i = 0; i < n; i++) /*�˴���ѭ��ʵ�ֽ�����ʽ��������b��*/
//	{
//		for (c = 0; c < n - 1; c++)
//		{
//			for (j = 0; j < n - 1; j++)
//			{
//				if (c < i)
//				{          /*����c�ж�ÿ�е��ƶ�����*/
//					p = 0; /*��p=0ʱ������ʽֻ�����ƣ�����ȥ��Ӧ�ĵ�һ�е���*/
//				}
//				else
//				{ /*��������ʽ���ƺ�������*/
//					p = 1;
//				}
//				b[c][j] = a[c + p][j + 1];
//			}
//		}
//		if (i % 2 == 0)
//		{ /*i+j����ʱj=0����ֻ����i��Ϊż�����ӷ�Ԥ��*/
//			x = 1;
//		}
//		else
//		{ /*i+jΪ��������������*/
//			x = (-1);
//		}
//		sum += a[i][0] * det_a(n - 1, b) * x; /*��������ʽ��ֵ*/
//
//	}
//	for (i = 0; i < n; i++)
//	{
//		delete[]b[i];
//	}
//	delete[]b;
//	return sum; /*��ֵ����*/
//}


void random_choose_two(int*org, int num, int**in)
{
	//从给定的org当中，任取两个，填充到in里面
	int i = 0, j = 0;
	//int total = num * (num - 1) / 2;	
	int count = 0;
	for (i = 0; i < num; i++)
	{		
		for (j = i+1; j < num; j++)
		{
			in[count][0] = org[i];
			in[count][1] = org[j];
			//cout << in[count][0] << "," << in[count][1] << endl;
			count++;
		}		
	}
	//cin.get();
	
}

int judge_falg_file(char* path)
{
	//成功返回1，失败返回0	
	ifstream fin;
	char real[200];
	strcpy(real, path);
	strcat(real, "/flag");
	fin.open(real,ios::in);
	string temp;
	fin >> temp;
	fin.close();
	if (temp == "Different")
		return 1;
	else if(temp.find("All")!=-1)
		return 0;
	else
	{
		cout << "wrong!flag file!" << endl;
		cin.get();
		return 0;
	}
}
