//
//  main.cpp
//  Compare_++
//
//  Created by 翁谋毅 on 2018/8/3.
//  Copyright © 2018年 翁谋毅. All rights reserved.
//
#define N_CEHNGSHU 1
#define DISTANCE 2.5
#include <fstream>
#include <stdio.h>
#include<string.h>
#include <string>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include <vector>
using namespace std;
const int cengshu = 3;
const int yanshen = cengshu * cengshu*cengshu;
int total = 0;
int count = 0;
const int meatal_xuhao[89] = { 3, 4, 11, 12, 13, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 37, 38, 39, 40, 41, 42, 43, 44, 45,46, 47, 48, 49, 50,  51,55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 85,87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111 };
const int main_group_element[16] = { 1, 5,6,7,8,9,14,15,16,17,33,34,35,52,53,84 };//at为金属元素
const int rare_gas[6] = { 2,10,18,36,54,86 };
const int other_elements = {};

const bool IF_IGNORE_TYPE =true; //是否完全忽略元素,也就是只判断灰色原子,这个如果为true，则下面的不起作用
const bool IF_ONLY_POSITIVE_AND_NEGATIVE =false; //是否只识别元素的正负，仅在判断同构时有用
const double ridus_plus_factor = 1.15;
const double val_radius_factor = 1.2;
const double metal_ridus_factor = 1.1;
const int metal_num = 89;
const int main_groupnum = 16;
const int rare_gasnum = 6;

const double H_MATAL_RULE = 2;
const double OF_RULE = 2;
double dist[120][120];//避免冲突规定了不同情况需要看的rule，这个默认是共价的rule
double dist_ri[120][120];//离子半径的rule
double dist_me[120][120];//金属键的rule


double dist_b[120][120];
int chengshu_num;
char name3[100];
int classify_metal_maingroup(int &atomic_number);


class element {
public:
	//element();
	int atomic_num;//atomic number
	char name[3];   // element name
	double vdw_radius_min;
	double vdw_radius_max;
	//van der waals radius,min and max
	double cov_radius;//covalence radius
	int num_metal_radius;  // number of metallic radius
	double * metal_radius; // the radius of metallic radius
	int num_common_val;  // number of common valence
	int * common_val;  // common valence
	int num_unusual_val; // number of unusual valence
	int * unusual_val;  // unusual valence
	double electron_negativity;   //electronic negativity
	double first_ionization_energy;
};

class cell
{
public:
	cell(char *jiegou_name, element* e,vector<vector<double>>&max_ion,int falg = 0);//0表示读取一般的结构，1表示读取带有半径信息的结构
	void judge_lonely(vector<vector<double>>max_ion);
	void change_to_pure_positiveandnegative();
	int num = 0;
	double **letice;
	double **p;
	double ***p_real;
	double ***real_position;
	double **ridus;
	int *type;
	int type_num;
	int *type_save;
	int *positive;


	int*if_positive;
	int* my_classify;
	int* if_lonely;//1孤独金属，2不孤独金属，3其他
};



class atom {
public:
	int xuhao;
	int atom_num;
	int l[3];           //cell place
	double p[3];        // portion place (without add l[3])
	double rp[3];      // real place
	int type;           // atom type number
	int near_num;        //
	atom *near;
	int level;
	void ini_atom(int *l_in, int n_in, cell a_in, int total);
	int find_near(cell in, int n,int flag);
	int find_near_iter(cell in, int n,int flag);
	int print_atom(cell in, int n,FILE*out);
	int print_atom(cell in, int n);
	double dis(double *p1, double *p2);
	int if_lonely;//1孤独金属，2不孤独金属，3其他
	int if_positive;//这个原子是不是正价
	int my_clsssif;//标记是金属还是非金属，1是金属，2是非金属，0是稀有气体

};
atom temp_aa[2000];
atom temp_bb[2000];

void read_atom(char *name_atom, double **d, double **p, int *type);
void read_atom_n(char *name_atom, int *n);
void make_matrix(int n, int **m, int *m_l, int *m_s, int *m_ss, int *m_p,int flag);
int compare_matrix(int n, int **m1, int **m2, int nn, int *m, int *m1_l, int *m2_l, int *m1_s, int *m2_s, int *m1_ss, int *m2_ss, int *m1_p, int *m2_p);
double dis(double *p1, double *p2);
void new_get_style(char *style);
void read_element(element *e, string& file_element_r, string& file_colvance, string& file_electronic_negativity, string& file_first_ionization_energy);
void read_dis();
void new_get_style(cell&cell_a, cell&cell_b, vector<vector<double>>max_ionic, element* e);
void Sort_temp(int total_num);
//进行相关元素的距离信息的调整
void showRidusInformation(element*e, vector<vector<double>> &max_ionic_riuds);
void showSingleInformation(element*e, vector<vector<double>> &max_ionic_riuds, int index);
void changeElementInformation(element*e, vector<vector<double>> &max_ionic_riuds);
const char atom_name[120][3] = { " ","H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lw", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg" };

int main(int argc, const char * argv[]) {
	string data_path = "C:\\Users\\王志\\Desktop\\硕士毕业\\审稿意见\\1T和2H\\";
	//string data_path ;
	string supporting_data_path= "C:\\Users\\王志\\Desktop\\graph_series\\graph_compare\\data\\";
	//string supporting_data_path = "./data/";
	//string supporting_data_path = "/share/home/wangz/projects/new_connection1231/bin/x64/Debug/";
	string file_name1 = supporting_data_path + "input_ionic";
	string file_name2 = supporting_data_path + "input_ionic_plus";
	string file_element_r = supporting_data_path + "ridus";
	string file_colvance = supporting_data_path + "colvance";
	string file_file_nagetivity = supporting_data_path + "negativity.txt";
	string file_first_ionization_energy = supporting_data_path + "first_ionazation_energy.txt";
	string max_ionic = supporting_data_path + "max_ionic";

	
	int l[3];
	int i, j, k, ii;
	int n1, n2;
	int **m1, **m2;
	int *m1_l, *m2_l;
	int *m1_s, *m2_s;
	int *m1_ss, *m2_ss;
	int *m1_p, *m2_p;
	int *squ;
	int *flag;
	FILE *out;
	int result = 1;
	double max1 = 0;
	double max2 = 0;
	double temp;
	char name1[100]="73711-0_1d.config", name2[100]="18208-0_1d.config", name4[5];
	for (i = 0; i < argc; i++)
	{
		if (i == 1)
		{
			for (j = 0; j < strlen(argv[i]); j++)
			{
				name1[j] = argv[i][j];
			}
			name1[strlen(argv[i])] = '\0';
		}
		if (i == 2)
		{
			for (j = 0; j < strlen(argv[i]); j++)
			{
				name2[j] = argv[i][j];
			}
			name2[strlen(argv[i])] = '\0';
		}
		if (i == 3)
		{
			for (j = 0; j < strlen(argv[i]); j++)
			{
				name3[j] = argv[i][j];
			}
			name3[strlen(argv[i])] = '\0';
		}
		if (i == 4)
		{
			for (j = 0; j < strlen(argv[i]); j++)
			{
				name4[j] = argv[i][j];
			}
			name4[strlen(argv[i])] = '\0';
		}
	}

	//先读取共价半径的信息
	
	element *e;
	e = new element[120];
	//开始依次读取相关信息
	//read_radius(ir, file_name1, file_name2);
	read_element(e, file_element_r, file_colvance, file_file_nagetivity, file_first_ionization_energy);
	ifstream fin;
	fin.open(max_ionic, ios::in);
	if (!fin.is_open()) {
		cout << "no file:" << max_ionic << endl;
		cin.get();
	}
	vector<vector<double>> max_ionic_riuds(120);
	for (i = 0; i < 120; i++)
	{
		max_ionic_riuds[i].resize(2);
		max_ionic_riuds[i][0] = -1;
		max_ionic_riuds[i][1] = -1;
		fin >> max_ionic_riuds[i][0];
		fin >> max_ionic_riuds[i][1];

	}
	fin.close();
	//showRidusInformation(e, max_ionic_riuds);
	changeElementInformation(e, max_ionic_riuds);
	cell cell_a(const_cast<char*>((data_path+name1).c_str()), e, max_ionic_riuds,2);
	cell cell_b(const_cast<char*>((data_path + name2).c_str()),e, max_ionic_riuds, 2);
	
	
	cell_a.judge_lonely(max_ionic_riuds);
	cell_b.judge_lonely(max_ionic_riuds);
	new_get_style(cell_a, cell_b, max_ionic_riuds, e);
	//元素只区分正负
	chengshu_num = 3;
	//打印一下距离阈值
	cout << dist[16][16] << "," << dist[52][52] << dist[33][52]<<endl;
	cout << dist_ri[16][40] << "," << dist_ri[16][50] << endl;
	cout << dist_me[40][40] << "," << dist_me[50][50] << endl;
	printf("%d\n", chengshu_num);

	atom *b;
	atom *a;
	a = (atom *)malloc(cell_a.num * sizeof(atom));
	b = (atom *)malloc(cell_b.num * sizeof(atom));

	for (i = 0; i < 3; i++)
	{
		l[i] = 0;
	}
	//在这里实现原子的初始化
	for (i = 0; i < cell_a.num; i++)//分别对两个晶胞中的原子进行距离、序数等信息的赋值
	{
		total++;
		a[i].ini_atom(l, i, cell_a, total);
		a[i].if_positive = cell_a.if_positive[i];
		a[i].my_clsssif = cell_a.my_classify[i];
		a[i].if_lonely = cell_a.if_lonely[i];

	}
	for (i = 0; i < cell_b.num; i++)
	{
		total++;
		b[i].ini_atom(l, i, cell_b, total);		
		b[i].if_positive = cell_b.if_positive[i];
		b[i].my_clsssif = cell_b.my_classify[i];
		b[i].if_lonely = cell_b.if_lonely[i];
	}
	////检查一下原子的赋值情况
	/*for (i = 0; i < cell_a.num; i++)
	{
		cout << atom_name[cell_a.type[i]]<<":"<<a[i].if_positive << "," << a[i].my_clsssif << "," << a[i].if_lonely << endl;
	}
	for (i = 0; i < cell_b.num; i++)
	{
		cout << atom_name[cell_b.type[i]] << ":"<< b[i].if_positive << "," << b[i].my_clsssif << "," << b[i].if_lonely << endl;
	}*/


	flag = (int *)malloc((cell_a.num * sizeof(int)));
	for (i = 0; i < cell_a.num; i++)
	{
		flag[i] = 0;
	}
	for (i = 0; i < cell_a.num; i++)
	{
		//printf("A:%d\n",i);
		a[i].level = 0;
		
		a[i].find_near_iter(cell_a, chengshu_num,0);
		out = fopen("tem1.xyz", "wb");
		fprintf(out, "           \nIteration 0\n");
		n1 = a[i].print_atom(cell_a, -1,out);
		Sort_temp(n1 + 1);//对得到的数据进行排序
		fseek(out, 0, SEEK_SET);
		fprintf(out, "   %d", n1 + 1);
		fclose(out);
		m1 = (int**)malloc((n1 + 2) * sizeof(int *));
		for (j = 0; j < n1 + 2; j++)
		{
			m1[j] = (int *)malloc((n1 + 2) * sizeof(int));
		}
		//这里加上初始化
		for (int i = 0; i < n1 + 2; i++)
		{
			for (int j = 0; j < n1 + 2; j++)
				m1[i][j] = 0;
		}
		m1_l = (int *)malloc((n1 + 2) * sizeof(int));
		m1_s = (int *)malloc((n1 + 2) * sizeof(int));
		m1_ss = (int *)malloc((n1 + 2) * sizeof(int));
		m1_p = (int *)malloc((n1 + 2) * sizeof(int));
		make_matrix(n1 + 1, m1, m1_l, m1_s, m1_ss, m1_p,1);
		out= fopen("m1.txt", "wb");
		for (j = 0; j < n1 + 2; j++)
		{
			for (k = 0; k < n1 + 2; k++)
			{
				fprintf(out, "%d\t", m1[j][k]);
			}
			fprintf(out, "\n");
		}
		fclose(out);//向m1文件中写入数据（矩阵）
		out = fopen("test1.xyz", "wb");
		fprintf(out, "   %d\n", n1 + 1);
		fprintf(out, "Iteration 0\n");
		for (j = 0; j < n1 + 1; j++)
		{
			fprintf(out, "%s\t%lf\t%lf\t%lf\t#%d\t%d\t%d\t%d\n", atom_name[temp_bb[j].type], temp_bb[j].rp[0], temp_bb[j].rp[1], temp_bb[j].rp[2], m1_l[j + 1], m1_s[j + 1], m1_ss[j + 1], m1_p[j + 1]);
			fprintf(out, "  %s  %lf  %lf  %lf\n",atom_name[temp_bb[j].type], temp_bb[j].rp[0], temp_bb[j].rp[1], temp_bb[j].rp[2]);
		}
		fclose(out);
		for (k = 0; k < cell_b.num; k++)
		{
			//printf("B:%d\n",k);
			/*if (i==4 && k==5)
			{
				printf("test\t atom1:%d, atom2:%d\n!!!!!\n\n", i+1,k+1);
			}
			if (i==4 && k==6)
			{
				printf("test\t atom1:%d, atom2:%d\n!!!!!\n\n", i+1,k+1);
			}*/
			/*if (cell_b.type[k]!=cell_a.type[i])
			{
				continue;
			}*/
			
			b[k].find_near_iter(cell_b, chengshu_num,0);
			out = fopen("tem2.xyz", "wb");
			fprintf(out, "           \niteration 0\n");
			n2 = b[k].print_atom(cell_b, -1,out);
			Sort_temp(n2 + 1);
			fseek(out, 0, SEEK_SET);
			fprintf(out, "   %d", n2 + 1);
			fclose(out);
			if (n2 == n1)
			{
				m2 = (int**)malloc((n1 + 2) * sizeof(int *));
				for (j = 0; j < n1 + 2; j++)
				{
					m2[j] = (int *)malloc((n1 + 2) * sizeof(int));
				}
				for (int i = 0; i < n1 + 2; i++)
				{
					for (int j = 0; j < n1 + 2; j++)
						m2[i][j] = 0;
				}
				m2_l = (int *)malloc((n1 + 2) * sizeof(int));
				m2_s = (int *)malloc((n1 + 2) * sizeof(int));
				m2_ss = (int *)malloc((n1 + 2) * sizeof(int));
				m2_p = (int *)malloc((n1 + 2) * sizeof(int));
				squ = (int *)malloc((n1 + 1) * sizeof(int));
				for (j = 0; j < n1 + 1; j++)
				{
					squ[j] = -2;
				}
				make_matrix(n2 + 1, m2, m2_l, m2_s, m2_ss, m2_p,1);
				out = fopen("m2.txt", "wb");
				for (ii = 0; ii < n2 + 2; ii++)
				{
					for (j = 0; j < n2 + 2; j++)
					{
						fprintf(out, "%d\t", m2[ii][j]);
					}
					fprintf(out, "\n");
				}
				fclose(out);
				out = fopen("test2.xyz", "wb");
				fprintf(out, "   %d\n", n1 + 1);
				fprintf(out, "Iteration 0\n");
				for (j = 0; j < n1 + 1; j++)
				{
					fprintf(out, "%s\t%lf\t%lf\t%lf\t#%d\t%d\t%d\t%d\n", atom_name[temp_bb[j].type], temp_bb[j].rp[0], temp_bb[j].rp[1], temp_bb[j].rp[2], m2_l[j + 1], m2_s[j + 1], m2_ss[j + 1], m2_p[j + 1]);
					fprintf(out, "  %s  %lf  %lf  %lf\n",atom_name[temp_bb[j].type], temp_bb[j].rp[0], temp_bb[j].rp[1], temp_bb[j].rp[2]);
				}
				fclose(out);
				::count = 0;//到这里同样也对第二个晶胞完成了构建矩阵的任务
				flag[i] = compare_matrix(n2 + 1, m1, m2, 0, squ, m1_l, m2_l, m1_s, m2_s, m1_ss, m2_ss, m1_p, m2_p);
				printf("flag=%d\n", flag[i]);
				for (j = 0; j < n1 + 2; j++)
				{
					free(m2[j]);
				}
				free(m2);
				free(m2_l);
				free(squ);
				free(m2_s);
				free(m2_ss);
				if (flag[i] == 1)
				{
					break;
				}
			}
		}
		if (flag[i] == 1)
		{
			//printf("atom1:%d\t==\tatom2:%d\n", i + 1, k + 1);
		}
		if (flag[i] == 0)
		{
			//printf("config1:atom%d\tdo not have same position in 2\n", i + 1);
			break;
		}
	}


	for (i = 0; i < cell_a.num; i++)
	{
		//printf("flag%d=%d\n", i, flag[i]);
		result = result * flag[i];
	}
	//printf("result=%d\n", result);


	for (i = 0; i < n1 + 2; i++)
	{
		free(m1[i]);
	}
	free(m1);
	free(m1_l);
	free(m1_s);
	free(m1_ss);

	if (result == 1)
	{
		cout << "1" << endl;
		system("pause\n");
		return 4;
	}
	else
	{
		 cout<< "0" << endl;
		system("pause\n");
		return 5;

	}
}
//修改特定元素的相关信息
void changeElementInformation(element*e, vector<vector<double>> &max_ionic_riuds) {
	e[14].cov_radius = 105;
	max_ionic_riuds[15][1] = 1.65;
	
}
//用来检查相关元素情况
void showRidusInformation(element*e,  vector<vector<double>> &max_ionic_riuds) {
	vector<int>index = { 8,17,25 };
	for( auto inn:index){
		showSingleInformation(e, max_ionic_riuds,inn);
	}
}
void showSingleInformation(element*e, vector<vector<double>> &max_ionic_riuds, int index) {
	cout << "element :" << atom_name[index]<<":";
	cout << "ionic ridus:" << max_ionic_riuds[index][0] << "," << max_ionic_riuds[index][1] << endl;
	cout << "metal ridus:" << e[index].metal_radius[0]/100.0 << endl;
	cout << "vdw :" << e[index].cov_radius/100.0 << endl;
	return;
}
void new_get_style(char *style)
{
	ifstream fin;
	fin.seekg(ios::beg);
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
				//cout << dist[i][j] << endl;
			}
		}
	}
	fin.close();
	return;
}

int compare_matrix(int n, int **m1, int **m2, int nn, int *m, int *m1_l, int *m2_l, int *m1_s, int *m2_s, int *m1_ss, int *m2_ss, int *m1_p, int *m2_p)
{//比较两个矩阵是否是同构的
	::count++;
	int flag = 0;
	int i, j;
	int sum1, sum2;
	int temp;
	//printf("n=%d\tnn=%d\n",n, nn);
	if (n == nn)
	{
		//printf("Yes\n");
		return 1;
	}
	if (nn == 0)
	{
		sum1 = 0;
		sum2 = 0;
		for (i = 0; i < n; i++)
		{
			sum1 = sum1 + m1[0][i + 1];
			sum2 = sum2 + m2[0][i + 1];
		}
		if (sum1 != sum2)
		{
			return 0;
		}
		sum1 = 0;
		sum2 = 0;
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
			{
				sum1 = sum1 + m1[i + 1][j + 1];
				sum2 = sum2 + m2[i + 1][j + 1];
			}
		}
		if (sum1 != sum2)
		{
			return 0;
		}
	}

	for (i = 0; i < n; i++)
	{
		/*if (count==3264 &&  i==34 && nn==22)
		{
			printf("test\n");
		}*/
		flag = 1;
		for (j = 0; j < nn; j++)
		{
			if (i == m[j])
			{
				flag = 0;
				break;
			}
			if (m1[i + 1][m[j] + 1] != m2[nn + 1][j + 1])
			{
				flag = 0;
				break;
			}
		}
		if (IF_IGNORE_TYPE == false)
		{
			if (m1_l[i + 1] == m2_l[nn + 1] && m1_s[i + 1] == m2_s[nn + 1] && m1_ss[i + 1] == m2_ss[nn + 1] && m1_p[i + 1] == m2_p[nn + 1] && m1[i + 1][0] == m2[nn + 1][0] && flag == 1)
			{
				flag = 1;
			}
			else
			{
				continue;
			}
		}		
		else
		{
			if (m1_l[i + 1] == m2_l[nn + 1]  && m1_p[i + 1] == m2_p[nn + 1] && m1[i + 1][0] == m2[nn + 1][0] && flag == 1)
			{
				flag = 1;
			}
			else
			{
				continue;
			}
		}
		/*
		sum1=0;
		for (j=0;j<n+1;j++)
		{
			sum1=sum1+m1[i+1][j];
		}
		sum2=0;
		for (j=0;j<n+1;j++)
		{
			sum2=sum2+m2[nn+1][j];
		}
		if (sum1==sum2)
		{
			flag=1;
		}
		else{
			continue;
		}*/


		if (flag == 1)
		{
			m[nn] = i;

			if (::count % 500000 == 0)
			{
				
				for (j = 0; j < nn + 1; j++)
				{
					//printf("%d\t", m[j]);
				}
				//printf("\n");
				printf("%d\t", ::count);
				return 0;
				
			}

			temp = compare_matrix(n, m1, m2, nn + 1, m, m1_l, m2_l, m1_s, m2_s, m1_ss, m2_ss, m1_p, m2_p);
			if (temp == 1)
			{
				return temp;
			}
		}
	}
	return 0;
}
void make_matrix(int n, int **m, int *m_l, int *m_s, int *m_ss, int *m_p,int flag)//生成矩阵，如果距离小于之前规定的距离，则认为两个原子之间有联系，用于后面的比较
{
	int i, j;
	double *r1, *r2;
	int sum;
	/*
	FILE *out;
	out=fopen("test.xyz", "wb");
	fprintf(out, "   %d\n", n);
	fprintf(out, "Iteration 0\n");
	for (i=0;i<n;i++)
	{
		fprintf(out, "  %s  %lf  %lf  %lf\n",a[temp_aa[i].type-1+temp_aa[i].level], temp_aa[i].rp[0], temp_aa[i].rp[1], temp_aa[i].rp[2]);
	}
	fclose(out);
	 */
	r1 = (double *)malloc(3 * sizeof(double));
	r2 = (double *)malloc(3 * sizeof(double));
	for (i = 0; i < n; i++)
	{
		m_l[i + 1] = temp_bb[i].level;
		//m[i+1][0]=temp_bb[i].type;
		//m[0][i+1]=temp_bb[i].type;
		m[i + 1][0] = 1;
		m[0][i + 1] = 1;
	}
	for (i = 0; i < n + 1; i++)
	{
		m[i][i] = -1;
	}
	for (i = 0; i < n; i++)
	{
		for (j = i; j < n; j++)
		{
			r1[0] = temp_bb[i].rp[0];
			r1[1] = temp_bb[i].rp[1];
			r1[2] = temp_bb[i].rp[2];
			r2[0] = temp_bb[j].rp[0];
			r2[1] = temp_bb[j].rp[1];
			r2[2] = temp_bb[j].rp[2];
			if (flag == 1)
			{
				if (temp_bb[i].if_lonely == 1)
				{
					if (temp_bb[j].my_clsssif == 1)
					{
						if (dis(r1, r2) < dist_me[temp_bb[i].type][temp_bb[j].type] && dis(r1, r2) > 0.1)
						{
							m[i + 1][j + 1] = 1;
							m[j + 1][i + 1] = 1;
						}
						else
						{
							m[i + 1][j + 1] = 0;
							m[j + 1][i + 1] = 0;
						}
					}
									

				}
				
				else if (temp_bb[i].if_lonely == 2)
				{
					if (temp_bb[i].if_positive * temp_bb[j].if_positive == -1)
					{

						if (dis(r1, r2) < dist_ri[temp_bb[i].type][temp_bb[j].type] && dis(r1, r2) > 0.1)
						{
							m[i + 1][j + 1] = 1;
							m[j + 1][i + 1] = 1;
						}
						else
						{
							m[i + 1][j + 1] = 0;
							m[j + 1][i + 1] = 0;
						}
					}
					
					
				}
				else if (temp_bb[i].if_lonely == 3)//会尝试连接金属或者非金属
				{
					//此时这是一个非金属
					int right_flag = 0;

					if (temp_bb[j].my_clsssif== 1 && temp_bb[i].if_positive * temp_bb[j].if_positive == -1)//连接金属
					{
						if (dis(r1, r2) < dist_ri[temp_bb[i].type][temp_bb[j].type] && dis(r1, r2) > 0.1)
						{
							m[i + 1][j + 1] = 1;
							m[j + 1][i + 1] = 1;
						}
						else
						{
							m[i + 1][j + 1] = 0;
							m[j + 1][i + 1] = 0;
						}
					}
					else if (temp_bb[j].my_clsssif == 2)//连接非金属
					{
						if (dis(r1, r2) < dist[temp_bb[i].type][temp_bb[j].type] && dis(r1, r2) > 0.1)
						{
							m[i + 1][j + 1] = 1;
							m[j + 1][i + 1] = 1;
						}
						else
						{
							m[i + 1][j + 1] = 0;
							m[j + 1][i + 1] = 0;
						}
					}
					else
					{
						//稀有气体不管
					}				

				}
				else
				{					
					cout << "connect the atom wrong!please check!" << endl;
					cin.get();
				}
			}			
			else if (flag == 0)
			{
				if (dis(r1, r2) > dist[temp_bb[i].type][temp_bb[j].type])
				{
					m[i + 1][j + 1] = 0;
					m[j + 1][i + 1] = 0;
				}
				else
				{
					m[i + 1][j + 1] = 1;
					m[j + 1][i + 1] = 1;
				}
			}
			else
			{
				cout << "wrong flag!please check!" << endl;
				cin.get();
			}
			
		}
	}
	for (i = 0; i < n + 1; i++)
	{
		sum = 0;
		for (j = 0; j < n + 1; j++)
		{
			sum = sum + m[i][j];
		}
		m_s[i] = sum;
	}
	for (i = 1; i < n + 1; i++)
	{
		sum = 0;
		for (j = 1; j < n + 1; j++)
		{
			if (m[i][j] == 1)
			{
				//为了只比较正负离子，在这里加上一条
				if (IF_ONLY_POSITIVE_AND_NEGATIVE == true)
				{
					cout << "yaosi" << endl;
					sum = sum + temp_bb[j - 1].level + temp_bb[j - 1].if_positive*10;
					cout << sum << endl;
				}
				else
					sum = sum + temp_bb[j - 1].level + temp_bb[j - 1].type;
			}
		}
		m_ss[i] = sum;
	}
	for (i = 1; i < n + 1; i++)
	{
		sum = 0;
		for (j = 1; j < n + 1; j++)
		{
			if (m[i][j] == 1)
			{
				sum = sum + m_s[j];
			}
		}
		m_p[i] = sum;
	}
	free(r1);
	free(r2);
	return;
}

double atom::dis(double *p1, double *p2)//同样是计算距离的函数
{
	return sqrt(pow(p1[0] - p2[0], 2) + pow(p1[1] - p2[1], 2) + pow(p1[2] - p2[2], 2));
}
double dis(double *p1, double *p2)//距离计算，返回一个距离值
{
	return sqrt(pow(p1[0] - p2[0], 2) + pow(p1[1] - p2[1], 2) + pow(p1[2] - p2[2], 2));
}

void atom::ini_atom(int *l_in, int n_in, cell a_in, int total)
{
	level = 0;
	atom_num = n_in;
	xuhao = total;
	l[0] = l_in[0];
	l[1] = l_in[1];
	l[2] = l_in[2];
	p[0] = a_in.p[n_in][0];
	p[1] = a_in.p[n_in][1];
	p[2] = a_in.p[n_in][2];
	rp[0] = p[0] * a_in.letice[0][0] + p[1] * a_in.letice[1][0] + p[2] * a_in.letice[2][0];
	rp[1] = p[0] * a_in.letice[0][1] + p[1] * a_in.letice[1][1] + p[2] * a_in.letice[2][1];
	rp[2] = p[0] * a_in.letice[0][2] + p[1] * a_in.letice[1][2] + p[2] * a_in.letice[2][2];
	type = a_in.type[n_in];
	near_num = 0;
	near = NULL;
	return;
}
int atom::find_near_iter(cell in, int n, int flag)
{
	find_near(in, n, flag);
	int i;
	if (n > 1)
	{
		for (i = 0; i < near_num; i++)
		{
			near[i].find_near_iter(in, n - 1, flag);
		}
		return 1;
	}
	if (n == 1)
	{
		return 0;
	}
	return 1;
}

int atom::find_near(cell in, int n, int flag)
{
	int i, j, k, ii;
	double temp_p[3];
	int num = 0;
	int pl[50][3];
	double pp[50][3];
	double prp[50][3];
	int ptype[50];
	int pxuhao[50];
	int patom_num[50];
	int plevel[50];

	int plonely[50];//标记新的原子的孤独标记
	int classify[50];//标记元素是主族还是金属
	int ppsotive[50];//标记是正还是负
	for (i = -1; i < 2; i++)
	{
		for (j = -1; j < 2; j++)
		{
			for (k = -1; k < 2; k++)
			{
				for (ii = 0; ii < in.num; ii++)
				{
					//printf("%d %d %d %d\n", i,j,k,ii);
					temp_p[0] = in.real_position[13][ii][0] + (i + l[0])*in.letice[0][0] + (j + l[1])*in.letice[1][0] + (k + l[2])*in.letice[2][0];
					temp_p[1] = in.real_position[13][ii][1] + (i + l[0])*in.letice[0][1] + (j + l[1])*in.letice[1][1] + (k + l[2])*in.letice[2][1];
					temp_p[2] = in.real_position[13][ii][2] + (i + l[0])*in.letice[0][2] + (j + l[1])*in.letice[1][2] + (k + l[2])*in.letice[2][2];
					/*if (i==0&&j==0&&k==0)
					{
						printf("%lf\n", in.rp[ii][2]);
						printf("%d\t%d\t%d\n", l[0], l[1],l[2]);
						printf("%lf\t%lf\t%lf\n", temp_p[0], temp_p[1], temp_p[2]);
						printf("%lf\t%lf\t%lf\n", rp[0], rp[1], rp[2]);
						printf("%d\t%d\t%d\t%d\t%lf\n", i,j,k,ii,dis(temp_p, rp));
					}*/
					//cout << type << endl;
					//这里我们要考虑金属键的问题
					if (flag == 0)
					{
						//cout << if_lonely << endl;
						if (if_lonely == 1)//针对孤独金属情况
						{
							if (in.my_classify[ii] == 1)
							{
								if (dis(temp_p, rp) < dist_me[type][in.type[ii]] && dis(temp_p, rp) > 0.1)
								{
									total++;
									plevel[num] = chengshu_num - n + 1;
									pxuhao[num] = total;
									patom_num[num] = ii;
									pl[num][0] = i + l[0];
									pl[num][1] = j + l[1];
									pl[num][2] = k + l[2];
									pp[num][0] = in.p[ii][0];
									pp[num][1] = in.p[ii][1];
									pp[num][2] = in.p[ii][2];
									prp[num][0] = temp_p[0];
									prp[num][1] = temp_p[1];
									prp[num][2] = temp_p[2];
									ptype[num] = in.type[ii];
									plonely[num] = in.if_lonely[ii];
									classify[num] = 1;
									ppsotive[num] = in.if_positive[ii];
									//printf("Pair%d:\n%lf\t%lf\t%lf\n%lf\t%lf\t%lf\n",num, rp[0],rp[1],rp[2], temp_p[0], temp_p[1], temp_p[2]);
									//printf("%d\t%d\t%d\n", i,j,k);
									//cout << "connect lonely metal!" << endl;
									num++;
								}
							}
						}
						else if (if_lonely == 2)
						{
							if ( in.if_positive[ii]*if_positive==-1)
							{
								
								if (dis(temp_p, rp) < dist_ri[type][in.type[ii]] && dis(temp_p, rp) > 0.1)
								{
									
									total++;
									plevel[num] = chengshu_num - n + 1;
									pxuhao[num] = total;
									patom_num[num] = ii;
									pl[num][0] = i + l[0];
									pl[num][1] = j + l[1];
									pl[num][2] = k + l[2];
									pp[num][0] = in.p[ii][0];
									pp[num][1] = in.p[ii][1];
									pp[num][2] = in.p[ii][2];
									prp[num][0] = temp_p[0];
									prp[num][1] = temp_p[1];
									prp[num][2] = temp_p[2];
									ptype[num] = in.type[ii];

									plonely[num] = in.if_lonely[ii];
									classify[num] = in.my_classify[ii];
									ppsotive[num] = in.if_positive[ii];
									//printf("Pair%d:\n%lf\t%lf\t%lf\n%lf\t%lf\t%lf\n",num, rp[0],rp[1],rp[2], temp_p[0], temp_p[1], temp_p[2]);
									//printf("%d\t%d\t%d\n", i,j,k);
									
									/*cout << "a" << endl;
									cout << plonely[num] << endl;
									cout << in.type[ii] << endl;
									cin.get();*/
									//cout << "connect unlonely metal!" << endl;
									num++;
								}
							}
						}
						else if (if_lonely == 3)//会尝试连接金属或者非金属
						{
							//此时这是一个非金属
							//cout << "connect unmetal!" << endl;
							int right_flag = 0;
							if (in.my_classify[ii] == 1 && in.if_positive[ii] * if_positive == -1)//连接金属
							{
								if (dis(temp_p, rp) < dist_ri[type][in.type[ii]] && dis(temp_p, rp) > 0.1)
								{
									//cout << dis(temp_p, rp) << ":metal:"<<dist_ri[type][in.type[ii]] << endl;
									right_flag = 1;
								}
								else {
									if(dis(temp_p, rp)<4)
										cout << dis(temp_p, rp) << ":metal:"<<dist_ri[type][in.type[ii]] << endl;
								}
							}
							else if (in.my_classify[ii] == 2)//非连接金属
							{
								if (dis(temp_p, rp) < dist[type][in.type[ii]] && dis(temp_p, rp) > 0.1)
								{
									//cout << dis(temp_p, rp) << ":nonmetal" << dist[type][in.type[ii]] << endl;
									right_flag = 1;

								}
								else {
									if (dis(temp_p, rp) < 4)
										cout << dis(temp_p, rp) << ":nonmetal" << dist[type][in.type[ii]] << endl;
								}
									
							}
							else
							{
								//稀有气体不管
							}
							if (right_flag == 1)
							{
								total++;
								plevel[num] = chengshu_num - n + 1;
								pxuhao[num] = total;
								patom_num[num] = ii;
								pl[num][0] = i + l[0];
								pl[num][1] = j + l[1];
								pl[num][2] = k + l[2];
								pp[num][0] = in.p[ii][0];
								pp[num][1] = in.p[ii][1];
								pp[num][2] = in.p[ii][2];
								prp[num][0] = temp_p[0];
								prp[num][1] = temp_p[1];
								prp[num][2] = temp_p[2];
								ptype[num] = in.type[ii];

								plonely[num] = in.if_lonely[ii];
								classify[num] = in.my_classify[ii];
								ppsotive[num] = in.if_positive[ii];
								//printf("Pair%d:\n%lf\t%lf\t%lf\n%lf\t%lf\t%lf\n",num, rp[0],rp[1],rp[2], temp_p[0], temp_p[1], temp_p[2]);
								//printf("%d\t%d\t%d\n", i,j,k);
								num++;
							}

						}
						else
						{
							cout << if_lonely <<","<<n<< endl;
							cout << this->if_positive << "," << this->atom_num << endl;
							cout << "connect the atom wrong!please check!" << endl;
							cin.get();
						}


						
					}
					else
					{
						if (dis(temp_p, rp) < dist_b[type][in.type[ii]] && dis(temp_p, rp) > 0.1)
						{
							total++;
							plevel[num] = chengshu_num - n + 1;
							pxuhao[num] = total;
							patom_num[num] = ii;
							pl[num][0] = i + l[0];
							pl[num][1] = j + l[1];
							pl[num][2] = k + l[2];
							pp[num][0] = in.p[ii][0];
							pp[num][1] = in.p[ii][1];
							pp[num][2] = in.p[ii][2];
							prp[num][0] = temp_p[0];
							prp[num][1] = temp_p[1];
							prp[num][2] = temp_p[2];
							ptype[num] = in.type[ii];
							//printf("Pair%d:\n%lf\t%lf\t%lf\n%lf\t%lf\t%lf\n",num, rp[0],rp[1],rp[2], temp_p[0], temp_p[1], temp_p[2]);
							//printf("%d\t%d\t%d\n", i,j,k);
							num++;
						}
					}

				}
			}
		}
	}
	near_num = num;
	near = (atom *)malloc(num * sizeof(atom));
	for (i = 0; i < num; i++)
	{
		//cout << i << ":" << num << endl;
		near[i].atom_num = patom_num[i];
		near[i].xuhao = pxuhao[i];
		near[i].type = ptype[i];
		near[i].near_num = 0;
		for (j = 0; j < 3; j++)
		{
			near[i].l[j] = pl[i][j];
			near[i].p[j] = pp[i][j];
			near[i].rp[j] = prp[i][j];
		}
		near[i].near = NULL;
		near[i].level = plevel[i];
		near[i].if_lonely = plonely[i];
		near[i].my_clsssif = classify[i];
		near[i].if_positive = ppsotive[i];
	}

	/*printf("Atom %d has %d neighbors\n", xuhao, num);
	 printf("%lf\t%lf\t%lf\n", rp[0], rp[1], rp[2]);
	 for (j=0;j<near_num;j++)
	 {
	 printf("%d\t%lf\t%lf\t%lf\n", j, near[j].rp[0], near[j].rp[1], near[j].rp[2]);
	 }*/
	return 0;
}
int atom::print_atom(cell in, int n)//生成一个文件tem1和tem2，用来输出表示找到的4层以内的原子种类及位置
{
	int i;
	char a[200][3] = { "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lw", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg" };
	for (i = 0; i < n + 1; i++)
	{
		if (atom_num == temp_aa[i].atom_num && l[0] == temp_aa[i].l[0] && l[1] == temp_aa[i].l[1] && l[2] == temp_aa[i].l[2])
		{
			if (near_num > 0 && level < temp_aa[i].level)
			{
				temp_aa[i].level = level;
				for (i = 0; i < near_num; i++)
				{
					n = near[i].print_atom(in, n);
				}
			}
			return n;
		}
	}
	n++;
	/*if (n > 2000)
	{
		printf("ERROR: too much atoms\n");
	}
	fprintf(out, "%s  %lf  %lf  %lf\n", a[type - 1], rp[0], rp[1], rp[2]);*/
	temp_aa[n].type = type;
	temp_aa[n].rp[0] = rp[0];
	temp_aa[n].rp[1] = rp[1];
	temp_aa[n].rp[2] = rp[2];
	temp_aa[n].l[0] = l[0];
	temp_aa[n].l[1] = l[1];
	temp_aa[n].l[2] = l[2];
	temp_aa[n].atom_num = atom_num;
	temp_aa[n].level = level;

	temp_aa[n].if_lonely = if_lonely;
	temp_aa[n].if_positive = if_positive;
	temp_aa[n].my_clsssif = my_clsssif;
	if (near_num > 0)
	{
		for (i = 0; i < near_num; i++)
		{
			n = near[i].print_atom(in, n);
		}
	}

	if (near_num == 0)
	{
		return n;
	}
	return n;
}

int atom::print_atom(cell in, int n,FILE* out)//生成一个文件tem1和tem2，用来输出表示找到的4层以内的原子种类及位置
{
	int i;
	char a[200][3] = { "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lw", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg" };
	for (i = 0; i < n + 1; i++)
	{
		if (atom_num == temp_aa[i].atom_num && l[0] == temp_aa[i].l[0] && l[1] == temp_aa[i].l[1] && l[2] == temp_aa[i].l[2])
		{
			if (near_num > 0 && level < temp_aa[i].level)
			{
				temp_aa[i].level = level;
				for (i = 0; i < near_num; i++)
				{
					n = near[i].print_atom(in, n,out);
				}
			}
			return n;
		}
	}
	n++;
	if (n > 2000)
	{
		printf("ERROR: too much atoms\n");
	}
	fprintf(out, "%s  %lf  %lf  %lf\n", a[type - 1], rp[0], rp[1], rp[2]);
	temp_aa[n].type = type;
	temp_aa[n].rp[0] = rp[0];
	temp_aa[n].rp[1] = rp[1];
	temp_aa[n].rp[2] = rp[2];
	temp_aa[n].l[0] = l[0];
	temp_aa[n].l[1] = l[1];
	temp_aa[n].l[2] = l[2];
	temp_aa[n].atom_num = atom_num;
	temp_aa[n].level = level;

	temp_aa[n].if_lonely = if_lonely;
	temp_aa[n].if_positive = if_positive;
	temp_aa[n].my_clsssif = my_clsssif;
	if (near_num > 0)
	{
		for (i = 0; i < near_num; i++)
		{
			n = near[i].print_atom(in, n,out);
		}
	}

	if (near_num == 0)
	{
		return n;
	}
	return n;
}


cell::cell(char *name, element*e, vector<vector<double>>&max_ion, int flag)
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
	if (flag == 1)
		positive = new int[num];
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
	if (flag == 1)
	{
		ridus = new double*[num];
		for (i = 0; i < num; i++)
		{
			ridus[i] = new double[2];
		}
	}

	while (fgets(temp, 300, in) != NULL)
	{
		if (strstr(temp, "VECTOR") != NULL || strstr(temp, "vector") != NULL || strstr(temp, "LATTICE") != NULL)
			break;
	}
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			fscanf(in, "%lf", &letice[i][j]);
		}
		fgets(temp, 300, in);
	}

	fgets(temp, 300, in);
	//cout << temp << endl;
	char line[10];
	for (i = 0; i < num; i++)
	{

		fscanf(in, "%d", &type[i]);
		fscanf(in, "%lf", &p[i][0]);
		fscanf(in, "%lf", &p[i][1]);
		fscanf(in, "%lf", &p[i][2]);
		if (flag == 1)
		{
			fscanf(in, "%s", line);
			if (line[0] == 'p')
			{
				positive[i] = 1;
				fscanf(in, "%lf", &ridus[i][0]);
			}
			else if (line[0] == 's')//读到s的时候说明是单质
			{
				fscanf(in, "%s", line);
				fscanf(in, "%lf", &ridus[i][0]);
			}
			else//如果是negative的话，就读两个ridus值
			{
				positive[i] = -1;
				fscanf(in, "%lf", &ridus[i][0]);//第一个是负价态对应的半径
				fscanf(in, "%lf", &ridus[i][1]);//对应正价态的半径
			}
		}


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

	if (flag == 2)
	{
		my_classify = new int[num];
		if_positive = new int[num];
		if_lonely = new int[num];
		//开始按照之前的结果进行标定
		for (i = 0; i < num; i++)
		{
			my_classify[i] = classify_metal_maingroup(type[i]);

		}
		for (i = 0; i < num; i++)
		{
			//针对H和O的特殊情况
			if (type[i] == 1)
			{
				if_positive[i] = 1;
				for (j = 0; j < yanshen; j++)
				{
					for (k = 0; k < num; k++)
					{
						if (my_classify[k] == 1 && e[type[k]].electron_negativity < 2.20)
						{
							double temp_dis = dis(real_position[13][i], real_position[j][k]);
							if (temp_dis < (max_ion[1][1] + max_ion[type[k]][0])*ridus_plus_factor)
							{
								if_positive[i] = -1;
								j = yanshen + 1;
								break;
							}
						}
					}
				}
			}
			/*else if (type[i] == 8)
			{
				if_positive[i] = -1;
				for (j = 0; j < num; j++)
				{
					if (type[j] == 9)
					{
						for (k = 0; k < yanshen; k++)
						{
							double temp_dis = dis(real_position[13][i],real_position[k][j]);
							if (temp_dis < OF_RULE)
							{
								if_positive[i] = 1;
								j = num + 1;
								break;
							}
						}
					}
				}
			}*/
			else
			{
				if (my_classify[i] == 0)
				{
					if_positive[i] = 0;
				}
				else if (my_classify[i] == 1)
				{
					if_positive[i] = 1;
				}
				else if (my_classify[i] == 2)
				{
					/*cout << "negative!" << endl;
					cin.get();*/
					if_positive[i] = -1;
				}
			}
			//cout << atom_name[type[i]] << ":" << my_classify[i] <<"positive:"<< if_positive[i] << endl;

		}
	}

	fclose(in);
}

//重载版本
void new_get_style(cell&cell_a, cell&cell_b, vector<vector<double>>max_ionic, element* e)
{
	//根据两个结构填充进去对应的判断距离
	//建立重载版本的构建连接关系
	int i = 0, j = 0;
	int m = 0, n = 0;
	double temp = 0.0;

	int num = cell_a.num;
	//首先尝试使用共价键连接原子
	for (i = 0; i < num; i++)
	{
		if (cell_a.my_classify[i%cell_a.num] == 2)
			for (j = i + 1; j < num; j++)
			{
				if (cell_a.my_classify[j%cell_a.num] == 2)
				{

					if (e[cell_a.type[i%cell_a.num]].cov_radius == -1 || e[cell_a.type[j%cell_a.num]].cov_radius == -1)
					{
						//cout << "unknown covridus!:atomic :" << cell_a.type[i%cell_a.num] << "," << cell_a.type[j%cell_a.num] << endl;
						//cin.get();
						continue;
					}
					double rule = (e[cell_a.type[i%cell_a.num]].cov_radius + e[cell_a.type[j%cell_a.num]].cov_radius) / 100.0*val_radius_factor;
					dist[cell_a.type[i]][cell_a.type[j]] = rule;
					dist[cell_a.type[j]][cell_a.type[i]] = rule;
				}
			}
	}

	//然后使用离子半径连接正负基团
	for (i = 0; i < num; i++)
	{
		if (cell_a.if_positive[i%cell_a.num] == 1)
		{
			for (j = 0; j < num; j++)
			{
				if (cell_a.if_positive[j%cell_a.num] == -1)
				{

					if (max_ionic[cell_a.type[i%cell_a.num]][0] < -100 || max_ionic[cell_a.type[j%cell_a.num]][1] < -100)
					{
						cout << "unknown ionic!:atomic :" << atom_name[cell_a.type[i%cell_a.num]] << "," << atom_name[cell_a.type[j%cell_a.num]] << endl;
						cin.get();
						continue;
					}
					double rule = (max_ionic[cell_a.type[i%cell_a.num]][0] + max_ionic[cell_a.type[j%cell_a.num]][1])*ridus_plus_factor;
					dist_ri[cell_a.type[i%cell_a.num]][cell_a.type[j%cell_a.num]] = rule;
					dist_ri[cell_a.type[j%cell_a.num]][cell_a.type[i%cell_a.num]] = rule;
				}
			}
		}
	}


	//然后检查金属的情况，如果金属周围是孤立的，则尝试用金属半径连接金属与金属
	//这里我们只是写上连接信息，判断原子由下步完成
	for (i = 0; i < num; i++)
	{
		if (cell_a.my_classify[i%cell_a.num] == 1)
		{
			for (j = 0; j < yanshen*num; j++)
			{
				if (cell_a.my_classify[j%cell_a.num] == 1)
				{
					double rule = (e[cell_a.type[i]].metal_radius[0] + e[cell_a.type[j%cell_a.num]].metal_radius[0]) / 100.0*metal_ridus_factor;
					dist_me[cell_a.type[i]][cell_a.type[j%cell_a.num]] = rule;
					dist_me[cell_a.type[j%cell_a.num]][cell_a.type[i]] = rule;

				}

			}
		}
	}



	//然后是针对结构B的情况
		//首先尝试使用共价键连接原子
	num = cell_b.num;
	for (i = 0; i < num; i++)
	{
		if (cell_b.my_classify[i%cell_b.num] == 2)
			for (j = i + 1; j < num; j++)
			{
				if (cell_b.my_classify[j%cell_b.num] == 2)
				{

					if (e[cell_b.type[i%cell_b.num]].cov_radius == -1 || e[cell_b.type[j%cell_b.num]].cov_radius == -1)
					{
						//cout << "unknown covridus!:atomic :" << cell_b.type[i%cell_b.num] << "," << cell_b.type[j%cell_b.num] << endl;
						//cin.get();
						continue;
					}
					double rule = (e[cell_b.type[i%cell_b.num]].cov_radius + e[cell_b.type[j%cell_b.num]].cov_radius) / 100.0*val_radius_factor;
					dist[cell_b.type[i]][cell_b.type[j]] = rule;
					dist[cell_b.type[j]][cell_b.type[i]] = rule;
				}
			}
	}

	//然后使用离子半径连接正负基团
	for (i = 0; i < num; i++)
	{
		if (cell_b.if_positive[i%cell_b.num] == 1)
		{
			for (j = 0; j < num; j++)
			{
				if (cell_b.if_positive[j%cell_b.num] == -1)
				{

					if (max_ionic[cell_b.type[i%cell_b.num]][0] < -100 || max_ionic[cell_b.type[j%cell_b.num]][1] < -100)
					{
						//cout << "unknown ionic!:atomic :" << atom_name[cell_b.type[i%cell_b.num]] << "," << atom_name[cell_b.type[j%cell_b.num]] << endl;
						//cin.get();
						continue;
					}
					double rule = (max_ionic[cell_b.type[i%cell_b.num]][0] + max_ionic[cell_b.type[j%cell_b.num]][1])*ridus_plus_factor;
					dist_ri[cell_b.type[i%cell_b.num]][cell_b.type[j%cell_b.num]] = rule;
					dist_ri[cell_b.type[j%cell_b.num]][cell_b.type[i%cell_b.num]] = rule;
				}
			}
		}
	}


	//然后检查金属的情况，如果金属周围是孤立的，则尝试用金属半径连接金属与金属
	//这里我们只是写上连接信息，判断原子由下步完成
	for (i = 0; i < num; i++)
	{
		if (cell_b.my_classify[i%cell_b.num] == 1)
		{
			for (j = 0; j < yanshen*num; j++)
			{
				if (cell_b.my_classify[j%cell_b.num] == 1)
				{
					double rule = (e[cell_b.type[i]].metal_radius[0] + e[cell_b.type[j%cell_b.num]].metal_radius[0]) / 100.0*metal_ridus_factor;
					dist_me[cell_b.type[i]][cell_b.type[j%cell_b.num]] = rule;
					dist_me[cell_b.type[j%cell_b.num]][cell_b.type[i]] = rule;

				}

			}
		}
	}


	return;
}
void cell::judge_lonely(vector<vector<double>>max_ions)
{
	int i = 0, j = 0;
	double temp;
	for (i = 0; i < this->num; i++)
	{
		if (my_classify[i] == 2 || my_classify[i] == 0)
		{
			if_lonely[i] = 3;
		}
		else
		{
			int iso_flag = 1;
			if (this->if_positive[i] == -1)
			{
				if_lonely[i] = 2;
			}
			else if (this->if_positive[i] == 1)
			{
				for (j = 0; j < yanshen*num; j++)
				{

					if (this->my_classify[j%this->num] == 2 && this->if_positive[j%this->num] == -1)
					{
						temp = dis(this->real_position[13][i], this->real_position[j / this->num][j % this->num]);
						if (abs(temp) < 1e-3)
							continue;
						double rule = (max_ions[type[i]][0] + max_ions[type[j%this->num]][1])*ridus_plus_factor;
						if (temp < rule)
						{
							iso_flag = 0;
							this->if_lonely[i] = 2;
							break;
						}
					}

				}
				if (iso_flag == 1)
				{
					this->if_lonely[i] = 1;
				}

			}
			else
				continue;

			
		}
		
	}
}

void Sort_temp(int total_num)//对tme1进行一定的排序
{
	int i, j, k;
	int ccn = 0;
	for (i = 0; i <= chengshu_num; i++)
	{
		for (j = 0; j < total_num; j++)
		{
			if (temp_aa[j].level == i)
			{
				temp_bb[ccn].type = temp_aa[j].type;
				temp_bb[ccn].atom_num = temp_aa[j].atom_num;
				temp_bb[ccn].xuhao = temp_aa[j].xuhao;
				temp_bb[ccn].near_num = temp_aa[j].near_num;
				temp_bb[ccn].level = temp_aa[j].level;
				//继续补充
				temp_bb[ccn].if_lonely = temp_aa[j].if_lonely;
				temp_bb[ccn].if_positive = temp_aa[j].if_positive;
				temp_bb[ccn].my_clsssif = temp_aa[j].my_clsssif;
				for (k = 0; k < 3; k++)
				{
					temp_bb[ccn].l[k] = temp_aa[j].l[k];
					temp_bb[ccn].p[k] = temp_aa[j].p[k];
					temp_bb[ccn].rp[k] = temp_aa[j].rp[k];
				}
				ccn++;
			}
		}
	}
}


void read_dis()//读取style。ini文件
{
	char a[120][3] = { " ","H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lw", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg" };
	int i, j;
	int num1=0, num2=0;
	FILE *in;
	char str[300];
	char name1[5];
	char name2[5];
	char temp[100];
	for (i = 0; i < 120; i++)
	{
		for (j = 0; j < 120; j++)
		{
			dist[i][j] = 0;
		}
	}
	in = fopen("style.ini", "r");
	if (in == NULL)
	{
		printf("ERROR Read style.ini file\n");
		return;
	}
	while (fgets(str, 300, in) != NULL)
	{
		if (strstr(str, "SBOND") != NULL)
		{
			break;
		}
	}
	for (i = 0; i < 914; i++)
	{
		fscanf(in, "%s", temp);
		fscanf(in, "%s", name1);
		fscanf(in, "%s", name2);
		fscanf(in, "%s", temp);
		for (j = 0; j < 120; j++)//循环读取了原子序号
		{
			if (a[j][0] == name1[0] && a[j][1] == name1[1])
			{
				num1 = j;
				break;
			}
		}
		for (j = 0; j < 120; j++)
		{
			if (a[j][0] == name2[0] && a[j][1] == name2[1])
			{
				num2 = j;
				break;
			}
		}
		if (j == 119)
		{
			printf("ERROR Read %s\n", name3);
			return;
		}
		fscanf(in, "%lf", &dist[num1][num2]);//读取两个原子之间的最小成键长度
		dist[num2][num1] = dist[num1][num2];
		fgets(temp, 100, in);
	}
	fclose(in);

}
void read_element(element *e, string &file_element_r)
{
	int i, j, k;
	char atom_name[120][3] = { "D", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg" };

	//先初始化信息
	for (i = 0; i < 120; i++)
	{
		//strcpy(e[i].name, atom_name[i]);
		e[i].atomic_num = i;
		//e[i].num_common_val = 0;
		e[i].num_metal_radius = 0;
		//e[i].num_unusual_val = 0;
	}
	i = 0, j = 0, k = 0;
	int line_in;
	ifstream fin;
	int atom_num_in;
	double vdw_r_min_in, vdw_r_max_in;
	double cov_r;
	int num_metal_r;
	double *temp;
	char str[5];
	char str_temp[500];
	int num_common_in;
	int num_unsual_in;
	//读取共价键范德华键和金属键

	fin.open(file_element_r.c_str(), ios::in);
	if (!fin.is_open())
	{
		cout << "i can not find the file !" << file_element_r.c_str() << endl;
		cin.get();
	}
	while (fin.peek() != EOF && fin.good())
	{
		fin >> atom_num_in;
		fin >> vdw_r_min_in;
		fin >> vdw_r_max_in;
		fin >> cov_r;
		if (cov_r == -1)
		{
			//cout << atom_name[atom_num_in]<<","<< atom_num_in << endl;
		}
		fin >> num_metal_r;

		e[atom_num_in].metal_radius = new double[num_metal_r];
		double *temp = new double[num_metal_r];
		for (i = 0; i < num_metal_r; i++)
		{
			fin >> temp[i];
		}
		e[atom_num_in].atomic_num = atom_num_in;
		e[atom_num_in].vdw_radius_max = vdw_r_max_in;
		e[atom_num_in].vdw_radius_min = vdw_r_min_in;
		e[atom_num_in].cov_radius = cov_r;

		e[atom_num_in].num_metal_radius = num_metal_r;
		for (i = 0; i < num_metal_r; i++)
		{
			e[atom_num_in].metal_radius[i] = temp[i];
		}
		free(temp);
	}
	fin.close();
	return;
}
int classify_metal_maingroup(int &atomic_number)//0 rare gas,1 metal ,2 main group
{

	int j = 0, k = 0, m = 0;
	while (j < main_groupnum || k < metal_num || m < rare_gasnum)
	{
		if (rare_gas[m] == atomic_number)
		{
			/*cout << "0 "<<atomic_number << endl;
			cin.get();*/
			return 0;
		}
		if (meatal_xuhao[k] == atomic_number)
		{
			
			return 1;
		}
		if (main_group_element[j] == atomic_number)
		{
			
			return 2;
		}


		if (meatal_xuhao[k] < atomic_number)
		{
			k++;
			continue;
		}
		if (rare_gas[m] < atomic_number)
		{
			m++;
			continue;
		}
		if (main_group_element[j] < atomic_number)
		{
			j++;
			continue;
		}
	}
	cout << "what's wrong!" << atomic_number << endl;
	cin.get();
	return -1;

}


void read_element(element *e, string& file_element_r, string& file_colvance, string& file_electronic_negativity, string& file_first_ionization_energy)
{
	int i, j, k;
	char atom_name[120][3] = { "D","H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg" };
	//cout << atom_name[87] << endl;
	//先初始化信息
	for (i = 0; i < 120; i++)
	{
		strcpy(e[i].name, atom_name[i]);
		e[i].atomic_num = i;
		e[i].num_common_val = 0;
		e[i].num_metal_radius = 0;
		e[i].num_unusual_val = 0;
	}
	i = 0, j = 0, k = 0;
	int line_in;
	ifstream fin;
	int atom_num_in;
	double vdw_r_min_in, vdw_r_max_in;
	double cov_r;
	int num_metal_r;
	double *temp;
	char str[5];
	char str_temp[500];
	int num_common_in;
	int num_unsual_in;
	//读取共价键范德华键和金属键

	fin.open(file_element_r.c_str(), ios::in);
	if (!fin.is_open())
	{
		cout << "i can not find the file !" << file_element_r.c_str() << endl;
		cin.get();
	}
	while (fin.peek() != EOF && fin.good())
	{
		fin >> atom_num_in;
		fin >> vdw_r_min_in;
		fin >> vdw_r_max_in;
		fin >> cov_r;
		fin >> num_metal_r;

		e[atom_num_in].metal_radius = new double[num_metal_r];
		double *temp = new double[num_metal_r];
		for (i = 0; i < num_metal_r; i++)
		{
			fin >> temp[i];
		}
		e[atom_num_in].atomic_num = atom_num_in;
		e[atom_num_in].vdw_radius_max = vdw_r_max_in;
		e[atom_num_in].vdw_radius_min = vdw_r_min_in;
		e[atom_num_in].cov_radius = cov_r;
		e[atom_num_in].num_metal_radius = num_metal_r;
		for (i = 0; i < num_metal_r; i++)
		{
			e[atom_num_in].metal_radius[i] = temp[i];
		}
		free(temp);
	}
	fin.close();


	////然后开始读取价态相关信息
	//fin.clear();
	//fin.open(file_colvance.c_str(), ios::in);
	//if (!fin.is_open())
	//{
	//	cout << "i can not find the file !" << file_colvance.c_str() << endl;
	//	cin.get();
	//}
	//while (fin.peek() != EOF && fin.good())
	//{
	//	fin >> atom_num_in;
	//	//printf("%d\n", atom_num_in);
	//	//fin >> num_common_in;
	//	fin >> e[atom_num_in].num_common_val;
	//	e[atom_num_in].common_val = new int[e[atom_num_in].num_common_val];
	//	for (i = 0; i < e[atom_num_in].num_common_val; i++)
	//	{
	//		fin >> e[atom_num_in].common_val[i];
	//	}
	//	fin >> e[atom_num_in].num_unusual_val;
	//	e[atom_num_in].unusual_val = new int[e[atom_num_in].num_unusual_val];
	//	for (i = 0; i < e[atom_num_in].num_unusual_val; i++)
	//	{
	//		fin >> e[atom_num_in].unusual_val[i];
	//	}
	//	//delete[]e[atom_num_in].common_val;
	//}
	//fin.close();

	//开始读取电负性相关信息
	string temp_test;
	fin.clear();
	fin.open(file_electronic_negativity.c_str(), ios::in);
	if (!fin.is_open())
	{
		cout << "i can not find the file !" << file_electronic_negativity.c_str() << endl;
		cin.get();
	}
	while (fin.peek() != EOF && fin.good())
	{
		fin >> atom_num_in;
		//cout << atom_num_in << endl;
		fin >> str;
		//fin >> temp_test;
		//cout <<temp_test << endl;
		/*if (str[0] == '\0')
		{
			break;
		}*/
		fin >> e[atom_num_in].electron_negativity;
		//cout << e[atom_num_in].electron_negativity << endl;
		if (str[0] != atom_name[atom_num_in][0] || str[1] != atom_name[atom_num_in][1])
		{
			printf("ERROR atom name %s %s %d\n", str, atom_name[atom_num_in], atom_num_in);
			cout << file_electronic_negativity.c_str() << endl;
			cin.get();
		}
	}
	fin.close();

	////最后开始读取第一电离能的相关信息
	//fin.clear();
	//fin.open(file_first_ionization_energy.c_str(), ios::in);
	//if (!fin.is_open())
	//{
	//	cout << "i can not find the file !" << file_first_ionization_energy << endl;
	//	cin.get();
	//}
	//while (fin.peek() != EOF && fin.good())
	//{
	//	fin >> atom_num_in;
	//	//cout << atom_num_in << endl;
	//	fin >> str;
	//	/*if (str[0] = '\0')
	//	{
	//		break;
	//	}*/
	//	fin >> e[atom_num_in].first_ionization_energy;
	//	//cout << e[atom_num_in].first_ionization_energy << endl;
	//	if (str[0] != atom[atom_num_in][0] || str[1] != atom[atom_num_in][1])
	//	{
	//		printf("ERROR atom name %s %s %d\n", str, atom[atom_num_in], atom_num_in);
	//		cout << file_first_ionization_energy << endl;
	//		cin.get();
	//	}
	//}
	//fin.close();
	return;
}


void cell::change_to_pure_positiveandnegative()
{
	//将相关元素信息完全转化为正负的信息
	for (int i = 0; i < this->num; i++)
	{
		this->type[i] = this->if_positive[i];
	}
	return;
}