//
//  main.cpp
//  Compare_++
//
//  Created by 翁谋毅 on 2018/8/3.
//  Copyright © 2018年 翁谋毅. All rights reserved.
//
#define N_CEHNGSHU 3
#define DISTANCE 2.5
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include <vector>
int total = 0;
int count = 0;
double dist[120][120];
int chengshu_num;
char name3[100];
class cell {
public:
	cell(char *name_atom);
	double **d;//晶格常数
	double **p;//相对坐标
	double **rp;//真实坐标
	int *type;
	int num;

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
	int find_near(cell in, int n);
	int find_near_iter(cell in, int n);
	int print_atom(cell in, int n, FILE *out);
	double dis(double *p1, double *p2);
};
atom temp_aa[2000];
atom temp_bb[2000];

void read_atom(char *name_atom, double **d, double **p, int *type);
void read_atom_n(char *name_atom, int *n);
void make_matrix(int n, int **m, int *m_l, int *m_s, int *m_ss, int *m_p);
int compare_matrix(int n, int **m1, int **m2, int nn, int *m, int *m1_l, int *m2_l, int *m1_s, int *m2_s, int *m1_ss, int *m2_ss, int *m1_p, int *m2_p);
double dis(double *p1, double *p2);
void read_dis();


void Sort_temp(int total_num);


int main(int argc, const char * argv[]) {
	char atom_name[120][3] = { " ","H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lw", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg" };
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
	char name1[100]="atom (1094).config", name2[100]="atom (700).config", name4[5];
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

	strcpy(name3, "style.ini");
	cell cell_a(name1);
	cell cell_b(name2);

	//cell cell_a(name1);
	//cell cell_b(name2);

	read_dis();
	for (i = 0; i < 3; i++)//开始计算晶胞中两个原子之间的距离
	{
		temp = 0;
		for (j = 0; j < 3; j++)
		{
			temp = temp + cell_a.d[i][j] * cell_a.d[i][j];
		}
		if (max1 < sqrt(temp))
		{
			max1 = sqrt(temp);
		}
	}
	for (i = 0; i < 3; i++)
	{
		temp = 0;
		for (j = 0; j < 3; j++)
		{
			temp = temp + cell_b.d[i][j] * cell_b.d[i][j];
		}
		if (max2 < sqrt(temp))
		{
			max2 = sqrt(temp);
		}
	}
	if (max1 >= max2)
	{
		chengshu_num = (int)(max1 / 4 + 1);
	}
	if (max1 < max2)
	{
		chengshu_num = (int)(max2 / 4 + 1);
	}
	chengshu_num = name4[0] - '0';
	chengshu_num = 3;
	printf("%d\n", chengshu_num);

	atom *b;
	atom *a;
	a = (atom *)malloc(cell_a.num * sizeof(atom));
	b = (atom *)malloc(cell_b.num * sizeof(atom));

	for (i = 0; i < 3; i++)
	{
		l[i] = 0;
	}
	for (i = 0; i < cell_a.num; i++)//分别对两个晶胞中的原子进行距离、序数等信息的赋值
	{
		total++;
		a[i].ini_atom(l, i, cell_a, total);
	}
	for (i = 0; i < cell_b.num; i++)
	{
		total++;
		b[i].ini_atom(l, i, cell_b, total);
	}
	/*
	for (i=0;i<cell_b.num;i++)
	{
		printf("%d\t%lf\t%lf\t%lf\n", b[i].type, b[i].rp[0], b[i].rp[1], b[i].rp[2]);
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
		a[i].find_near_iter(cell_a, chengshu_num);
		out = fopen("tem1.xyz", "wb");
		fprintf(out, "           \nIteration 0\n");
		n1 = a[i].print_atom(cell_a, -1, out);
		Sort_temp(n1 + 1);//对得到的数据进行排序
		fseek(out, 0, SEEK_SET);
		fprintf(out, "   %d", n1 + 1);
		fclose(out);
		m1 = (int**)malloc((n1 + 2) * sizeof(int *));
		for (j = 0; j < n1 + 2; j++)
		{
			m1[j] = (int *)malloc((n1 + 2) * sizeof(int));
		}
		m1_l = (int *)malloc((n1 + 2) * sizeof(int));
		m1_s = (int *)malloc((n1 + 2) * sizeof(int));
		m1_ss = (int *)malloc((n1 + 2) * sizeof(int));
		m1_p = (int *)malloc((n1 + 2) * sizeof(int));
		make_matrix(n1 + 1, m1, m1_l, m1_s, m1_ss, m1_p);
		out = fopen("m1.txt", "wb");
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
			//fprintf(out, "  %s  %lf  %lf  %lf\n",atom_name[temp_bb[j].type], temp_bb[j].rp[0], temp_bb[j].rp[1], temp_bb[j].rp[2]);
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
			b[k].find_near_iter(cell_b, chengshu_num);
			out = fopen("tem2.xyz", "wb");
			fprintf(out, "           \nIteration 0\n");
			n2 = b[k].print_atom(cell_b, -1, out);
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
				m2_l = (int *)malloc((n1 + 2) * sizeof(int));
				m2_s = (int *)malloc((n1 + 2) * sizeof(int));
				m2_ss = (int *)malloc((n1 + 2) * sizeof(int));
				m2_p = (int *)malloc((n1 + 2) * sizeof(int));
				squ = (int *)malloc((n1 + 1) * sizeof(int));
				for (j = 0; j < n1 + 1; j++)
				{
					squ[j] = -2;
				}
				make_matrix(n2 + 1, m2, m2_l, m2_s, m2_ss, m2_p);
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
					//fprintf(out, "  %s  %lf  %lf  %lf\n",atom_name[temp_bb[j].type], temp_bb[j].rp[0], temp_bb[j].rp[1], temp_bb[j].rp[2]);
				}
				fclose(out);
				count = 0;//到这里同样也对第二个晶胞完成了构建矩阵的任务
				flag[i] = compare_matrix(n2 + 1, m1, m2, 0, squ, m1_l, m2_l, m1_s, m2_s, m1_ss, m2_ss, m1_p, m2_p);
				//printf("flag=%d\n", flag[i]);
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
			printf("atom1:%d\t==\tatom2:%d\n", i + 1, k + 1);
		}
		if (flag[i] == 0)
		{
			printf("config1:atom%d\tdo not have same position in 2\n", i + 1);
			break;
		}
	}


	for (i = 0; i < cell_a.num; i++)
	{
		//printf("flag%d=%d\n", i, flag[i]);
		result = result * flag[i];
	}
	printf("result=%d\n", result);
	system("pause");
	//scanf("%d", a);
	for (i = 0; i < n1 + 2; i++)
	{
		free(m1[i]);
	}
	free(m1);
	free(m1_l);
	free(m1_s);
	free(m1_ss);
	if (result == 1)
		return 25;
	else
		return 0;

}

void read_dis()//读取style。ini文件
{
	char a[120][3] = { " ","H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lw", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg" };
	int i, j;
	int num1, num2;
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
	in = fopen(name3, "r");
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

int compare_matrix(int n, int **m1, int **m2, int nn, int *m, int *m1_l, int *m2_l, int *m1_s, int *m2_s, int *m1_ss, int *m2_ss, int *m1_p, int *m2_p)
{//比较两个矩阵是否是同构的
	count++;
	printf("aaa\n");
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
		printf("1\n");
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
	printf("2\n");
	//给程序一个向下运行的标志
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
		if (m1_l[i + 1] == m2_l[nn + 1] && m1_p[i + 1] == m2_p[nn + 1] && m1[i + 1][0] == m2[nn + 1][0] && flag == 1)
		{
			flag = 1;
		}
		else
		{
			continue;
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

			if (count % 500000 == 0)//这里这个count的作用是防止产生无限循环，每过一段时间出现一次
			{
				printf("%d\t", count);
				for (j = 0; j < nn + 1; j++)
				{
					printf("%d\t", m[j]);
				}
				printf("\n");
			}

			temp = compare_matrix(n, m1, m2, nn + 1, m, m1_l, m2_l, m1_s, m2_s, m1_ss, m2_ss, m1_p, m2_p);
			if (temp == 1)
			{
				return temp;//这里的temp作用是，作为递归函数总要返回一个数值
			}
		}
	}
	return 0;
}
void make_matrix(int n, int **m, int *m_l, int *m_s, int *m_ss, int *m_p)//生成矩阵，如果距离小于之前规定的距离，则认为两个原子之间有联系，用于后面的比较
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

void atom::ini_atom(int *l_in, int n_in, cell a_in, int total)//这里将晶胞中的位置、相对位置、原子序数赋给每一个原子
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
	rp[0] = p[0] * a_in.d[0][0] + p[1] * a_in.d[1][0] + p[2] * a_in.d[2][0];
	rp[1] = p[0] * a_in.d[0][1] + p[1] * a_in.d[1][1] + p[2] * a_in.d[2][1];
	rp[2] = p[0] * a_in.d[0][2] + p[1] * a_in.d[1][2] + p[2] * a_in.d[2][2];
	type = a_in.type[n_in];
	near_num = 0;
	near = NULL;
	return;
}
int atom::find_near_iter(cell in, int n)//根据设置完成找几层的原子
{
	find_near(in, n);
	int i;
	if (n > 1)
	{
		for (i = 0; i < near_num; i++)
		{
			near[i].find_near_iter(in, n - 1);
		}
		return 1;
	}
	if (n == 1)
	{
		return 0;
	}
	return 1;
}

int atom::find_near(cell in, int n)//在一层中完成了找相邻原子的任务
{
	int i, j, k, ii;
	double temp_p[3];
	int num = 0;
	int pl[20][3];
	double pp[20][3];
	double prp[20][3];
	int ptype[20];
	int pxuhao[20];
	int patom_num[20];
	int plevel[20];
	for (i = -1; i < 2; i++)
	{
		for (j = -1; j < 2; j++)
		{
			for (k = -1; k < 2; k++)
			{
				for (ii = 0; ii < in.num; ii++)
				{
					//printf("%d %d %d %d\n", i,j,k,ii);
					temp_p[0] = in.rp[ii][0] + (i + l[0])*in.d[0][0] + (j + l[1])*in.d[1][0] + (k + l[2])*in.d[2][0];
					temp_p[1] = in.rp[ii][1] + (i + l[0])*in.d[0][1] + (j + l[1])*in.d[1][1] + (k + l[2])*in.d[2][1];
					temp_p[2] = in.rp[ii][2] + (i + l[0])*in.d[0][2] + (j + l[1])*in.d[1][2] + (k + l[2])*in.d[2][2];
					/*if (i==0&&j==0&&k==0)
					{
						printf("%lf\n", in.rp[ii][2]);
						printf("%d\t%d\t%d\n", l[0], l[1],l[2]);
						printf("%lf\t%lf\t%lf\n", temp_p[0], temp_p[1], temp_p[2]);
						printf("%lf\t%lf\t%lf\n", rp[0], rp[1], rp[2]);
						printf("%d\t%d\t%d\t%d\t%lf\n", i,j,k,ii,dis(temp_p, rp));
					}*/

					if (dis(temp_p, rp) < dist[type][in.type[ii]] && dis(temp_p, rp) > 0.1)
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
	near_num = num;
	near = (atom *)malloc(num * sizeof(atom));
	for (i = 0; i < num; i++)
	{
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
	}

	/*printf("Atom %d has %d neighbors\n", xuhao, num);
	 printf("%lf\t%lf\t%lf\n", rp[0], rp[1], rp[2]);
	 for (j=0;j<near_num;j++)
	 {
	 printf("%d\t%lf\t%lf\t%lf\n", j, near[j].rp[0], near[j].rp[1], near[j].rp[2]);
	 }*/
	return 0;
}

int atom::print_atom(cell in, int n, FILE *out)//生成一个文件tem1和tem2，用来输出表示找到的4层以内的原子种类及位置
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
					n = near[i].print_atom(in, n, out);
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
	if (near_num > 0)
	{
		for (i = 0; i < near_num; i++)
		{
			n = near[i].print_atom(in, n, out);
		}
	}

	if (near_num == 0)
	{
		return n;
	}
	return n;
}


cell::cell(char *name_atom) {//生成cell变量后进行初始化，读取相关信息
	FILE *in;
	in = fopen(name_atom, "r");
	if (in == NULL)
	{
		printf("No inputfile: %s\n", name_atom);
		return;
	}
	char str[500];
	fscanf(in, "%d", &num);
	fgets(str, 500, in);
	fgets(str, 100, in);
	int i, j;
	d = (double **)malloc(3 * sizeof(double *));
	for (i = 0; i < 3; i++)
	{
		d[i] = (double *)malloc(3 * sizeof(double));
	}
	p = (double **)malloc(num * sizeof(double *));
	for (i = 0; i < num; i++)
	{
		p[i] = (double *)malloc(3 * sizeof(double));
	}
	rp = (double **)malloc(num * sizeof(double *));
	for (i = 0; i < num; i++)
	{
		rp[i] = (double *)malloc(3 * sizeof(double));
	}
	type = (int *)malloc(num * sizeof(int));
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			fscanf(in, "%lf", &d[i][j]);
		}
	}
	fgets(str, 100, in);
	fgets(str, 100, in);
	for (i = 0; i < num; i++)
	{
		fscanf(in, "%d", &type[i]);
		fscanf(in, "%lf", &p[i][0]);
		fscanf(in, "%lf", &p[i][1]);
		fscanf(in, "%lf", &p[i][2]);
		fgets(str, 100, in);
	}
	for (i = 0; i < num; i++)
	{
		rp[i][0] = p[i][0] * d[0][0] + p[i][1] * d[1][0] + p[i][2] * d[2][0];
		rp[i][1] = p[i][0] * d[0][1] + p[i][1] * d[1][1] + p[i][2] * d[2][1];
		rp[i][2] = p[i][0] * d[0][2] + p[i][1] * d[1][2] + p[i][2] * d[2][2];
	}
	fclose(in);
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


