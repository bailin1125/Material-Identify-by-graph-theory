//
// compare_multifiles.cpp
//  compare structures and delete the same files
//
//  Created by 王志 on 2019/6/27.
//  Copyright © 2019年 王志. All rights reserved.
//
#include <omp.h>
#include <cstdio>
#include <cstring>
//#include "read_atom.hpp"
#include <unistd.h>
#include <time.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sstream>
#pragma warning(disable : 4996)
using namespace std;
int main(int argc, const char *argv[])
{
	string data_path = "/share/home/wangz/2d_search/ridus_cut/analyse_1020/all/bulk/config/";
	pid_t status;
	const int NUM = 5743;
	string* name = new string[NUM];
	int stru_jump[NUM] = { 0 };//用来控制把谁扔出去,0表示未见重复	
	int i = 0, j = 0, k = 0;
	ifstream fin;
	fin.open("bulk_name", ios::in);
	if (!fin.is_open())
	{
		cout << "i can not find the name txt!" << endl;
		cin.get();
	}
	while (fin.good() && fin.peek() != EOF)
	{
		fin >> name[i];
		name[i] = data_path + name[i];
		//name[i] = "./data/" + name[i];
		i++;
	}
	fin.close();
	string atom[25][2];
	//开始进行命令拼接
	char command_temp[200] = "./compare.out ";
	int **compare = new int *[NUM];
	for (i = 0; i < NUM; i++)
	{
		compare[i] = new int[NUM];
	}
	for (i = 0; i < NUM; i++)
	{
		for (j = 0; j < NUM; j++)
		{
			compare[i][j] = 0;
			if (i == j)
				compare[i][j] = 1;
		}
	}
	//int compare[NUM][NUM] = { {0} };
	int n = 0;
	char command[25][300];
	//从这里开始并行的进行图同构去重

#pragma omp parallel for  private(n,j) shared(stru_jump) schedule(dynamic)
	for (i = 0; i < 244; i++)
	{
		//cout << "now is the " << omp_get_thread_num() << " fork" << endl;
		//if (name[i] == "" || name[i] == "\0" || name[i] == " ")
		//	break;
		if (i % 100 == 0)
		{
			cout << "has gone " << i << " files" << endl;
		}


#pragma omp critical
		{
			n = omp_get_thread_num();
			//fscanf(name, "%s", now_name[n]);
			//cout << "the threads xuhao  is:" << n << endl;
			//cout << "the name is:" << now_name[n] << endl;
			//cin.get();
		}

		for (j = i + 1; j < NUM; j++)
		{
			if (name[j] == "")
				break;
			atom[n][0] = name[i];
			atom[n][1] = name[j];
			if (stru_jump[j] == 0)
			{
				strcpy(command[n], command_temp);
				strcat(command[n], atom[n][0].c_str());
				strcat(command[n], " ");
				//strcat(command, "./files/");
				strcat(command[n], atom[n][1].c_str());
				//cout << command << endl;
				//cin.get();
				status = system(command[n]);
				if (WEXITSTATUS(status) != 0)
				{
					atom[n][0] = name[j];
					atom[n][1] = name[i];
					strcpy(command[n], command_temp);
					strcat(command[n], atom[n][0].c_str());
					strcat(command[n], " ");
					//strcat(command, "./files/");
					strcat(command[n], atom[n][1].c_str());
					status = system(command[n]);
					if (WEXITSTATUS(status) != 0)
					{
						stru_jump[j] = 1;
						compare[i][j] = compare[j][i] = 1;
					}
				}
				else
				{
					compare[i][j] = compare[j][i] = 0;
				}
				//cout << "result is:" << compare[i][j] << endl;
			}
		}
	}

	//这样的话每个节点的任务只是出输出要丢的任务名称
	ofstream fout;
	fout.open("jump_str_name2.txt", ios::out);
	for (i = 0; i < NUM; i++)
	{
		if (stru_jump[i] == 1)
		{
			fout << name[i];
			fout << endl;
		}
	}
	fout.close();



	//cout << "stop3" << endl;
	//到这里完成了全部的比对工作
	//对每个节点输出自己的文件
	/*ofstream fout;
	fout.open("compare_NUM.txt", ios::out);
	for (i = 0; i < NUM; i++)
	{
		for (j = 0; j < NUM; j++)
		{
			fout << compare[i][j];
			fout << "\t";
		}
		fout << endl;
	}
	fout.close();*/

	////下面需要根据比对结果，确立重复文件额名单，并且删除重复文件
	//int num = i;//这个就是文件个数
	////cout << "the file geshu is:" << num << endl;
	//int *flag = new int[NUM];
	//for (i = 0; i < NUM; i++)
	//{
	//	flag[i] = 0;
	//}
	////int flag[NUM] = { 0 };//这个用来标记说是不是每个被用到了
	///*for (i = 0; i < NUM; i++)
	//{
	//	cout << flag[i] << endl;
	//}*/
	//string ** compare_result = new string *[NUM];
	//for (i = 0; i < NUM; i++)
	//{
	//	compare_result[i] = new string[30];
	//}
	////string compare_result[NUM][100];
	//for (i = 0; i < num; i++)
	//{
	//	if (flag[i] != 0)
	//		continue;
	//	for (j = 0; j < num; j++)
	//	{
	//		if (compare[i][j] != 0 && compare[j][i] != 0)
	//		{
	//			flag[j] = 1;
	//			compare_result[i][k] = name[j];
	//			k++;
	//			if (k >= 30)
	//			{
	//				cout << "oo,maybe 30 is too small!" << endl;
	//				cin.get();
	//			}
	//		}
	//	}
	//	k = 0;
	//}
	////cout << "stop4" << endl;
	//ofstream fout;
	//fout.open("result_compare.txt", ios::out);
	//for (i = 0; i < num; i++)
	//{
	//	if (compare_result[i][0] == "")
	//	{
	//		continue;
	//	}
	//	
	//	for (j = 0; j < num; j++)
	//	{
	//		if (compare_result[i][j] == "" )
	//			break;
	//		fout << compare_result[i][j];
	//		fout << "\t";
	//	}
	//	fout << endl;
	//	
	//}
	//fout.close();


	//char rm_command[100] = "rm -f ";
	//
	//for (i = 0; i < num; i++)
	//{
	//	if (compare_result[i][0] == "" || compare_result[i][1] == "")
	//		continue;
	//	for (j = 1; j < num; j++)
	//	{
	//		if (compare_result[i][j] == "")
	//			break;
	//		else
	//		{
	//			char rm_temp[100];
	//			strcpy(rm_temp, rm_command);
	//			strcat(rm_temp, compare_result[i][j].c_str());
	//			cout << rm_temp << endl;
	//			system(rm_temp);
	//		}
	//	}
	//}

	cout << "all total work done!" << endl;
	for (i = 0; i < NUM; i++)
	{
		delete[]compare[i];
		//delete[]compare_result[i];
	}
	//delete[]compare_result;
	delete[]compare;
	delete[]name;
	//delete[]flag;
	//cin.get();
	cout << "all total work done!" << endl;
	return 0;



}