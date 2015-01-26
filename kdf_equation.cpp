//Author: Yavorskiy Alexandr

#include <vcl.h>
#pragma hdrstop
#include <iostream>
#include <math.h>
#include "conio.h"
#include <fstream>

#include "ololo.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma link "iComponent"
#pragma link "iCustomComponent"
#pragma link "iPlotComponent"
#pragma link "iVCLComponent"
#pragma link "iXYPlot"
#pragma link "iAnalogOutput"
#pragma link "iEditCustom"
#pragma resource "*.dfm"
using namespace std;
TForm1 *Form1;
//---------------------------------------------------------------------------
__fastcall TForm1::TForm1(TComponent* Owner)
        : TForm(Owner)
{
}
//---------------------------------------------------------------------------


void __fastcall TForm1::Button1Click(TObject *Sender)
{
        double t1;
        t1= iAnalogOutput1->Value;
        double c=0.00005,dt=0.0004; //параметр скоростного уравнения
						//КдВ (3 солитона)
	int N=200;		//число узлов пространственной
						//сетки
	//!шаг интегрирования по времени
	//!вспомогательные переменные
	int N1; //N1=N-1
			//число Пи
		//постоянные разностной схемы
	double dtdx1;	//dtdx1=dt/(2*dx)
	double dtdx2;	//dtdx2=c*dt/(2*dx**3)
	double dx;		//шаг по координате dx=1/N
		//таблица сдвига индекса
	int dj[5]={-2,-1,0,1,2};
		//массив ближайших индексов
	int jn[5];
		//массивы для хранения решения уравнения КдВ
	double *v;
	v=(double *)malloc((N+3)*sizeof(double));//текущее решение
	int i,j;	//индекс пространственного узла
	double x;	//пространственная координата
		//параметры начального решения уравнения:
	double xm=0.5;	//центр горба
	double A=0.5;	//амплитуда горба
	double s=0.08;	//ширина горба
		//вычисление постоянных разностной схемы интегрирования
		//и вспомогательных переменных
	dx = 1./N;	//шаг интегрирования по координате
	dtdx1 = dt/(2*dx);
	dtdx2 = c*dt/(2*pow(dx,3));
	N1    = N-1;
		//располагаем массивы для хранения решения
		//в оперативной памяти
		//определяем начальное значение для скорости v(x,0)
	for (j=0;j<=N;j++)
	{
		x=j*dx;
		v[j]=A*exp(-pow(x-xm,2)/(2*pow(s,2)));
	}
//далее численное интегрирование ур-я КдФ:
	double t, vn1,vd;
	//нелинейный и дисперсионный вклад
	//в скорость
	for (i=0;;i++)	//цикл интегрирования по времени
	{
		t=i*dt;
                if (t>t1) break;
		for (j=0;j<=N;j++)	//цикл интегрирования по координате
		{
			//вычисляем значения ближайших индексов
			//пространственной сетки
			if(j==0)
			{
				jn[0]=N1; jn[1]=N;
				for (int iter=2;iter<=4;iter++)
				{
					jn[iter]=j+dj[iter];
				}
			}
			else if(j==1)
			{
				jn[0]=N;
				for (int iter=1;iter<=4;iter++)
				{
					jn[iter]=j+dj[iter];
				}
			}
			else if(j==N-1)
			{
				for(int iter=0;iter<=3;iter++)
				{
					jn[iter]=j+dj[iter];
				}
				jn[4]=0;
			}
			else if(j==N)
			{
				for(int iter=0;iter<=2;iter++)
				{
					jn[iter]=j+dj[iter];
				}
			        jn[3]=0; jn[4]=1;
			}
			else
			{
				for(int iter=0;iter<=4;iter++)
				{
					jn[iter]=j+dj[iter];
				}
			}
			//вычисляем значения скоростей в узлах
			//пространственной сетки
			vn1 =-dtdx1*v[jn[2]]*(v[jn[3]]-v[jn[1]]);
			vd  =-dtdx2*(v[jn[4]]-2.0*v[jn[3]]+2.0*v[jn[1]]-v[jn[0]]);
			v[j]=v[jn[2]]+vn1+vd;   //переписываем текущее решение в старое
		}
        }
         for (j=0;j<=N;j++)
                {
                        x=j*dx;
                        iXYPlot1->Channel[0]->AddXY(x,v[j]);
                }
}
//---------------------------------------------------------------------------

void __fastcall TForm1::Gohome1Click(TObject *Sender)
{
exit(0);        
}
//---------------------------------------------------------------------------

void __fastcall TForm1::Button2Click(TObject *Sender)
{
        iXYPlot1->ClearAllData();
}
//---------------------------------------------------------------------------
