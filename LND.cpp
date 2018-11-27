#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>
#include <string>

using namespace std;

void Init(int Size,int countMatr,int density,int **Matr,int **MatrS,int *traceP,int *traceS )
{
	for (int i = 0; i < Size; i++)
	{
		traceP[i] = i;
		traceS[i]=i;
	}		
	for (int i = 0; i < Size; i++)
	{
		for (int j = i + 1; j < Size; j++)
		{
			if (i != j)
			{
				MatrS[i][j] = MatrS[j][i] =	Matr[i][j] = Matr[j][i] = (rand() % 100 < density ? (1) : 0);// в 75% - дуги
			}
			else
			{
				MatrS[i][j] = MatrS[j][i]= Matr[i][j] = Matr[j][i] = 1;
			}
		}
	}
	for (int i = 0; i < Size; i++)
	{
		MatrS[i][i] = Matr[i][i] = 0;
	}
}
void PrintLastMatr(int **Matr, int size)
{
	for (int i = 0; i < size; i++) 
	{
		for (int j = 0; j < size; j++)
		{
			cout << Matr[i][j] << " ";
		}
		cout << endl;

	}
	cout << endl;	
}
void Cut(int **Matr, int * trace, int size, int countMatr, bool last,int rank)
{
	if(size<60)//вывод матриц + запись данных для проверки корректности
	{
		PrintLastMatr(Matr,size);


		if(rank==0)
		{

			FILE *f;
			string Buf;
			f = fopen("check.txt","w");
			ofstream fout1 ("check.txt"); 		
			for (int i = 0; i <size; i++) 
			{  		
				for (int j = 0; j < size; j++)
				{		
					fout1 << Matr[i][j];
					fclose(f);
				}		
			}				
		}
	}	


	int *mas = new int[size];
	int *part = new int[size];

	part[0] = 0;
	mas[0] = 0;
	for (int i = 1; i < size; i++)
	{
		mas[i] = 1;
		part[i]= -1;
	}

	int k = 1;//номер первой используемой вершины
	int l = 0;

	//первое нумерование вершин (поиск граничной вершины)

	for (int i=1; i < size;) {
		for (int j = 0; j < size; j++) {
			if (l == i)
				break;
			if (Matr[part[l]][j] == 1 && (mas[j] > mas[part[l]] + 1))
			{
				part[i] = j;
				mas[j] = mas[part[l]] + 1;
				i++;
			}
		}
		if (l == i)
		{
			break;
		}
		else
		{
			l++;
		}
		k=i;
	}

	int start = part[k - 1];

	for (int i = 0; i < size; i++)
	{
		part[i] = i;
	}

	int temp = part[0];
	part[start] = part[0];
	part[0] = temp;


	mas[start] = 0; 
	for (int i = 0; i < size; i++)
	{
		mas[i] = 1;
	}

	l = 0;

	int SizeSubMatr = size * (countMatr / 2) / countMatr;

	//второе нумерование вершин (оптимальное, начинающееся от граничной вершины)
	for (int i=1; i < SizeSubMatr;) 
	{
		for (int j = 0; j < size; j++)
			if (Matr[part[l]][part[j]] == 1 && (mas[part[j]] > mas[part[l]] + 1))
			{
				mas[part[j]] = mas[part[l]] + 1;

				int tmp = part[i];
				part[j] = part[i];
				part[i] = tmp;

				i++;
			}
			if (l == i)
			{
				break;
			}
			else
			{
				l++;
			}
	}

	int tmp;
	//перезапись наличия дуг(значений матрицы смежности)
	for (int i = 0; i < SizeSubMatr; i++)
	{
		//одну половину матрицы смежности
		for (int j = 0; j < size; j++)
		{
			tmp = Matr[i][j];
			Matr[i][j] = Matr[part[i]][j];
			Matr[part[i]][j] = tmp;
		}
		//вторую (симметричную)
		for (int j = 0; j < size; j++)
		{
			tmp = Matr[j][i];
			Matr[j][i] = Matr[j][part[i]];
			Matr[j][part[i]] = tmp;
		}
	}

	//конечная операция на итерации алгоритма

	//записываем вектор вершин для перестановки в порядке формирования двух графов	
	for (int i = 0; i < size; i++)
	{
		part[i] = trace[part[i]];
	}
	//перезаписываем вектор вершин для перестановки
	for (int i = 0; i < size; i++)
	{
		trace[i] = part[i];
	}		

	//в итоге получим один вектор в котором содержится новый порядок создания 2 матриц смежности

	//если нужно делить до конца без сторонних перераспределений (последовательно на одном процессе)
	//функция работает рекурсивно		 
	if (last) //если сработал сигнал, что делим до конца 
	{ 
		


		int **Submatr = 0;//создаем вторую парную матрицу
		//её парой будет часть исходной 

		if (countMatr > 2)//если надо получить более 2х матриц смежности
		{
			Submatr = new int*[size - SizeSubMatr];
			for (int i = SizeSubMatr; i < size; i++)
			{
				Submatr[i - SizeSubMatr] = &Matr[i][SizeSubMatr];//записываем нужную часть от исходной матрицы
			}
		}


		if (countMatr == 3)//если нужно произвести ровно 1 деление (в итоге получим 2 новых графа)
		{
			Cut(Submatr, &trace[SizeSubMatr], size - SizeSubMatr, 2, last,rank);//делим только одну матрицу (вторая уже разделена)

		}

		if (countMatr > 3)
		{
			//продолжаем делить обе
			Cut(Matr, trace, SizeSubMatr, countMatr/2, last,rank);
			Cut(Submatr, &trace[SizeSubMatr], size - SizeSubMatr, countMatr - (countMatr/2), last,rank);
		}
	}
}
void Tree(int **  a, int n,int k,int*trace,  int root,  int h) {
	
	int rank; 
	int ProcNum; 
	int size = (int)(powl(2, h)); 
	int unnatureSize = (int)(powl(2, h)) + root; 

	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int unnatureRank = rank, unnatureProcNum = ProcNum + root;
	unnatureRank = rank - root;

	if ((unnatureRank >= 0) && (unnatureRank < size)) {
		int unnatureReceiverRank = unnatureRank + unnatureSize , receiverRank = rank + size;
	
		if (unnatureReceiverRank < unnatureProcNum) {

			if (receiverRank < ProcNum) {
							
	// отправляем часть исходной матрицы на обработку другим процессом
		    int n1 = n * (k / 2) / k;
			int k1 = k / 2;
			
			MPI_Send(&n1, 1, MPI_INT,receiverRank, 0, MPI_COMM_WORLD);
			MPI_Send(&k1, 1, MPI_INT, receiverRank, 0, MPI_COMM_WORLD);			
			for (int i = 0; i < n; i++)
			{		
			MPI_Send(a[i], n1, MPI_INT, receiverRank, 0, MPI_COMM_WORLD);
			}			
			MPI_Send(trace, n1, MPI_INT,receiverRank, 0, MPI_COMM_WORLD);
			

			//перезапись текущей матрицы в процессе
			//заменяем её на остаточную половину
			int n2 = n - n1;
			int k2 = k - k / 2;
			int **b = new int*[n - n1];
			int* buf=new int[n2];
		for (int i = n * (k / 2) / k; i < n; i++)
		{
			buf[i - n * (k / 2) / k]=trace[i - n1];
			b[i - n * (k / 2) / k] = &a[i][n1];
		}
	
trace = new int[n2];
  a = new int*[n2];
  for (int i = 0; i < n2; i++)
{ 
	trace[i]=buf[i];
a[i] = new int[n2];
  a[i]=b[i];
  }
		n=n2;
		k=k2;

		//если распределение завершено

				if(ProcNum==size*2)
			{	
				cout << "Procrank = " << rank << " Size = " << n << " -> " << k << " matrices "<<endl;
				Cut(a, trace, n, k,true,rank);// обработка самостоятельно имеющегося куска
				return;
			}
			
			}
			

		}
		else {
			return;
		}
	}
	else {
			if ((unnatureRank >= size) && (unnatureRank < 2 * size))
			{
				int senderRank = rank - size;// -3
				
				if (senderRank >= 0 && ( rank > root))
				{
					
			MPI_Recv(&n, 1, MPI_INT,  senderRank, 0, MPI_COMM_WORLD,MPI_STATUSES_IGNORE);
			MPI_Recv(&k, 1, MPI_INT,  senderRank, 0, MPI_COMM_WORLD,MPI_STATUSES_IGNORE);
			
			a = new int*[n];
			for (int i = 0; i < n; i++)
			{a[i] = new int[n];}
			trace = new int[n];
			for (int i = 0; i < n; i++)
			{		
			MPI_Recv(a[i], n, MPI_INT,  senderRank, 0, MPI_COMM_WORLD,MPI_STATUSES_IGNORE);
			}
			
			MPI_Recv(trace, n, MPI_INT, senderRank, 0, MPI_COMM_WORLD,MPI_STATUSES_IGNORE);
		
		
						if(ProcNum==size*2)						
			{				
				cout << "Procrank = " << rank << " Size = " << n << " -> " << k << " matrices "<<endl;	
				Cut(a, trace, n, k,true,rank);			
				return;
			}						
				}
			}
	}
	
	h++;//высота дерева распределения

	if (size > sqrtl(ProcNum))
	{				
		return;
	}
	
	Tree(a, n,k,trace , root, h);
	
}

int main(int argc, char* argv[])
{
	int Procrank;
	int ProcNum;

	int Size;//число вершин графа

	int countMatr;//число разбиений

	//вектора перестановок
	int *trace = NULL;
	int *traceS= NULL;

	double t1=0, t2=0;
	double t1s=0, t2s=0;

	int density;//плотность дуг

	//матрицы смежности
	int **Matr = NULL;
	int **MatrS= NULL;

	int sizeCheck=0;
	int sizeCheckS=0;

	string CheckBuf;
	string CheckBufS;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &Procrank);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);

		if (argc > 1)
		{	Size = atoi(argv[1]);}
		else
		{	Size = 12000;}
		if (argc > 2)
		{	countMatr = atoi(argv[2]);}
		else
		{	countMatr = 255;}
		if (argc > 3)
		{density = atoi(argv[3]);}
		else
			{density = 75;}

		if (Procrank == 0)
	{ 
		trace = new int[Size];
		traceS = new int[Size];

		Matr = new int*[Size];
		MatrS = new int*[Size];
		for (int i = 0; i <Size; i++)
			{
		Matr[i] = new int[Size];
	    MatrS[i] = new int[Size];
		}
		Init(Size,countMatr,density,Matr,MatrS,trace,traceS);
		
	t1s = MPI_Wtime();

	Cut(MatrS, traceS, Size, countMatr, true,0);

	t2s = MPI_Wtime();
if (Size<60)
{
FILE *t;
string bufS;
t = fopen("check.txt","r");
ifstream fin1 ("check.txt");
fin1 >> CheckBufS;
fclose(t);
}
	



	
		
	t1 = MPI_Wtime();
	}			

Tree(Matr,Size,countMatr,trace , 0, 0);	


	MPI_Barrier(MPI_COMM_WORLD);
	if (Procrank == 0)
	{
		t2 = MPI_Wtime();

	cout<< endl << "Size of adjacency matrix: " <<Size << endl;
	cout << "Number of dividual matrices: " << countMatr << endl; 
	cout << "Number of processes: " << ProcNum << endl;
	cout<<"Density of graph: "<<density<<"%"<< endl;

	cout<< endl << "Time Serial: "<<t2s-t1s<<endl;	
	cout << endl << "Time Parallel: " << t2 - t1 << endl;
	cout << endl <<"Acceleration: "<<(t2s-t1s)/(t2-t1)<<endl;

		delete[] *Matr;
		delete[]*MatrS;

if(Size<60)
{
FILE *t;
string bufS;
t = fopen("check.txt","r");
ifstream fin1 ("check.txt");
fin1 >> CheckBuf;

fclose(t);

cout << endl <<"Last graph of Serial LND:  "<< CheckBufS << endl;
cout << endl <<"Last graph of Parallel LND:" << CheckBuf << endl;

if(CheckBufS==CheckBuf){cout<< endl<<"Algorithm is correct";}
}

	
	}
	MPI_Finalize();
	return 0;
}

