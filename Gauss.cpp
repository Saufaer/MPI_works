
///////////////////////////////////////Extra.h/////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <time.h>
#include <math.h> 
#include <mpi.h>
#include <iostream>

using namespace std;

void ShowResultP(double*& Result, int Size,int*KeyPosP)
{
	for(int i = 0; i < Size; i++)
	{
		cout<<endl << "x" <<i+1<<"="<< Result[KeyPosP[i]] ;
	}
}

void ShowResultS(double*& Result, int Size)
{
	for(int i = 0; i < Size; i++)
	{
		cout<<endl << "x" <<i+1<<"="<< Result[i] ;
	}
}
void ShowElMatrixS(double* Matrix, double* Vector,int Size,int*KeyPosS)
{
	for(int i = 0; i < Size; i++)
	{
		for(int j = 0; j < Size; j++)
		{			
			double intpartM;
			modf(Matrix[Size * KeyPosS[i] + j]*10000, &intpartM);//intpart = целая часть от первого аргумента
			cout << intpartM/10000 << " ";			
		}
		double intpartV;
		modf(Vector[KeyPosS[i]]*10000, &intpartV);
		cout <<" | "<< intpartV/10000 <<endl;		
	}
}

void ShowMatrix(double*& Matrix,double*& Vector, int RowNum,int Size)
{
	for(int i = 0; i < RowNum; i++)
	{
		for(int j = 0; j < Size; j++)
		{
			double intpartM;
			modf(Matrix[RowNum * i + j]*10000, &intpartM);

			cout << intpartM/10000 << " ";	
		}
		double intpartV;
		modf(Vector[i]*10000,&intpartV);
		cout << "| "<< intpartV/10000 <<endl ;
	}
}


void RandomInit(double* Matrix, double* Vector, int Size) 
{    
	//srand(time(0));
	for (int i=0; i<Size; i++)
	{    
		Vector[i] = rand()%10;     
		for (int j=0; j<Size; j++) 
		{         			 
			Matrix[i*Size+j] = rand()%10;     			
		}  
	} 
}



/////////////////////////////////////Serial.h///////////////////////////////////////////

using namespace std;
int* KeyPosS;  
int* KeyIterS; 

void EliminateGaussS(double* Matrix,double* Vector,int Size)
{     
	KeyPosS  = new int [Size];   
	KeyIterS = new int [Size];   
	for (int i=0; i<Size; i++)
	{     
		KeyIterS[i] = -1;  
	}   
	// поиск ведущей строки

	int KeyRow=-1;   
	for (int i=0; i<Size; i++)
	{     

		double MaxValue =  0;   

		for (int j=0; j<Size; j++) 
		{   
			if ((KeyIterS[j] == -1) && (fabs(Matrix[j*Size+i]) > MaxValue)) 
			{      
				KeyRow = j;       
				MaxValue = fabs(Matrix[j*Size+i]);   
			}  
		}   
		KeyPosS[i] = KeyRow;     
		KeyIterS[KeyRow] = i; 

		//преобразования
		double KeyValue, KeyFactor;   
		KeyValue = Matrix[KeyRow*Size+i]; 
		for (int k=0; k<Size; k++) 
		{     
			if (KeyIterS[k] == -1) 
			{     
				KeyFactor = Matrix[k*Size+i] / KeyValue; 
				for (int f=i; f<Size; f++)
				{         
					Matrix[k*Size + f] -= KeyFactor * Matrix[KeyRow*Size+f];  
				}       
				Vector[k] -= KeyFactor * Vector[KeyRow];   
			}   		
		} 
	} 
} 

void BackGaussS (double* Matrix, double* Vector,    double* Result, int Size) 
{  
	int RowIndex;  
	for (int i=Size-1; i>=0; i--)
	{    
		RowIndex = KeyPosS[i];   
		Result[i] = Vector[RowIndex]/Matrix[Size*RowIndex+i];
		for (int j=0; j<Size; j++) 
		{   
			if(KeyIterS[j]<i)
			{
				Vector[j] -= Matrix[j*Size+i]*Result[i];   
				Matrix[j*Size+i] = 0;   
			}
		}   
	} 
}


///////////////////////////////////////////Parallel.h/////////////////////////////////////////////////
using namespace std;
int ProcNum;            
int ProcRank;    

int* ProcInd; 
int* RowProcNum; 

int *KeyPosP; 
int *ProcKeyIter;                       

void Init (double* &Matrix, double* &Vector, double* &Result, double* &ProcRows, double* &ProcVector,double* &ProcResult, int &Size, int &RowNum)
{		

	if (ProcRank == 0) 
	{    
		cout<<"Size=";
		cin>>Size;
		Matrix = new double [Size*Size];    
		Vector = new double [Size];    
		Result = new double [Size];      
		RandomInit(Matrix, Vector, Size); 
	}   
	MPI_Bcast(&Size, 1, MPI_INT, 0, MPI_COMM_WORLD);   
	int WaitRows = Size;  

	for (int i=0; i<ProcRank; i++)  
	{
		WaitRows = WaitRows-WaitRows/(ProcNum-i);  
	}
	RowNum = WaitRows/(ProcNum-ProcRank);
	ProcRows = new double [RowNum*Size];   
	ProcVector = new double [RowNum];  
	ProcResult = new double [RowNum];  
	KeyPosP = new int [Size];      

	ProcKeyIter = new int [RowNum];  

	for (int i=0; i<RowNum; i++)   
	{
		ProcKeyIter[i] = -1;   
	}

	ProcInd = new int [ProcNum];   
	RowProcNum = new int [ProcNum];   
}


void Locate(double* Matrix, double* ProcRows, double* Vector,   double* ProcVector, int Size, int RowNum) 
{  
	int *Recvcount;   
	int *DisplsM;                    	
	int WaitRows=Size;
	DisplsM = new int [ProcNum];   
	Recvcount = new int [ProcNum];   

	RowNum = Size/ProcNum;  
	Recvcount[0] = RowNum*Size;  
	DisplsM[0] = 0;   

	int WaitR = Size;  
	ProcInd[0] = 0;   
	RowProcNum[0] = Size/ProcNum;   

	for (int i=1; i<ProcNum; i++)
	{    

		WaitRows -= RowNum;   
		RowNum = WaitRows/(ProcNum-i);    
		Recvcount[i] = RowNum*Size;  
		DisplsM[i] = DisplsM[i-1]+Recvcount[i-1];   

		WaitR -= RowProcNum[i-1];    
		RowProcNum[i] = WaitR/(ProcNum-i);   
		ProcInd[i] = ProcInd[i-1]+RowProcNum[i-1]; 
	}   

	MPI_Scatterv(Matrix, Recvcount, DisplsM, MPI_DOUBLE, ProcRows, Recvcount[ProcRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);  
	MPI_Scatterv(Vector, RowProcNum, ProcInd, MPI_DOUBLE, ProcVector, RowProcNum[ProcRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);    

}


void EliminateGaussP (double* ProcRows, double* ProcVector,   int Size, int RowNum)
{  
	int    KeyPos;   

	struct { double val; int rank; } ProcKey;
	struct { double MaxVal; int rankMax; }  Key;  


	double* KeyRow = new double [Size];
	double KeyCoeff;
	for (int i=0; i<Size; i++) 
	{    
		double Max = 0;     
		for (int j=0; j<RowNum; j++) 
		{      
			if ((ProcKeyIter[j] < 0 ) &&  (Max < fabs(ProcRows[j*Size+i]))) 
			{       
				Max = fabs(ProcRows[j*Size+i]);   
				KeyPos = j;        
			}   
		}    

		ProcKey.val = Max;    
		ProcKey.rank = ProcRank;   

		MPI_Allreduce(&ProcKey, &Key, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);  

		if ( ProcRank == Key.rankMax )
		{     
			ProcKeyIter[KeyPos]= i;   
			KeyPosP[i]= ProcInd[ProcRank] + KeyPos; 
		}  
		MPI_Bcast(&KeyPosP[i], 1, MPI_INT, Key.rankMax,  MPI_COMM_WORLD);     
		if ( ProcRank == Key.rankMax )
		{     

			for (int j=0; j<Size; j++) 
			{        
				KeyRow[j] = ProcRows[KeyPos*Size + j];  
			}      
			KeyCoeff = ProcVector[KeyPos];   
		}  
		MPI_Bcast(KeyRow, Size+1, MPI_DOUBLE, Key.rankMax,  MPI_COMM_WORLD); 
		MPI_Bcast(&KeyCoeff, 1, MPI_DOUBLE, Key.rankMax,  MPI_COMM_WORLD); 


		double coeff=0;    
		for (int e=0; e<RowNum; e++)
		{     
			if (ProcKeyIter[e] == -1)
			{     
				coeff = ProcRows[e*Size+i] / KeyRow[i];  
				for (int j=i; j<Size; j++) 
				{       
					ProcRows[e*Size + j] -= coeff*KeyRow[j];       
				}       
				ProcVector[e] -= coeff*KeyCoeff;     
			}   
		}     
	}
}



void BackGaussP (double* ProcRows, double* ProcVector,   double* ProcResult, int Size, int RowNum) 
{  
	int IterProcRank;   
	int IterKeyPos;    
	double IterResult;   
	for (int i=Size-1; i>=0; i--) 
	{   


		for (int f=0; f<ProcNum-1; f++)
		{    
			if ((ProcInd[f]<=KeyPosP[i]) && (KeyPosP[i]<ProcInd[f+1]))  
			{
				IterProcRank = f;  

			}
		}  
		if (KeyPosP[i] >= ProcInd[ProcNum-1]) 
		{
			IterProcRank = ProcNum-1; 
		}
		IterKeyPos = KeyPosP[i] - ProcInd[IterProcRank];


		if (ProcRank == IterProcRank)
		{      
			IterResult =   ProcVector[IterKeyPos]/ProcRows[IterKeyPos*Size+i]; 
			ProcResult[IterKeyPos] = IterResult;     
		}    


		MPI_Bcast(&IterResult, 1, MPI_DOUBLE, IterProcRank, MPI_COMM_WORLD);    

		for (int j=0; j<RowNum; j++) 
		{
			if ( ProcKeyIter[j] < i )
			{         
				ProcVector[j]-=ProcRows[j*Size + i] * IterResult;  			
			}  
		}
	}
}

//////////////////////////////////////////////////main.cpp////////////////////////////////////////////////


//#include "Extra.h"
//#include "Serial.h"
//#include "Parallel.h"

void main(int argc, char* argv[]) 
{  

	double* MatrixP,* MatrixS;        
	double* VectorP,* VectorS;          
	double* ResultP,* ResultS;        
	double *ProcRows;                 
	double *ProcVector;               
	double *ProcResult;               
	int Size;                       
	int RowNum;                    
	double startP, finishP, TimeP;  
	double startS, finishS,TimeS;   


	MPI_Init ( &argc, &argv );   
	MPI_Comm_rank ( MPI_COMM_WORLD, &ProcRank ); 
	MPI_Comm_size ( MPI_COMM_WORLD, &ProcNum );      

	Init(MatrixP, VectorP, ResultP, ProcRows, ProcVector, ProcResult, Size, RowNum);  

	if(ProcRank==0)
	{
		MatrixS = new double [Size*Size];
		VectorS = new double [Size]; 
		ResultS = new double [Size];   
		for (int i=0; i<Size; i++) 
		{   
			for (int j=0; j<Size; j++) 
			{ 
				MatrixS[i*Size+j]= MatrixP[i*Size+j];
			}
			VectorS[i]=VectorP[i];
			ResultS[i]=0;
		}  	
		if(Size<8)
		{
			cout<<"\n----------START MATRIX------------\n\n";

			ShowMatrix(MatrixS,VectorS, Size, Size);	
		}

		startS = MPI_Wtime(); 

		EliminateGaussS (MatrixS, VectorS, Size); 	

		BackGaussS (MatrixS, VectorS, ResultS, Size); 

		finishS = MPI_Wtime();

		TimeS = finishS-startS;

		if(Size<8)
		{
			cout<<"\n--------AFTER ELIMINATION---------\n\n";
			ShowElMatrixS(MatrixS,VectorS, Size, KeyPosS);
		}
	}

	startP = MPI_Wtime(); 

	Locate(MatrixP, ProcRows, VectorP, ProcVector, Size, RowNum); 

	EliminateGaussP (ProcRows, ProcVector, Size, RowNum);

	BackGaussP (ProcRows, ProcVector, ProcResult, Size, RowNum);

	MPI_Gatherv(ProcResult, RowProcNum[ProcRank], MPI_DOUBLE, ResultP,  RowProcNum, ProcInd, MPI_DOUBLE, 0, MPI_COMM_WORLD);    

	finishP = MPI_Wtime(); 

	TimeP = finishP-startP;  

	if (ProcRank == 0)
	{    

		if(Size<8)
		{
			cout<<"\n-------------RESULT--------------\n";
			cout<<"\nSerial:";	  
			ShowResultS(ResultS,Size);
			cout<<"\n\nParallel:";	
			ShowResultP(ResultP,Size,KeyPosP);
			cout<<"\n";
		}
		cout<<"\n-----------TEST RESULTS-----------\n";
		double check=0;
		for(int i=0;i<Size;i++)
		{
			check+=ResultP[KeyPosP[i]]-ResultS[i];	
		}
		cout<<"\nDifference in results: "<<check<<"\n";

		cout<<"\n----------TIME RESULTS------------\n";
		cout<<"\nSerial: "<< TimeS<<"\n";
		cout<<"\nParallel: "<< TimeP<<"\n";
		double jump=TimeS/TimeP;
		cout<<"\nJump: "<< jump<<"\n";
	} 
	if (ProcRank == 0) 
	{    
		delete [] MatrixP; 
		delete [] VectorP; 
		delete [] ResultP; 
		delete [] MatrixS; 
		delete [] VectorS; 
		delete [] ResultS; 
	}  
	delete [] ProcRows;   
	delete [] ProcVector;   
	delete [] ProcResult; 
	delete [] KeyPosP; 
	delete [] ProcKeyIter;   
	delete [] ProcInd;  
	delete [] RowProcNum;

	MPI_Finalize();
}