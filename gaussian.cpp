// reading a text file
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
using namespace std;

int readtxt(float** a,float* b){  // parameters are one multidimensional and one one-dimensional dynamically allocatoed matrices
	
   string line;
   int n=0;
   ifstream myfile ("A.txt");    // opening the first .txt file
   if (myfile.is_open())
   {
     while ( getline (myfile,line) ) 
     {
       n++;      // indicating the row and the column(which are equal since it is square matrix) number of the matrix in .txt file
     }
     myfile.close();
   }

   else cout << "Unable to open file";

   ifstream second ("A.txt");

   for( int k=0;k<n;k++){
   	for(int j=0;j<n;j++){
	    second>>a[k][j];   // copying the values in .txt file to matrix a
	  }
   }

    ifstream B ("B.txt");
    for(int i=0;i<n;i++){
    	B>>b[i];            // copying the values in .txt file to matrix b
    }
	return n;   // row number is returned, since it will be frequently used by other functions in the main
}

void get_augmented(float** augm,float** mat1,float* mat2,int row){  // dynamically allocated matrices used for parameters
	
   for(int i=0;i<row;i++){
   	  for(int j=0;j<row+1;j++){
   	  	if(j<row)
   	  	augm[i][j]=mat1[i][j];   // first n columns of augmented matrix consists of the columns of a 
   	  	else
   	  	augm[i][j]=mat2[i]; 	  // last column of augmmented matrix is taken from the b	
		}
   }
}

void gaussian_elimination(float** augm,int row){
    float rate;
     int holder;
    
    for(int m=0;m<row;m++){
    
       for(int l=m;l<m;l++){
       	if(abs(augm[l][m])>abs(augm[m][m])){  // for partial pivoting magnitudes of rows are compared
       		for(int p=0;p<m;p++){             
       			holder=augm[l][p];            // interchange of rows, if needed
       			augm[l][p]=augm[m][p];        
       			augm[m][p]= holder;
			   }
		   }
	   }
    	
    	for(int k=m;k<row-1;k++){            // loop for gaussian elimination
    		rate=augm[k+1][m]/augm[m][m];    // calculating ratio of the pivot elements and the entries above them
    		for(int t=0;t<row+1;t++){
      		augm[k+1][t]=augm[k+1][t]-rate*augm[m][t];   // making the entries above pivot 0
      	}	 	
	}
	}
}

int singularity_test(float** aug,int row){
	float e=10;
	float mach;             // calculating the machine epsilon for singularity detection
    while ((1+e) != 1) 
    { 
        // copying the value of e to machine epsilon
        mach= e; 
  
        // dividing epsilon by 2 to get more accurate result step by step, as in Zenon's paradox.
        e/=2; 
    } 
    
	float control=0;
    
    		if(abs(aug[row-1][row-1])>mach){  
    			control++;    // if matrix is singular, last pivot at the right corner would be less than machine epsilon in magnitude
			}                 // if last pivot is less than machine e, control variable will stay zero indicating singularity for the main function
		
	return control;
}

void makepivots1 (float** aug,int row){

    for(int i=0;i<row;i++){
    	float divider=aug[i][i];
    	for(int j=0;j<row+1;j++){
    		aug[i][j]=aug[i][j]/divider;   // dividing every column to its pivot to make all pivots 1
		}                                  // for the ease of back substitution
	}
}

void seperate_augmented (float** aug,float** a,float* b,int row){

	for(int i=0;i<row;i++){
		for(int j=0;j<row;j++){ 
			a[i][j]=aug[i][j];      // first n columnns of augmented matrix will be new version of matrix a
    }
	}

	for(int i=0;i<row;i++){
		for(int j=row;j<row+1;j++){
			b[i]=aug[i][j];          // last column of augmented will be held by matrix b
	}
}
}

void back_substitution(float* x,float** a,float* b, int row){
									
	for(int k=row-1;k>=0;k--){	
	x[k]=b[k];                     // since pivots are 1, last entry of x matrix directly will be equal to b
	for(int i=k-1;i>=0;i--){
		b[i]=b[i]-x[k]*a[i][k];    // portions of x calculated above(with its appropriate coefficients) is substracted from other matrix b elements
	}	
}
}

float condnum_calculator(float** a){
	float condnum,coeff,maxcol_a,maxrow_a;

	maxcol_a=abs(a[0][0])+abs(a[1][0]);      
	if(abs(a[0][1])+abs(a[1][1])>maxcol_a){    // comparing absolute column sums of the matrix a
		maxcol_a=abs(a[0][1])+abs(a[1][1]);    // holding the value of greatest absolute column sum of matrix 
	}
	
	maxrow_a=abs(a[0][0])+abs(a[0][1]);         // maximum absolute column sum of inverse of a will actually be equal to max absolute row sums of matrix a
	if(abs(a[1][1])+abs(a[1][0])>maxrow_a){
		maxrow_a=abs(a[1][1])+abs(a[1][0]);     
	}

	coeff=1/(a[0][0]*a[1][1]-a[0][1]*a[1][0]);    // calculating determinant 
	 
	condnum= coeff*maxcol_a*maxrow_a;  // using condition number formula (cond(A) = ||A|| · ||A^-1||)
	
	return condnum;    // returning cond num for 2x2 cases
}

int main () {
    int n;
    
    float **a=new float*[n];      // defining dynamically allocated multidimensional array for firrst .txt file
    for(int i=0;i<n;i++){
   		a[i]=new float[n];
    }
   
    float*b;               // defining dynamically allocated array for second .txt file
    b=new float[n];

	n=readtxt(a,b);   // function for taking the values of .txt file and holding them in matrices and returning the number of rows
	
	if(n==2){         // showing cond num for 2x2 matrices
	float condition;
		condition=condnum_calculator(a);
		cout<<"The condition numbers using 1-norm and infinity-norm are same and is equal to "<< condition << endl;
    }
    
	float **aug=new float*[n];    // defining dynamically allocated multidimensional array for augmented matrix
    for(int i=0;i<n+1;i++){
   		aug[i]=new float[n];
    }

    get_augmented(aug,a,b,n);     
    
    gaussian_elimination(aug,n);	

		float control=0;
		control=singularity_test(aug,n);   // unchanged control variable indicates singulariy

		if (control==0){
		cout<<"Error! Matrix is singular";
		exit;
     	}
     	
    makepivots1(aug,n); 	
     	
	seperate_augmented(aug,a,b,n);
	
    float*x;                       // defining dynamically allocated array for solution matrix x
    x=new float[n];      
	
	back_substitution(x,a,b,n);
	
	 ofstream mytext;               // opening .txt file to write in the solution
     mytext.open ("x.txt");
   
     for(int i=0;i<n;i++){
    	mytext<<"x"<<i<<": "<<x[i]<<endl;    // writing values of matrix is one by one to x.txt file
	 }
      mytext.close();

			
   return 0;
}

