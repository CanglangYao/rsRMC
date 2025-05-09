#include <iostream>
#include <fstream>
#include <math.h>  
#include <string>
#include <iomanip>
using namespace std;

void allocatetwoD(double **&arr, int m, int n ){
    arr = new double*[m];
    for (int r = 0; r<m; r++){
        arr[r]=new double[n];
	for (int c = 0; c <n; c++) {
            arr[r][c]=0.0;
	}
    }
}

void allocateoneD(double *&arr, int m){
    arr = new double[m];
}

void deallocatetwoD(double **&arr, int m){
    for (int r = 0; r<m; r++){
        delete[] arr[r];
    }
    delete[] arr;
    arr = nullptr;
}

void deallocateoneD(double *&arr) {
    delete[] arr;
    arr = nullptr;
}


int main(int argc, char *argv[]) {
	  double Qdamp=0.03;
          std::string s(argv[1]);
	  int atom_num = 384; 
	  std::ifstream file( argv[1] );
	  double **geom = nullptr;
	  allocatetwoD(geom, atom_num, 4);
    for (int r = 0; r < atom_num; r++) //Outer loop for rows
    {
        for (int c = 0; c < 4; c++) //inner loop for columns
        {
          file >> geom[r][c];  //Take input from file and put into myArray
        }
    }

    double table[4][9]={
	            {2.31000, 20.8439, 1.02000, 10.2075, 1.58860, 0.568700, 0.865000, 51.6512, 0.215600},
		    {12.2126, 0.005700, 3.13220, 9.89330, 2.01250, 28.9975, 1.16630, 0.582600, -11.529},
	            {20.1472, 4.34700, 18.9949, 0.381400, 7.51380, 27.7660, 2.27350, 66.8776, 4.07120},
                    {31.0617, 0.690200, 13.0637, 2.35760, 18.4420, 8.61800, 5.96960, 47.2579, 13.4118}
                    };
    
    double la=25.24;
    double lb=25.24;
    double lc=25.24;
    double rou = atom_num/la/lb/lc;
    double pi = 3.1415926;
    int ix, iy, iz, nd, w_id;
    double x, y, z;
    int N=280;
    int type = 4;
    double step = 0.02;
    int cut = 630;
    double d;
    double **Fq=nullptr;
    double **sf=nullptr;
    double **reduced=nullptr;
    double **spec=nullptr;
    allocatetwoD(spec, 3000, type*type);
    allocatetwoD(Fq, N, 2);
    allocatetwoD(sf, N, type);
    allocatetwoD(reduced, cut, 2);
    double ave_sf;
    for (int ii=0; ii < atom_num-1; ii=ii+1){
	for (int jj =ii+1; jj < atom_num; jj=jj+1){
	    ix=(int) ((geom[ii][1]-geom[jj][1])/la*2);
	    iy=(int) ((geom[ii][2]-geom[jj][2])/lb*2);
	    iz=(int) ((geom[ii][3]-geom[jj][3])/lc*2);
	    ix=ix+(int) (-0.5*ix);
	    iy=iy+(int) (-0.5*iy);
	    iz=iz+(int) (-0.5*iz);
	    x=(geom[ii][1]-geom[jj][1])-ix*la;
	    y=(geom[ii][2]-geom[jj][2])-iy*lb;
	    z=(geom[ii][3]-geom[jj][3])-iz*lc;
	    d=sqrt(x*x+y*y+z*z);
	    nd = (int) (d/step);
	    w_id=(int) geom[ii][0]*type+(int) geom[jj][0];
	    spec[nd][w_id]=spec[nd][w_id]+1;
	 }
    }

    for (int nq = 0; nq < N; nq++){
 	double **Gr = nullptr;
	allocatetwoD(Gr, 3000, 2);
	double *weight= new double[type*type];
	for(int r = 0; r<3000; r++){
	    Gr[r][0] = r*step+0.01;
	} 
        Fq[nq][0]=(nq)/((double) 10);
	double q =(nq)/((double) (10*4*pi));
        for (int k = 0; k < type; k++){
            for (int l = 0; l < 4; l++){
                sf[nq][k]=sf[nq][k]+table[k][l*2]*exp(-table[k][l*2+1]*q*q);
            }
            sf[nq][k]=sf[nq][k]+table[k][8];
        }
	sf[nq][0]=sf[nq][0]*1.5;
	sf[nq][1]=sf[nq][1]*1.4286;
	ave_sf=(sf[nq][0]+sf[nq][1]+sf[nq][2]*3+sf[nq][3])/6.0;
        for (int i = 0; i< type; i++){
	    for (int j=0; j <type; j++){
		weight[i*type+j]=sf[nq][i]*sf[nq][j]/ave_sf/ave_sf;
	    }
	}	
        for (int r=0; r< 3000; r++){
	    for (int c = 0; c<type*type; c++){
	        Gr[r][1]=Gr[r][1]+spec[r][c]*weight[c]*2;
	    }
            Gr[r][1]=Gr[r][1]/atom_num/step;
	    Gr[r][1]=(Gr[r][1]/Gr[r][0]-4*pi*Gr[r][0]*rou);
	}
        delete[] weight;
	for (int i=0; i < cut; i++){
	    Fq[nq][1]=Fq[nq][1]+Gr[i][1]*sin(Fq[nq][0]*Gr[i][0])*step;
	} 
	deallocatetwoD(Gr, 3000);
    }
    double dq = Fq[1][0]-Fq[0][0];
    for(int i =0; i< cut; i++){
        reduced[i][0]=i*step;
	for(int j = 0; j <N; j++){
	    reduced[i][1]=reduced[i][1]+Fq[j][1]*sin(reduced[i][0]*Fq[j][0])*dq;
	}
	reduced[i][1]=reduced[i][1]*2/pi*exp(-Qdamp*reduced[i][0]*Qdamp*reduced[i][0]/2.0);
    }

    ofstream myfile (s + "." + "rpdf");
    if (myfile.is_open()){
	for(int count = 0; count < cut; count ++){
	    myfile << reduced[count][0] << "   " << reduced[count][1] << '\n';
        }
	myfile.close();
    } 
/*    ofstream myfile (s + "." + "sq");
    if (myfile.is_open()){
        for(int count = 0; count < N; count ++){
            myfile << Fq[count][0] << "   " << Fq[count][1] << '\n';
        }
        myfile.close();
    } */
    deallocatetwoD(geom, atom_num);
    deallocatetwoD(reduced, cut);
    deallocatetwoD(spec, 3000);
    deallocatetwoD(Fq,N);
    deallocatetwoD(sf,N);
}

