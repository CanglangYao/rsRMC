#include <iostream>
#include <fstream>
#include <math.h>  
#include <string>
#include <iomanip>
#include <mpi.h>
#include <Eigen/Dense>

using namespace std;
using Eigen::MatrixXd;
int atom_num = 384;
int block_num = 320;
double Qdamp=0.03;
double la=25.24;
double lb=25.24;
double lc=25.24;
double rou = atom_num/la/lb/lc;
double pi = 3.1415926;
double step=0.02;
int N = 280;
int type = 4;
int cut = 630;
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

void pdf_calc(double **g, double **&pdf){
    int ix, iy, iz, nd, w_id;
    double x, y, z, d;
    for (int ii=0; ii < atom_num-1; ii=ii+1){
        for (int jj =ii+1; jj < atom_num; jj=jj+1){
	    ix=(int) ((g[ii][1]-g[jj][1])/la*2);
            iy=(int) ((g[ii][2]-g[jj][2])/lb*2);
	    iz=(int) ((g[ii][3]-g[jj][3])/lc*2);
	    ix=ix+(int) (-0.5*ix);
	    iy=iy+(int) (-0.5*iy);
	    iz=iz+(int) (-0.5*iz);
            x=(g[ii][1]-g[jj][1])-ix*la;
            y=(g[ii][2]-g[jj][2])-iy*lb;
            z=(g[ii][3]-g[jj][3])-iz*lc;
            d=sqrt(x*x+y*y+z*z);
	    nd = (int) (d/step);
	    w_id=(int) g[ii][0]*type+(int) g[jj][0];
	    pdf[nd][w_id]=pdf[nd][w_id]+1;
	}
    }
}
void calc_chi2(double **&reduced, double **fq){
    double dq = fq[1][0]-fq[0][0];
    for(int i =0; i< cut; i++){
	reduced[i][0]=i*step;
	for(int j = 0; j <N; j++){
	    reduced[i][1]=reduced[i][1]+fq[j][1]*sin(reduced[i][0]*fq[j][0])*dq;
    	}
	reduced[i][1]=reduced[i][1]*2.0/pi*exp(-Qdamp*reduced[i][0]*Qdamp*reduced[i][0]/2.0);
    }
}
void sort(double *&allmoves, int num){
    double minchi2=allmoves[6];
    for (int i=1; i<num; i++){
        if (allmoves[i*7+6]<minchi2){
	    minchi2=allmoves[i*7+6];
	    allmoves[0]=allmoves[i*7+0];
	    allmoves[1]=allmoves[i*7+1];
	    allmoves[2]=allmoves[i*7+2];
	    allmoves[3]=allmoves[i*7+3];
	    allmoves[4]=allmoves[i*7+4];
	    allmoves[5]=allmoves[i*7+5];
	    allmoves[6]=allmoves[i*7+6];
	}
    }
}

int calc_dist(double **g, double coord_x, double coord_y, double coord_z, int id){
    int iix, iiy, iiz;
    double xx, yy, zz, d_2;
    int too_short=1;
    for (int i_i=0; i_i < atom_num; i_i++){
        iix=(int) ((g[i_i][1]-coord_x)/la*2);
        iiy=(int) ((g[i_i][2]-coord_y)/lb*2);
        iiz=(int) ((g[i_i][3]-coord_z)/lc*2);
	iix=iix+(int) (-0.5*iix);
	iiy=iiy+(int) (-0.5*iiy);
	iiz=iiz+(int) (-0.5*iiz);
        xx=(g[i_i][1]-coord_x)-iix*la;
        yy=(g[i_i][2]-coord_y)-iiy*lb;
        zz=(g[i_i][3]-coord_z)-iiz*lc;
        d_2=xx*xx+yy*yy+zz*zz; 
	if (d_2<7.3 and i_i!=id){
	    if(id>127){
	        too_short=0;
	        break;
	    }
	    else if(id <64){
	        if (i_i != id+64){
		    too_short=0;
		    break;
		}
            }
	    else{
	        if(i_i != id-64){
		    too_short=0;
		    break;
		}
	    }
	}
    }
    return too_short;
}

void rd_str_generate(double **&str, double *sub_list, double *&collec, int id, int rd_id){
    if (id<64){
	double m_center[3];
        double theta_r = sub_list[0]*2*pi;
	double theta_p = sub_list[1]*2*pi;
	double theta_y = sub_list[2]*2*pi;
	double translation = 2.0;
	double vibration = 0.4;
	double C_array[3]={-0.17858, -0.724015, 0.002285};
	double N_array[3]={0.17858, 0.724015, -0.002285};
	double rol[3][3]={{1,0,0}, {0, cos(theta_r), -sin(theta_r)}, {0, sin(theta_r), cos(theta_r)}};
	double pit[3][3]={{cos(theta_p), 0, sin(theta_p)}, {0,1,0}, {-sin(theta_p), 0, cos(theta_p)}};
	double yaw[3][3]={{cos(theta_y), -sin(theta_y), 0}, {sin(theta_y), cos(theta_y), 0}, {0,0,1}};
	MatrixXd C(3,1), N(3,1), C_orig(3,1), N_orig(3,1);
        MatrixXd mtx_rol(3,3), mtx_pit(3,3), mtx_yaw(3,3), mtx_rot(3,3);
	for (int i = 0; i<3; i++){
	    m_center[i]=(str[id][i+1]+str[id+64][i+1])*0.5+translation*(sub_list[i+3]-0.5);
	    C_orig(i,0)=C_array[i];
	    N_orig(i,0)=N_array[i];
	    for (int j = 0; j<3; j++){
	    	mtx_rol(i,j)=rol[i][j];
		mtx_pit(i,j)=pit[i][j];
		mtx_yaw(i,j)=yaw[i][j];
	    }
	}
        mtx_rot=mtx_rol*mtx_pit*mtx_yaw;
	C = mtx_rot*C_orig;
	N = mtx_rot*N_orig;
	for (int k=0; k<3; k++){
	    str[id][k+1]=C(k,0)+m_center[k]+(sub_list[k+6]-0.5)*vibration;
	    str[id+64][k+1]=N(k,0)+m_center[k]+(sub_list[k+9]-0.5)*vibration;
	    collec[rd_id*7+k]=str[id][k+1];
	    collec[rd_id*7+3+k]=str[id+64][k+1];
	}
//	cout << collec[rd_id*7] << "  "<< collec[rd_id*7+1] <<"  " << collec[rd_id*7+2] << endl;
    }
    else if (id > 127){
	for (int i=0; i<3; i++){
	    str[id][i+1]=str[id][i+1]+(sub_list[i]-0.5)*3.0;
            collec[rd_id*7+i]=str[id][i+1];
	}
    }
    else {
        cout << "you have inputed a N atom, which however does not need to be computed twice." << endl;
    }
}


int main(int argc, char *argv[]) {
          std::string s(argv[1]);
	  std::ifstream file( argv[1] );
	  double frames = stod(argv[2]);
	  int rank, size;
	  double **geom = nullptr;
	  double **backup = nullptr;
	  allocatetwoD(geom, atom_num, 4);
	  allocatetwoD(backup, atom_num, 4);
	  MPI_Init(&argc,&argv);
	  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	  MPI_Comm_size(MPI_COMM_WORLD,&size);
    for (int r = 0; r < atom_num; r++) //Outer loop for rows
    {
        for (int c = 0; c < 4; c++) //inner loop for columns
        {
          file >> geom[r][c];  //Take input from file and put into myArray
	  backup[r][c]=geom[r][c];
        }
    }
    double **exp_pdf;
    allocatetwoD(exp_pdf, cut, 2);
    std::ifstream exp_p;
    exp_p.open("./exp-pdf");
    for (int r = 0; r < cut; r++){
        for(int c = 0; c < 2; c++){
            exp_p >> exp_pdf[r][c];
        }
    }
    double **prev_pdf;
    allocatetwoD(prev_pdf, cut, 2);
    std::ifstream prev;
    prev.open("./previous_pdf");
    for (int r = 0; r < cut; r++){
        for(int c = 0; c < 2; c++){
	    prev >> prev_pdf[r][c];
	}
    }
    double factor;
    std::ifstream scale;
    scale.open("./factor");
    scale >> factor;

    int *order = new int[block_num];
    float od_temp;
    std::ifstream at_od;
    at_od.open("./order");
    for (int r =0; r<block_num; r++){
	at_od >> od_temp;
	order[r]=(int)od_temp-1;
    }
    double table[4][9]={ 
	            {2.31000, 20.8439, 1.02000, 10.2075, 1.58860, 0.568700, 0.865000, 51.6512, 0.215600},
	            {12.2126, 0.005700, 3.13220, 9.89330, 2.01250, 28.9975, 1.16630, 0.582600, -11.529},
	            {20.1472, 4.34700, 18.9949, 0.381400, 7.51380, 27.7660, 2.27350, 66.8776, 4.07120},
                    {31.0617, 0.690200, 13.0637, 2.35760, 18.4420, 8.61800, 5.96960, 47.2579, 13.4118}
                    };
    
    int rd_num = 20;
    double chi2_last_run=10000000.0;
    int idx;
    srand(time(NULL)+rank);
    double ave_sf, delta;
    for(int aod = 0; aod < block_num; aod++){
	idx=order[aod];
        double *selec=new double[rd_num*7];
	double *best=new double[size*rd_num*7];
        for (int nrd = 0; nrd< rd_num; nrd++){
	    double *rd_list=new double[size*12];
            double **Fq=nullptr;
	    double **sf=nullptr;
	    double **rpdf = nullptr;
	    double **spec = nullptr;
	    allocatetwoD(Fq, N, 2);
	    allocatetwoD(sf, N, type);
	    allocatetwoD(rpdf, 3000, 2);
	    allocatetwoD(spec, 3000, type*type);
	    double *atom=new double[12];
	    double chi2=0.0;
/*	    if (rank==0){
	        for (int i=0; i< size*12; i++){
		    rd_list[i]=(rand()%1000)/((double)1000);
		}
	    }
	    MPI_Barrier(MPI_COMM_WORLD);
	    MPI_Scatter(rd_list, 12, MPI_DOUBLE, atom, 12, MPI_DOUBLE, 0, MPI_COMM_WORLD); */
            int shortest_1, shortest_2;
	    for (int rd_t=0; rd_t<100; rd_t++){
		for (int i=0; i<12; i++){
		    atom[i]=(rand()%600)/((double)600);
		}
	        rd_str_generate(geom, atom, selec, idx, nrd);
	        shortest_1=calc_dist(geom, selec[nrd*7+0], selec[nrd*7+1], selec[nrd*7+2], idx);
    		if (shortest_1 >0){
		    if(idx>127){
		        break;
		    } 
		    if(idx<64){
		        shortest_2=calc_dist(geom, selec[nrd*7+3], selec[nrd*7+4], selec[nrd*7+5], idx+64);
			if (shortest_2 >0){
			    break;
			}
		    }
                }
		else{
		    geom[idx][1]=backup[idx][1];
		    geom[idx][2]=backup[idx][2];
		    geom[idx][3]=backup[idx][3];
		    if (idx<64){
		        geom[idx+64][1]=backup[idx+64][1];
			geom[idx+64][2]=backup[idx+64][2];
			geom[idx+64][3]=backup[idx+64][3];
		    }
		}
	    }
	    pdf_calc(geom, spec);
            for (int nq = 0; nq < N; nq++){
 	        double **Gr = nullptr;
		double *weight = new double[type*type];
	        allocatetwoD(Gr, 3000, 2);
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
	        ave_sf=(sf[nq][0]+sf[nq][1]+sf[nq][2]*3+sf[nq][3])/6;  
		for (int i=0; i<type; i++){
		    for (int j=0; j < type; j++){
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
            calc_chi2(rpdf, Fq);
	    for(int aid=0; aid <3; aid++){
	        geom[idx][aid+1]=backup[idx][aid+1];
		if (idx <64) {
		    geom[idx+64][aid+1]=backup[idx+64][aid+1];
		}
	    }
	    for (int i=0; i<cut; i++){
	        delta=(rpdf[i][1]+prev_pdf[i][1])*factor/frames-exp_pdf[i][1];
	        chi2=chi2+delta*delta;
	    }
	    selec[nrd*7+6]=chi2;
            deallocatetwoD(rpdf, 3000);	
	    deallocatetwoD(Fq,N);
	    deallocatetwoD(sf,N);
	    deallocatetwoD(spec, 3000);
	    delete[] rd_list;
	    delete[] atom;
        }
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Gather(selec, rd_num*7, MPI_DOUBLE, best, rd_num*7, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
	if (rank==0){
            sort(best, rd_num*size);
	    cout << aod << "  " << best[6] << endl;
	}
	MPI_Bcast(best, 7, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	if (best[6]<chi2_last_run){
	    chi2_last_run=best[6];
	    geom[idx][1]=best[0];
	    geom[idx][2]=best[1];
	    geom[idx][3]=best[2];
	    if (idx <64){
	        geom[idx+64][1]=best[3];
	        geom[idx+64][2]=best[4];
	        geom[idx+64][3]=best[5];
	    }
	}
	MPI_Barrier(MPI_COMM_WORLD);
        delete[] selec;
	delete[] best;
    }
    if (rank==0){
//	cout << atom_num<< endl;
        ofstream myfile (s + "-new");
        if (myfile.is_open()){
            for(int r = 0; r < atom_num; r ++){
	        for(int c=0; c<4; c++){
                    myfile << geom[r][c] << "   ";
    	        }
	        myfile << endl;
	    }
            myfile.close();
        } 
    }
    MPI_Barrier(MPI_COMM_WORLD);
    deallocatetwoD(exp_pdf, cut);
    deallocatetwoD(geom, atom_num);
    deallocatetwoD(backup, atom_num);
    delete[] order;
    delete[] prev_pdf;
}

