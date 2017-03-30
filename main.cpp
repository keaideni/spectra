#include <iostream>
#include <fstream>
#include "Mat.h"
#include "mpi.h"

int OP::Max;
using namespace std;

int main(int argc, char* argv[])
{
        ifstream inpara("./data/QNosave.txt");
        if(!inpara)
        {
                cerr<<"the file QNosave.txt doesn't exit!"<<endl;
        }

        JC_Parameter para(inpara);
        OP::Max=para.ParticleNo();

        inpara.close();
        //para.show();

        /*MPI_Status status;

        int myid, numprocess;

        int groupn(48);
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);
        MPI_Comm_size(MPI_COMM_WORLD, &numprocess);

        int everygroup(48/numprocess);//cout<<everygroup;
        if(myid==0)
        {
                double gl(0.003), gr(0.027);
                Mat H0(para, gl, gr);


                ofstream spectra("./result/spectra");

                spectra<<gl<<"\t";
                for(int i=0; i<5; ++i)
                        spectra<<H0.Eigenvalues()(i)<<"\t";

                spectra<<endl;

                for(int i=0; i<everygroup; ++i)
                {
                        gl=0.003+0.0005*(i+1);
                        gr=0.03-gl;
                        Mat H(para, gl, gr);

                        spectra<<gl<<"\t";
                        for(int i=0; i<5; ++i)
                                spectra<<H.Eigenvalues()(i)<<"\t";

                        spectra<<endl;

                }
                for(int id=1; id<numprocess; ++id)
                {
                        VectorXd gl1(everygroup);
                        //VectorXd overlap1(everygroup);
                        MatrixXd spectra1(everygroup,5);//here the 5 is for the number of eigenvalues we obtained.
                        
                        MPI_Recv(&gl1(0), everygroup, MPI_DOUBLE, id, id,
                                MPI_COMM_WORLD, &status);
                        
                        //cout<<gl1<<endl;

                        
                        MPI_Recv(&spectra1(0,0), everygroup*5, MPI_DOUBLE, id,
                                id+numprocess, MPI_COMM_WORLD, &status);
                        //cout<<"from 0:"<<endl;
                        //cout<<overlap1<<endl;



                        
                        for(int i=0; i<everygroup;++i)
                        {
                                spectra<<gl1(i)<<"\t";
                                for(int j=0; j<5; ++j)
                                        spectra<<spectra1(i,j)<<"\t";
                                spectra<<endl;
                        }

                        //delete gl1;
                        //delete energy1;
                        
                }
                spectra.close();

        }
        else
        {
                
                
                VectorXd gl(everygroup);
                
                MatrixXd spectra(everygroup, 5);


                for(int i=everygroup*myid; i<everygroup*(myid+1);++i)
                {
                        //cout<<i-myid*everygroup<<endl;
                        
                        gl((i-myid*everygroup))=0.003+0.0005*(i+1);
                        double gr=0.03-gl[i-myid*everygroup];
                        Mat H(para, gl[i-myid*everygroup], gr);

                        spectra.row(i-myid*everygroup)=H.Eigenvalues().transpose();

                        
                        
                }
                //cout<<"from "<<myid<<endl;
                //cout<<overlap<<endl;
                
                

                MPI_Send(&gl(0), everygroup, MPI_DOUBLE, 0, 
                        myid, MPI_COMM_WORLD);
                MPI_Send(&spectra(0,0), everygroup*5, MPI_DOUBLE, 0, 
                        myid+numprocess, MPI_COMM_WORLD);
                //cout<<myid+numprocess<<endl;
                


                //delete gl;
                //delete energy;
                
                
        }

        


        
        MPI_Finalize();//this is so important. to prove the early exit.

        

        //com(para);*/

        cout.precision(10);

        Mat H1(para, 0.00, 0.03), H2(para, 0.03, 0.00);
        cout<<H1.Eigenvalues()<<endl<<endl<<H2.Eigenvalues()<<endl;

        
}