#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <sstream>
#include <string>
#include <unistd.h>
#include <array>
#include <vector>
#include <random>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <errno.h>

#include "src/Eigen/Core" 
#include "src/Eigen/Dense"
#include "src/functions.hh"
#include "src/update_parameters.hh"
#include "src/kalman.hh"

using namespace std;
using namespace Eigen;

int Nlin; //Number of (super)lineages 
int Ndeme; //Number of demes 
int T; //Number of timepoints
int mcmc_max=10000;//Total MCMC steps
double frac_burnin=0.0;//Fraction of burn in 
int numprint=10000;//Number of MCMC steps at which the state is recorded

int noise_mode=0;//noise ansatz
int C_mode=0;//0: infer deviation from uniform sampling

double Neff_ini=1000;//Initial population size
double C_ini = 1;//Initial value of deviation from uniform sampling

string Q_DB="nonDB"; 

std::random_device rd;
std::mt19937 gen(rd());
normal_distribution<double> norm_dist(0,1);
uniform_real_distribution<double> uniform(0.0, 1.0);

int main(int argc, char * argv[])
{
    //clock_t c_init = clock(); // Initial time; used to output run time
    using namespace std;
    time_t start_time, end_time;
    start_time = time(NULL);

    string infilename;
    string outfilename;
    //string IO_dir="";
    

    int c = 0;

    srand(time(NULL));//initialize the seed for rand function

    while ((c = getopt (argc, argv, "m:f:b:n:N:D:p:C:")) != -1)
    {
        if (c == 'm')
            mcmc_max = atoi(optarg);
        else if(c=='f')
            infilename = optarg;//filename of the input
        // else if(c=='g')
        //     outfilename = optarg;//filename of the input
        else if(c=='b')
            frac_burnin = atof(optarg);
        else if(c=='n')
            noise_mode = atoi(optarg);
        else if(c=='N')
            Neff_ini = atof(optarg);
        else if(c=='D')
            Q_DB= optarg;
        else if(c=='p')
            numprint= atoi(optarg);
        else if(c=='C')
            C_mode= atoi(optarg);
    }
    outfilename = infilename;

    mkpath("output/");
    //files for output
    //cout<<"OUTPUT files @"<<"output/"+IO_dir<<endl;
    cout<<"Filename: "<< infilename<<endl;
    cout<<C_mode<<": IF 0, C is inferred:"<<endl;

    //open output files
    string filename;
    ofstream f_A;
    filename="output/A_"+outfilename+".csv";
    f_A.open(filename);

    ofstream f_Ne;
    filename="output/Ne_"+outfilename+".csv";
    f_Ne.open(filename);

    ofstream f_C;
    filename="output/C_"+outfilename+".csv";
    f_C.open(filename);

    ofstream f_logLH;
    filename="output/logLH_"+outfilename+".csv";
    f_logLH.open(filename);

    ofstream f_log;
    filename="output/logfile_"+outfilename+".csv";
    f_log.open(filename);

    //read input files
    vector<vector<double> >  shape;
    filename="input/shape_"+infilename+".csv";
    Getmatrix(filename, shape);
    T=shape[0][0];
    Nlin=shape[1][0];
    Ndeme=shape[2][0];
    
    cout<<"T="<<T<<", Nlin="<<Nlin<<", Ndeme="<<Ndeme<<endl<<"noisemode="<<noise_mode<<endl<<Q_DB<<endl;
    f_log<<"T="<<T<<", Nlin="<<Nlin<<", Ndeme="<<Ndeme<<endl<<"noisemode="<<noise_mode<<endl<<Q_DB<<endl;

    filename="input/counts_"+infilename+".csv";
    vector<vector<double> > countsaux;
    Getmatrix(filename, countsaux);
    vector<vector<vector<double> > > counts(Nlin, vector<vector<double> >(T, vector<double>(Ndeme)));
    for (int i=0; i<countsaux.size(); i++) {
        int lin_label;
        int t;
        lin_label = (int)(i/T);
        t = i%T;
        for(int d=0;d<Ndeme;d++){
            counts[lin_label][t][d] =countsaux[i][d];
        }
    }

    filename="input/totcounts_"+infilename+".csv";
    vector<vector<double> > totcounts_aux;
    Getmatrix(filename, totcounts_aux);
    // cout<<totcounts_aux.size()<<endl;
    // cout<<totcounts_aux[0].size()<<endl;

    vector<vector<double> > totcounts( Ndeme , vector<double> (T));
    for (int t=0; t<T;t++) {
        for(int i=0;i<Ndeme;i++){
            totcounts[i][t] = totcounts_aux[t][i]+1;//+1 is added. In computing B below, pseudo counts are used. To make B<1, +1 here is necessary.
        }
    }


    vector<vector<vector<double> > > B(Nlin, vector<vector<double> >(T, vector<double>(Ndeme)));
    for(int t=0;t<T;t++){
        for (int i = 0; i < Ndeme; i++){
            for (int  l = 0; l < Nlin; l++){
                B[l][t][i] = (counts[l][t][i]+1)/totcounts[i][t];
    }}}  


    cout<<"Counts"<<endl;
    for(int i=0;i<Ndeme;i++){
        cout<<"@"<<i<<" :";
        for (int j=0; j<T; j++) {
            cout<<totcounts[i][j]-1<<" ";
        }
        cout<<endl;
    }
    cout<<endl;
        
        
    MatrixXd Imat = MatrixXd::Zero(Ndeme,Ndeme);
    for (int i=0; i<Ndeme; i++) {
            Imat(i,i)=1.0;
    }

/* Initialize  the state of a random walker*/
    MatrixXd A_ini(Ndeme,Ndeme);
    calc_A_ini(Ndeme, T, Nlin, B, A_ini);// A symmetric A close to the least square estimation.(A_ini satisfies DB).
    
    MatrixXd A_old(Ndeme,Ndeme);

    A_old =A_ini;
    vector<double> Ne_old (Ndeme);
    for (int i=0; i<Ndeme; i++) {
        Ne_old[i]=Neff_ini;
    }

    vector<double> C_old (Ndeme);
    for (int i=0; i<Ndeme; i++) {
        C_old[i]=C_ini;
    }

    vector<vector<double> > totcounts_eff( Ndeme , vector<double> (T));
    for (int t=0; t<T;t++) {
        for(int i=0;i<Ndeme;i++){
            totcounts_eff[i][t] = totcounts[i][t]/C_old[i];
        }
    }

    
    vector<double> Pi_old (Ndeme);
    calc_Pi(Ndeme, A_old,  Pi_old,"print_n");

    double logLH_old;

    cout<<"Aold"<<endl<<A_old<<endl;


    logLH_old = calc_LH(Nlin, Ndeme, T, B, Ne_old, totcounts_eff,A_old, noise_mode, Imat);
    cout<<"logLH_old"<< std::fixed << std::setprecision(8)<<logLH_old<<endl;
    int accept=0;

    //Declare variables used for proposed state
    MatrixXd A_new(Ndeme,Ndeme);
    vector<double> Pi_new (Ndeme);
    vector<double> Ne_new(Ndeme);//Controls noise in transition
    vector<double> C_new(Ndeme);

/* Run MCMC */
    int dstep=(int)(mcmc_max*(1-frac_burnin)/(double)numprint);// Steps between write output to files 
    if (dstep<1) { dstep=1;}
    int dstep_checkpoint=(int)((double)mcmc_max/20.0);// Steps between cout
    if(dstep_checkpoint<1){dstep_checkpoint=1;}


    for (int mcmc_step=0; mcmc_step<mcmc_max; mcmc_step++) {

        /* Cout "the state" */
        if(mcmc_step%dstep_checkpoint==0 &&mcmc_step>0){
            calc_Pi(Ndeme, A_old,  Pi_old,"print_n");
            cout<<"Step = "<<mcmc_step;
            cout<<", "<<" acc = "<<round((double)100.0*accept/dstep_checkpoint)<<"%, ";
            cout<<"lnLH="<<std::fixed << std::setprecision(6)<<logLH_old<<", ";
            check_DB(Ndeme, Pi_old,A_old);

            f_log<<mcmc_step<<", "<<" acc (%)= "<<round((double)100.0*accept/dstep_checkpoint)<<", ";
            f_log<<"lnLH="<<logLH_old<<endl;

            accept=0;//reset the counter
        }

        double p_proposal;
        double logLH_new;
        double p_acc;

        /* Imposed DB or not */
        if(Q_DB=="DB"){
            // To improve numerical stability, make sure the normalization of A, and compute the stationary distribution, and then impose DB on  A.  
            // Note: The ordering of the three operations is important, because the stationary distribution Pi exists only for the normalized A and because the DB can be imposed on A straightforwardly with Pi being fixed.
            normalize_A(Ndeme,A_old);
            calc_Pi(Ndeme, A_old,  Pi_old,"print_n");
            recover_DB(Ndeme, Pi_old,A_old);

            double r1; 
            r1 = uniform(gen);
            if(r1<0.5){//With probability 0.5, perform the reversible update that maintain the starionary distribution 
                update_Ne(Ndeme, Ne_old,Ne_new);
                if(C_mode==0){
                    update_C(Ndeme, C_old,C_new);
                }
                else{
                    C_new = C_old;
                }
                p_proposal=DB_reversible_update(Ndeme,Pi_old,Pi_new, A_old,A_new);
            }
            else{//Else, rescale a row that affects the starionary distribution 
                update_Ne(Ndeme, Ne_old,Ne_new);
                if(C_mode==0){
                    update_C(Ndeme, C_old,C_new);
                }
                else{
                    C_new = C_old;
                }
                p_proposal=DB_row_update(Ndeme,Pi_old,Pi_new, A_old,A_new);
            }
        }

        else if(Q_DB=="nonDB"){
            // To improve numerical stability, make sure the normalization of A 
            normalize_A(Ndeme,A_old);
            update_Ne(Ndeme, Ne_old,Ne_new);

            if(C_mode==0){
                update_C(Ndeme, C_old,C_new);
            }
            else{
                C_new = C_old;
            }


            p_proposal=nonrev_update(Ndeme, A_old,A_new);   
        }
        else{
            cout<<"Specify DB or nonDB"<<endl;
            return 0;
        }

        /*Compute the acceptance probability*/
        for (int t=0; t<T;t++) {
            for(int i=0;i<Ndeme;i++){
                totcounts_eff[i][t] = totcounts[i][t]/C_new[i];
            }
        }

        logLH_new = calc_LH(Nlin, Ndeme, T, B, Ne_new, totcounts_eff,A_new, noise_mode, Imat);
        p_acc =  p_proposal*exp(logLH_new-logLH_old);  

        /*Accept the proposal with Pr = p_acc*/
        double r;
        r = uniform(gen);
        if(r<p_acc){
            logLH_old=logLH_new;
            A_old =A_new;
            Ne_old=Ne_new;
            C_old=C_new;
            Pi_old = Pi_new;
            accept+=1;
        }

         /* Write "the state" to files*/
        if(mcmc_step%dstep==0 && mcmc_step>=frac_burnin*mcmc_max) {
            //output A
            for (int i=0; i<Ndeme; i++) {
                for (int j=0; j<Ndeme;j++) {
                    if (i==Ndeme-1 and j==Ndeme-1) {
                        f_A<<A_old(i,j)<<endl;
                    }else{
                        f_A<<A_old(i,j)<<",";
                    }
                }
            }
            //output Ne
            for (int i=0; i<Ndeme; i++) {
                if(i<Ndeme-1){
                    f_Ne<<Ne_old[i]<<",";
                }else{
                    f_Ne<<Ne_old[i]<<endl;
                }
            }
            //output C
            for (int i=0; i<Ndeme; i++) {
                if(i<Ndeme-1){
                    f_C<<C_old[i]<<",";
                }else{
                    f_C<<C_old[i]<<endl;
                }
            }
            //output logLH
            f_logLH<<mcmc_step <<","<<logLH_old<<endl;
        }

    }//End of MCMC


    end_time = time(NULL);

    cout<<"Afinal"<<endl<<A_old<<endl;
    f_log << "Run time: " << double(end_time-start_time)<< endl;

    return 0;

};
