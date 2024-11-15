using namespace std;
using namespace Eigen;

template <typename MatrixType>
inline typename MatrixType::Scalar logdet(const MatrixType& M, bool use_cholesky = false) {
  using namespace Eigen;
  using std::log;
  typedef typename MatrixType::Scalar Scalar;
  Scalar ld = 0;
  if (use_cholesky) {
    LLT<Matrix<Scalar,Dynamic,Dynamic> > chol(M);
    auto& U = chol.matrixL();
    for (unsigned i = 0; i < M.rows(); ++i)
      ld += log(U(i,i));
    ld *= 2;
  } else {
    PartialPivLU<Matrix<Scalar,Dynamic,Dynamic> > lu(M);
    auto& LU = lu.matrixLU();
    Scalar c = lu.permutationP().determinant(); // -1 or 1
    for (unsigned i = 0; i < LU.rows(); ++i) {
      const auto& lii = LU(i,i);
      if (lii < Scalar(0)) c *= -1;
      ld += log(abs(lii));
    }
    ld += log(c);
  }
  return ld;
};


double logdet_2dvec(Eigen::MatrixXd& mat){
    return logdet(mat);   
};


double logpdf_gauss(int n, vector<double> & x,vector<double> & mu, Eigen::MatrixXd& cov){

    double res=0;
    MatrixXd covinv(n, n);
    matinv(cov, covinv);
    
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            res+=-0.5*(x[i]-mu[i])*covinv(i,j)*(x[j]-mu[j]);
        }
    }

    res +=-0.5*n*log(2*3.14159265359) - 0.5*logdet_2dvec(cov);
    
    return res; 
};


void calc_mu_V_logLH_0(int n, vector<double> & x,Eigen::MatrixXd& Imat,  vector<double> & totcounts_t, vector<double> & Ne,  vector<double> & mu,Eigen::MatrixXd& V,double &logc_0, int noise_mode, double freq_scale, vector<double> & freq_scale_vec){

    MatrixXd Sigma=MatrixXd::Zero(n,n);//Noise in emission
   // MatrixXd Gamma=MatrixXd::Zero(n,n);//Noise in transition

    if(noise_mode==0){
        for (int i=0; i<n; i++){
            Sigma(i,i)=x[i]*(1-x[i])/totcounts_t[i];
        }
    }
    else if (noise_mode==1) {
        for (int i=0; i<n; i++){
            Sigma(i,i)=freq_scale*(1-freq_scale)/totcounts_t[i];
        }
    }
    else if (noise_mode==2) {
        for (int i=0; i<n; i++){
            Sigma(i,i)=freq_scale_vec[i]*(1-freq_scale_vec[i])/totcounts_t[i];
        }
    }

    //calculate K_0 from V_pre
    MatrixXd K=MatrixXd::Zero(n,n);
    MatrixXd PsumSigma=MatrixXd::Zero(n,n);
    for (int i=0; i<n; i++)for (int j=0; j<n; j++)PsumSigma(i,j)=V(i,j)+ Sigma(i,j);

    MatrixXd invPsumSigma(n,n);
    matinv(PsumSigma, invPsumSigma);
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            K(i,j) = 0;
            for (int k=0; k<n; k++) {
                K(i,j) +=V(i,k)*invPsumSigma(k,j);
            };
    }}
    
   //calculate mu_0
    vector<double> mu_next(n,0);
    for (int i=0; i<n; i++) {
        mu_next[i]+=mu[i];
        for (int k=0; k<n; k++)mu_next[i]+=K(i,k)*(x[k] -mu[k]);// This line can be commented out because mu is assumed to be x at t=0
    }

    //calculate V_0
    MatrixXd V_next=MatrixXd::Zero(n,n);
    for (int i=0; i<n; i++) {
        for (int j=0;j<n; j++) {
            for (int k=0; k<n; k++) {
                V_next(i,j)+=(Imat(i,k) - K(i,k))*V(k,j);
            }
        }
    }
    
    //calculate c_0
    MatrixXd cov(n,n);
    for (int i=0; i<n; i++) {
        for (int j=0;j<n; j++) {
            cov(i,j) = V(i,j) + Sigma(i,j);
        };
    };
    
    logc_0=logpdf_gauss(n, x, mu, cov);
    
   //update mu&V
    mu=mu_next;
    V=V_next;
};


void calc_mu_V_logLH_t(int n,vector<double> & x,vector<double> & x_pre, Eigen::MatrixXd& Imat,Eigen::MatrixXd& A, vector<double> & totcounts_t, vector<double> & Ne,  vector<double> & mu,Eigen::MatrixXd& V,double &logc_t,int noise_mode, double freq_scale,  vector<double> & freq_scale_vec){
  
    //Function that calculates mu_t, V_t, c_t from mu_t-1, V_t-1, c_t-1
    vector<double> Amu(n,0);
    for (int i=0; i<n; i++)for (int k=0; k<n; k++){Amu[i]+=A(i,k)*mu[k];}
    
    // vector<double> Ax(n,0);
    // for (int i=0; i<n; i++)for (int k=0; k<n; k++){Ax[i]+=A(i,k)*x[k];}
  
    MatrixXd Sigma=MatrixXd::Zero(n,n);//Noise in emission
    MatrixXd Gamma=MatrixXd::Zero(n,n);//Noise in transition

    if(noise_mode==0){
            for (int i=0; i<n; i++){
                Sigma(i,i)=x[i]*(1-x[i])/totcounts_t[i];
                Gamma(i,i)=x[i]*(1-x[i])/Ne[i];
                //Gamma(i,i)=x_pre[i]*(1-x_pre[i])/Ne[i];
            }
    }
    else if(noise_mode==1){
        for (int i=0; i<n; i++){
                Sigma(i,i)=freq_scale*(1-freq_scale)/totcounts_t[i];
        }
        for (int i=0; i<n; i++){
                Gamma(i,i)=freq_scale*(1-freq_scale)/Ne[i];
        }
    }
    else if(noise_mode==2){
        for (int i=0; i<n; i++){
                Sigma(i,i)=freq_scale_vec[i]*(1-freq_scale_vec[i])/totcounts_t[i];
        }
        for (int i=0; i<n; i++){
                Gamma(i,i)=freq_scale_vec[i]*(1-freq_scale_vec[i])/Ne[i];
        }
    }

    //calculate P_t-1
    MatrixXd P=MatrixXd::Zero(n,n);
    for (int i=0; i<n; i++)for (int j=0; j<n; j++)P(i,j)+=Gamma(i,j);
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            for (int k=0; k<n; k++) {
                for (int l=0; l<n; l++) {
                    P(i,j)+=A(i,k)*V(k,l)*A(j,l);
                }
            }
        }
    }
    
    //calculate K_t from P_t-1
    MatrixXd K=MatrixXd::Zero(n,n);
    MatrixXd PsumSigma=MatrixXd::Zero(n,n);

    for (int i=0; i<n; i++)for (int j=0; j<n; j++)PsumSigma(i,j)=P(i,j)+ Sigma(i,j);

    MatrixXd invPsumSigma=MatrixXd::Zero(n,n);
    matinv(PsumSigma, invPsumSigma);
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            K(i,j) = 0;
            for (int k=0; k<n; k++) {
                K(i,j) +=P(i,k)*invPsumSigma(k,j);
            };
    }}
    
    //calculate mu_t
    vector<double> mu_next(n,0);
    for (int i=0; i<n; i++) {
        mu_next[i]=Amu[i];
        for (int k=0; k<n; k++)mu_next[i]+=K(i,k)*(x[k] -Amu[k]);
    }

    //calculate V_t
    MatrixXd V_next=MatrixXd::Zero(n,n);
    for (int i=0; i<n; i++) {
        for (int j=0;j<n; j++) {
            for (int k=0; k<n; k++) {
                V_next(i,j)+=(Imat(i,k) - K(i,k))*P(k,j);
            }
        }
    }
    
    //calculate c_t
    MatrixXd cov(n,n);
    for (int i=0; i<n; i++) {
        for (int j=0;j<n; j++) {
            cov(i,j) = P(i,j) + Sigma(i,j);
        };
    };

    logc_t=logpdf_gauss(n,x,Amu, cov);

    //update mu&V
    mu=mu_next;
    V=V_next;
}


double calc_LH(int Nlin, int n, int T, vector<vector<vector<double> > >& B, vector<double>& Ne, vector<vector<double> > & totcounts, Eigen::MatrixXd& A, int noise_mode, Eigen::MatrixXd& Imat){

       size_t lin;
       double logLH=0;
       
//#pragma omp parallel for reduction(+:logLH)
       for (lin =0; lin<Nlin; lin++) {

           double freq_scale=0;// for noise_mode=1
           for(int d=0;d<n;d++){
                for(int t=0;t<T;t++){
                    freq_scale+=B[lin][t][d];
                }}
           freq_scale*=1./((double)n*(double)T);

           vector<double> freq_scale_vec(n,0);// for noise_mode=2
           for(int d=0;d<n;d++){
                for(int t=0;t<T;t++){
                    freq_scale_vec[d]+=B[lin][t][d];
                }
                freq_scale_vec[d]*=1./((double)T);
            }

           vector<double> x_0(n);
           x_0=B[lin][0];

           //Compute mu(0),V(0),c(0)
           vector<double> mu(n,0);
           for(int i=0;i<n;i++){
                mu[i] = x_0[i];
           }
           
           MatrixXd V=MatrixXd::Zero(n,n);

           for (int i=0; i<n;i++)V(i,i) = mu[i]*(1-mu[i])/totcounts[i][0];

           vector<double> x_t(n);
           vector<double> totcounts_t(n);

           double logc_t;
           
           //Compute mu(t),V(t),c(t) from mu(t-1),V(t-1),c(t-1)
           for(int t=0;t<T;t++){
               x_t=B[lin][t];
            
               for(int i=0;i<n;i++){totcounts_t[i]= totcounts[i][t];}

               if(t>0){
                   vector<double> x_pre(n);
                   x_pre=B[lin][t-1];
                   calc_mu_V_logLH_t(n,x_t,x_pre, Imat, A, totcounts_t,Ne, mu, V, logc_t, noise_mode,freq_scale,freq_scale_vec);}
               else{
                   calc_mu_V_logLH_0(n,x_t,Imat, totcounts_t, Ne, mu, V, logc_t,noise_mode,freq_scale,freq_scale_vec);
               }

               if(t>0){//add logLH only for t>0. This is because we assumed mu=x for t=0
               logLH+=logc_t;
               }

           }
       
       }

       return logLH;
}