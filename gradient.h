/*
 class grad
 {
 public:
 grad(
 const Eigen::MatrixXcd& H0, const int& N1, const int& N2, const int& N3
 )
 {
 h0 = H0;
 n1 = N1;
 n2 = N2;
 n3 = N3;
 }
 
 column_vector operator() (const column_vector& spinconfiguration) const
 {
 Eigen::MatrixXcd h = constructspins(h0, flattospins(spinconfiguration,n1,n2,n3),n1,n2,n3);
 Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(h);
 Eigen::VectorXcd eigvals = es.eigenvalues();
 Eigen::MatrixXcd eigvects = es.eigenvectors();
 
 map<int, Eigen::Vector3d> Aspins = flattospins(spinconfiguration,n1,n2,n3);
 map<int, Eigen::Vector3d> gradient;
 Eigen::Vector3d partial;
 double J = 0.00035; double Jk = 0.007;
 map<int, int> nearAsites; map<int, int> nearBsites;
 map<int, int> mu;
 complex<double> chem = chempot(eigvals);
 Eigen::MatrixXcd eigvectsadjoint = eigvects.adjoint();
 int c=0;
 for (int m = 0;m<eigvals.size();m++){if (real(eigvals[m])<real(chem)){mu.insert(make_pair(c,m));} c++;}
 for (int n = 0;n<4*n1*n2*n3;n++){
 partial<<0,0,0;
 nearAsites = nearto(n,"A","A",n1,n2,n3);
 nearBsites = nearto(n,"A","B",n1,n2,n3);
 
 for (int o = 0; o<nearBsites.size();o++){
 for (int p = 0; p<mu.size();p++){
 partial = partial + Jk*(eigvects(nearBsites[o],mu[p])*tau(nearBsites[o],mu[p],n1,n2,n3)*eigvectsadjoint(nearBsites[o],mu[p])).real();
 }
 }
 
 for (int q = 0; q<nearAsites.size();q++){
 partial = partial + J*Aspins[nearAsites[q]];
 }
 gradient.insert(make_pair(n,partial));
 }
 
 return spinstoflat(gradient,n1,n2,n3);
 }
 
 private:
 Eigen::MatrixXcd h0; int n1; int n2; int n3;
 };
 *//*
    class bgrad
    {
    public:
    bgrad(
    const MatrixXcd& H0, const int& N1, const int& N2, const int& N3
    )
    {
    h0 = H0;
    n1 = N1;
    n2 = N2;
    n3 = N3;
    }
    
    column_vector operator() (const column_vector& spinconfiguration) const
    {
    MatrixXcd h = constructspins(h0, flattospins(spinconfiguration,n1,n2,n3),n1,n2,n3);
    ComplexEigenSolver<MatrixXcd> es(h);
    VectorXcd eigvals = es.eigenvalues();
    MatrixXcd eigvects = es.eigenvectors();
    
    map<int, Vector3d> Bspins = flattospins(spinconfiguration,n1,n2,n3);
    map<int,Vector3d> gradient;
    Vector3d partial;
    
    double U = 0.0;
    
    map<int, int> mu;
    complex<double> chem = chempot(eigvals);
    MatrixXcd eigvectsadjoint = eigvects.adjoint();
    int c=0;
    for (int m = 0;m<eigvals.size();m++){if (real(eigvals[m])<real(chem)){mu.insert(make_pair(c,m));} c++;}
    
    for (int i=0; i<4*n1*n2*n3; i++){
    partial = U*Bspins[i];
    for (int p = 0; p<mu.size();p++){
    partial = partial - U*(eigvects(i,mu[p])*tau(i,mu[p],n1,n2,n3)*eigvectsadjoint(i,mu[p])).real();
    }
    gradient.insert(make_pair(i,partial));
    }
    
    return spinstoflat(gradient,n1,n2,n3);
    }
    
    private:
    MatrixXcd h0; int n1; int n2; int n3;
    };*/
/*
 class ngrad
 {
 public:
 ngrad(
 const MatrixXcd& H0, const int& N1, const int& N2, const int& N3
 )
 {
 h0 = H0;
 n1 = N1;
 n2 = N2;
 n3 = N3;
 }
 
 column_vector operator() (const column_vector& spinconfiguration) const
 {
 double delta = 0.1;
 MatrixXcd h = constructspins(h0, flattospins(spinconfiguration,n1,n2,n3),n1,n2,n3);
 MatrixXcd hf; MatrixXcd hb; MatrixXcd hfd; MatrixXcd hbd;
 ComplexEigenSolver<MatrixXcd> es(h);
 VectorXcd eigvals = es.eigenvalues();
 VectorXcd eigf;VectorXcd eigb;
 double energyf; double energyb;
 
 MatrixXcd eigvects = es.eigenvectors();
 column_vector gradient;
 column_vector stepforward = spinconfiguration;
 column_vector stepback = spinconfiguration;
 MatrixXcd eigvectsadjoint = eigvects.adjoint();
 for (int n = 0;n<spinconfiguration.size();n++){
 stepforward = spinconfiguration;
 stepback = spinconfiguration;
 stepforward(n) = stepforward(n) + delta;
 stepback(n) = stepback(n) - delta;
 hf = constructspins(h0, flattospins(stepforward,n1,n2,n3),n1,n2,n3);
 hb = constructspins(h0, flattospins(stepback,n1,n2,n3),n1,n2,n3);
 hfd = eigvects*hf*eigvectsadjoint;
 hbd = eigvects*hb*eigvectsadjoint;
 for (int i = 0;i<2*4*n1*n2*n3;i++){
 eigf[i] = hfd(i,i);
 eigb[i] = hbd(i,i);
 }
 energyf = real(groundenergy(eigf));
 energyb = real(groundenergy(eigb));
 gradient(n) = (energyf - energyb)/(2.0*delta);
 }
 
 
 return gradient;
 }
 
 private:
 MatrixXcd h0; int n1; int n2; int n3;
 };*/