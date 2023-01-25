Eigen::VectorXcd hk8(Eigen::MatrixXcd h, Eigen::Vector3d k, int N1, int N2, int N3){
    Eigen::VectorXcd eigvals;
    Eigen::MatrixXcd hk(8,8);
    hk<<
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0;

    double normalization = pow(2*4*N1*N2*N3,-0.5);
    Eigen::Vector3d rij; Eigen::Vector3d Rij; rij<<0,0,0; Rij<<0,0,0;
    Eigen::Vector3d a1; Eigen::Vector3d a2; Eigen::Vector3d a3;double A = 1.0;
    a1<<0,0.5*A,0.5*A; a2<<0.5*A,0,0.5*A; a3<<0.5*A,0.5*A,0;
    double a, b, c;

    for (int i = 0;i<2*4*N1*N2*N3;i++){
        for (int j = 0;j<2*4*N1*N2*N3;j++){
            if (h(i,j) != 0.0){
                a=0.0;b=0.0;c=0.0;
                int mui = iton(i,N1,N2,N3)[3];
                int muj = iton(j,N1,N2,N3)[3];
                rij = ntor(iton(j%(4*N1*N2*N3),N1,N2,N3))-ntor(iton(i%(4*N1*N2*N3),N1,N2,N3));
                if (pow((rij+N1*a1).dot((rij+N1*a1)),0.5) < 0.5){a=1.0;}
                if (pow((rij+N1*a1+N2*a2).dot((rij+N1*a1+N2*a2)),0.5) < 0.5){a=1.0;b=1.0;}
                if (pow((rij+N1*a1-N2*a2).dot((rij+N1*a1-N2*a2)),0.5) < 0.5){a=1.0;b=-1.0;}
                if (pow((rij+N1*a1+N3*a3).dot((rij+N1*a1+N3*a3)),0.5) < 0.5){a=1.0;c=1.0;}
                if (pow((rij+N1*a1-N3*a3).dot((rij+N1*a1-N3*a3)),0.5) < 0.5){a=1.0;c=-1.0;}
                if (pow((rij+N1*a1+N2*a2+N3*a3).dot((rij+N1*a1+N2*a2+N3*a3)),0.5) < 0.5){a=1.0;b=1.0;c=1.0;}
                if (pow((rij+N1*a1+N2*a2-N3*a3).dot((rij+N1*a1+N2*a2-N3*a3)),0.5) < 0.5){a=1.0;b=1.0;c=-1.0;}
                if (pow((rij+N1*a1-N2*a2+N3*a3).dot((rij+N1*a1-N2*a2+N3*a3)),0.5) < 0.5){a=1.0;b=-1.0;c=1.0;}
                if (pow((rij+N1*a1-N2*a2-N3*a3).dot((rij+N1*a1-N2*a2-N3*a3)),0.5) < 0.5){a=1.0;b=-1.0;c=-1.0;}
                if (pow((rij-N1*a1).dot((rij-N1*a1)),0.5) < 0.5){a=-1.0;}
                if (pow((rij-N1*a1+N2*a2).dot((rij-N1*a1+N2*a2)),0.5) < 0.5){a=-1.0;b=1.0;}
                if (pow((rij-N1*a1-N2*a2).dot((rij-N1*a1-N2*a2)),0.5) < 0.5){a=-1.0;b=-1.0;}
                if (pow((rij-N1*a1+N3*a3).dot((rij-N1*a1+N3*a3)),0.5) < 0.5){a=-1.0;c=1.0;}
                if (pow((rij-N1*a1-N3*a3).dot((rij-N1*a1-N3*a3)),0.5) < 0.5){a=-1.0;c=-1.0;}
                if (pow((rij-N1*a1+N2*a2+N3*a3).dot((rij-N1*a1+N2*a2+N3*a3)),0.5) < 0.5){a=-1.0;b=1.0;c=1.0;}
                if (pow((rij-N1*a1+N2*a2-N3*a3).dot((rij-N1*a1+N2*a2-N3*a3)),0.5) < 0.5){a=-1.0;b=1.0;c=-1.0;}
                if (pow((rij-N1*a1-N2*a2+N3*a3).dot((rij-N1*a1-N2*a2+N3*a3)),0.5) < 0.5){a=-1.0;b=-1.0;c=1.0;}
                if (pow((rij-N1*a1-N2*a2-N3*a3).dot((rij-N1*a1-N2*a2-N3*a3)),0.5) < 0.5){a=-1.0;b=-1.0;c=-1.0;}
                if (pow((rij+N2*a2).dot((rij+N2*a2)),0.5) < 0.5){b=1.0;}
                if (pow((rij+N2*a2+N3*a3).dot((rij+N2*a2+N3*a3)),0.5) < 0.5){b=1.0;c=1.0;}
                if (pow((rij+N2*a2-N3*a3).dot((rij+N2*a2-N3*a3)),0.5) < 0.5){b=1.0;c=-1.0;}
                if (pow((rij-N2*a2).dot((rij-N2*a2)),0.5) < 0.5){b=-1.0;}
                if (pow((rij-N2*a2+N3*a3).dot((rij-N2*a2+N3*a3)),0.5) < 0.5){b=-1.0;c=1.0;}
                if (pow((rij-N2*a2-N3*a3).dot((rij-N2*a2-N3*a3)),0.5) < 0.5){b=-1.0;c=-1.0;}
                if (pow((rij+N3*a3).dot((rij+N3*a3)),0.5) < 0.5){c=1.0;}
                if (pow((rij-N3*a3).dot((rij-N3*a3)),0.5) < 0.5){c=-1.0;}
                Rij = rij +a*N1*a1 +b*N2*a2 +c*N3*a3;
                //Is this okay?
                if (i<4*N1*N2*N3 && j<4*N1*N2*N3){hk(mui,muj) += normalization*h(i,j)*exp(-I*k.dot(Rij));}
                if (i>=4*N1*N2*N3 && j<4*N1*N2*N3){hk(mui,muj) += normalization*h(i,j)*exp(-I*k.dot(Rij));}
                if (i<4*N1*N2*N3 && j>=4*N1*N2*N3){hk(mui,muj) += normalization*h(i,j)*exp(-I*k.dot(Rij));}
                if (i>=4*N1*N2*N3 && j>=4*N1*N2*N3){hk(mui,muj) += normalization*h(i,j)*exp(-I*k.dot(Rij));}
            }
        }
    }
    //cout<<hk<<endl;
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(hk);
    eigvals = pow(N1+N2+N3,-1)*es.eigenvalues();
    //double chem = chempot(eigvals).real();
    //for (int i=0;i<eigvals.size();i++){
    //    eigvals[i] = eigvals[i] - chem;
    //}

    return order(eigvals);
}
    //Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(hk);
  //  eigvals = pow(N1+N2+N3,-1)*es.eigenvalues();
    //
    //return order(eigvals);
//}

map<int, Eigen::VectorXcd> Wbands(Eigen::MatrixXcd h, int N1, int N2, int N3){
    map<int, Eigen::VectorXcd> patheigvals;
    Eigen::Vector3d k;k<<0.0,0.0,0.0;
    Eigen::Vector3d G; Eigen::Vector3d L; Eigen::Vector3d X; Eigen::Vector3d W; Eigen::Vector3d K; Eigen::Vector3d U;
    double pi = 3.1415926535;
    G<<0.0,0.0,0.0;X<<0.0,2.0*pi,0.0;L<<pi,pi,pi;W<<pi,2.0*pi,0.0;U<<pi/2.0,2.0*pi,pi/2.0;K<<3.0*pi/2.0,3.0*pi/2.0,0.0;

    for (double s = 0;s<1000.0;s++){
        if (s<333.0){k = (X+(s/333.0)*(W-X));}
        if (s>=333.0 && s<666.0){k = (W+((s-333.0)/333.0)*(G-W));}
        if (s>=666.0 && s<1000.0){k = (G+((s-666.0)/333.0)*(L-G));}
        patheigvals.insert(make_pair(s,hk8(h,k,N1,N2,N3)));
    }
    double eF = chempot(patheigvals);

    for (int s=0;s<1000;s++){
        for (int i=0;i<patheigvals[s].size();i++){
            patheigvals[s][i] = patheigvals[s][i] - eF;
        }
    }

    return patheigvals;
}

map<int, Eigen::VectorXcd> KYIbands(Eigen::MatrixXcd h, int N1, int N2, int N3){
    map<int, Eigen::VectorXcd> patheigvals;
    Eigen::Vector3d k;k<<0.0,0.0,0.0;
    Eigen::Vector3d G; Eigen::Vector3d L; Eigen::Vector3d X; Eigen::Vector3d W; Eigen::Vector3d K; Eigen::Vector3d U;
    double pi = 3.1415926535;
    G<<0.0,0.0,0.0;X<<0.0,2.0*pi,0.0;L<<pi,pi,pi;W<<pi,2.0*pi,0.0;U<<pi/2.0,2.0*pi,pi/2.0;K<<3.0*pi/2.0,3.0*pi/2.0,0.0;

    for (double s = 0.0;s<1000.0;s++){
        if (s<200.0){k = (G+((s)/200.0)*(X-G));}
        if (s>=200.0 && s<400.0){k = (X+((s-200.0)/200.0)*(W-X));}
        if (s>=400.0 && s<600.0){k = (W+((s-400.0)/200.0)*(K-W));}
        if (s>=600.0 && s<800.0){k = (K+((s-600.0)/200.0)*(L-K));}
        if (s>=800.0 && s<1000.0){k = (L+((s-800.0)/200.0)*(G-L));}
        Eigen::VectorXcd eigvals = hk8(h,k,N1,N2,N3);
        patheigvals.insert(make_pair(s,eigvals));
    }
    double eF = chempot(patheigvals);

    for (int s=0;s<1000;s++){
        for (int i=0;i<patheigvals[s].size();i++){
            patheigvals[s][i] = (patheigvals[s][i] - eF).real();
        }
    }

    return patheigvals;
}

map<int, Eigen::VectorXcd> WKKbands(Eigen::MatrixXcd h, int N1, int N2, int N3){
    map<int, Eigen::VectorXcd> patheigvals;
    Eigen::Vector3d k;k<<0.0,0.0,0.0;
    Eigen::Vector3d G; Eigen::Vector3d L; Eigen::Vector3d X; Eigen::Vector3d W; Eigen::Vector3d K; Eigen::Vector3d U;
    double pi = 3.1415926535;
    G<<0.0,0.0,0.0;X<<0.0,2.0*pi,0.0;L<<pi,pi,pi;W<<pi,2.0*pi,0.0;U<<pi/2.0,2.0*pi,pi/2.0;K<<3.0*pi/2.0,3.0*pi/2.0,0.0;

    for (double s = 0;s<1000.0;s++){
        if (s<167.0){k = (G+((s)/167.0)*(X-G));}
        if (s>=167.0 && s<333.0){k = (X+((s-167.0)/167.0)*(W-X));}
        if (s>=333.0 && s<500.0){k = (W+((s-333.0)/167.0)*(L-W));}
        if (s>=500.0 && s<667.0){k = (L+((s-500.0)/167.0)*(G-L));}
        if (s>=667.0 && s<833.0){k = (G+((s-667.0)/167.0)*(K-G));}
        if (s>=833.0 && s<1000.0){k = (K+((s-833.0)/167.0)*(X-K));}
        Eigen::VectorXcd eigvals = hk8(h,k,N1,N2,N3);
        patheigvals.insert(make_pair(s,eigvals));
    }
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(h);
    Eigen::VectorXcd eigvals = es.eigenvalues();
    double eF = chempot(patheigvals);

    for (int s=0;s<1000;s++){
        for (int i=0;i<patheigvals[s].size();i++){
            patheigvals[s][i] = (patheigvals[s][i] - eF).real();
        }
    }
    return patheigvals;
}

map<int, Eigen::VectorXcd> DFTbands(Eigen::MatrixXcd h, int N1, int N2, int N3){
    map<int, Eigen::VectorXcd> patheigvals;
    Eigen::Vector3d k;k<<0.0,0.0,0.0;
    Eigen::Vector3d G; Eigen::Vector3d L; Eigen::Vector3d X; Eigen::Vector3d W; Eigen::Vector3d K; Eigen::Vector3d U;
    double pi = 3.1415926535;
    G<<0.0,0.0,0.0;X<<0.0,2.0*pi,0.0;L<<pi,pi,pi;W<<pi,2.0*pi,0.0;U<<pi/2.0,2.0*pi,pi/2.0;K<<3.0*pi/2.0,3.0*pi/2.0,0.0;

    for (double s = 0;s<1000.0;s++){
        if (s<157.0){k = (X+((s)/157.0)*(W-X));}
        if (s>=157.0 && s<507.0){k = (W+((s-157.0)/350.0)*(G-W));}
        if (s>=507.0 && s<778.0){k = (G+((s-507.0)/271.0)*(L-G));}
        if (s>=778.0 && s<1000.0){k = (L+((s-778.0)/222.0)*(W-L));}
        Eigen::VectorXcd eigvals = hk8(h,k,N1,N2,N3);
        patheigvals.insert(make_pair(s,eigvals));
    }
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(h);
    Eigen::VectorXcd eigvals = es.eigenvalues();
    double eF = chempot(patheigvals);

    for (int s=0;s<1000;s++){
        for (int i=0;i<patheigvals[s].size();i++){
            patheigvals[s][i] = (patheigvals[s][i] - eF).real();
        }
    }
    return patheigvals;
}
