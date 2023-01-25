double tG = -1.0;
double lambdaG = 0.67;
double JG = 0.05;
double JkG = 0.01;
double UG = 2.0;

//Set of hopping terms
map<int, TermElement> HoppingTerm(int N1, int N2, int N3, double t = tG){
    map<int, pair<int,int> > pairs = nnpairs(N1,N2,N3);
    map<int, TermElement> hoppings;
    Eigen::Matrix2cd m;
    //m<<
    //t,0,
    //0,t;

    m<<
    t,0,
    0,t;

    TermElement Term;

    int c=0;
    pair<int,int> pr;
    for (int i=0;i<pairs.size();i++){
        pr = pairs[i];
        Term.Bond = pr;
        Term.m = m;
        hoppings.insert(make_pair(c,Term));
        c++;
    }
    return hoppings;
}

//Set of spin-spin terms
map<int, TermElement> SpinSpinTerm(map<int, Eigen::Vector3d> Aspins, int N1, int N2, int N3, double J = JG){
    map<int, pair<int,int> > pairs = nnpairs(N1,N2,N3,"A","A");
    map<int, ALatticeSite> sites = ASites(N1,N2,N3,Aspins);
    map<int, TermElement> spinspinterms;
    Eigen::Matrix2cd m;

    TermElement Term;

    int c=0;
    pair<int,int> pr;
    for (int i=0;i<pairs.size();i++){
        pr = pairs[i];
        Term.Bond = pr;
        int a = pr.first; int b = pr.second;
        ALatticeSite site1 = sites[a];
        ALatticeSite site2 = sites[b];
        Eigen::Vector3d spin1 = site1.spin; Eigen::Vector3d spin2 = site2.spin;
        m<<
        J*spin1.dot(spin2),0.0,
        0.0,J*spin1.dot(spin2);

        Term.m = m;
        spinspinterms.insert(make_pair(c,Term));
        c++;
    }
    return spinspinterms;
}


//Set of spin field terms
map<int, TermElement> SpinFieldTerm(map<int, Eigen::Vector3d> Aspins, int N1, int N2, int N3, double Jk = JkG){
    map<int, pair<int,int> > pairs = nnpairs(N1,N2,N3,"A","B");
    map<int, TermElement> spinfieldterms;
    map<int, ALatticeSite> Asites = ASites(N1,N2,N3,Aspins);
    Eigen::Matrix2cd m;
    Eigen::Matrix2cd sigmax; Eigen::Matrix2cd sigmay; Eigen::Matrix2cd sigmaz;

    sigmax<<
    0.0,1.0,
    1.0,0.0;
    sigmay<<
    0.0,-I,
    I,0.0;
    sigmaz<<
    1.0,0.0,
    0.0,-1.0;

    TermElement Term;

    int c=0;
    pair<int,int> pr;
    for (int i=0;i<pairs.size();i++){
        pr = pairs[i];
        int asite = pairs[i].first;
        int bsite = pairs[i].second;
        Term.Bond = make_pair(bsite,bsite);
        Eigen::Vector3d asitespin = Asites[asite].spin;
        //cout<<asitespin[0]<<", "<<asitespin[1]<<", "<<asitespin[2]<<endl;
        m = Jk*(asitespin[0]*sigmax + asitespin[1]*sigmay + asitespin[2]*sigmaz);
        //0,0,
        //0,0;
        Term.m = m;
        spinfieldterms.insert(make_pair(c,Term));
        c++;
    }
    return spinfieldterms;
}

//Set of spin-orbit terms
map<int, TermElement> SpinOrbitTerm(int N1, int N2, int N3, double lambda = lambdaG){
    map<int, pair<int,int> > pairs = nnpairs(N1,N2,N3);
    map<int, TermElement> sohoppings;
    Eigen::Matrix2cd mat;
    Eigen::Matrix2cd sigmax; Eigen::Matrix2cd sigmay; Eigen::Matrix2cd sigmaz;

    sigmax<<
    0.0,1.0,
    1.0,0.0;
    sigmay<<
    0.0,-I,
    I,0.0;
    sigmaz<<
    1.0,0.0,
    0.0,-1.0;

    TermElement Term;

    int c=0;
    pair<int,int> pr;
    Eigen::Vector3d v;
    for (int i=0;i<pairs.size();i++){
        pr = pairs[i];
        v = nu(pairs[i],N1,N2,N3);
        Term.Bond = pr;
        mat = I*lambda*(v[0]*sigmax + v[1]*sigmay + v[2]*sigmaz);
        Term.m = mat;
        sohoppings.insert(make_pair(c,Term));
        c++;
    }
    return sohoppings;
}

//Set of repulsion terms
map<int, TermElement> RepulsionTerm(map<int, Eigen::Vector3d> Bspins, int N1, int N2, int N3, double U = UG){
    map<int, TermElement> repulsionterms;
    Eigen::Matrix2cd m;
    TermElement Term;

    map<int, BLatticeSite> Bsites = BSites(N1,N2,N3,Bspins);

    Eigen::Matrix2cd identity; Eigen::Matrix2cd sigmax; Eigen::Matrix2cd sigmay; Eigen::Matrix2cd sigmaz;

    identity<<
    1.0,0.0,
    0.0,1.0;
    sigmax<<
    0.0,1.0,
    1.0,0.0;
    sigmay<<
    0.0,-I,
    I,0.0;
    sigmaz<<
    1.0,0.0,
    0.0,-1.0;

    Eigen::Vector3d bsitespin;

    m<<
    U*0.25,0.0,
    0.0,U*0.25;
    int c = 0;

    for (int i=0;i<4*N1*N2*N3;i++){
        Term.Bond = make_pair(i, i);
        bsitespin = Bsites[i].spin;

        m = m - U*( (bsitespin[0]*sigmax + bsitespin[1]*sigmay + bsitespin[2]*sigmaz) );
        m = m + U*identity*( (bsitespin[0]*bsitespin[0] + bsitespin[1]*bsitespin[1] + bsitespin[2]*bsitespin[2]) );
        Term.m = m;
        repulsionterms.insert(make_pair(c,Term));
        c++;
    }

    return repulsionterms;
}

/*
map<int, TermElement> BSpinSpinTerm(map<int, Eigen::Vector3d> Bspins, int N1, int N2, int N3, double Jb = JbG, double Db = 0.0){
    map<int, pair<int,int> > pairs = nnpairs(N1,N2,N3,"B","B");
    map<int, BLatticeSite> sites = BSites(N1,N2,N3,Bspins);
    map<int, TermElement> spinspinterms;
    Eigen::Matrix2cd m;
    Eigen::Vector3d D; D<<0,0,0;

    TermElement Term;

    int c=0;
    pair<int,int> pr;
    for (int i=0;i<pairs.size();i++){
        pr = pairs[i];
        Term.Bond = pr;
        int a = pr.first; int b = pr.second;
        BLatticeSite site1 = sites[a];
        BLatticeSite site2 = sites[b];
        Eigen::Vector3d spin1 = site1.spin; Eigen::Vector3d spin2 = site2.spin;
        D = nu(pr,N1,N2,N3);
        m<<
        Jb*spin1.dot(spin2)+Db*D.dot(spin1.cross(spin2)),0,
        0,Jb*spin1.dot(spin2)+Db*D.dot(spin1.cross(spin2));

        Term.m = m;
        spinspinterms.insert(make_pair(c,Term));
        c++;
    }
    return spinspinterms;
}*/

//Hamiltonian builder
Eigen::MatrixXcd construct(Eigen::MatrixXcd h, map<int, TermElement> TermMap, int N1, int N2, int N3){
    int pr1;
    int pr2;
    complex<double> g;

    for (int i=0;i<TermMap.size();i++){
        Eigen::Matrix2cd mat = TermMap[i].m;

        pr1 = TermMap[i].Bond.first;
        pr2 = TermMap[i].Bond.second;
        g = mat(0,0);
        h(pr1,pr2) = h(pr1,pr2) + g;

        pr1 = TermMap[i].Bond.first;
        pr2 = TermMap[i].Bond.second;
        g = mat(1,0);
        h(pr1+4*N1*N2*N3,pr2) = h(pr1+4*N1*N2*N3,pr2) + g;

        pr1 = TermMap[i].Bond.first;
        pr2 = TermMap[i].Bond.second;
        g = mat(0,1);
        h(pr1,pr2+4*N1*N2*N3) = h(pr1,pr2+4*N1*N2*N3) + g;

        pr1 = TermMap[i].Bond.first;
        pr2 = TermMap[i].Bond.second;
        g = mat(1,1);
        h(pr1+4*N1*N2*N3,pr2+4*N1*N2*N3) = h(pr1+4*N1*N2*N3,pr2+4*N1*N2*N3) + g;
    }
    return h;
}

Eigen::MatrixXcd constructall(map<int, Eigen::Vector3d> Aspins, map<int, Eigen::Vector3d> Bspins, int N1, int N2, int N3, double t = tG, double lambda = lambdaG, double J = JG, double Jk = JG, double U = UG){
    Eigen::MatrixXcd h(2*4*N1*N2*N3,2*4*N1*N2*N3);
    for (int j=0;j<2*4*N1*N2*N3;j++){
        for (int i=0;i<2*4*N1*N2*N3;i++){
            h(i,j) = 0.0;
        }
    }
    h = construct(h, HoppingTerm(N1,N2,N3,t),N1,N2,N3);
    h = construct(h, SpinOrbitTerm(N1,N2,N3,lambda),N1,N2,N3);
    h = construct(h, SpinSpinTerm(Aspins, N1,N2,N3,J),N1,N2,N3);
    h = construct(h, SpinFieldTerm(Aspins, N1,N2,N3,Jk),N1,N2,N3);
    h = construct(h, RepulsionTerm(Bspins, N1,N2,N3,U),N1,N2,N3);
    return h;
}

Eigen::MatrixXcd constructraw(int N1, int N2, int N3, double t = tG, double lambda = lambdaG, double U = UG){
    Eigen::MatrixXcd h(2*4*N1*N2*N3,2*4*N1*N2*N3);
    for (int j=0;j<2*4*N1*N2*N3;j++){
        for (int i=0;i<2*4*N1*N2*N3;i++){
            h(i,j) = 0.0;
        }
    }
    h = construct(h, HoppingTerm(N1,N2,N3,t),N1,N2,N3);
    h = construct(h, SpinOrbitTerm(N1,N2,N3,lambda),N1,N2,N3);
    return h;
}

Eigen::MatrixXcd constructspins(Eigen::MatrixXcd h, map<int, Eigen::Vector3d> Aspins, map<int, Eigen::Vector3d> Bspins, int N1, int N2, int N3, double J = JG, double Jk = JG, double U = UG){
    h = construct(h, SpinSpinTerm(Aspins, N1,N2,N3,J),N1,N2,N3);
    h = construct(h, SpinFieldTerm(Aspins, N1,N2,N3,Jk),N1,N2,N3);
    h = construct(h, RepulsionTerm(Bspins, N1,N2,N3,U),N1,N2,N3);

    return h;
}
