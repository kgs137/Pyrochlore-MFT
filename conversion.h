double Rho = 1.0/2.0;

//Transcribe from (n1,n2,n3,mu) to i.
int ntoi (Eigen::Vector4d v,int N1,int N2,int N3) {
    int n1=v(0); int n2=v(1); int n3=v(2); int mu=v(3);
    int i = n1 + n2*N1 + n3*N1*N2 + mu*N1*N2*N3;
    return i;
}

//Transcribe from i to (n1,n2,n3,mu).
Eigen::Vector4d iton (int i,int N1,int N2,int N3) {
    int n1, n2, n3, mu;
    mu = i/(N1*N2*N3);
    i = i - mu*(N1*N2*N3);
    n3 = i/(N1*N2);
    i = i - n3*(N1*N2);
    n2 = i/(N1);
    i = i - n2*(N1);
    n1 = i;

    Eigen::Vector4d v;
    v<<n1,n2,n3,mu;

    return v;
}

//Transcribe from (n1,n2,n3,mu) to real space.
Eigen::Vector3d ntor (Eigen::Vector4d v, string sitetype = "B"){
    double a = 1.0;
    int n1=v(0); int n2=v(1); int n3=v(2); int mu=v(3);
    Eigen::Vector3d r; Eigen::Vector3d a1;Eigen::Vector3d a2;Eigen::Vector3d a3; Eigen::Vector3d d1;Eigen::Vector3d d2;Eigen::Vector3d d3;Eigen::Vector3d d4;

    if (sitetype=="A"){a1<<0,0.5*a,0.5*a; a2<<0.5*a,0.0,0.5*a; a3<<0.5*a,0.5*a,0.0; d1=0.5*(a2+a3); d2=0.5*(a3+2.0*a2-a1); d3=0.5*(2.0*a3+a2-a1); d4=0.5*(a3+a2-a1);}
    if (sitetype=="B"){a1<<0.0,0.5*a,0.5*a; a2<<0.5*a,0.0,0.5*a; a3<<0.5*a,0.5*a,0.0; d1=0.5*a1; d2=0.5*a2; d3=0.5*a3; d4<<0,0,0;}

    Eigen::Vector3d dmu;
    if (mu==0){dmu = d4;}
    if (mu==1){dmu = d1;}
    if (mu==2){dmu = d2;}
    if (mu==3){dmu = d3;}

    r = n1*a1 + n2*a2 + n3*a3 + dmu;

    return r;
}

//Transcribe from real space to (n1,n2,n3,mu).
Eigen::Vector4d rton (Eigen::Vector3d r, string sitetype = "B"){
    double a = 1.0;

    Eigen::Matrix3d basisinverse;
    basisinverse<<
    -1/a,1/a,1/a,
    1/a,-1/a,1/a,
    1/a,1/a,-1/a;

    Eigen::Vector3d basisprojection;

    basisprojection = basisinverse*r;
    //This is the projection onto the (a1,a2,a3) basis.
    //The dmu needs to be separated out. We'll use the fact
    //that dmu is half integer of a1,a2,a3.

    int mu = 1;
    double p1 = round(basisprojection(0)-0.5);
    double p2 = round(basisprojection(1)-0.5);
    double p3 = round(basisprojection(2)-0.5);

    if (round(basisprojection(0)-0.5)==basisprojection(0)){mu=1;}
    if (round(basisprojection(2)-0.5)==basisprojection(2)){mu=2;}
    if (round(basisprojection(1)-0.5)==basisprojection(1)){mu=3;}
    if (p1==basisprojection(0) && p2==basisprojection(1) && p3==basisprojection(2)){mu=0;}

    //Now the dmu portion needs to be subtracted.
    Eigen::Vector3d a1;Eigen::Vector3d a2;Eigen::Vector3d a3; Eigen::Vector3d d1;Eigen::Vector3d d2;Eigen::Vector3d d3;Eigen::Vector3d d4;

    if (sitetype=="A"){a1<<0.0,0.5*a,0.5*a; a2<<0.5*a,0.0,0.5*a; a3<<0.5*a,0.5*a,0.0; d1=0.5*(a2+a3); d2=0.5*(a3+2.0*a2-a1); d3=0.5*(2.0*a3+a2-a1); d4=0.5*(a3+a2-a1);}
    if (sitetype=="B"){a1<<0.0,0.5*a,0.5*a; a2<<0.5*a,0.0,0.5*a; a3<<0.5*a,0.5*a,0.0; d1=0.5*a1; d2=0.5*a2; d3=0.5*a3; d4<<0.0,0.0,0.0;}

    Eigen::Vector3d dmu;
    if (mu==0){dmu = d4;}
    if (mu==1){dmu = d1;}
    if (mu==2){dmu = d2;}
    if (mu==3){dmu = d3;}

    Eigen::Vector3d integerprojection;
    integerprojection = basisprojection;// - dmu;
    //Now we can read off the n1, n2, n3.

    int n1, n2, n3;

    n1 = integerprojection(0);
    n2 = integerprojection(1);
    n3 = integerprojection(2);

    Eigen::Vector4d v;
    v<<n1,n2,n3,mu;

    return v;
}

map<int, Eigen::Vector3d> iniSpins(int N1, int N2, int N3){
    map<int, Eigen::Vector3d> spins;
    Eigen::Vector3d vect;
    for (int i=0;i<4*N1*N2*N3;i++){
        //vect<<pow(27.0,-0.5)*((rand()%200)*0.01-1),pow(27.0,-0.5)*((rand()%200)*0.01-1),pow(27.0,-0.5)*((rand()%200)*0.01-1);
        double y = pow(27.0,-0.5);
        //vect<<y*1.0,y*1.0,y*1.0;
        //int mu = iton(i,N1,N2,N3)[3];
        //if (mu == 0){vect<<1*y,1*y,1*y;};//1*y,1*y,1*y;};
        //if (mu == 1){vect<<1*y,-1*y,-1*y;};//-1*y,-1*y,1*y;};
        //if (mu == 2){vect<<-1*y,1*y,-1*y;};//-1*y,1*y,-1*y;};
        //if (mu == 3){vect<<-1*y,-1*y,1*y;};//1*y,-1*y,-1*y;};
        vect<<0.0,0.0,0.0;
        spins.insert(make_pair(i,vect));
    }
    return spins;
}

Eigen::Vector3d AngulartoCartesian(double Theta, double Phi){
    Eigen::Vector3d X;
    double x = 0.0; double y = 0.0; double z = 0.0;
    double r = Rho;
    x = r*sin(Theta)*cos(Phi); y = r*sin(Theta)*sin(Phi); z = r*cos(Theta);
    X<<x,y,z;
    return X;
}

Eigen::Vector2d CartesiantoAngular(double x, double y, double z){
    Eigen::Vector2d A;
    double Theta; double Phi;
    Theta = atan2(pow(pow(x,2.0)+pow(y,2.0),0.5),z);
    Phi = atan2(y,x);
    A<<Theta,Phi;
    return A;
}

column_vector spinstoflat(map<int, Eigen::Vector3d> spins, int N1, int N2, int N3, string type = "Cartesian"){
    column_vector flat(3* 4*N1*N2*N3);
    if (type == "Cartesian"){
        for (int i=0;i<4*N1*N2*N3;i++){
            flat(3*i+0) = spins[i][0];
            flat(3*i+1) = spins[i][1];
            flat(3*i+2) = spins[i][2];
        }
    }
    if (type == "Angular"){
        flat.set_size(2*4*N1*N2*N3);
        Eigen::Vector2d temp;
        for (int i=0;i<4*N1*N2*N3;i++){
            temp = CartesiantoAngular(spins[i][0],spins[i][1],spins[i][2]);
            flat(2*i+0) = temp[0];
            flat(2*i+1) = temp[1];
        }
    }
    return flat;
}

column_vector combine(column_vector flat_a, column_vector flat_b, int N1, int N2, int N3){
    column_vector flat_ab(5* 4*N1*N2*N3);
    for (int i = 0; i < 5* 4*N1*N2*N3;i++){
      if (i < 2*4*N1*N2*N3){
        flat_ab(i) = flat_a(i);
      }
      if (i >= 2*4*N1*N2*N3){
        flat_ab(i) = flat_b(i - 2*4*N1*N2*N3);
      }
    }
    return flat_ab;
}

column_vector separate(column_vector flat_ab, int N1, int N2, int N3, int aorb){
    column_vector flat_a(2* 4*N1*N2*N3);
    column_vector flat_b(3* 4*N1*N2*N3);
    for (int i = 0; i < 5* 4*N1*N2*N3;i++){
      if (i < 2*4*N1*N2*N3){
        flat_a(i) = flat_ab(i);
      }
      if (i >= 2*4*N1*N2*N3){
        flat_b(i - 2*4*N1*N2*N3) = flat_ab(i);
      }
    }
    if (aorb==0){return flat_a;}
    else if (aorb==1){return flat_b;}
    else{cout<<"Failed separation"; return flat_ab;}
}

map<int, Eigen::Vector3d> flattospins(column_vector flat, int N1, int N2, int N3, string type = "Cartesian"){
    map<int, Eigen::Vector3d> spins;
    Eigen::Vector3d spin;
    if (type == "Cartesian"){
        for (int i=0;i<4*N1*N2*N3;i++){
            spin<<flat(3*i+0),flat(3*i+1),flat(3*i+2);
            spins.insert(make_pair(i,spin));
        }
    }
    if (type == "Angular"){
        for (int i=0;i<4*N1*N2*N3;i++){
            spin<<AngulartoCartesian(flat(2*i+0),flat(2*i+1));
            spins.insert(make_pair(i,spin));
        }
    }
    return spins;
}

//
map<int, ALatticeSite> ASites(int N1, int N2, int N3, map<int,Eigen::Vector3d> Aspins){
    map<int, ALatticeSite> sites;
    ALatticeSite site;

    for (int i=0;i<4*N1*N2*N3;i++){
        site.i = i;
        Eigen::Vector4d nmu = iton(i,N1,N2,N3);
        site.nmu = nmu;
        site.r = ntor(nmu,"A");
        site.spin<<Aspins[i];
        sites.insert(make_pair(i,site));
    }
    return sites;
}

//
map<int, BLatticeSite> BSites(int N1, int N2, int N3, map<int, Eigen::Vector3d> Bspins){
    map<int, BLatticeSite> sites;
    BLatticeSite site;

    for (int i=0;i<4*N1*N2*N3;i++){
        site.i = i;
        Eigen::Vector4d nmu = iton(i,N1,N2,N3);
        site.nmu = nmu;
        site.r = ntor(nmu,"B");
        site.spin<<Bspins[i];
        sites.insert(make_pair(i,site));
    }
    return sites;
}
