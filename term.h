//Finding nearest neighbor pairs.
map<int, pair<int,int> > nnpairs(int N1, int N2, int N3, string sitetype1 = "B", string sitetype2 = "B"){
    Eigen::MatrixXd positions(4*N1*N2*N3,3);
    Eigen::MatrixXd positionsA(4*N1*N2*N3,3);
    Eigen::MatrixXd positionsB(4*N1*N2*N3,3);
    
    Eigen::Vector3d p;
    
    //Constructing a list of positions.
    for (int i=0;i<4*N1*N2*N3;i++){
        if (sitetype1==sitetype2){
            if (sitetype1=="A"){
                p = ntor(iton(i,N1,N2,N3),"A");
                positions(i,0) = p(0);
                positions(i,1) = p(1);
                positions(i,2) = p(2);
            }
            if (sitetype1=="B"){
                p = ntor(iton(i,N1,N2,N3),"B");
                positions(i,0) = p(0);
                positions(i,1) = p(1);
                positions(i,2) = p(2);
            }
        }else{
            p = ntor(iton(i,N1,N2,N3),"A");
            positionsA(i,0) = p(0);
            positionsA(i,1) = p(1);
            positionsA(i,2) = p(2);
            p = ntor(iton(i,N1,N2,N3),"B");
            positionsB(i,0) = p(0);
            positionsB(i,1) = p(1);
            positionsB(i,2) = p(2);
        }
    }
    
    Eigen::Vector3d ri;
    Eigen::Vector3d rj;
    map<int, pair<int,int> >  pairs;
    pair<int, int> pr;
    int c = 0;
    
    Eigen::Vector4d ni;
    double a = 1.0;
    Eigen::Vector3d a1; Eigen::Vector3d a2; Eigen::Vector3d a3;
    a1<<0,0.5*a,0.5*a; a2<<0.5*a,0,0.5*a; a3<<0.5*a,0.5*a,0;
    
    //Comparing pairwise separation.
    
    if (sitetype1==sitetype2){
        for (int i=0;i<(4*N1*N2*N3);i++){
            ri(0) = positions(i,0);
            ri(1) = positions(i,1);
            ri(2) = positions(i,2);
            for (int j=0;j<4*N1*N2*N3;j++){
                rj(0) = positions(j,0);
                rj(1) = positions(j,1);
                rj(2) = positions(j,2);
                if (pow(((rj-ri).dot(rj-ri)),0.5)<0.5 && pow(((rj-ri).dot(rj-ri)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri+N1*a1).dot(rj-ri+N1*a1)),0.5)<0.5 && pow(((rj-ri+N1*a1).dot(rj-ri+N1*a1)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri+N1*a1+N2*a2).dot(rj-ri+N1*a1+N2*a2)),0.5)<0.5 && pow(((rj-ri+N1*a1+N2*a2).dot(rj-ri+N1*a1+N2*a2)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri+N1*a1-N2*a2).dot(rj-ri+N1*a1-N2*a2)),0.5)<0.5 && pow(((rj-ri+N1*a1-N2*a2).dot(rj-ri+N1*a1-N2*a2)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri+N1*a1+N3*a3).dot(rj-ri+N1*a1+N3*a3)),0.5)<0.5 && pow(((rj-ri+N1*a1+N3*a3).dot(rj-ri+N1*a1+N3*a3)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri+N1*a1-N3*a3).dot(rj-ri+N1*a1-N3*a3)),0.5)<0.5 && pow(((rj-ri+N1*a1-N3*a3).dot(rj-ri+N1*a1-N3*a3)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri+N1*a1+N2*a2+N3*a3).dot(rj-ri+N1*a1+N2*a2+N3*a3)),0.5)<0.5 && pow(((rj-ri+N1*a1+N2*a2+N3*a3).dot(rj-ri+N1*a1+N2*a2+N3*a3)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri+N1*a1+N2*a2-N3*a3).dot(rj-ri+N1*a1+N2*a2-N3*a3)),0.5)<0.5 && pow(((rj-ri+N1*a1+N2*a2-N3*a3).dot(rj-ri+N1*a1+N2*a2-N3*a3)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri+N1*a1-N2*a2+N3*a3).dot(rj-ri+N1*a1-N2*a2+N3*a3)),0.5)<0.5 && pow(((rj-ri+N1*a1-N2*a2+N3*a3).dot(rj-ri+N1*a1-N2*a2+N3*a3)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri+N1*a1-N2*a2-N3*a3).dot(rj-ri+N1*a1-N2*a2-N3*a3)),0.5)<0.5 && pow(((rj-ri+N1*a1-N2*a2-N3*a3).dot(rj-ri+N1*a1-N2*a2-N3*a3)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri-N1*a1).dot(rj-ri-N1*a1)),0.5)<0.5 && pow(((rj-ri-N1*a1).dot(rj-ri-N1*a1)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri-N1*a1+N2*a2).dot(rj-ri-N1*a1+N2*a2)),0.5)<0.5 && pow(((rj-ri-N1*a1+N2*a2).dot(rj-ri-N1*a1+N2*a2)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri-N1*a1-N2*a2).dot(rj-ri-N1*a1-N2*a2)),0.5)<0.5 && pow(((rj-ri-N1*a1-N2*a2).dot(rj-ri-N1*a1-N2*a2)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri-N1*a1+N3*a3).dot(rj-ri-N1*a1+N3*a3)),0.5)<0.5 && pow(((rj-ri-N1*a1+N3*a3).dot(rj-ri-N1*a1+N3*a3)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri-N1*a1-N3*a3).dot(rj-ri-N1*a1-N3*a3)),0.5)<0.5 && pow(((rj-ri-N1*a1-N3*a3).dot(rj-ri-N1*a1-N3*a3)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri-N1*a1+N2*a2+N3*a3).dot(rj-ri-N1*a1+N2*a2+N3*a3)),0.5)<0.5 && pow(((rj-ri-N1*a1+N2*a2+N3*a3).dot(rj-ri-N1*a1+N2*a2+N3*a3)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri-N1*a1+N2*a2-N3*a3).dot(rj-ri-N1*a1+N2*a2-N3*a3)),0.5)<0.5 && pow(((rj-ri-N1*a1+N2*a2-N3*a3).dot(rj-ri-N1*a1+N2*a2-N3*a3)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri-N1*a1-N2*a2+N3*a3).dot(rj-ri-N1*a1-N2*a2+N3*a3)),0.5)<0.5 && pow(((rj-ri-N1*a1-N2*a2+N3*a3).dot(rj-ri-N1*a1-N2*a2+N3*a3)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri-N1*a1-N2*a2-N3*a3).dot(rj-ri-N1*a1-N2*a2-N3*a3)),0.5)<0.5 && pow(((rj-ri-N1*a1-N2*a2-N3*a3).dot(rj-ri-N1*a1-N2*a2-N3*a3)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri+N2*a2).dot(rj-ri+N2*a2)),0.5)<0.5 && pow(((rj-ri+N2*a2).dot(rj-ri+N2*a2)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri-N2*a2).dot(rj-ri-N2*a2)),0.5)<0.5 && pow(((rj-ri-N2*a2).dot(rj-ri-N2*a2)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri+N2*a2+N3*a3).dot(rj-ri+N2*a2+N3*a3)),0.5)<0.5 && pow(((rj-ri+N2*a2+N3*a3).dot(rj-ri+N2*a2+N3*a3)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri+N2*a2-N3*a3).dot(rj-ri+N2*a2-N3*a3)),0.5)<0.5 && pow(((rj-ri+N2*a2-N3*a3).dot(rj-ri+N2*a2-N3*a3)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri+N3*a3).dot(rj-ri+N3*a3)),0.5)<0.5 && pow(((rj-ri+N3*a3).dot(rj-ri+N3*a3)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri-N3*a3).dot(rj-ri-N3*a3)),0.5)<0.5 && pow(((rj-ri-N3*a3).dot(rj-ri-N3*a3)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                //These were missed earlier:
                if (pow(((rj-ri-N2*a2+N3*a3).dot(rj-ri-N2*a2+N3*a3)),0.5)<0.5 && pow(((rj-ri-N2*a2+N3*a3).dot(rj-ri-N2*a2+N3*a3)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri-N2*a2-N3*a3).dot(rj-ri-N2*a2-N3*a3)),0.5)<0.5 && pow(((rj-ri-N2*a2-N3*a3).dot(rj-ri-N2*a2-N3*a3)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
            }
        }
        //Nearest neighbors amongst A and B
    }else{
        for (int i=0;i<(4*N1*N2*N3);i++){
            ri(0) = positionsA(i,0);
            ri(1) = positionsA(i,1);
            ri(2) = positionsA(i,2);
            for (int j=0;j<4*N1*N2*N3;j++){
                rj(0) = positionsB(j,0);
                rj(1) = positionsB(j,1);
                rj(2) = positionsB(j,2);
                if (pow(((rj-ri).dot(rj-ri)),0.5)<0.5 && pow(((rj-ri).dot(rj-ri)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri+N1*a1).dot(rj-ri+N1*a1)),0.5)<0.5 && pow(((rj-ri+N1*a1).dot(rj-ri+N1*a1)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri+N1*a1+N2*a2).dot(rj-ri+N1*a1+N2*a2)),0.5)<0.5 && pow(((rj-ri+N1*a1+N2*a2).dot(rj-ri+N1*a1+N2*a2)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri+N1*a1-N2*a2).dot(rj-ri+N1*a1-N2*a2)),0.5)<0.5 && pow(((rj-ri+N1*a1-N2*a2).dot(rj-ri+N1*a1-N2*a2)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri+N1*a1+N3*a3).dot(rj-ri+N1*a1+N3*a3)),0.5)<0.5 && pow(((rj-ri+N1*a1+N3*a3).dot(rj-ri+N1*a1+N3*a3)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri+N1*a1-N3*a3).dot(rj-ri+N1*a1-N3*a3)),0.5)<0.5 && pow(((rj-ri+N1*a1-N3*a3).dot(rj-ri+N1*a1-N3*a3)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri+N1*a1+N2*a2+N3*a3).dot(rj-ri+N1*a1+N2*a2+N3*a3)),0.5)<0.5 && pow(((rj-ri+N1*a1+N2*a2+N3*a3).dot(rj-ri+N1*a1+N2*a2+N3*a3)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri+N1*a1+N2*a2-N3*a3).dot(rj-ri+N1*a1+N2*a2-N3*a3)),0.5)<0.5 && pow(((rj-ri+N1*a1+N2*a2-N3*a3).dot(rj-ri+N1*a1+N2*a2-N3*a3)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri+N1*a1-N2*a2+N3*a3).dot(rj-ri+N1*a1-N2*a2+N3*a3)),0.5)<0.5 && pow(((rj-ri+N1*a1-N2*a2+N3*a3).dot(rj-ri+N1*a1-N2*a2+N3*a3)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri+N1*a1-N2*a2-N3*a3).dot(rj-ri+N1*a1-N2*a2-N3*a3)),0.5)<0.5 && pow(((rj-ri+N1*a1-N2*a2-N3*a3).dot(rj-ri+N1*a1-N2*a2-N3*a3)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri-N1*a1).dot(rj-ri-N1*a1)),0.5)<0.5 && pow(((rj-ri-N1*a1).dot(rj-ri-N1*a1)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri-N1*a1+N2*a2).dot(rj-ri-N1*a1+N2*a2)),0.5)<0.5 && pow(((rj-ri-N1*a1+N2*a2).dot(rj-ri-N1*a1+N2*a2)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri-N1*a1-N2*a2).dot(rj-ri-N1*a1-N2*a2)),0.5)<0.5 && pow(((rj-ri-N1*a1-N2*a2).dot(rj-ri-N1*a1-N2*a2)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri-N1*a1+N3*a3).dot(rj-ri-N1*a1+N3*a3)),0.5)<0.5 && pow(((rj-ri-N1*a1+N3*a3).dot(rj-ri-N1*a1+N3*a3)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri-N1*a1-N3*a3).dot(rj-ri-N1*a1-N3*a3)),0.5)<0.5 && pow(((rj-ri-N1*a1-N3*a3).dot(rj-ri-N1*a1-N3*a3)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri-N1*a1+N2*a2+N3*a3).dot(rj-ri-N1*a1+N2*a2+N3*a3)),0.5)<0.5 && pow(((rj-ri-N1*a1+N2*a2+N3*a3).dot(rj-ri-N1*a1+N2*a2+N3*a3)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri-N1*a1+N2*a2-N3*a3).dot(rj-ri-N1*a1+N2*a2-N3*a3)),0.5)<0.5 && pow(((rj-ri-N1*a1+N2*a2-N3*a3).dot(rj-ri-N1*a1+N2*a2-N3*a3)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri-N1*a1-N2*a2+N3*a3).dot(rj-ri-N1*a1-N2*a2+N3*a3)),0.5)<0.5 && pow(((rj-ri-N1*a1-N2*a2+N3*a3).dot(rj-ri-N1*a1-N2*a2+N3*a3)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri-N1*a1-N2*a2-N3*a3).dot(rj-ri-N1*a1-N2*a2-N3*a3)),0.5)<0.5 && pow(((rj-ri-N1*a1-N2*a2-N3*a3).dot(rj-ri-N1*a1-N2*a2-N3*a3)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri+N2*a2).dot(rj-ri+N2*a2)),0.5)<0.5 && pow(((rj-ri+N2*a2).dot(rj-ri+N2*a2)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri-N2*a2).dot(rj-ri-N2*a2)),0.5)<0.5 && pow(((rj-ri-N2*a2).dot(rj-ri-N2*a2)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri+N2*a2+N3*a3).dot(rj-ri+N2*a2+N3*a3)),0.5)<0.5 && pow(((rj-ri+N2*a2+N3*a3).dot(rj-ri+N2*a2+N3*a3)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri+N2*a2-N3*a3).dot(rj-ri+N2*a2-N3*a3)),0.5)<0.5 && pow(((rj-ri+N2*a2-N3*a3).dot(rj-ri+N2*a2-N3*a3)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri+N3*a3).dot(rj-ri+N3*a3)),0.5)<0.5 && pow(((rj-ri+N3*a3).dot(rj-ri+N3*a3)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri-N3*a3).dot(rj-ri-N3*a3)),0.5)<0.5 && pow(((rj-ri-N3*a3).dot(rj-ri-N3*a3)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri-N2*a2+N3*a3).dot(rj-ri-N2*a2+N3*a3)),0.5)<0.5 && pow(((rj-ri-N2*a2+N3*a3).dot(rj-ri-N2*a2+N3*a3)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
                if (pow(((rj-ri-N2*a2-N3*a3).dot(rj-ri-N2*a2-N3*a3)),0.5)<0.5 && pow(((rj-ri-N2*a2-N3*a3).dot(rj-ri-N2*a2-N3*a3)),0.5)!=0){
                    pr.first = i;
                    pr.second = j;
                    pairs.insert(make_pair(c,pr));
                    c++;
                }
            }
        }
    }
    
    return pairs;
}

//Spin-Orbit Nu vector given nearest neighbor pair
Eigen::VectorXd nu(pair<int, int> nnpair, int N1, int N2, int N3){
    Eigen::Vector3d nu;
    int i = nnpair.first;
    int j = nnpair.second;
    
    Eigen::Vector4d nmi = iton(i,N1,N2,N3);
    Eigen::Vector4d nmj = iton(j,N1,N2,N3);
    Eigen::Vector3d ri = ntor(nmi);
    Eigen::Vector3d rj = ntor(nmj);
    
    double a = 1.0;
    Eigen::Vector3d R1; Eigen::Vector3d R2; Eigen::Vector3d R3; Eigen::Vector3d R4; Eigen::Vector3d R5; Eigen::Vector3d R6;
    R1<<0.0,0.25,0.25;R2<<0.25,0.0,0.25;R3<<0.25,0.25,0.0;R4<<0.0,-0.25,0.25;R5<<0.25,0.0,-0.25;R6<<-0.25,0.25,0.0;
    Eigen::Vector3d a1; Eigen::Vector3d a2; Eigen::Vector3d a3; a1<<0,0.5*a,0.5*a; a2<<0.5*a,0,0.5*a; a3<<0.5*a,0.5*a,0;
    
    double A = pow(2.0,-0.5);   //Ensures unit-vector result
    Eigen::Vector3d rij; rij<<0,0,0;
    
        if (pow(((rj-ri).dot(rj-ri)),0.5)<0.5 && pow(((rj-ri).dot(rj-ri)),0.5)!=0){
            rij = rj-ri;
        }
        if (pow(((rj-ri+N1*a1).dot(rj-ri+N1*a1)),0.5)<0.5 && pow(((rj-ri+N1*a1).dot(rj-ri+N1*a1)),0.5)!=0){
            rij = rj-ri+N1*a1;
        }
        if (pow(((rj-ri+N1*a1+N2*a2).dot(rj-ri+N1*a1+N2*a2)),0.5)<0.5 && pow(((rj-ri+N1*a1+N2*a2).dot(rj-ri+N1*a1+N2*a2)),0.5)!=0){
            rij = rj-ri+N1*a1+N2*a2;
        }
        if (pow(((rj-ri+N1*a1-N2*a2).dot(rj-ri+N1*a1-N2*a2)),0.5)<0.5 && pow(((rj-ri+N1*a1-N2*a2).dot(rj-ri+N1*a1-N2*a2)),0.5)!=0){
            rij = rj-ri+N1*a1-N2*a2;
        }
        if (pow(((rj-ri+N1*a1+N3*a3).dot(rj-ri+N1*a1+N3*a3)),0.5)<0.5 && pow(((rj-ri+N1*a1+N3*a3).dot(rj-ri+N1*a1+N3*a3)),0.5)!=0){
            rij = rj-ri+N1*a1+N3*a3;
        }
        if (pow(((rj-ri+N1*a1-N3*a3).dot(rj-ri+N1*a1-N3*a3)),0.5)<0.5 && pow(((rj-ri+N1*a1-N3*a3).dot(rj-ri+N1*a1-N3*a3)),0.5)!=0){
            rij = rj-ri+N1*a1-N3*a3;
        }
        if (pow(((rj-ri+N1*a1+N2*a2+N3*a3).dot(rj-ri+N1*a1+N2*a2+N3*a3)),0.5)<0.5 && pow(((rj-ri+N1*a1+N2*a2+N3*a3).dot(rj-ri+N1*a1+N2*a2+N3*a3)),0.5)!=0){
            rij = rj-ri+N1*a1+N2*a2+N3*a3;
        }
        if (pow(((rj-ri+N1*a1+N2*a2-N3*a3).dot(rj-ri+N1*a1+N2*a2-N3*a3)),0.5)<0.5 && pow(((rj-ri+N1*a1+N2*a2-N3*a3).dot(rj-ri+N1*a1+N2*a2-N3*a3)),0.5)!=0){
            rij = rj-ri+N1*a1+N2*a2-N3*a3;
        }
        if (pow(((rj-ri+N1*a1-N2*a2+N3*a3).dot(rj-ri+N1*a1-N2*a2+N3*a3)),0.5)<0.5 && pow(((rj-ri+N1*a1-N2*a2+N3*a3).dot(rj-ri+N1*a1-N2*a2+N3*a3)),0.5)!=0){
            rij = rj-ri+N1*a1-N2*a2+N3*a3;
        }
        if (pow(((rj-ri+N1*a1-N2*a2-N3*a3).dot(rj-ri+N1*a1-N2*a2-N3*a3)),0.5)<0.5 && pow(((rj-ri+N1*a1-N2*a2-N3*a3).dot(rj-ri+N1*a1-N2*a2-N3*a3)),0.5)!=0){
            rij = rj-ri+N1*a1-N2*a2-N3*a3;
        }
        if (pow(((rj-ri-N1*a1).dot(rj-ri-N1*a1)),0.5)<0.5 && pow(((rj-ri-N1*a1).dot(rj-ri-N1*a1)),0.5)!=0){
            rij = rj-ri-N1*a1;
        }
        if (pow(((rj-ri-N1*a1+N2*a2).dot(rj-ri-N1*a1+N2*a2)),0.5)<0.5 && pow(((rj-ri-N1*a1+N2*a2).dot(rj-ri-N1*a1+N2*a2)),0.5)!=0){
            rij = rj-ri-N1*a1+N2*a2;
        }
        if (pow(((rj-ri-N1*a1-N2*a2).dot(rj-ri-N1*a1-N2*a2)),0.5)<0.5 && pow(((rj-ri-N1*a1-N2*a2).dot(rj-ri-N1*a1-N2*a2)),0.5)!=0){
            rij = rj-ri-N1*a1-N2*a2;
        }
        if (pow(((rj-ri-N1*a1+N3*a3).dot(rj-ri-N1*a1+N3*a3)),0.5)<0.5 && pow(((rj-ri-N1*a1+N3*a3).dot(rj-ri-N1*a1+N3*a3)),0.5)!=0){
            rij = rj-ri-N1*a1+N3*a3;
        }
        if (pow(((rj-ri-N1*a1-N3*a3).dot(rj-ri-N1*a1-N3*a3)),0.5)<0.5 && pow(((rj-ri-N1*a1-N3*a3).dot(rj-ri-N1*a1-N3*a3)),0.5)!=0){
            rij = rj-ri-N1*a1-N3*a3;
        }
        if (pow(((rj-ri-N1*a1+N2*a2+N3*a3).dot(rj-ri-N1*a1+N2*a2+N3*a3)),0.5)<0.5 && pow(((rj-ri-N1*a1+N2*a2+N3*a3).dot(rj-ri-N1*a1+N2*a2+N3*a3)),0.5)!=0){
            rij = rj-ri-N1*a1+N2*a2+N3*a3;
        }
        if (pow(((rj-ri-N1*a1+N2*a2-N3*a3).dot(rj-ri-N1*a1+N2*a2-N3*a3)),0.5)<0.5 && pow(((rj-ri-N1*a1+N2*a2-N3*a3).dot(rj-ri-N1*a1+N2*a2-N3*a3)),0.5)!=0){
            rij = rj-ri-N1*a1+N2*a2-N3*a3;
        }
        if (pow(((rj-ri-N1*a1-N2*a2+N3*a3).dot(rj-ri-N1*a1-N2*a2+N3*a3)),0.5)<0.5 && pow(((rj-ri-N1*a1-N2*a2+N3*a3).dot(rj-ri-N1*a1-N2*a2+N3*a3)),0.5)!=0){
            rij = rj-ri-N1*a1-N2*a2+N3*a3;
        }
        if (pow(((rj-ri-N1*a1-N2*a2-N3*a3).dot(rj-ri-N1*a1-N2*a2-N3*a3)),0.5)<0.5 && pow(((rj-ri-N1*a1-N2*a2-N3*a3).dot(rj-ri-N1*a1-N2*a2-N3*a3)),0.5)!=0){
            rij = rj-ri-N1*a1-N2*a2-N3*a3;
        }
        if (pow(((rj-ri+N2*a2).dot(rj-ri+N2*a2)),0.5)<0.5 && pow(((rj-ri+N2*a2).dot(rj-ri+N2*a2)),0.5)!=0){
            rij = rj-ri+N2*a2;
        }
        if (pow(((rj-ri-N2*a2).dot(rj-ri-N2*a2)),0.5)<0.5 && pow(((rj-ri-N2*a2).dot(rj-ri-N2*a2)),0.5)!=0){
            rij = rj-ri-N2*a2;
        }
        if (pow(((rj-ri+N2*a2+N3*a3).dot(rj-ri+N2*a2+N3*a3)),0.5)<0.5 && pow(((rj-ri+N2*a2+N3*a3).dot(rj-ri+N2*a2+N3*a3)),0.5)!=0){
            rij = rj-ri+N2*a2+N3*a3;
        }
        if (pow(((rj-ri+N2*a2-N3*a3).dot(rj-ri+N2*a2-N3*a3)),0.5)<0.5 && pow(((rj-ri+N2*a2-N3*a3).dot(rj-ri+N2*a2-N3*a3)),0.5)!=0){
            rij = rj-ri+N2*a2-N3*a3;
        }
        if (pow(((rj-ri+N3*a3).dot(rj-ri+N3*a3)),0.5)<0.5 && pow(((rj-ri+N3*a3).dot(rj-ri+N3*a3)),0.5)!=0){
            rij = rj-ri+N3*a3;
        }
        if (pow(((rj-ri-N3*a3).dot(rj-ri-N3*a3)),0.5)<0.5 && pow(((rj-ri-N3*a3).dot(rj-ri-N3*a3)),0.5)!=0){
            rij = rj-ri-N3*a3;
        }
        if (pow(((rj-ri-N2*a2+N3*a3).dot(rj-ri-N2*a2+N3*a3)),0.5)<0.5 && pow(((rj-ri-N2*a2+N3*a3).dot(rj-ri-N2*a2+N3*a3)),0.5)!=0){
            rij = rj-ri-N2*a2+N3*a3;
        }
        if (pow(((rj-ri-N2*a2-N3*a3).dot(rj-ri-N2*a2-N3*a3)),0.5)<0.5 && pow(((rj-ri-N2*a2-N3*a3).dot(rj-ri-N2*a2-N3*a3)),0.5)!=0){
            rij = rj-ri-N2*a2-N3*a3;
        }
    
    //Right Side
    if(nmi(0)==nmj(0) && nmi(1)==nmj(1) && nmi(2)==nmj(2)){     //if(same cell)
        A = A;}         //A = A;
    //Left Side
    if(nmi(0)!=nmj(0) || nmi(1)!=nmj(1) || nmi(2)!=nmj(2)){     //if(not)
        A = -A;}        //A = -A;
    
    if(rij == R1){nu=A*R4;
    }else if(rij == R2){nu=A*R5;
    }else if(rij == R3){nu=A*R6;
    }else if(rij == R4){nu=A*R1;
    }else if(rij == R5){nu=A*R2;
    }else if(rij == R6){nu=A*R3;
    }else if(rij == -R1){nu=-A*R4;
    }else if(rij == -R2){nu=-A*R5;
    }else if(rij == -R3){nu=-A*R6;
    }else if(rij == -R4){nu=-A*R1;
    }else if(rij == -R5){nu=-A*R2;
    }else if(rij == -R6){nu=-A*R3;}
    
    return pow(2.0, 0.5)*4.0*nu;
}

/*
//Spin-Orbit Nu vector given nearest neighbor pair
Eigen::VectorXd nu(pair<int, int> nnpair, int N1, int N2, int N3){
    Eigen::Vector3d nu;
    int i = nnpair.first;
    int j = nnpair.second;
 
    Eigen::Vector4d nmi = iton(i,N1,N2,N3);
    Eigen::Vector4d nmj = iton(j,N1,N2,N3);
    Eigen::Vector3d ri = ntor(nmi);
    Eigen::Vector3d rj = ntor(nmj);
    Eigen::Vector3d r = (rj - ri);      //Preliminary
    //cout<<pow(r.dot(r),0.5)<<endl;
    
    Eigen::Vector3d R1; Eigen::Vector3d R2; Eigen::Vector3d R3; Eigen::Vector3d R4; Eigen::Vector3d R5; Eigen::Vector3d R6;
    R1<<0.0,0.25,0.25;R2<<0.25,0.0,0.25;R3<<0.25,0.25,0.0;R4<<0.0,-0.25,0.25;R5<<0.25,0.0,-0.25;R6<<-0.25,0.25,0.0;
    
    double A = pow(2.0,-0.5);   //Ensures unit-vector result
    
    if(r == R1){nu=A*R4;
    }else if(r == R2){nu=A*R5;
    }else if(r == R3){nu=A*R6;
    }else if(r == R4){nu=A*R1;
    }else if(r == R5){nu=A*R2;
    }else if(r == R6){nu=A*R3;
    }else if(r == -R1){nu=-A*R4;
    }else if(r == -R2){nu=-A*R5;
    }else if(r == -R3){nu=-A*R6;
    }else if(r == -R4){nu=-A*R1;
    }else if(r == -R5){nu=-A*R2;
    }else if(r == -R6){nu=-A*R3;
    }else if(pow(r.dot(r),0.5)>=0.5){               //Shift accounting for boundary pairs
        
        ri = ntor(nmi);
        rj = ntor(nmj);
        
        int shift = 0;
        int s = 0;
        
        if(nmi(0)==N1-1){ri = ri - 2.0*N1*R1;shift=1;}
        if(nmi(1)==N2-1){ri = ri - 2.0*N2*R2;shift=1;}
        if(nmi(2)==N3-1){ri = ri - 2.0*N3*R3;shift=1;}
        //if(nmj(0)==N1-1){rj = rj - 2.0*N1*R1;shift=1;}
        //if(nmj(1)==N2-1){rj = rj - 2.0*N2*R2;shift=1;}
        //if(nmj(2)==N3-1){rj = rj - 2.0*N3*R3;shift=1;}
        
        if(shift==0){cout<<"Unshifted ("<<i<<", "<<j<<") "<<endl<<r<<endl;}
        r<<0,0,0;
        r = (rj - ri);      //Adjusted
    
        if(r == R1){nu=A*R4;
        }else if(r == R2){nu=A*R5;s=1;
        }else if(r == R3){nu=A*R6;s=1;
        }else if(r == R4){nu=A*R1;s=1;
        }else if(r == R5){nu=A*R2;s=1;
        }else if(r == R6){nu=A*R3;s=1;
        }else if(r == -R1){nu=-A*R4;s=1;
        }else if(r == -R2){nu=-A*R5;s=1;
        }else if(r == -R3){nu=-A*R6;s=1;
        }else if(r == -R4){nu=-A*R1;s=1;
        }else if(r == -R5){nu=-A*R2;s=1;
        }else if(r == -R6){nu=-A*R3;s=1;}
        
        ri = ntor(nmi);
        rj = ntor(nmj);
        
        //if(nmi(0)==N1-1){ri = ri - 2.0*N1*R1;shift=1;}
        //if(nmi(1)==N2-1){ri = ri - 2.0*N2*R2;shift=1;}
        //if(nmi(2)==N3-1){ri = ri - 2.0*N3*R3;shift=1;}
        if(nmj(0)==N1-1){rj = rj - 2.0*N1*R1;shift=1;}
        if(nmj(1)==N2-1){rj = rj - 2.0*N2*R2;shift=1;}
        if(nmj(2)==N3-1){rj = rj - 2.0*N3*R3;shift=1;}
        
        if(shift==0){cout<<"Unshifted 1 ("<<i<<", "<<j<<") "<<endl<<r<<endl;}
        r<<0,0,0;
        r = (rj - ri);      //Adjusted
        
        if(r == R1){nu=A*R4;
        }else if(r == R2){nu=A*R5;s=1;
        }else if(r == R3){nu=A*R6;s=1;
        }else if(r == R4){nu=A*R1;s=1;
        }else if(r == R5){nu=A*R2;s=1;
        }else if(r == R6){nu=A*R3;s=1;
        }else if(r == -R1){nu=-A*R4;s=1;
        }else if(r == -R2){nu=-A*R5;s=1;
        }else if(r == -R3){nu=-A*R6;s=1;
        }else if(r == -R4){nu=-A*R1;s=1;
        }else if(r == -R5){nu=-A*R2;s=1;
        }else if(r == -R6){nu=-A*R3;s=1;}
        
        ri = ntor(nmi);
        rj = ntor(nmj);
        
        if(nmi(0)==N1-1){ri = ri - 2.0*N1*R1;shift=1;}
        if(nmi(1)==N2-1){ri = ri - 2.0*N2*R2;shift=1;}
        if(nmi(2)==N3-1){ri = ri - 2.0*N3*R3;shift=1;}
        if(nmj(0)==N1-1){rj = rj - 2.0*N1*R1;shift=1;}
        if(nmj(1)==N2-1){rj = rj - 2.0*N2*R2;shift=1;}
        if(nmj(2)==N3-1){rj = rj - 2.0*N3*R3;shift=1;}
        
        if(shift==0){cout<<"Unshifted 2 ("<<i<<", "<<j<<") "<<endl<<r<<endl;}
        r<<0,0,0;
        r = (rj - ri);      //Adjusted
        
        if(r == R1){nu=A*R4;
        }else if(r == R2){nu=A*R5;s=1;
        }else if(r == R3){nu=A*R6;s=1;
        }else if(r == R4){nu=A*R1;s=1;
        }else if(r == R5){nu=A*R2;s=1;
        }else if(r == R6){nu=A*R3;s=1;
        }else if(r == -R1){nu=-A*R4;s=1;
        }else if(r == -R2){nu=-A*R5;s=1;
        }else if(r == -R3){nu=-A*R6;s=1;
        }else if(r == -R4){nu=-A*R1;s=1;
        }else if(r == -R5){nu=-A*R2;s=1;
        }else if(r == -R6){nu=-A*R3;s=1;}
        
        if(s==0){cout<<"Missed shifted nu ("<<i<<", "<<j<<") "<<endl<<r<<endl;}
    }else{cout<<"Missed one."<<endl;}
    
    return 2.0*4.0*nu;
}*/

Eigen::Vector3cd tau(int i, int mu, int N1, int N2, int N3){
    Eigen::Vector3cd tauexpectation;
    Eigen::Vector3cd tuu, tud, tdu, tdd;
    
    //spinor.tau.spinor
    tuu<< 0.0,0.0,1.0;
    tud<< 1.0,-I,0.0;
    tdu<< 1.0,I,0.0;
    tdd<< 0.0,0.0,-1.0;
    if (i<4*N1*N2*N3 && mu<4*N1*N2*N3){tauexpectation = tuu;};
    if (i<4*N1*N2*N3 && mu>=4*N1*N2*N3){tauexpectation = tud;};
    if (i>=4*N1*N2*N3 && mu<4*N1*N2*N3){tauexpectation = tdu;};
    if (i>=4*N1*N2*N3 && mu>=4*N1*N2*N3){tauexpectation = tdd;};
    return tauexpectation;
}

map<int, int> nearto(int origin, string origintype,string neighbortype,int N1,int N2,int N3){
    map<int, pair<int,int> > neighbors;
    map<int, int> nearest;
    int c = 0;
    
    if (origintype == "A"){
        if (neighbortype =="A"){neighbors = nnpairs(N1,N2,N3,"A","A");};
        if (neighbortype =="B"){neighbors = nnpairs(N1,N2,N3,"A","B");};}
    if (origintype =="B"){
        if (neighbortype =="A"){neighbors = nnpairs(N1,N2,N3,"B","A");};
        if (neighbortype =="B"){neighbors = nnpairs(N1,N2,N3,"B","B");};}
    
    for (int n=0;n<neighbors.size();n++){
        if (neighbors[n].first == origin){
            nearest.insert(make_pair(c,neighbors[n].second));
            c++;
        }
    }
    return nearest;
}

Eigen::VectorXcd order(Eigen::VectorXcd vals){
    Eigen::VectorXcd orderedvals(vals.size());
    int place; int same;
    for(int i=0;i<vals.size();i++){
        place = 0; same = 0;
        for(int j=0;j<vals.size();j++){
            if(vals[j].real()<vals[i].real()){place++;}
            if(vals[j].real()==vals[i].real()){same++;}
        }
        for(int s=0;s<same;s++){
            orderedvals[place+s] = vals[i];
        }
    }
    return orderedvals;
}

double chempot(map<int, Eigen::VectorXcd> orderedpath){
    //This is half-filling. Hard coded for 8 bands. Allow for filling by adjusting valence and conduct.
    double mingap = orderedpath[0][4].real()-orderedpath[0][3].real(); int smin = 0;
    //cout<<mingap;
    for(int s=0;s<orderedpath.size();s++){
        //cout<<mingap;
        if((orderedpath[s][4].real()-orderedpath[s][3].real())<mingap){// && orderedpath[s][3].real() > orderedpath[smin][3].real()){
            mingap=orderedpath[s][4].real()-orderedpath[s][3].real();
            smin = s;
        }
    }
    return orderedpath[smin][3].real();
}

double chempot2(Eigen::VectorXcd eigvals){
    //This is half-filling. Hard coded for 8 bands. Allow for filling by adjusting valence and conduct.
    double eF = 0.0;
    eF = order(eigvals)[31].real();
    return eF;
}

double gap(map<int, Eigen::VectorXcd> orderedpath){
    //This is half-filling. Hard coded for 8 bands. Allow for filling by adjusting valence and conduct.
    double mingap = orderedpath[0][4].real()-orderedpath[0][3].real(); int smin = 0;
    //cout<<mingap;
    for(int s=0;s<orderedpath.size();s++){
        //cout<<mingap;
        if((orderedpath[s][4].real()-orderedpath[s][3].real())<mingap){// && orderedpath[s][3].real() > orderedpath[smin][3].real()){
            mingap=orderedpath[s][4].real()-orderedpath[s][3].real();
            smin = s;
        }
    }
    return mingap;
}

//double groundenergy(Eigen::VectorXcd eigvals){
//    return
//}
//double groundenergy(map<int, Eigen::VectorXcd> orderedpath){
//    double energy = 0.0;
//    for(int s=0;s<orderedpath.size();s++){
//        for(int i=0;i<4;i++){
//            energy += 0.001*orderedpath[s][i].real();
//        }
//    }
//    return energy;
//}