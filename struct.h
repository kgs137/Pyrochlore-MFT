//
struct ALatticeSite{
    int i;
    Eigen::Vector4d nmu;
    Eigen::Vector3d r;
    Eigen::Vector3d spin;
};

//
struct BLatticeSite{
    int i;
    Eigen::Vector4d nmu;
    Eigen::Vector3d r;
    Eigen::Vector3d spin;
};

//
struct TermElement{
    pair<int,int> Bond;
    Eigen::Matrix2cd m;
};