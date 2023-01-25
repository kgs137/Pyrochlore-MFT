//Band Analysis

//Band Comparison between data sets
//Associate Bands (Automated later)
//Path-Point-wise Comparison of distance
//Average of distances ofver each band
//Average of distances over all bands
//Output as double to be minimized via dlib
//Try to make this self-contained
//Input two band structures output optimal parameters

//Likeness comparison of two 8-band energy-ordered band structures
double band_comparison(map<int,Eigen::VectorXcd> bandstructure1, map<int,Eigen::VectorXcd> bandstructure2){
    //Assume bandstructure1 to be longer than bandstructure2 for now
    int pathlength1 = bandstructure1.size(); int pathlength2 = bandstructure2.size();

    double sqdif = 0.0; //Square of energetic distance between bands at particular k-point

    double sum0,sum1,sum2,sum3,sum4,sum5,sum6,sum7;  //Sum of sqdif for each band
    sum0=sum1=sum2=sum3=sum4=sum5=sum6=sum7=0.0;

    double allsum = 0.0;  //Total of band sums

    for(int s=0; s<pathlength2; s++){
        //Could accomodate arbitrary number of bands here
        sqdif = pow((bandstructure2[s][0]-bandstructure1[(round(pathlength1/pathlength2-0.5))*s][0]).real(),2.0);
        if(bandstructure2[s][0]!=0.0){sum0 += pow(2.71828,-17.3287*pow((bandstructure2[s][0]).real(),2.0))*sqdif;}
        sqdif = pow((bandstructure2[s][1]-bandstructure1[(round(pathlength1/pathlength2-0.5))*s][1]).real(),2.0);
        if(bandstructure2[s][1]!=0.0){sum1 += pow(2.71828,-17.3287*pow((bandstructure2[s][1]).real(),2.0))*sqdif;}
        sqdif = pow((bandstructure2[s][2]-bandstructure1[(round(pathlength1/pathlength2-0.5))*s][2]).real(),2.0);
        if(bandstructure2[s][2]!=0.0){sum2 += pow(2.71828,-17.3287*pow((bandstructure2[s][2]).real(),2.0))*sqdif;}
        sqdif = pow((bandstructure2[s][3]-bandstructure1[(round(pathlength1/pathlength2-0.5))*s][3]).real(),2.0);
        if(bandstructure2[s][3]!=0.0){sum3 += pow(2.71828,-17.3287*pow((bandstructure2[s][3]).real(),2.0))*sqdif;}
        sqdif = pow((bandstructure2[s][4]-bandstructure1[(round(pathlength1/pathlength2-0.5))*s][4]).real(),2.0);
        if(bandstructure2[s][4]!=0.0){sum4 += pow(2.71828,-17.3287*pow((bandstructure2[s][4]).real(),2.0))*sqdif;}
        sqdif = pow((bandstructure2[s][5]-bandstructure1[(round(pathlength1/pathlength2-0.5))*s][5]).real(),2.0);
        if(bandstructure2[s][5]!=0.0){sum5 += pow(2.71828,-17.3287*pow((bandstructure2[s][5]).real(),2.0))*sqdif;}
        sqdif = pow((bandstructure2[s][6]-bandstructure1[(round(pathlength1/pathlength2-0.5))*s][6]).real(),2.0);
        if(bandstructure2[s][6]!=0.0){sum6 += pow(2.71828,-17.3287*pow((bandstructure2[s][6]).real(),2.0))*sqdif;}
        sqdif = pow((bandstructure2[s][7]-bandstructure1[(round(pathlength1/pathlength2-0.5))*s][7]).real(),2.0);
        if(bandstructure2[s][7]!=0.0){sum7 += pow(2.71828,-17.3287*pow((bandstructure2[s][7]).real(),2.0))*sqdif;}
    }
    allsum = sum0+sum1+sum2+sum3+sum4+sum5+sum6+sum7;
    //cout<<allsum/(8*pathlength2);

    return allsum/(8*pathlength2);
}

class par_func
{
public:
    par_func(
             column_vector starting_point, const int& N1, const int& N2, const int& N3, const map<int,Eigen::VectorXcd>& bandstructure_data
         )
    {
        initial_state = starting_point;
        data_set = bandstructure_data;
        n1 = N1;
        n2 = N2;
        n3 = N3;
    }

    double operator() (const column_vector& par) const
    {
        double allsum;

        //Exchange par for constants here
        Eigen::MatrixXcd h = constructall(flattospins(initial_state,n1,n2,n3,"Cartesian"),flattospins(initial_state,n1,n2,n3,"Cartesian"),n1,n2,n3,par(0),par(1),0.0,0.0,par(2));
        map<int,Eigen::VectorXcd> bandstructure = DFTbands(h,n1,n2,n3);

        allsum = band_comparison(bandstructure,data_set);

        return allsum;
    }

private:
    Eigen::MatrixXcd h0; column_vector initial_state; int n1; int n2; int n3; map<int,Eigen::VectorXcd> data_set;
};


//Parameters are to be optimized such that 'band_comparison' is minimized
//Input hamiltonian as functor of variables, initial parameters, data set to compare against
//Output Optimal parameters and minimum 'band_comparison'

//Band Fitting
column_vector band_fitting(column_vector inipar, map<int,Eigen::VectorXcd> bandstructure_data){
    column_vector par = inipar; int N1,N2,N3; N1=N2=N3=2;
    column_vector starting_point;
    starting_point = spinstoflat(iniSpins(N1,N2,N3),N1,N2,N3);
    double delta = 10e-14;

    find_min_box_constrained(dlib::cg_search_strategy(),dlib::objective_delta_stop_strategy(delta),par_func(starting_point, N1, N2, N3, bandstructure_data),dlib::derivative(par_func(starting_point, N1, N2, N3, bandstructure_data)), par, -10.0, 10.0);

    return par;
}
