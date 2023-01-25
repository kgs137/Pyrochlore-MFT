//  main.cpp
//  Pyrochlore Iridate Hamiltonian Minimization
//
//  Created by Kyle G Sherman on 9/15.
//  Copyright © 2016 Kyle G Sherman. All rights reserved.

#include <iostream>
#include <fstream>
#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "dlib/optimization.h"
#include "dlib/statistics.h"
#include <cmath>
#include <map>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

complex<double> I(0.0,1.0);
typedef dlib::matrix<double,0,1> column_vector;

#include "struct.h"
#include "conversion.h"
#include "term.h"
#include "construct.h"
#include "ft.h"
#include "gradient.h"
#include "data_import.h"
#include "analysis.h"
//

class func
{
public:
    func(
         const Eigen::MatrixXcd& H0, column_vector initial_state, const int& N1, const int& N2, const int& N3, const double J, const double JK, const double U
         )
    {
        h0 = H0;
        n1 = N1;
        n2 = N2;
        n3 = N3;
        j = J;
        jk = JK;
        u = U;
    }

    double operator() (const column_vector& spinconfiguration) const
    {
        map<int, Eigen::Vector3d> ConfigurationPass_a;
        map<int, Eigen::Vector3d> ConfigurationPass_b;

        ConfigurationPass_a = flattospins(separate(spinconfiguration, n1,n2,n3, 0),n1, n2, n3,"Angular");
        ConfigurationPass_b = flattospins(separate(spinconfiguration, n1,n2,n3, 1),n1, n2, n3,"Cartesian");

        Eigen::MatrixXcd h = constructspins(h0, ConfigurationPass_a, ConfigurationPass_b, n1,n2,n3,j,jk,u);
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(h);
        Eigen::VectorXcd eigvals = es.eigenvalues();
        Eigen::VectorXcd orderedeigvals = order(eigvals);
        double ground_state_energy = 0.0;
        for(int i=0;i<0.5*orderedeigvals.size();i++){
            ground_state_energy=+orderedeigvals[i].real();
        }
        return ground_state_energy;
    }

private:
    Eigen::MatrixXcd h0; column_vector initial_state; int n1; int n2; int n3; double j; double jk; double u;
};

//
int main(){
    int N1,N2,N3;
    N1=N2=N3=2;

    //column_vector spin_state;
    //spin_state = spinstoflat(iniSpins(N1,N2,N3),N1,N2,N3);

    //Random seed
    int a; for (int seed=1000;seed<15936;seed++){a = rand();}

    double delta=0.0000001;
    map<int,double> energies;

    for (int seed=1000;seed<15936;seed++){
    cout<<"seed: "+to_string(seed)<<endl;
    for (double j = 0.0002; j <= 0.0002; j = j + 0.0002){
        cout << "j: " + to_string(j) << endl;
        for (double jk = 0.004; jk <= 0.004; jk = jk + 0.004){
          cout<<"jk: "+to_string(jk)<<endl;
          for (double u = 0.0; u <= 2.0; u = u + 1.0){
            cout<<"u: "+to_string(u)<<endl;
            for (double lambda = 0.15; lambda >= 0.15; lambda = lambda - 0.1){
              for (double t = -0.45; t <= -0.45; t = t + 0.1){
            column_vector spin_state_a; column_vector spin_state_b;
            column_vector lower_bound_a; column_vector lower_bound_b;
            column_vector upper_bound_a; column_vector upper_bound_b;
            spin_state_a = spinstoflat(iniSpins(N1,N2,N3),N1,N2,N3,"Angular");
            spin_state_b = spinstoflat(iniSpins(N1,N2,N3),N1,N2,N3,"Cartesian");
            lower_bound_a = spinstoflat(iniSpins(N1,N2,N3),N1,N2,N3,"Angular");
            lower_bound_b = spinstoflat(iniSpins(N1,N2,N3),N1,N2,N3,"Cartesian");
            upper_bound_a = spinstoflat(iniSpins(N1,N2,N3),N1,N2,N3,"Angular");
            upper_bound_b = spinstoflat(iniSpins(N1,N2,N3),N1,N2,N3,"Cartesian");

            for (int i=0;i<spin_state_a.size();i++){
              spin_state_a(i) += 3.141592*((1-fmod(i,2))+2*fmod(i,2))*(0.01*fmod(rand(),100));
              lower_bound_a(i) += 0;
              upper_bound_a(i) += 3.141592*((1-fmod(i,2))+2*fmod(i,2));
            }
            for (int i=0;i<spin_state_b.size();i++){
              spin_state_b(i) += pow(12,-0.5)*(0.01*(fmod(rand(),200)+1-100));
              lower_bound_b(i) += -0.5;
              upper_bound_b(i) += 0.5;
            }
            column_vector spin_state_ab, lower_bound_ab, upper_bound_ab;
            spin_state_ab = combine(spin_state_a, spin_state_b, N1, N2, N3);
            lower_bound_ab = combine(lower_bound_a, lower_bound_b, N1, N2, N3);
            upper_bound_ab = combine(upper_bound_a, upper_bound_b, N1, N2, N3);
            string folder = "new_ones";
            // Check EnergyandParameters is empty or openable(exists), if so then run. If not, skip.
            ofstream h_real ("phase diagram/"+folder+"/h_real_"+to_string(N1)+to_string(N2)+to_string(N3)+"_"+to_string(t)+"_"+to_string(lambda)+"_"+to_string(j)+"_"+to_string(jk)+"_"+to_string(u)+"_"+to_string(seed)+".csv");
            ofstream h_imag ("phase diagram/"+folder+"/h_imag_"+to_string(N1)+to_string(N2)+to_string(N3)+"_"+to_string(t)+"_"+to_string(lambda)+"_"+to_string(j)+"_"+to_string(jk)+"_"+to_string(u)+"_"+to_string(seed)+".csv");
            ofstream h_real_ini ("phase diagram/"+folder+"/h_real_ini_"+to_string(N1)+to_string(N2)+to_string(N3)+"_"+to_string(t)+"_"+to_string(lambda)+"_"+to_string(j)+"_"+to_string(jk)+"_"+to_string(u)+"_"+to_string(seed)+".csv");
            ofstream h_imag_ini ("phase diagram/"+folder+"/h_imag_ini_"+to_string(N1)+to_string(N2)+to_string(N3)+"_"+to_string(t)+"_"+to_string(lambda)+"_"+to_string(j)+"_"+to_string(jk)+"_"+to_string(u)+"_"+to_string(seed)+".csv");
            ofstream eig ("phase diagram/"+folder+"/eig_"+to_string(N1)+to_string(N2)+to_string(N3)+"_"+to_string(t)+"_"+to_string(lambda)+"_"+to_string(j)+"_"+to_string(jk)+"_"+to_string(u)+"_"+to_string(seed)+".csv");
            ofstream eig_ini ("phase diagram/"+folder+"/eig_ini_"+to_string(N1)+to_string(N2)+to_string(N3)+"_"+to_string(t)+"_"+to_string(lambda)+"_"+to_string(j)+"_"+to_string(jk)+"_"+to_string(u)+"_"+to_string(seed)+".csv");
            ofstream initial_state_a ("phase diagram/"+folder+"/initial_state_a_"+to_string(N1)+to_string(N2)+to_string(N3)+"_"+to_string(t)+"_"+to_string(lambda)+"_"+to_string(j)+"_"+to_string(jk)+"_"+to_string(u)+"_"+to_string(seed)+".csv");
            ofstream initial_state_b ("phase diagram/"+folder+"/initial_state_b_"+to_string(N1)+to_string(N2)+to_string(N3)+"_"+to_string(t)+"_"+to_string(lambda)+"_"+to_string(j)+"_"+to_string(jk)+"_"+to_string(u)+"_"+to_string(seed)+".csv");
            ofstream final_state_a ("phase diagram/"+folder+"/final_state_a_"+to_string(N1)+to_string(N2)+to_string(N3)+"_"+to_string(t)+"_"+to_string(lambda)+"_"+to_string(j)+"_"+to_string(jk)+"_"+to_string(u)+"_"+to_string(seed)+".csv");
            ofstream final_state_b ("phase diagram/"+folder+"/final_state_b_"+to_string(N1)+to_string(N2)+to_string(N3)+"_"+to_string(t)+"_"+to_string(lambda)+"_"+to_string(j)+"_"+to_string(jk)+"_"+to_string(u)+"_"+to_string(seed)+".csv");
            ofstream EnergyandParameters ("phase diagram/"+folder+"/EnergyandParameters_"+to_string(N1)+to_string(N2)+to_string(N3)+"_"+to_string(t)+"_"+to_string(lambda)+"_"+to_string(j)+"_"+to_string(jk)+"_"+to_string(u)+"_"+to_string(seed)+".csv");
            ofstream bands ("phase diagram/"+folder+"/bands_"+to_string(N1)+to_string(N2)+to_string(N3)+"_"+to_string(t)+"_"+to_string(lambda)+"_"+to_string(j)+"_"+to_string(jk)+"_"+to_string(u)+"_"+to_string(seed)+".csv");
            ofstream bands_ini ("phase diagram/"+folder+"/bands_ini_"+to_string(N1)+to_string(N2)+to_string(N3)+"_"+to_string(t)+"_"+to_string(lambda)+"_"+to_string(j)+"_"+to_string(jk)+"_"+to_string(u)+"_"+to_string(seed)+".csv");
  //###########################################################################This is next
            map<int, Eigen::Vector3d> spin_vectors_a = flattospins(spin_state_a,N1,N2,N3,"Angular");
            map<int, Eigen::Vector3d> spin_vectors_b = flattospins(spin_state_b,N1,N2,N3,"Cartesian");
            for (int m = 0; m < 4*N1*N2*N3; m++){
              for (int n = 0; n < 3; n++){
                initial_state_a << spin_vectors_a[m][n];
                initial_state_a << ",";
                initial_state_b << spin_vectors_b[m][n];
                initial_state_b << ",";
              }
              initial_state_a << endl;
              initial_state_b << endl;
            }
              initial_state_a.close();
              initial_state_b.close();

              Eigen::MatrixXcd h_ini(2*4*N1*N2*N3,2*4*N1*N2*N3);
              h_ini = constructall(flattospins(spin_state_a,N1,N2,N3,"Angular"),flattospins(spin_state_b,N1,N2,N3,"Cartesian"),N1,N2,N3, t, lambda, j, jk, u);

              for (int m = 0; m < 2*4*N1*N2*N3; m++){
                for (int n = 0; n < 2*4*N1*N2*N3; n++){
                  h_real_ini << real(h_ini(m,n));
                  h_real_ini << ",";
                  h_imag_ini << imag(h_ini(m,n));
                  h_imag_ini << ",";
                }
                h_real_ini << endl;
                h_imag_ini << endl;
              }
              h_real_ini.close();
              h_imag_ini.close();

              Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es_ini(h_ini);
              Eigen::VectorXcd eigvals_ini = es_ini.eigenvalues();
              eigvals_ini = order(eigvals_ini);

              for (int e = 0; e < eigvals_ini.size(); e++){
                eig_ini<<eigvals_ini[e]<<", ";
                eig_ini << endl;
              }
              eig_ini.close();

              double initial_energy = 0.0;
              for(int i=0;i<0.5*eigvals_ini.size();i++){
                  initial_energy=+eigvals_ini[i].real();
              }

              map<int,Eigen::VectorXcd> bandstructure_ini = DFTbands(h_ini,N1,N2,N3);

              for (int i = 0; i < 8; i++){
                for (int s=0;s<1000;s++){
                    bands_ini << real(bandstructure_ini[s][i]);
                    bands_ini << ",";
                }
                bands_ini << endl;
              }
              bands_ini.close();

            energies.insert(make_pair(0,find_min_box_constrained(dlib::cg_search_strategy(),dlib::objective_delta_stop_strategy(delta),func(constructraw(N1,N2,N3,t,lambda),spin_state_ab,N1,N2,N3,j,jk,u),dlib::derivative(func(constructraw(N1,N2,N3,t,lambda),spin_state_ab,N1,N2,N3,j,jk,u)), spin_state_ab, lower_bound_ab, upper_bound_ab)));

            spin_state_a = separate(spin_state_ab, N1,N2,N3,0);
            spin_state_b = separate(spin_state_ab, N1,N2,N3,1);

            Eigen::MatrixXcd h(2*4*N1*N2*N3,2*4*N1*N2*N3);
            h = constructall(flattospins(spin_state_a,N1,N2,N3,"Angular"),flattospins(spin_state_b,N1,N2,N3,"Cartesian"),N1,N2,N3, t, lambda, j, jk, u);

            spin_vectors_a = flattospins(spin_state_a,N1,N2,N3,"Angular");
            spin_vectors_b = flattospins(spin_state_b,N1,N2,N3,"Cartesian");

            for (int m = 0; m < 4*N1*N2*N3; m++){
              for (int n = 0; n < 3; n++){
                final_state_a << spin_vectors_a[m][n];
                final_state_a << ",";
                final_state_b << spin_vectors_b[m][n];
                final_state_b << ",";
              }
              for (int n = 0; n < 2*4*N1*N2*N3; n++){
                h_real << real(h(m,n));
                h_real << ",";
                h_imag << imag(h(m,n));
                h_imag << ",";
              }
              h_real << endl;
              h_imag << endl;
              final_state_a << endl;
              final_state_b << endl;
            }

            h_real.close();
            h_imag.close();
            final_state_a.close();
            final_state_b.close();
            //cout << h;
            Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(h);
            Eigen::VectorXcd eigvals = es.eigenvalues();
            eigvals = order(eigvals);

            for (int e = 0; e < eigvals.size(); e++){
              eig<<eigvals[e]<<", ";
              eig << endl;
            }
            eig.close();

            double final_energy = 0.0;
            for(int i=0;i<0.5*eigvals.size();i++){
                final_energy=+eigvals[i].real();
            }

            map<int,Eigen::VectorXcd> bandstructure = DFTbands(h,N1,N2,N3);

            for (int i = 0; i < 8; i++){
              for (int s=0;s<1000;s++){
                  bands << real(bandstructure[s][i]);
                  bands << ",";
              }
              bands << endl;
            }
            bands.close();
            EnergyandParameters<<"t = "<<to_string(t)<<","<<endl<<"lambda = "<<to_string(lambda)<<","<<endl<<" U = "<<to_string(u)<<","<<endl<<" Jk = "<<to_string(jk)<<","<<endl<<" J = "<<to_string(j)<<","<<endl<<"Initial_Energy = "<<initial_energy<<","<<endl<<"Final_Energy = "<<final_energy<<","<<endl;
            EnergyandParameters.close();
                  }
                }
              }
            }
          }
        }

    return 0;
}
