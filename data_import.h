map<int,Eigen::VectorXcd> import_band_structure(string dir){
    map<int,Eigen::VectorXcd> band_structure;
    
    Eigen::VectorXcd blank_vector(8);
    blank_vector<<0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0;
    for (int s = 0; s<1000;s++){band_structure.insert(make_pair(s,blank_vector));}
    
	string line0 = "";
    string line1 = "";
    string line2 = "";
    string line3 = "";
    string line4 = "";
    string line5 = "";
    string line6 = "";
    string line7 = "";
    
    int s = 0;
    
	int n = 0;
    
    double a = 0.0;
    double b = 0.0;
    
    double energy;
    
    //int count = 0;
	
	ifstream data0 (dir + "band0.txt");
    ifstream data1 (dir + "band1.txt");
    ifstream data2 (dir + "band2.txt");
    ifstream data3 (dir + "band3.txt");
    ifstream data4 (dir + "band4.txt");
    ifstream data5 (dir + "band5.txt");
    ifstream data6 (dir + "band6.txt");
    ifstream data7 (dir + "band7.txt");
    
	if (data0.is_open()){
        while(getline(data0,line0)){
            n = line0.find_last_of(".");
            s = stof(line0.substr(0,n-1));
            energy = stof(line0.substr(n-2));
            band_structure[s][0] = energy;
        }
		data0.close();
	}
    if (data1.is_open()){
        while(getline(data1,line1)){
            n = line1.find_last_of(".");
            s = stof(line1.substr(0,n-1));
            energy = stof(line1.substr(n-2));
            band_structure[s][1] = energy;
        }
        data1.close();
    }
    if (data2.is_open()){
        while(getline(data2,line2)){
            n = line2.find_last_of(".");
            s = stof(line2.substr(0,n-1));
            energy = stof(line2.substr(n-2));
            band_structure[s][2] = energy;
        }
        data2.close();
    }
    if (data3.is_open()){
        while(getline(data3,line3)){
            n = line3.find_last_of(".");
            s = stof(line3.substr(0,n-1));
            energy = stof(line3.substr(n-2));
            band_structure[s][3] = energy;
        }
        data3.close();
    }
    if (data4.is_open()){
        while(getline(data4,line4)){
            n = line4.find_last_of(".");
            s = stof(line4.substr(0,n-1));
            energy = stof(line4.substr(n-2));
            band_structure[s][4] = energy;
        }
        data4.close();
    }
    if (data5.is_open()){
        while(getline(data5,line5)){
            n = line5.find_last_of(".");
            s = stof(line5.substr(0,n-1));
            energy = stof(line5.substr(n-2));
            band_structure[s][5] = energy;
        }
        data5.close();
    }
    if (data6.is_open()){
        while(getline(data6,line6)){
            n = line6.find_last_of(".");
            s = stof(line6.substr(0,n-1));
            energy = stof(line6.substr(n-2));
            band_structure[s][6] = energy;
        }
        data6.close();
    }
    if (data7.is_open()){
        while(getline(data7,line7)){
            n = line7.find_last_of(".");
            s = stof(line7.substr(0,n-1));
            energy = stof(line7.substr(n-2));
            band_structure[s][7] = energy;
        }
        data7.close();
    }
	
	return band_structure;
}

map<int,Eigen::VectorXcd> import_dft_band_structure(string dir){
    map<int,Eigen::VectorXcd> band_structure;
    
    Eigen::VectorXcd blank_vector(8);
    blank_vector<<0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0;
    for (int s = 0; s<1000;s++){band_structure.insert(make_pair(s,blank_vector));}
    
    string line = "";
    
    double e0 = 0.0;double e1 = 0.0;double e2 = 0.0;double e3 = 0.0;
    double e4 = 0.0;double e5 = 0.0;double e6 = 0.0;double e7 = 0.0;
    string rem = "";
    
    double ef = 5.2504;
    
    int n = 0;
    int s = 0;
    
    double a = 0.0;
    double b = 0.0;
    
    double energy;
    
    //int count = 0;
    
    ifstream data (dir);
    int m = 0;
    if (data.is_open()){
        while(getline(data,line) && m<=99){
            n = line.find_first_of(".");
            e0 = stod(line.substr(n-1,6));
            rem = line.substr(n+5);
            n = rem.find_first_of(".");
            e1 = stod(rem.substr(n-1,6));
            rem = rem.substr(n+5);
            n = rem.find_first_of(".");
            e2 = stod(rem.substr(n-1,6));
            rem = rem.substr(n+5);
            n = rem.find_first_of(".");
            e3 = stod(rem.substr(n-1,6));
            rem = rem.substr(n+5);
            n = rem.find_first_of(".");
            e4 = stod(rem.substr(n-1,6));
            rem = rem.substr(n+5);
            n = rem.find_first_of(".");
            e5 = stod(rem.substr(n-1,6));
            rem = rem.substr(n+5);
            n = rem.find_first_of(".");
            e6 = stod(rem.substr(n-1,6));
            rem = rem.substr(n+5);
            n = rem.find_first_of(".");
            e7 = stod(rem.substr(n-1,6));
            
            
            band_structure[s]<<e0-ef, e1-ef, e2-ef, e3-ef, e4-ef, e5-ef, e6-ef, e7-ef;
            s++;
            m++;
        }
        data.close();
    }
    return band_structure;
}


int compare_filed_energies(){
    
    string diagram_array [50][40];
    string dir = "/Users/Kyle/Desktop/Y-227/j-jk/";//Code/pyrochlore-iridates/PyIr/la227+f_diagram_q4/";
    int n = 0;
    int m = 0; string g = ""; double gap_aiao,gap_2i2o,gap_3i1o;
    gap_aiao=gap_2i2o=gap_3i1o=0.0; double gap = 0.0;
    string s = ""; string line = ""; string order = "";
    double energy_aiao = 0.0; double energy_2i2o = 0.0; double energy_3i1o = 0.0;
    
    ofstream diagram (dir + "diagram.txt");
    ofstream diagram_gap (dir + "diagram_gap.txt");
    
    for (double j = 0.0; j <= 0.4; j = j + 0.004){
        for (double jk = 0.0;jk<=0.4;jk = jk + 0.004){//(double jk = 0.0; jk <= 0.1; jk = jk + 0.05){
            
            ifstream file_aiao (dir + "aiao/J"+to_string(j)+"Jk"+to_string(jk)+ "EnergyandParameters.txt");// + "Jk" + to_string(jk)+ "EnergyandParameters.txt";
            ifstream file_2i2o (dir + "2i2o/J"+to_string(j)+"Jk"+to_string(jk)+ "EnergyandParameters.txt");// + "Jk" + to_string(jk)+ "EnergyandParameters.txt";
            ifstream file_3i1o (dir + "3i1o/J"+to_string(j)+"Jk"+to_string(jk)+ "EnergyandParameters.txt");// + "Jk" + to_string(jk)+ "EnergyandParameters.txt";
            
            if (file_aiao.is_open()){
                while(getline(file_aiao,line)){
                    n = line.find_last_of("Energy");
                    if (n!=-1){s = line.substr(5+4,7);energy_aiao = stof(s);}
                    n=-1;
                    m = line.find_last_of("Gap");
                    if (m!=-1){g = line.substr(6);}//gap_aiao = stof(g);}
                    m=-1;
                }
                file_aiao.close();
            }
            if (file_2i2o.is_open()){
                while(getline(file_2i2o,line)){
                    n = line.find_last_of("Energy");
                    if (n!=-1){s = line.substr(5+4,7);energy_2i2o = stof(s);}
                    n=-1;
                    m = line.find_last_of("Gap");
                    if (m!=-1){g = line.substr(6);}//gap_2i2o = stof(g);}
                    m=-1;
                }
                file_2i2o.close();
            }
            if (file_3i1o.is_open()){
                while(getline(file_3i1o,line)){
                    n = line.find_last_of("Energy");
                    if (n!=-1){s = line.substr(5+4,7);energy_3i1o = stof(s);cout<<line;}
                    n=-1;
                    m = line.find_last_of("Gap");
                    if (m!=-1){g = line.substr(6);}//gap_3i1o = stof(g);}
                    m=-1;
                }
                file_3i1o.close();
            }
            if(energy_aiao<energy_2i2o){
                if(energy_aiao<energy_3i1o){order="AiAo"; gap = gap_aiao;}
                if(energy_aiao>energy_3i1o){order="3i1o"; gap = gap_3i1o;}
            }
            if(energy_aiao>energy_2i2o){
                if(energy_2i2o<energy_3i1o){order="2i2o"; gap = gap_2i2o;}
                if(energy_2i2o>energy_3i1o){order="3i1o"; gap = gap_3i1o;}
            }
            //diagram_array[static_cast<int>(10*u+0.5)][static_cast<int>(100*jk+0.5)] = order;
            diagram<<order<<", ";
            diagram_gap<<g<<", ";
            cout<<j<<", "<<jk<<"->"<<energy_aiao<<" "<<energy_3i1o<<" "<<energy_2i2o<<endl;
            order = "";
            gap = 0.0;
        }
        diagram<<endl;
        diagram_gap<<endl;
    }
    return 0;
}