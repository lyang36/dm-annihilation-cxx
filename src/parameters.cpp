#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cmath>
#include <rpc/types.h>
#include <rpc/xdr.h>
#include "datatypes.h"
#include "parameters.h"

using namespace std;


string getstring_double( double * d, long numx);
string getstring_int( int * d, int numx);

void print_out_natconst( Natconst * natconst);
void print_out_halo( Halo * halo);
void print_out_dm( Dm * dm);

void print_out_codeunits( Codeunits * codeunits);

void print_out_params( Params * params);

void print_out_master( Parameters * master);

void print_out_map( Map * map);


void Parameters::setupUnits(){
    //Natconst
    natconst.h100 = info.h100;
	natconst.rho_crit_in_cgs = 1.8788309e-29 * pow(info.h100, 2);
	natconst.Omega_m = info.Omega_M;
	natconst.Rvir_MW_in_Mpc = 0.275;
    
    // definition of the virial overdensity (Bryan & Norman 1998)  (rho/rho_0)
  	double Omega_m_at_z = natconst.Omega_m * pow((1.0 + redshift), 3) / (natconst.Omega_m * pow((1.0 + redshift), 3) + (1.0 - natconst.Omega_m));
  	double x = Omega_m_at_z - 1.0;
  	natconst.Delta_vir = (18.0 * PI * PI + 82.0 * x - 39.0 * x * x) / (1.0 + x);
  	natconst.rho_crit_in_cgs = 3.0 * pow((100.0 * natconst.h100/(units.Mpc_in_cgs * 1e-5)), 2) / (8 * PI * natconst.G_in_cgs);

    
    //setup code units
    codeunits.length_to_Mpc = info.Lbox_in_Mpc;
    codeunits.length_to_cgs = info.Lbox_in_Mpc*units.Mpc_in_cgs;
    codeunits.time_to_cgs = pow((natconst.rho_crit_in_cgs * natconst.G_in_cgs), (-0.5));
    codeunits.density_to_cgs = natconst.rho_crit_in_cgs;
    
  	codeunits.mass_to_cgs = codeunits.density_to_cgs * pow(codeunits.length_to_cgs, 3);
  	codeunits.mass_to_Msun = codeunits.mass_to_cgs / units.Msun_in_cgs;
  	codeunits.velocity_to_cgs = codeunits.length_to_cgs /codeunits.time_to_cgs;
	codeunits.annihilation_flux_to_cgs = (codeunits.density_to_cgs * codeunits.mass_to_cgs / pow(codeunits.length_to_cgs, 2));

}

void Parameters::setupHalo(){
    //setup halo
    double shape_radii = 0.0;
	double halo_shape =0.0;
    
    halo.Mvir_in_Msun = info.Mvir_in_Msun;
    halo.M200_in_Msun = info.M200_in_Msun;
    halo.M200crit_in_Msun = info.M200crit_in_Msun;
    halo.Rvir_in_Mpc = info.Rvir_in_Mpc;
    halo.R200_in_Mpc = info.R200_in_Mpc;
    halo.R200crit_in_Mpc = info.R200crit_in_Mpc;
    halo.Vmax_in_kms = info.Vmax_in_kms;
    halo.RVmax_in_kpc = info.RVmax_in_kpc;
    halo.rconverg_in_kpc = info.rconverg_in_kpc;
    
    for(int i = 0; i < 3; i++){
        halo.params_NFW[i] = info.params_NFW[i];
        halo.params_GNFW[i] = info.params_GNFW[i];
        halo.params_Einasto[i] = info.params_Einasto[i];
    }
    
    halo.shape_r = shape_radii;
    halo.shape = halo_shape;
    
}


void Parameters::setupParams(){
    
    // setting for parameters
	params.z = redshift;
    
    for(int i = 0; i < 3; i ++){
        params.cpos[i] = info.centerpos[i];
        params.cvel[i] = info.centervel[i];
    }
	params.otheta = 0.0;
	params.ophi = 0.0;
	params.Lbox_in_Mpc = info.Lbox_in_Mpc;
	for( int i = 0, j = 0; i < 10; i++){
		params.particle_numbers[i] = 0;
		params.particle_masses[i] = 0;
		params.particle_masses_in_Msun[i] = 0;
		if (info.particle_masses[i] > 0){
		    params.particle_numbers[j] = info.particle_numbers[i];
    		params.particle_masses[j] = info.particle_masses[i];
    		params.particle_masses_in_Msun[j] = info.particle_masses[i] * codeunits.mass_to_cgs / units.Msun_in_cgs;
			j++;
		}
	}
}

void Parameters::setupOthers(){
    //dark matter properties
    dm.M_in_GeV = 46.0;
    dm.M_in_cgs = 0.0;
    dm.sigma_v_in_cgs = 5.0e-26;
    dm.M_in_cgs = dm.M_in_GeV * 1e9 * units.eV_in_cgs / pow(natconst.c_in_cgs, 2);
}

void Parameters::setupRotation(){
    
}

void Parameters::initialize(){
    //nSide = 512;
    map.setNside(512);
    memParts = 10000000;
    testNum = -1;
}

Parameters::Parameters(string ifname){
    initialize();
    
    const int MAXLENGTH = 300;
    redshift = 0;
    for( int i = 0; i<10; i++){
		info.particle_masses[i] = 0;
		info.particle_numbers[i] = 0;
	}
    
    string line;
    ifstream infofile(ifname.c_str(), ios_base::in);
    
    if (infofile.is_open()){
		while (infofile.good()){
			getline (infofile, line);
			char linestr[MAXLENGTH];
			strcpy (linestr, line.c_str());
			if (line.length() == 0)
				continue;
			char * toks = strtok(linestr, " ");
			if (toks == NULL)
				continue;
			string temp = toks;
			if (toks[0] == '#')
				continue;
			char * val = strtok(NULL, " ");
			val = strtok(NULL, " ");
			if (temp == "Omega_M"){
				info.Omega_M = atof(val);
			}else if (temp == "h100"){
				info.h100 = atof(val);
			}else if (temp == "n_s"){
                info.n_s = atof(val);
            }else if (temp == "sigma_8"){
                info.sigma_8 = atof(val);
            }else if (temp == "Lbox_in_Mpc"){
                info.Lbox_in_Mpc = atof(val);
            }else if (temp == "particle_masses"){
				int i = 0;
				do{
					info.particle_masses[i] = atof(val);
					i++;
					val = strtok(NULL, " ");
				}while (val != NULL && i<10);
            }else if (temp == "particle_numbers"){
                int i = 0;
			    do{
                 	info.particle_numbers[i] = atol(val);
					i++;
					val = strtok(NULL, " ");
				}while (val != NULL && i<10);
            }else if (temp == "centerpos"){
                info.centerpos[0] = atof(val);
				val = strtok(NULL, " ");
				info.centerpos[1] = atof(val);
				val = strtok(NULL, " ");
				info.centerpos[2] = atof(val);
			}else if (temp == "centervel"){
                info.centervel[0] = atof(val);
				val = strtok(NULL, " ");
				info.centervel[1] = atof(val);
				val = strtok(NULL, " ");
				info.centervel[2] = atof(val);
			}else if (temp == "Mvir_in_Msun"){
				info.Mvir_in_Msun = atof(val);
			}else if (temp == "M200_in_Msun"){
                info.M200_in_Msun = atof(val);
			}else if (temp == "M200crit_in_Msun"){
                info.M200crit_in_Msun = atof(val);
			}else if (temp == "Rvir_in_Mpc"){
                info.Rvir_in_Mpc = atof(val);
			}else if (temp == "R200_in_Mpc"){
                info.R200_in_Mpc = atof(val);
			}else if (temp == "R200crit_in_Mpc"){
                info.R200crit_in_Mpc = atof(val);
			}else if (temp == "Vmax_in_kms"){
                info.Vmax_in_kms = atof(val);
			}else if (temp == "RVmax_in_kpc"){
                info.RVmax_in_kpc = atof(val);
			}else if (temp == "rconverg_in_kpc"){
                info.rconverg_in_kpc = atof(val);
			}else if (temp == "params_NFW"){
                info.params_NFW[0] = atof(val);
                val = strtok(NULL, " ");
                info.params_NFW[1] = atof(val);
                val = strtok(NULL, " ");
                info.params_NFW[2] = atof(val);
			}else if (temp == "params_GNFW"){
                info.params_GNFW[0] = atof(val);
                val = strtok(NULL, " ");
                info.params_GNFW[1] = atof(val);
                val = strtok(NULL, " ");
                info.params_GNFW[2] = atof(val);
			}else if (temp == "params_Einasto"){
                info.params_Einasto[0] = atof(val);
                val = strtok(NULL, " ");
                info.params_Einasto[1] = atof(val);
                val = strtok(NULL, " ");
                info.params_Einasto[2] = atof(val);
			}else if (temp == "align_vec"){
                align_vector[0] = atof(val);
                val = strtok(NULL, " ");
                align_vector[1] = atof(val);
                val = strtok(NULL, " ");
                align_vector[2] = atof(val);
			}else if (temp == "obs_pos"){
                obs_pos[0] = atof(val);
                val = strtok(NULL, " ");
                obs_pos[1] = atof(val);
                val = strtok(NULL, " ");
                obs_pos[2] = atof(val);
			}else if (temp == "nside"){
                int nSide = atoi(val);
                map.setNside(nSide);
			}else if (temp == "test"){
                testNum = atoi(val);
			}
            
		}
	}else{
        printf("Info file incorrect!!!\n");
        exit(1);
    }
    infofile.close();
    
    setupUnits();
    setupOthers();
    setupHalo();
    setupParams();
    setupRotation();
    
}

void Parameters::printParameters(){
    print_out_master(this);
}


string getstring_double( double * d, long numx){
    //stringstream ss(stringstream::in | stringstream::out);
    string oputs;
    string temps;
    for(long i = 0; i< numx; i++){
        stringstream ss;
        ss << d[i];
        //cout << d[i] << endl;
        ss >> temps;
        oputs += temps + " ";
    }
    //ss >> oputs;
    return oputs;
}

string getstring_int( int * d, int numx){
    //stringstream ss(stringstream::in | stringstream::out);
    string oputs;
    string temps;
    for(int i = 0; i< numx; i++){
        stringstream ss;
        ss << d[i];
        //cout << d[i] << endl;
        ss >> temps;
        oputs += temps + " ";
    }
    //ss >> oputs;
    return oputs;
}


void print_out_natconst( Natconst * natconst){
	cout << "Out put Natconst:-----------------------------------------------------" << endl;
    cout << "h100 " << natconst -> h100 << endl;
    cout << "G_in_cgs " << natconst -> G_in_cgs << endl;
    cout << "rho_crit_in_cgs " << natconst -> rho_crit_in_cgs << endl;
    cout << "Omega_m " << natconst -> Omega_m << endl;
    cout << "Delta_vir " << natconst -> Delta_vir << endl;
    cout << "Rvir_MW_in_Mpc " << natconst -> Rvir_MW_in_Mpc << endl;
    cout << "c_in_cgs " << natconst -> c_in_cgs << endl;
}

void print_out_halo( Halo * halo){
	cout << "Out put Halo:-----------------------------------------------------" << endl;
	cout << "Mvir_in_Msun " << halo -> Mvir_in_Msun << endl;
	cout << "M200_in_Msun " << halo -> M200_in_Msun << endl;
	cout << "M200crit_in_Msun " << halo -> M200crit_in_Msun << endl;
	cout << "Rvir_in_Mpc " << halo -> Rvir_in_Mpc << endl;
	cout << "R200_in_Mpc " << halo -> R200_in_Mpc << endl;
	cout << "R200crit_in_Mpc " << halo -> R200crit_in_Mpc << endl;
	cout << "Vmax_in_kms " << halo -> Vmax_in_kms << endl;
	cout << "RVmax_in_kpc " << halo -> RVmax_in_kpc << endl;
	cout << "rconverg_in_kpc " << halo -> rconverg_in_kpc << endl;
	cout << "params_NFW " << getstring_double( halo -> params_NFW, 3) << endl;
	cout << "params_GNFW " << getstring_double(  halo -> params_GNFW, 3) << endl;
	cout << "params_Einasto " << getstring_double( halo -> params_Einasto, 3) << endl;
	cout << "shape_r " << halo -> shape_r << endl;
	cout << "shape " << halo -> shape << endl;
}

void print_out_dm( Dm * dm){
	cout << "DM:-----------------------------------------------------" << endl;
	cout << "M_in_GeV " << dm -> M_in_GeV << endl;// = 46.0;
	cout << "M_in_cgs " << dm -> M_in_cgs << endl;// = 0.0;
	cout << "sigma_v_in_cgs "<< dm -> sigma_v_in_cgs << endl;// = 5.0e-26;
    
}

void print_out_codeunits( Codeunits * codeunits){
	cout << "Codeunits:-----------------------------------------------------" << endl;
	cout << "mass_to_cgs " <<  codeunits -> mass_to_cgs << endl;
	cout << "mass_to_Msun " <<  codeunits -> mass_to_Msun << endl;
	cout << "length_to_Mpc " << codeunits -> length_to_Mpc << endl;
	cout << "length_to_cgs " << codeunits -> length_to_cgs << endl;
	cout << "time_to_cgs " << codeunits -> time_to_cgs << endl;
	cout << "velocity_to_cgs " << codeunits ->  velocity_to_cgs << endl;
	cout << "density_to_cgs " << codeunits -> density_to_cgs << endl;
	cout << "annihilation_flus_to_cgs " << codeunits -> annihilation_flux_to_cgs << endl;
}


void print_out_params( Params * params){
	cout << "Params:-----------------------------------------------------" << endl;
	cout << "z " <<  params -> z << endl;
	cout << "cpos " << getstring_double(params -> cpos, 3) << endl;
	cout << "cvel " << getstring_double(params -> cvel, 3) << endl;//cvel[3];
	cout << "opos " << getstring_double(params -> opos, 3) << endl; // opos[3];
	cout << "otheta " << params -> otheta << endl;
	cout << "ophi " << params -> ophi << endl;//;
	cout << "Lbox_in_Mpc " << params -> Lbox_in_Mpc << endl;//;
	
	cout << "particle_numbers " <<  getstring_int( (params -> particle_numbers), 10) << endl;//[10];
	cout << "particle_masses " <<  getstring_double( (params -> particle_masses), 10) << endl;//[10];
	cout << "particle_masses_in_Msun " << getstring_double( (params -> particle_masses_in_Msun), 10) << endl;// [10];
}


void print_out_master( Parameters * master){
	cout << "Parameters:-----------------------------------------------------" << endl;
	print_out_natconst( &(master -> natconst) );
	print_out_dm( &(master -> dm) );
	print_out_halo( &(master -> halo) );
	print_out_params( &(master -> params) );
	print_out_codeunits( &(master -> codeunits));
    print_out_map( &(master -> map));

	cout << "rotmatrix " << getstring_double( master -> rotmatrix[0], 3) << endl;
	cout << "rotmatrix " << getstring_double( master -> rotmatrix[1], 3) << endl;
	cout << "rotmatrix " << getstring_double( master -> rotmatrix[2], 3) << endl;
}



void print_out_map( Map * map){
	cout << "map:-----------------------------------------------------" << endl;
	cout << "projection " << map -> projection << endl;//;
	cout << "Nside " << map -> Nside << endl;//;
	cout << "Npix " << map -> Npix << endl;//;
	cout << "dOmega " << map -> dOmega << endl;//;
	cout << "theta0 " << map -> theta0 << endl;//;
}
