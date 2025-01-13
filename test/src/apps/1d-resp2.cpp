/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sstream>
#include "deom.hpp"

int main () {

    ifstream jsonFile("input.json");
    stringstream strStream;
    strStream << jsonFile.rdbuf();
    string jsonStr = strStream.str();
    string err;

    const Json json = Json::parse(jsonStr,err);
    if (!err.empty()) {
        printf ("Error in parsing input file: %s\n", err.c_str());
        return 0;
    }

    syst s(json["deom"]["syst"]);
    bath b(json["deom"]["bath"]);
    hidx h(json["deom"]["hidx"]);

    const int    nt = json["spec"]["nt"].int_value();
    const double dt = json["spec"]["dt"].number_value();
    const double staticErr = json["spec"]["staticErr"].number_value();
    const double steadyErr = json["spec"]["steadyErr"].number_value();
    const int    nk = json["spec"]["nk"].int_value();
    const char sch_hei = json["spec"]["sch_hei"].string_value().c_str()[0];
	const int    xp = json["spec"]["xpflag"].int_value();
    const string sdipFile = json["spec"]["sdipFile"].string_value();
    const string pdipFile = json["spec"]["pdipFile"].string_value();
    const string bdip1File = json["spec"]["bdip1File"].string_value();
    const string bdip2File = json["spec"]["bdip2File"].string_value();
    const string bdipFile = json["spec"]["bdipFile"].string_value();
    const string sdip0File = json["spec"]["sdip0File"].string_value();
    const string pdip0File = json["spec"]["pdip0File"].string_value();
    const string bdip0File = json["spec"]["bdip0File"].string_value();
    mat  sdip;
    cube pdip;
    cx_vec  bdip1;
    cx_vec  bdip2;
    cx_vec  bdip;
    mat  sdip0;
    cube pdip0;
    cx_vec  bdip0;
    if (sdip.load(sdipFile, arma_ascii)) {
        sdip.print(sdipFile);
    } else {
        printf("Fail to load sdip");
    }

    if (pdip.load(pdipFile, arma_ascii)) {
        pdip.print(pdipFile);
    } else {
        printf("Fail to load pdip");
    }

    if (bdip1.load(bdip1File, arma_ascii)) {
        bdip1.print(bdip1File);
    } else {
        printf("Fail to load bdip1");
    }

    if (bdip2.load(bdip2File, arma_ascii)) {
        bdip2.print(bdip2File);
    } else {
        printf("Fail to load bdip2");
    }

    if (bdip.load (bdipFile, arma_ascii)) {
        bdip.print(bdipFile);
    } else {
        printf("Fail to load bdip!\n");
    }

    if (sdip0.load(sdip0File, arma_ascii)) {
        sdip0.print(sdip0File);
    } else {
        printf("Fail to load sdip0");
    }

    if (pdip0.load(pdip0File, arma_ascii)) {
        pdip0.print(pdip0File);
    } else {
        printf("Fail to load pdip0");
    }

    if (bdip0.load(bdip0File, arma_ascii)) {
        bdip0.print(bdip0File);
    } else {
        printf("Fail to load bdip0");
    }


    deom d1(s,b,h);

    cx_cube rho_t0 = zeros<cx_cube>(d1.nsys,d1.nsys,d1.nmax);

    const mat& exph= expmat(-real(d1.ham1)/d1.temperature);
    rho_t0.slice(0).set_real(exph/trace(exph));
    //d1.equilibrium (rho_t0,dt,staticErr,nk);
    d1.equilibrium (rho_t0,dt,staticErr,nk);
    //cx_cube rhoconv (cx_cube& rho0, deom d1, const double dt, const double err);
   // rho_t0 = rhoconv(rho_t0, d1, dt, steadyErr);
    
    cx_cube rho_t1 = zeros<cx_cube>(d1.nsys,d1.nsys,d1.nmax);
    printf("##############################################\n");
    d1.oprAct(rho_t1,sdip0,pdip0,bdip0,rho_t0,'l',xp);
    //rho_t2(0,0,0) = rho_t0(1,0,0) - rho_t0(0,1,0); rho_t2(0,1,0) = rho_t0(1,1,0) - rho_t0(0,0,0);
    //rho_t2(1,0,0) = rho_t0(0,0,0) - rho_t0(1,1,0); rho_t2(1,1,0) = rho_t0(0,1,0) - rho_t0(1,0,0);
    //printf("##############################################");
    //rho_t2.slice(0).print();
    //exit(1);


    cx_vec ft1 = zeros<cx_vec>(nt); 
    cx_vec ft2 = zeros<cx_vec>(nt); 
    
    if (sch_hei == 's') { // sch-picture

        FILE *log_t1 = fopen("f1qt.dat","w");
	    FILE *log_t2 = fopen("f2qt.dat","w");
        for (int it=0; it<nt; ++it) {
            double t = it*dt;

            ft1(it) = d1.Trace(sdip,pdip,bdip1,rho_t1);
            ft2(it) = d1.Trace(sdip,pdip,bdip2,rho_t1);

            printf ("In 1d-correlation: it=%d, nddo=%d, lddo=%d\n", it, d1.nddo, d1.lddo);
            fprintf(log_t1, "%16.6e%16.6e%16.6e\n", t/deom_fs2unit, real(ft1(it)), imag(ft1(it)));
		    fprintf(log_t2, "%16.6e%16.6e%16.6e\n", t/deom_fs2unit, real(ft2(it)), imag(ft2(it)));
            d1.rk4 (rho_t1,t,dt);
        }
        fclose(log_t1);
	    fclose(log_t2);

    }

    // 1D FFT
    ft1 = ft1-ft1(nt-1);
    ft1(0) *= 0.5;
    const double dw1 = 2.0*deom_pi/(nt*dt);
    const cx_vec& fw1 = ifft(ft1)*nt*dt;

    
    FILE *log_w1 = fopen("f1qw.dat","w");
    for (int iw=nt/2; iw<nt; ++iw) {
        double w = (iw-nt)*dw1/deom_cm2unit;
        fprintf(log_w1, "%16.6e%16.6e%16.6e\n", w, real(fw1(iw)), imag(fw1(iw)));
    }
    for (int iw=0; iw<nt/2; ++iw) {
        double w = iw*dw1/deom_cm2unit;
        fprintf(log_w1, "%16.6e%16.6e%16.6e\n", w, real(fw1(iw)), imag(fw1(iw)));
    }
    fclose(log_w1);

    //printf("fw(0) = %16.6e  +%16.6eI\n", real(fw(0)), imag(fw(0)));
    // 1D FFT
    ft2 = ft2-ft2(nt-1);
    ft2(0) *= 0.5;
    const double dw2 = 2.0*deom_pi/(nt*dt);
    const cx_vec& fw2 = ifft(ft2)*nt*dt;

    
    FILE *log_w2 = fopen("f2qw.dat","w");
    for (int iw=nt/2; iw<nt; ++iw) {
        double w = (iw-nt)*dw2/deom_cm2unit;
        fprintf(log_w2, "%16.6e%16.6e%16.6e\n", w, real(fw2(iw)), imag(fw2(iw)));
    }
    for (int iw=0; iw<nt/2; ++iw) {
        double w = iw*dw2/deom_cm2unit;
        fprintf(log_w2, "%16.6e%16.6e%16.6e\n", w, real(fw2(iw)), imag(fw2(iw)));
    }
    fclose(log_w2);

    //printf("fw(0) = %16.6e  +%16.6eI\n", real(fw(0)), imag(fw(0)));

    return 0;
}

