/** @file fisher.c 
 * Hao Wang, 25.1.2013
 * Fisher driver of CLASS 
 */
 
#include "class.h"

#define NVARY 4 /* 4, dim of fisher matrix = NVARY x NVARY */
#define NREDSHIFT 3 /* 8, number of time slices */
#define ZMAX 2.0

int set_fiducial(
		            struct file_content * pfc,
	                ErrorMsg errmsg);

int class_assuming_bessels_computed(
				    struct file_content *pfc,
				    struct precision * ppr,
				    struct background * pba,
				    struct thermo * pth,
				    struct perturbs * ppt,
				    struct bessels * pbs,
				    struct transfers * ptr,
				    struct primordial * ppm,
				    struct spectra * psp,
				    struct nonlinear * pnl,
				    struct lensing * ple,
				    struct output * pop,

                    double z,
				    double * psCl,
                    int lmax,

				    ErrorMsg errmsg);

int vary_parameter(
				    struct file_content *pfc,
				    struct precision * ppr,
				    struct background * pba,
				    struct thermo * pth,
				    struct perturbs * ppt,
				    struct bessels * pbs,
				    struct transfers * ptr,
				    struct primordial * ppm,
				    struct spectra * psp,
				    struct nonlinear * pnl,
				    struct lensing * ple,
				    struct output * pop,

				    double z,
                    double * psCl1,
				    double * psCl0,
				    int lmax,
                    int ipara,
                    double * idCldPara,

				    ErrorMsg errmsg);
 
int count;
int main() {
count = 0;
  struct precision pr;        /* for precision parameters */
  struct background ba;       /* for cosmological background */
  struct thermo th;           /* for thermodynamics */
  struct perturbs pt;         /* for source functions */
  struct bessels bs;          /* for bessel functions */
  struct transfers tr;        /* for transfer functions */
  struct primordial pm;       /* for primordial spectra */
  struct spectra sp;          /* for output spectra */
  struct nonlinear nl;        /* for non-linear spectra */
  struct lensing le;          /* for lensed spectra */
  struct output op;           /* for output files */
  ErrorMsg errmsg;            /* for error messages */

  int i,j; 
  int num_ct_max=7;

  int lmax;
  int l; /* l starts from 2 (quadrupole) */

  double * psCl;  /*power spectrum in 2-D */
  double * psCl_fid;
  double * dCldPara[NVARY]; /* first derivatives, crucial for Fisher matrix */
  double * var_cl; /* uncertainty of Cl's using a_lm as observables???? Might differ from CMB */

  double z;
  int iz;

  double zs[NREDSHIFT];
  double zmin=0.5;

  double fisher[NVARY][NVARY]={0};

  /* parameter of the telescope */
  double tele_fwhm0 = 2.9; /* angular resolution in minute degree, here for FAST */ 
  double tele_fwhm; /* angular resolution at z: tele_fwhm = tele_fwhm0*(1+z) */
  double tele_temp_sys = 30.; /* system temperature in K; always mK in calculation */
  double tele_time_integ = 10.; /* integration time */
  double tele_freq = 1.e6; /* delta_nu in Hz; noise = temp_sys / sqrt(freq*time_integ) */
  double tele_fsky = 0.585;

  double noise_ins;
  double noise_fg = 0.;
  double theta_res, deltaN; 
  double tmp_cl_noise;
  double delta_N; 

  struct file_content fc;

  double deg2rad = M_PI/180.;
  tele_fwhm0 = tele_fwhm0/60.*deg2rad;

  for(iz=0;iz<NREDSHIFT;iz++) zs[iz] = zmin + iz*(ZMAX-zmin)/NREDSHIFT;

  /* radio equation in mK*/
  noise_ins = tele_temp_sys*1000. / sqrt(tele_freq*tele_time_integ);  

  if (set_fiducial(&fc,errmsg) == _FAILURE_) {
    printf("\n\nError set fiducial parameters \n=>%s\n",errmsg); 
    return _FAILURE_;
  }

  /* read input parameters and compute bessel functions once and for all */
  if (input_init(&fc,&pr,&ba,&th,&pt,&bs,&tr,&pm,&sp,&nl,&le,&op,errmsg) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg); 
    return _FAILURE_;
  }

  if (bessel_init(&pr,&bs) == _FAILURE_) {
    printf("\n\nError in bessel_init \n =>%s\n",bs.error_message);
    return _FAILURE_;
  }

  for (iz = 0; iz<NREDSHIFT; iz++) {
    z = zs[iz];
    tele_fwhm = tele_fwhm0 * (1.+z); /* [minute deg] tele_fwhm[minute degree] grows with z */
    lmax = (int) sqrt(4*M_PI / pow(tele_fwhm,2))-1; //(lmax+1)^2 = num_pixel
    psCl_fid = malloc((lmax+1)*sizeof(double));
    var_cl = malloc((lmax+1)*sizeof(double));
    for(i=0; i<NVARY; i++) dCldPara[i] = malloc((lmax+1)*sizeof(double));
    
    /* calls class (except the bessel module which has been called
     once and for all) and return the psCl[l]'s */
    class_assuming_bessels_computed(&fc,&pr,&ba,&th,&pt,&bs,&tr,&pm,&sp,&nl,&le,&op, \
                    z,psCl_fid,lmax,errmsg);

    /* now, loop over parameters */
    for (i=0; i<NVARY; i++) {
        /* vary parameter and get dCldPara[ipara][il]'s at redshift z; allocate
         * an array psCl first; psCl is used only inside the loop, so allocate
         * and free inside */
        psCl = malloc((lmax+1)*sizeof(double));
        vary_parameter(&fc,&pr,&ba,&th,&pt,&bs,&tr,&pm,&sp,&nl,&le,&op, \
                    z,psCl,psCl_fid,lmax,i,dCldPara[i],errmsg);  
        free(psCl);
    }

    /* get the value for instrument */
    delta_N = noise_ins + noise_fg;
    theta_res = tele_fwhm / sqrt(8*log(2));
    for (l=2; l<=lmax; l++) {
        tmp_cl_noise = pow(tele_fwhm*delta_N, 2) * exp(l*(l+1)*theta_res*theta_res);
        var_cl[l] = 2./((2*l+1)*tele_fsky) * pow((10*psCl_fid[l]+tmp_cl_noise), 2);
    }
        
    for (i = 0; i<NVARY; i++) 
      for (j = 0; j<NVARY; j++) 
        for (l = 2; l <= lmax; l++) 
            fisher[i][j] += dCldPara[i][l]*dCldPara[j][l]/var_cl[l]; 

    free(psCl_fid);
    free(var_cl);
    for(i=0; i<NVARY; i++) free(dCldPara[i]);
  }

  for (i=0; i<NVARY; i++){
    printf("\n"); 
    for (j=0; j<NVARY; j++)
      printf("%20.10e ", fisher[i][j]);
  }

  /* now free the bessel structure */
  if (bessel_free(&bs) == _FAILURE_)  {
    printf("\n\nError in bessel_free \n=>%s\n",bs.error_message);
    return _FAILURE_;
  }

  return _SUCCESS_;
  
}

int class_assuming_bessels_computed(
				    struct file_content *pfc,
				    struct precision * ppr,
				    struct background * pba,
				    struct thermo * pth,
				    struct perturbs * ppt,
				    struct bessels * pbs,
				    struct transfers * ptr,
				    struct primordial * ppm,
				    struct spectra * psp,
				    struct nonlinear * pnl,
				    struct lensing * ple,
				    struct output * pop,
                    double z,
				    double * psCl,
				    int lmax,
				    ErrorMsg errmsg) {

  /*local variables*/  
  double T21;
  int l;
  double tau; /* conformal age, in unit of Mpc */
  double pk_tmp;
  double k;
  double * pk_ic;

  if (input_init(pfc,ppr,pba,pth,ppt,pbs,ptr,ppm,psp,pnl,ple,pop,errmsg) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg); 
    return _FAILURE_;
  }

  if (background_init(ppr,pba) == _FAILURE_) {
    printf("\n\nError running background_init \n=>%s\n",pba->error_message);
    return _FAILURE_;
  }
    
  if (thermodynamics_init(ppr,pba,pth) == _FAILURE_) {
    printf("\n\nError in thermodynamics_init \n=>%s\n",pth->error_message);
    return _FAILURE_;
  }

  if (perturb_init(ppr,pba,pth,ppt) == _FAILURE_) {
    printf("\n\nError in perturb_init \n=>%s\n",ppt->error_message);
    return _FAILURE_;
  }

  if (transfer_init(ppr,pba,pth,ppt,pbs,ptr) == _FAILURE_) {
    printf("\n\nError in transfer_init \n=>%s\n",ptr->error_message);
    return _FAILURE_;
  }

  if (primordial_init(ppr,ppt,ppm) == _FAILURE_) {
    printf("\n\nError in primordial_init \n=>%s\n",ppm->error_message);
    return _FAILURE_;
  }

  if (spectra_init(ppr,pba,ppt,ptr,ppm,psp) == _FAILURE_) {
    printf("\n\nError in spectra_init \n=>%s\n",psp->error_message);
    return _FAILURE_;
  }

  if (nonlinear_init(ppr,pba,pth,ppt,pbs,ptr,ppm,psp,pnl) == _FAILURE_) {
    printf("\n\nError in nonlinear_init \n=>%s\n",pnl->error_message);
    return _FAILURE_;
  }

  if (lensing_init(ppr,ppt,psp,pnl,ple) == _FAILURE_) {
    printf("\n\nError in lensing_init \n=>%s\n",ple->error_message);
    return _FAILURE_;
  }

  if (output_init(pba,ppt,psp,pnl,ple,pop) == _FAILURE_) {
    printf("\n\nError in output_init \n=>%s\n",pop->error_message);
    return _FAILURE_;
  }
  double c_light  = 299792458.0;
  double c_planck = 6.626068e-34;
  double c_boltzmann = 1.3806503e-23;
  double c_nu21 = 1420.4e6;
  double c_emis = 2.85e-15;
  double c_newton = 6.67300e-11;
  double c_lyr2m = 9.4605284e15;
  double c_kmsmpc2s;
  double c_mH = 1.67e-27;
  double aver_Tau;
  double aver_x_H = 0.01;
  double bias_x_H = 1.;
  c_kmsmpc2s = 1000./(3.26e6*c_lyr2m);
  /* average optical depth */
  aver_Tau = 3*pow(c_light,3)*c_planck/(2*M_PI)*c_emis* \
       3*pow(100*pba->h*c_kmsmpc2s, 2)/(8*M_PI*c_newton)*pba->Omega0_b*pow(1.+z,3)/c_mH * aver_x_H/ \
       (16.*c_boltzmann*pow(c_nu21,2) *100*pba->h*c_kmsmpc2s *\
       sqrt((pba->Omega0_cdm+pba->Omega0_b)*pow(1.+z,3)+pba->Omega0_lambda));
  /* Tau = 3c^3*h*A_10*n_H / 16 
   * T21 = 23*(1+delta)*x_H * (T_s-T_cmb)/T_s * omegab/0.02 * 
   * sqrt(0.15/omegam * (1+z)/10) mK 
   * ==> (minus the mean value) ==>
   * T21 = 23*Y*delta * omegab/0.02 * sqrt(0.15/omega_m) mK */
  T21 = aver_Tau/(1.+z)*1000.; /* there is a factor 1/T_S missing in aver_Tau;
                                * T21 = (T_S - T_CMB)/aver_tau * (1+z), so
                                * as T_S>>T_CMB, we have T21 = aver_Tau/(1+z), 
                                * only note the unit is better to be mK (*1000) 
                                */

  background_tau_of_z(pba, z, &tau);

  /* THE FOLLOWING SEEMS INEFFICIENT */
  for (l=2;l<=lmax;l++){
    k = l/(pba->conformal_age - tau); /* tau in Mpc; k in spectra_pk_at_k_and_z is in [1/Mpc] */
    spectra_pk_at_k_and_z(pba, ppm, psp, k, z, &pk_tmp, pk_ic);
    /* 3-D vs 2-D:  l(l+1)Cl/2PI = k^3P21(k)/2PI^2, P21(k) = (T21*Y)^2 * P(k) */
    /* psCl[l]'s are in mK^2 */
    if(l%500==0) printf("z is %e, l is %d,  tau is %e, k is %e,  pk is %e\n", z, l, tau, k, pk_tmp);
    psCl[l] = 2*M_PI/(l*1.0*(l+1)) * T21*T21 * pow(1.+bias_x_H, 2) * pow(k, 3)*pk_tmp/(2*M_PI*M_PI);
  }


  /****** all calculations done, now free the structures ******/
  if (lensing_free(ple) == _FAILURE_) {
    printf("\n\nError in spectra_free \n=>%s\n",ple->error_message);
    return _FAILURE_;
  }

  if (nonlinear_free(pnl) == _FAILURE_) {
    printf("\n\nError in nonlinear_free \n=>%s\n",pnl->error_message);
    return _FAILURE_;
  }
    
  if (spectra_free(psp) == _FAILURE_) {
    printf("\n\nError in spectra_free \n=>%s\n",psp->error_message);
    return _FAILURE_;
  }
    
  if (primordial_free(ppm) == _FAILURE_) {
    printf("\n\nError in primordial_free \n=>%s\n",ppm->error_message);
    return _FAILURE_;
  }
  
  if (transfer_free(ptr) == _FAILURE_) {
    printf("\n\nError in transfer_free \n=>%s\n",ptr->error_message);
    return _FAILURE_;
  }

  if (perturb_free(ppt) == _FAILURE_) {
    printf("\n\nError in perturb_free \n=>%s\n",ppt->error_message);
    return _FAILURE_;
  }

  if (thermodynamics_free(pth) == _FAILURE_) {
    printf("\n\nError in thermodynamics_free \n=>%s\n",pth->error_message);
    return _FAILURE_;
  }

  if (background_free(pba) == _FAILURE_) {
    printf("\n\nError in background_free \n=>%s\n",pba->error_message);
    return _FAILURE_;
  }

  return _SUCCESS_;

}

int set_fiducial(
		struct file_content * pfc,
	        ErrorMsg errmsg
		){
	
  int i;
  char format[100];
  /* all parameters for which we don't want to keep default values
     should be passed to the code through a file_content
     structure. Create such a structure with the size you need: 9 in
     this exemple */
  parser_init(pfc,16,errmsg);

  /* assign values to these 9 parameters. Some will be fixed, some
     will be varied in the loop. */


/* the following is the parameters varied, indices help for-looping, otherwise,
 * use names, eg. pfc->h or pfc->omega_b */
  strcpy(pfc->name[0],"omega_b");
  sprintf(pfc->value[0],"%e",0.0226);

  strcpy(pfc->name[1],"n_s");
  sprintf(pfc->value[1],"%e",0.963);

  strcpy(pfc->name[2],"h");
  sprintf(pfc->value[2],"%e",0.704);

  strcpy(pfc->name[3],"omega_cdm");
  sprintf(pfc->value[3],"%e",0.113);
 /*end of varying parameter-------------------*/ 

/*
  strcpy(pfc->name[5],"m_ncdm");
  sprintf(pfc->value[5],"%e",0.04);

  strcpy(pfc->name[6],"N_ncdm");
  sprintf(pfc->value[6],"%d",1);
*/
  strcpy(pfc->name[7],"output");
  strcpy(pfc->value[7],"tCl,pCl,lCl,mPk");
/*
  strcpy(pfc->name[7],"selection_num"); 
  sprintf(pfc->value[7],"%d",NREDSHIFT);

  format = " "
  for(i=0;i<NREDSHIFT;i++) strcat(format, "%f,")
  strcpy(pfc->name[8],"selection_mean");
  sprintf(pfc->value[8],format, 0.5, 0.7, 0.9,  \
                                1.1, 1.3, 1.5,  \
                                1.7, 1.9);

  strcpy(pfc->name[9],"selection_width");
  sprintf(pfc->value[9],"%e",0.01);
*/
  strcpy(pfc->name[10],"z_max_pk");
  sprintf(pfc->value[10],"%e", ZMAX);

  strcpy(pfc->name[11],"l_max_scalars"); 
  sprintf(pfc->value[11],"%d", 3000);

  strcpy(pfc->name[12],"modes");
  sprintf(pfc->value[12],"%s","s");

  strcpy(pfc->name[13],"ic");
  sprintf(pfc->value[13],"%s","ad");
 
  strcpy(pfc->name[14], "Omega_k");
  sprintf(pfc->value[14], "%e", 0.);
 
//  strcpy(pfc->name[12], "background_verbose");
//  sprintf(pfc->value[12], "%d", 1);
}

int vary_parameter(
				    struct file_content *pfc,
				    struct precision * ppr,
				    struct background * pba,
				    struct thermo * pth,
				    struct perturbs * ppt,
				    struct bessels * pbs,
				    struct transfers * ptr,
				    struct primordial * ppm,
				    struct spectra * psp,
				    struct nonlinear * pnl,
				    struct lensing * ple,
				    struct output * pop,

                    double z,
				    double * psCl1,
				    double * psCl0,
				    int lmax,
                    int ipara,
                    double * idCldPara,

				    ErrorMsg errmsg) {
    int l;
    /* used in read_double, so must be double */
    int flag1;
    double param1, param2;
    /* proportional change of each parameter; so the fiducial parameter value should not be 0 */
    double delta = 0.1;   
    parser_read_double(pfc,pfc->name[ipara],&param1,&flag1,errmsg);
    param2 = param1 * (1.+delta);
    sprintf(pfc->value[ipara],"%e", param2);
printf("<=====================================================>");
    printf("%s from %e to %e\n", pfc->name[ipara], param1, param2);

    /* calls class again for *renewed* fc and return the P(k)'s*/
    class_assuming_bessels_computed(pfc,ppr,pba,pth,ppt,pbs,ptr,ppm,psp,pnl,ple,pop,z,psCl1,lmax,errmsg);

    for(l=2; l<=lmax; l++) idCldPara[l] = (psCl1[l]-psCl0[l])/(param1*delta); 

    /* reset fc to fiducial values */
    if (set_fiducial(pfc,errmsg) == _FAILURE_) {
      printf("\n\nError set fiducial parameters \n=>%s\n",errmsg); 
      return _FAILURE_;
    }
}
