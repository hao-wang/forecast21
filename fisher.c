/** @file fisher.c
 * Hao Wang, 23.5.2013
 *  use Omega_HI instead of Omega_b*x_HI
 *  add l-interval;
 *  remove theta_res = tele_fwhm / sqrt(8*log(2)) and consider
 *      theta_res=2.9';
 * Hao Wang, 25.1.2013
 * Fisher driver of CLASS
 * Telescope parameters are for FAST
 */

#include "class.h"

#define NVARY 6 /* dim of fisher matrix = NVARY x NVARY */
#define NREDSHIFT 21 /* (21-1 = ) 20 frequency intervals */
#define ZMAX 2.0
#define ZMIN 0.5
#define LMAX 2000 /* global l max, used in array dimensions */

/*char * trim(char *str);*/
    
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
    int lmin,
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
    int lmin,
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

    FILE * fp;
    FILE * fp2;

    char f_tmp[50];
    double fbest[NREDSHIFT][LMAX];

    int i,j;
    int num_ct_max=7;

    int lmax, lmin;
    int l; /* l starts from 2 (quadrupole) in spectrum calculation, but
              lmin for Fisher matrix */

    double * psCl;  /*power spectrum in 2-D */
    double * psCl_fid;
    double * dCldPara[NVARY]; /* first derivatives, crucial for Fisher matrix */
    double * var_cl; /* uncertainty of Cl's using a_lm as observables */

    double z;
    int iz;
    double khmax; /* use fitting formula to get smallest nonlinear kh, see 
                   mathematica notebook for fast21 */

    double zs[NREDSHIFT];
    double pvecback[100];

    double fisher[NVARY][NVARY]= {0};

    /* parameter of the telescope */
    double tele_res0 = 2.9; /* also FWHM, angular resolution in minute degree, will be converted to radian */
    double tele_res; /* angular resolution at z: tele_res = tele_res0*(1+z) */
    double tele_temp_sys = 20.; /* system temperature in K; always mK in calculation */
    double tele_time_pix; /* observation time per pixel, in [sec] */
    double tele_freq = 0.2e6; /* delta_nu in Hz; noise = temp_sys / sqrt(freq*time_pix) */
/*
    double tele_time_tot = 8. ;//survey time in [yr] 
    int    tele_num_beam=100;
    double tele_fsky = 0.575;
*/
    double tele_time_tot = 1. ;// survey time in [yr] 
    int    tele_num_beam=1;
    double tele_fsky = 0.0005;

    double noise_ins;
    double noise_fg = 0.; /* RMS of the residual foreground noise, in mK */
    double deltaN;
    double tmp_cl_noise;
    double delta_N;
    double freqH = 1420.4;
    double freq_step, center_f;

    struct file_content fc;

    double tau;
    double deg2rad = M_PI/180.;
    double yr2sec = 31556926;
    tele_res0 = tele_res0/60.*deg2rad;
    freq_step = freqH*(1./(1.+ZMIN) - 1./(1.+ZMAX))/(NREDSHIFT-1.);

    center_f=freqH*(1./(1.+ZMIN));
    for(iz=0; iz<NREDSHIFT; iz++,center_f-=freq_step) zs[iz] = freqH/center_f-1;

for(iz=0;iz<NREDSHIFT; iz++) printf("one z is %f \n",zs[iz]);
    if (set_fiducial(&fc,errmsg) == _FAILURE_) {
        printf("\n\nError set fiducial parameters \n=>%s\n",errmsg);
        return _FAILURE_;
    }

    /* read input parameters and compute bessel functions once and for all */
    if (input_init(&fc,&pr,&ba,&th,&pt,&bs,&tr,&pm,&sp,&nl,&le,&op,errmsg) == _FAILURE_) {
        printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg);
        return _FAILURE_;
    }

    if (background_init(&pr,&ba) == _FAILURE_) {
        printf("\n\nError running background_init \n=>%s\n",ba.error_message);
        return _FAILURE_;
    }

    if (bessel_init(&pr,&bs) == _FAILURE_) {
        printf("\n\nError in bessel_init \n =>%s\n",bs.error_message);
        return _FAILURE_;
    }

fp2=fopen("fsky_best.dat", "w");
    for (iz = 0; iz<NREDSHIFT; iz++) {
        z = zs[iz];
        tele_res = tele_res0 * (1.+z); /* [minute deg] tele_res[minute degree] grows with z */
        khmax=0.170255158514253+0.13037369068045573*z-0.10832440336426165* pow(z,2)+0.046760120277751734*pow(z,3) - 0.005948516308501343* pow(z,4); /* fitting formula for minimum NL k-modes */
        background_tau_of_z(&ba, z, &tau);
        lmax = (int) (khmax*ba.h * (ba.conformal_age - tau));
        printf("lmax = %d\n", lmax);
        background_functions(&ba, 1/(1.+z), 1, pvecback);
        lmin = (int) (ba.conformal_age - tau) * 2*M_PI/(pow(1+z,2)/pvecback[ba.index_bg_H]*freq_step/freqH);
        printf("lmin = %d\n", lmin);
        psCl_fid = malloc((lmax+1)*sizeof(double));
        var_cl = malloc((lmax+1)*sizeof(double));
        for(i=0; i<NVARY; i++) dCldPara[i] = malloc((lmax+1)*sizeof(double));

        /* calls class (except the bessel module which has been called
         once and for all) and return the psCl[l]'s */
        class_assuming_bessels_computed(&fc,&pr,&ba,&th,&pt,&bs,&tr,&pm,&sp,&nl,&le,&op, \
                                        z,psCl_fid,lmin,lmax,errmsg);
	for(i=lmin; i<=lmax; i+=10) printf("psCl_fid[%d] = %e\n", i, psCl_fid[i]); 

        /* now, loop over parameters */
        for (i=0; i<NVARY; i++) {
            /* vary parameter and get dCldPara[ipara][il]'s at redshift z; allocate
             * an array psCl first; psCl is used only inside the loop, so allocate
             * and free inside */
            psCl = malloc((lmax+1)*sizeof(double));
            vary_parameter(&fc,&pr,&ba,&th,&pt,&bs,&tr,&pm,&sp,&nl,&le,&op, \
                           z,psCl,psCl_fid,lmin,lmax,i,dCldPara[i],errmsg);
            free(psCl);
        }

        /* get the value for instrument */
        tele_res=tele_res0 * (1+z);
        tele_time_pix=tele_time_tot * yr2sec * tele_num_beam * pow(tele_res,2) / (4*M_PI*tele_fsky);
        /* radio equation in mK, in consistence with signal
         * dimension */
        noise_ins = tele_temp_sys*1000. / sqrt(tele_freq*tele_time_pix);

        delta_N = noise_ins + noise_fg;
        for (l=lmin; l<=lmax; l++) {
            tmp_cl_noise = pow(tele_res,2) * pow(delta_N,2) * exp(l*(l+1) * pow(tele_res,2)/(8*log(2)));
            var_cl[l] = 2./((2*l+1)*tele_fsky) * pow(psCl_fid[l]+tmp_cl_noise, 2); //variance of cl
            printf("noise vs signal: %f, %f\n", tmp_cl_noise,psCl_fid[l]);
fbest[iz][l]=(tele_freq*tele_num_beam*tele_time_tot*yr2sec*psCl_fid[l])/(4*M_PI*pow(tele_temp_sys*1000.,2)*exp(l*(l+1) * pow(tele_res,2)/(8*log(2))));
        }

     if(NVARY != 0)
        for (i = 0; i<NVARY; i++)
            for (j = 0; j<NVARY; j++)
                for (l = lmin; l <= lmax; l++)
                    fisher[i][j] += dCldPara[i][l]*dCldPara[j][l]/var_cl[l];

        free(psCl_fid);
        free(var_cl);
        for(i=0; i<NVARY; i++) free(dCldPara[i]);
    }
    for (iz=0; iz<NREDSHIFT; iz++)
    {
      for (l=0; l<LMAX; l++) fprintf(fp2, "%f ", fbest[iz][l]);
      fprintf(fp2, "\n");
    } 
    fclose(fp2);

    fp = fopen("fisher.mat","w"); 

    for (i=0; i<NVARY; i++) {
        for (j=0; j<NVARY; j++) {
            printf("%20.10e ", fisher[i][j]);
	    fprintf(fp, "%20.10e ", fisher[i][j]);
	}
    	fprintf(fp, "\n");
    }
    for (i=0; i<NVARY; i++) fprintf(fp, "%s ", fc.name[i]);
    fprintf(fp, "\n");
    for (i=0; i<NVARY; i++) fprintf(fp, "%s ", fc.value[i]);
    fprintf(fp, "\n");

    fclose(fp);

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
    int lmin,
    int lmax,
    ErrorMsg errmsg) {

    /*local variables*/
    double pvecback[100];
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

    background_functions(pba, 1/(1.+z), 1, pvecback);
    /* T21 := 7.59*10^-2 h (1+\delta) (1+z)^2 / E(z);
     * E[z]:= Sqrt[Omega_m(1+z)^3+Omega_l exp(-3 Integrate[(1+w)/a, a])] mK */
    T21 = 7.59 * 0.01 * pba->h * pow(1.+z,2) / (pvecback[pba->index_bg_H]/pba->H0); 

    background_tau_of_z(pba, z, &tau);

    for (l=lmin; l<=lmax; l+=1) {
        k = l/(pba->conformal_age - tau); /* tau in Mpc; k in spectra_pk_at_k_and_z is in [1/Mpc] */
        spectra_pk_at_k_and_z(pba, ppm, psp, k, z, &pk_tmp, pk_ic);
        /* 3-D vs 2-D:  l(l+1)Cl/2PI = k^3 P21(k)/2PI^2, P21(k) = (T21*Y)^2 * P(k) */
        /* psCl[l]'s are in mK^2 */
        if(l%100==0) printf("z is %e, l is %d,  tau is %e, k is %e,  pk is %e\n", z, l, tau, k, pk_tmp);
        psCl[l] = 2*M_PI/(l*1.0*(l+1)) * T21*T21 * pow(k, 3) * pk_tmp/(2*M_PI*M_PI);
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
    ErrorMsg errmsg) {

    int i;
    char format[100];
    /* all parameters for which we don't want to keep default values
       should be passed to the code through a file_content
       structure. Create such a structure with the size you need: 9 in
       this exemple */
    parser_init(pfc,16,errmsg);
    i=0;

    /* assign values to these 9 parameters. Some will be fixed, some
       will be varied in the loop. */


    /* the following is the parameters varied, indices help for-looping, otherwise,
     * use names, eg. pfc->h or pfc->omega_b */
    strcpy(pfc->name[i],"omega_b");
    sprintf(pfc->value[i],"%e",0.0220);
    i++;

    strcpy(pfc->name[i],"n_s");
    sprintf(pfc->value[i],"%e",0.9617);
    i++;

    strcpy(pfc->name[i],"h");
    sprintf(pfc->value[i],"%e",0.6927);
    i++;

    strcpy(pfc->name[i],"omega_cdm");
    sprintf(pfc->value[i],"%e",0.1207);
    i++;
/*
    strcpy(pfc->name[i], "Omega_k");
    sprintf(pfc->value[i], "%e", 0.);
    i++;
*/
    strcpy(pfc->name[i], "w0_fld");
    sprintf(pfc->value[i], "%e", -1.007);
    //here is the only place original class need to 
    //be modified: you need to comment out the range 
    //check for w0_fld in source/input.c.
    i++;

    strcpy(pfc->name[i], "wa_fld");
    sprintf(pfc->value[i], "%e", -0.29);
    i++;

    /*end of varying parameter-------------------*/

    strcpy(pfc->name[i], "Omega_fld"); //note that curved space is not supported
    sprintf(pfc->value[i], "%e", 0.6825);
    i++;

    /*
      strcpy(pfc->name[i],"m_ncdm");
      sprintf(pfc->value[i],"%e",0.04);
    i++;

      strcpy(pfc->name[i],"N_ncdm");
      sprintf(pfc->value[i],"%d",1);
    i++;
    */
    strcpy(pfc->name[i],"output");
    strcpy(pfc->value[i],"tCl,pCl,lCl,mPk");
    i++;
    /*
      strcpy(pfc->name[i],"selection_num");
      sprintf(pfc->value[i],"%d",NREDSHIFT);

      format = " "
      for(i=0;i<NREDSHIFT;i++) strcat(format, "%f,")
      strcpy(pfc->name[8],"selection_mean");
      sprintf(pfc->value[8],format, 0.5, 0.7, 0.9,  \
                                    1.1, 1.3, 1.5,  \
                                    1.7, 1.9);

      strcpy(pfc->name[i],"selection_width");
      sprintf(pfc->value[i],"%e",0.01);
    */
    strcpy(pfc->name[i],"z_max_pk");
    sprintf(pfc->value[i],"%e", ZMAX);
    i++;

    strcpy(pfc->name[i],"l_max_scalars");
    sprintf(pfc->value[i],"%d", 3000);
    i++;

    strcpy(pfc->name[i],"modes");
    sprintf(pfc->value[i],"%s","s");
    i++;

    strcpy(pfc->name[i],"ic");
    sprintf(pfc->value[i],"%s","ad");
    i++;

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
    int lmin,
    int lmax,
    int ipara,
    double * idCldPara,

    ErrorMsg errmsg) {
    int l;
    int flag1;
    /* used in read_double, so must be double */
    double param1, param2;
    /* proportional change of each parameter; so the fiducial parameter value should not be 0 */
    double delta = 0.1;
    parser_read_double(pfc,pfc->name[ipara],&param1,&flag1,errmsg);
    param2 = param1 * (1.+delta);
    sprintf(pfc->value[ipara],"%e", param2);
    printf("<=================================>");
    printf("%s: %e --> %e; z: %f\n", pfc->name[ipara], param1, param2, z);

    /* calls class again for *renewed* fc and return the P(k)'s*/
    class_assuming_bessels_computed(pfc,ppr,pba,pth,ppt,pbs,ptr,ppm,psp,pnl,ple,pop,z,psCl1,lmin,lmax,errmsg);

    for(l=lmin; l<=lmax; l++) idCldPara[l] = (psCl1[l]-psCl0[l])/(param1*delta);

    /* reset fc to fiducial values */
    if (set_fiducial(pfc,errmsg) == _FAILURE_) {
        printf("\n\nError set fiducial parameters \n=>%s\n",errmsg);
        return _FAILURE_;
    }
}
/*
char * trim(char *str)
{
  char *str_last, *str_cur;
  if(str==NULL) return;
  for(; *str==0x20 || *str=='\t'; ++str);
  for(str_last=str_cur=str;*str_cur!='\0';++str_cur)
    if(*str_cur!=0x20 && *str_cur!='\t')
      str_last=str_cur;
  *++str_last=0;
  return str;
}
*/
