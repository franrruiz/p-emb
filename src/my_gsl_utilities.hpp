#ifndef MY_GSL_UTILITIES_HPP
#define MY_GSL_UTILITIES_HPP

const double myINFINITY=1.0e100;

/************************ Misc. Utility functions ************************/

inline double my_min(double a, double b) {
    return((a<b)?a:b);
}

inline double my_max(double a, double b) {
    return((a<b)?b:a);
}

inline double my_exp(double val) {
    if(val>700.0) {
        return(gsl_sf_exp(700.0));
    } else if(val<-690.775528) {
        return(1e-300);
    } else {
        return(gsl_sf_exp(val));
    }
}

inline double my_sigmoid(double val) {
    if(val>0) {
        if(val>23.0259) {
            val = 23.0259;
        }
        return(1.0/(1.0+my_exp(-val)));
    } else {
        if(val<-23.0259) {
            val = -23.0259;
        }
        double aux = my_exp(val);
        return(aux/(1.0+aux));
    }
}

inline double my_pow2(double val) {
    return(val*val);
}

inline double my_log(double val) {
    return (val<=0?-myINFINITY:gsl_sf_log(val));
}

inline double my_logit(double val) {
    if(val>=1.0-1.0e-6) {
        val = 1.0-1.0e-6;
    } else if(val<1.0e-6) {
        val = 1.0e-6;
    }
    return(my_log(val)-my_log(1-val));
}

inline double logsumexp(double *x, int N1) {
    double maxim = -myINFINITY;
    double suma = 0;
    for(int n1=0; n1<N1; n1++) {
        maxim = (x[n1]>maxim)?x[n1]:maxim;
    }
    for(int n1=0; n1<N1; n1++) {
        suma += my_exp(x[n1]-maxim);
    }
    return(maxim+my_log(suma));
}

inline double my_gsl_ran_gamma(const gsl_rng *r, double shape, double scale) {
    double aux = gsl_ran_gamma(r,shape,scale);
    return((aux<1e-300)?1e-300:aux);
}

inline double gsl_ran_invgamma(const gsl_rng *r, double shape, double scale) {
    return(1.0/my_gsl_ran_gamma(r,shape,1.0/scale));
}

inline double my_gsl_ran_beta(const gsl_rng *r, double p1, double p2) {
    double aux = gsl_ran_beta(r,p1,p2);
    aux = (aux<1e-300)?1e-300:aux;
    return((aux>1.0-1e-300)?1.0-1e-300:aux);
}

inline double my_gsl_ran_gammaE(const gsl_rng *r, double shape, double mu) {
    return(my_gsl_ran_gamma(r,shape,mu/shape));
}

inline void my_gsl_ran_dirichlet(const gsl_rng *r, size_t K, const double alpha[], double theta[]) {
    gsl_ran_dirichlet(r,K,alpha,theta);
    for(unsigned int k=0; k<K; k++) {
        if(theta[k]<1e-300) {
            theta[k] = 1e-300;
        } else if(theta[k]>1.0-1e-300) {
            theta[k] = 1.0-1e-300;
        }
    }
}

inline double my_gsl_sf_lngamma(double val) {
    return((val<1e-300)?gsl_sf_lngamma(1e-300):gsl_sf_lngamma(val));
}

inline double my_logfactorial(double val) {
    return((val==0.0)?0.0:my_gsl_sf_lngamma(val+1.0));
}

inline double my_gsl_sf_psi(double val) {
    return((val<1e-300)?gsl_sf_psi(1e-300):gsl_sf_psi(val));
}

/*************************** Get/Set functions ***************************/

inline double getdouble_2D(const double *x, int n1, int n2, int N1, int N2) {
    return x[static_cast<unsigned long long>(N1*n2)+n1];
}

inline double getdouble_3D(const double *x, int n1, int n2, int n3, int N1, int N2, int N3) {
    return x[static_cast<unsigned long long>(N1*N2*n3)+static_cast<unsigned long long>(N1*n2)+n1];
}

inline double getdouble_4D(const double *x, int n1, int n2, int n3, int n4, int N1, int N2, int N3, int N4) {
    return x[static_cast<unsigned long long>(N1*N2*N3*n4)+static_cast<unsigned long long>(N1*N2*n3)+static_cast<unsigned long long>(N1*n2)+n1];
}

inline void setdouble_2D(double val, double *x, int n1, int n2, int N1, int N2) {
    x[static_cast<unsigned long long>(N1*n2)+n1] = val;
}

inline void setdouble_3D(double val, double *x, int n1, int n2, int n3, int N1, int N2, int N3) {
    x[static_cast<unsigned long long>(N1*N2*n3)+static_cast<unsigned long long>(N1*n2)+n1] = val;
}

inline void setdouble_4D(double val, double *x, int n1, int n2, int n3, int n4, int N1, int N2, int N3, int N4) {
    x[static_cast<unsigned long long>(N1*N2*N3*n4)+static_cast<unsigned long long>(N1*N2*n3)+static_cast<unsigned long long>(N1*n2)+n1] = val;
}

inline int getint_2D(const double *x, int n1, int n2, int N1, int N2) {
    return static_cast<int>(x[static_cast<unsigned long long>(N1*n2)+n1]);
}

inline int getint_2D(const int *x, int n1, int n2, int N1, int N2) {
    return x[static_cast<unsigned long long>(N1*n2)+n1];
}

inline int getint_3D(const double *x, int n1, int n2, int n3, int N1, int N2, int N3) {
    return static_cast<int>(x[static_cast<unsigned long long>(N1*N2*n3)+static_cast<unsigned long long>(N1*n2)+n1]);
}

inline int getint_3D(const int *x, int n1, int n2, int n3, int N1, int N2, int N3) {
    return x[static_cast<unsigned long long>(N1*N2*n3)+static_cast<unsigned long long>(N1*n2)+n1];
}

inline void setint_2D(int val, int *x, int n1, int n2, int N1, int N2) {
    x[static_cast<unsigned long long>(N1*n2)+n1] = val;
}

inline void setint_3D(int val, int *x, int n1, int n2, int n3, int N1, int N2, int N3) {
    x[static_cast<unsigned long long>(N1*N2*n3)+static_cast<unsigned long long>(N1*n2)+n1] = val;
}

inline bool getbool_2D(const bool *x, int n1, int n2, int N1, int N2) {
    return x[static_cast<unsigned long long>(N1*n2)+n1];
}

inline void setbool_2D(int val, bool *x, int n1, int n2, int N1, int N2) {
    x[static_cast<unsigned long long>(N1*n2)+n1] = val;
}

#endif

