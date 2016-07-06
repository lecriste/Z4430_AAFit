

//================ Decay Momentum ====================
// Momentum in 2-particle decay : m0->m1+m2

double RooMyPdf::sq_calc(double x,double y, double z) const
{
    return pow(x,2)+pow(y,2)+pow(z,2)-2.0*x*y-2.0*x*z+2.0*y*z;
}

double RooMyPdf::dec2mm (double m0, double m1, double m2) const
{
    double temp = sq_calc(m0*m0,m1*m1,m2*m2);
    return sqrt(temp)/(2.0*m0);
}
//================ Decay Momentum ====================


//================ Blatt-Weisskopf Form Factors ======
// l = spin
// q = momentum from "dec2mm"
// q0 = momentum from "dec2mm" with PDG mass
// r = meson radial parameter (hadron scale)
double RooMyPdf::bwff(int l, double q, double q0, double r) const
{
    double z = r*r*q*q;
    double z0 = r*r*q0*q0;
    double f;
    //########### spin 0 ###############
    if (l == 0) {
        f = 0;
    }
    //########### spin 1 ###############
    if (l == 1) {
        f = sqrt((1+z0)/(1+z));
    }
    //########### spin 2 ###############
    if (l == 2) {
        f = sqrt((z0*z0+3.0*z0+9.0)/(z*z+3.0*z+9.0));
    }
    //########### spin 3 ###############
    if (l == 3) {
        f = sqrt((z0*z0*z0+6.0*z0*z0+45.0*z0+225.0)/(z*z*z+6.0*z*z+45.0*z+225.0));
    }
    
    return f;
}
//================ Blatt-Weisskopf Form Factors ======

//================ Breit-Wigner Amplitude ============
// m0 = resonance mass (pdg)
// w0 = width (pdg)
// m = invariant mass of two daughters of the resonance
// m_d1, m_d2 = daughter masses
// l = relative angular momentum
// f = BW form factor
// q = momentum from "dec2mm"
// q0 = momentum from "dec2mm" with PDG mass

complex<double> RooMyPdf::bwamp(double m0,double w0,double m,double m_d1,double m_d2,int l,double f,double q0,double q) const
{
    double width = w0*pow((q/q0),2*l+1)*(m0/m)*f*f;
    double deno = (m0*m0 - m*m)*(m0*m0 - m*m) + m0*m0*width*width;
    double rl = f*pow((q/m),l)*(m0*m0 - m*m)/deno;
    double imag = f*pow((q/m),l)*m0*width/deno;
    complex<double> val(rl,imag);
    
    return val;
}

//================ Jacobi Polynomial =================
//Jacobi polynomial - order n
double RooMyPdf::jacobi_Pn (int n, double a, double b, double x) const
{
    
    if (n==0){
        return 1.0;
    }
    else if (n==1){
        return  0.5 * (a - b + (a + b + 2.0)*x);
    }
    else {
        
        double p0, p1, a1, a2, a3, a4, p2=0.0;
        int i;
        p0 = 1.0;
        p1 = 0.5 * (a - b + (a + b + 2)*x);
        
        for(i=1; i<n; ++i){
            a1 = 2.0*(i+1.0)*(i+a+b+1.0)*(2.0*i+a+b);
            a2 = (2.0*i+a+b+1.0)*(a*a-b*b);
            a3 = (2.0*i+a+b)*(2.0*i+a+b+1.0)*(2.0*i+a+b+2.0);
            a4 = 2.0*(i+a)*(i+b)*(2.0*i+a+b+2.0);
            p2 = 1.0/a1*( (a2 + a3*x)*p1 - a4*p0);
            
            p0 = p1;
            p1 = p2;
        }
        
        return p2;
    }
    
}
//================ Jacobi Polynomial =================


//================ factorial =========================
int RooMyPdf::Factorial(int x) const
{
    if (x==0) { return 1; }
    return (x == 1 ? x : x * Factorial(x - 1));
}
//================ factorial =========================


//================ combination =======================
int RooMyPdf::Combination(int n, int r) const
{
  return (Factorial(n)) / ((Factorial(n-r)) * Factorial(r));
}
//================ combination =======================



//================ wigner d calculations =============

double RooMyPdf::wigner_d (int j, int m1, int m2, double theta ) const
{
    int array[] = {j+m1, j-m1, j+m2, j-m2};
    int k = *min_element(array,array+4) ;
    int a = abs(m1-m2);//fabs?
    double lambda; //int - not working due to pow overload resolution
    if (k == j+m1) { lambda = 0;}
    else if (k == j-m1) { lambda = m1-m2;}
    else if (k == j+m2) { lambda = m1-m2;}
    else if (k == j-m2) { lambda = 0;}
    
    int b = 2*j-2*k-a;
    
    double value = pow(-1,lambda) * pow(Combination(2*j-k,k+a),0.5) * pow(Combination(k+b,b),-0.5) * pow(sin(0.5*theta),a) * pow(cos(0.5*theta),b) * jacobi_Pn(k,a,b,cos(theta));

    return value;
 
}
//================ wigner d calculations =============

//================ Signal Density Calculation ========
//pB = B0 3-momentum
double RooMyPdf::get_signal_density (double mBcalc, double mKPicalc, double mJpsicalc, double pB, double theta_k, double phi, double theta_jpsi ) const
{ // signal density begin

    double qB = pB/mB;
    double qB2 = qB*qB;
    double qB3 = qB2*qB;
    double qB4 = qB3*qB;
    double qB5 = qB4*qB;
    double q = dec2mm(mBcalc,mKPicalc,mJpsicalc); //dec2mm(mB,mKPi,mJpsi)
    //================ Amplitudes for the different K resonances===========
    // mKPi = K-Pi inv mass calculated from data
    double qK = dec2mm(mKPicalc,mK,mPi);
    //############## K*(892) ###################
    double qK892 = dec2mm(mK892,mK,mPi);
    double fK892 = bwff(1,qK,qK892,rR);
    complex<double> a_K_892 = bwamp(mK892,wK892,mKPicalc,mK,mPi,1,fK892,qK892,qK);

    //############## K0*(1430) ###################
    double q0_1430 = dec2mm(mB,mK0_1430,mJpsi);
    double qK0_1430 = dec2mm(mK0_1430,mK,mPi);
    double fK0_1430 = bwff(1,qK,qK0_1430,rR);
    complex<double> a_K0_1430 = bwamp(mK0_1430,wK0_1430,mKPicalc,mK,mPi,1,fK0_1430,qK0_1430,qK);
    a_K0_1430 = a_K0_1430 * qB * bwff(1,q,q0_1430,rB);

    //############## K2*(1430) ###################
    double q2_1430 = dec2mm(mB,mK2_1430,mJpsi);
    double qK2_1430 = dec2mm(mK2_1430,mK,mPi);
    double fK2_1430 = bwff(1,qK,qK2_1430,rR);
    complex<double> a_K2_1430 = bwamp(mK2_1430,wK2_1430,mKPicalc,mK,mPi,1,fK2_1430,qK2_1430,qK);
    a_K2_1430 = a_K2_1430 * qB * bwff(1,q,q2_1430,rB);

//******************** lepton pair helicity minus 1**************
    complex<double> index_minus1_m1(0.0,-1*phi);
    complex<double> a_K_892_minus1_m1 = a_K_892*wigner_d(1,-1,0,theta_k)*exp(index_minus1_m1)*wigner_d(1,-1,-1,theta_jpsi);
    
    complex<double> index_zero_m1(0.0,0.0);
    complex<double> a_K_892_zero_m1 = a_K_892*wigner_d(1,0,0,theta_k)*exp(index_zero_m1)*wigner_d(1,0,-1,theta_jpsi);
    
    complex<double> index_plus1_m1(0.0,phi);
    complex<double> a_K_892_plus1_m1 = a_K_892*wigner_d(1,1,0,theta_k)*exp(index_plus1_m1)*wigner_d(1,1,-1,theta_jpsi);
    
//******************** lepton pair helicity plus 1**************
    complex<double> index_minus1_p1(0.0,-1*phi);
    complex<double> a_K_892_minus1_p1 = a_K_892*wigner_d(1,-1,0,theta_k)*exp(index_minus1_p1)*wigner_d(1,-1,1,theta_jpsi);
    
    complex<double> index_zero_p1(0.0,0.0);
    complex<double> a_K_892_zero_p1 = a_K_892*wigner_d(1,0,0,theta_k)*exp(index_zero_p1)*wigner_d(1,0,1,theta_jpsi);
    
    complex<double> index_plus1_p1(0.0,phi);
    complex<double> a_K_892_plus1_p1 = a_K_892*wigner_d(1,1,0,theta_k)*exp(index_plus1_p1)*wigner_d(1,1,1,theta_jpsi);
    
    
//******************* helicity phase ***************************


    double val = pow(abs(a_K_892_minus1_m1 //* exp(helphase_index_K_892_minus1)
                        +a_K_892_zero_m1   //* exp(helphase_index_K_892_zero)
                        +a_K_892_plus1_m1  //* exp(helphase_index_K_892_plus1)
                        
                        +a_K_892_minus1_p1 //* exp(helphase_index_K_892_minus1)
                        +a_K_892_zero_p1   //* exp(helphase_index_K_892_zero)
                        +a_K_892_plus1_p1 ),2); //* exp(helphase_index_K_892_plus1)),2);
    return val ;

    

//================ Amplitudes for the different K resonances===========

}// signal density end
//================ Signal Density Calculation ========



