// Halo Boost Factor module
#define print_debug_info 0
#define Use_Conde_Concentration 1
#define Boost_nmh 1000
#define Boost_Interp_Table_Size 200
#include "Tables.h"

double Halo_Concentration_Conde(double m, double z)
{
    /*
    Halo concentration at z=0 following Eq.1 of MNRAS 442, 2271–2277 (2014)
    */
   double x, c0, c1, c2, c3, c4, c5, h, r;
   h = 0.6766;
   x = log(m*h);
   c0 = 37.5153;
   c1 = -1.5093;
   c2 = 1.636E-2;
   c3 = 3.66E-4;
   c4 = -2.89237E-5;
   c5 = 5.32E-7;
   
   r = c0 + c1 * x + c2 * pow(x, 2.0) + c3 * pow(x, 3.0) + c4 * pow(x, 4.0) + c5 * pow(x, 5.0);
   r = r / (1.0+z);
   return r;
}

double HaloProfile_Kernel(double z, double mh, double r, int ProfileType)
{
    /* Halo density profile
    ---- inputs ----
    z : redshift
    mh : halo mass in msun, dm+baryon
    r : distance to center in pc
    ProfileType : result type
                  0 - viral radius in pc
                  1 - Rho_DM in msun/pc^3
    */

    // printf("Check that the integrated mass converges to m\n");
    double OmM, OmC, m, OmR, OmL, h, pi, rho_cr0, zp, OmMz, d, Delta_C, log10_c, c, delta_c, rv1, rv2, rv3, r_vir, x, cx, rho_cr, RhoDM;

    // Some settings
    OmM = 0.30964168161;
    OmR = 9.1e-5;
    OmC = 0.260667;
    OmL = 0.69026731839;
    h = 0.6766;
    pi = 3.141592653589793;
    rho_cr0 = 2.775e-7 * pow(h, 2.); // critical density in msun/pc^3

    // Pre-requisites
    m = mh * OmC / OmM; // DM mass
    zp = 1. + z;
    OmMz = OmM * pow(zp, 3.) / (OmM * pow(zp, 3.) + OmL);
    d = OmMz - 1.;
    Delta_C = 18. * pow(pi, 2.) + 82. * d - 39. * pow(d, 2.);
    log10_c = 1.071 - 0.098 * (log10(m) - 12.);
    c = pow(10., log10_c) / zp; // concentration, see appdx.A of Zip.et for the additional (1+z) factor
    delta_c = Delta_C * pow(c, 3.) / (3. * (log(1. + c) - c / (1. + c)));

    rv1 = 0.784 * pow(m * h / 1.0e8, 1. / 3.);
    rv2 = pow(OmM * Delta_C / (OmMz * 18. * pow(pi, 2.)), -1. / 3.);
    rv3 = (10. / zp) / h * 1000.;
    r_vir = rv1 * rv2 * rv3;

    x = r / r_vir;
    if (x > 1.)
    {
        RhoDM = 0.;
        // printf("%E  %E\n", r, r_vir);
    }
    else
    {
        cx = c * x;
        // rho_cr = rho_cr0 * (OmL + OmM * zp**3 + OmR * zp**4);
        rho_cr = rho_cr0 * (OmL + OmM * pow(zp, 3.) + OmR * pow(zp, 4.));
        RhoDM = rho_cr * delta_c / (cx * pow(1. + cx, 2.));
    }
    if (ProfileType == 0)
    {
        return r_vir;
    }
    else
    {
        return RhoDM;
    }
}

double HaloProfile_Integrator(double z, double mh)
{
    /* Do the following HaloProfile integration analytically:
    \int dr r^2 \rho^2_{dm}
    ---- inputs ----
    z : redshift
    mh : halo mass in msun, dm+baryon
    ---- output unit ----
    msun^2/Mpc^3
    */

    // printf("Check that the integrated mass converges to m\n");
    double OmM, OmC, m, OmR, OmL, h, pi, rho_cr0, zp, OmMz, d, Delta_C, log10_c, c, delta_c, rv1, rv2, rv3, r_vir, rho_cr, res;

    // Some settings
    OmM = 0.30964168161;
    OmR = 9.1e-5;
    OmC = 0.260667;
    OmL = 0.69026731839;
    h = 0.6766;
    pi = 3.141592653589793;
    rho_cr0 = 2.775e-7 * pow(h, 2.); // critical density in msun/pc^3

    // Pre-requisites
    m = mh * OmC / OmM; // DM mass
    zp = 1. + z;
    OmMz = OmM * pow(zp, 3.) / (OmM * pow(zp, 3.) + OmL);
    d = OmMz - 1.;
    Delta_C = 18. * pow(pi, 2.) + 82. * d - 39. * pow(d, 2.);
    log10_c = 1.071 - 0.098 * (log10(m) - 12.);
    c = pow(10., log10_c) / zp; // concentration, see appdx.A of Zip.et for the additional (1+z) factor
    if (Use_Conde_Concentration)
    {
        c = Halo_Concentration_Conde(mh, z);
        // printf("Using Conde concentration model\n");
    }
    delta_c = Delta_C * pow(c, 3.) / (3. * (log(1. + c) - c / (1. + c)));

    rv1 = 0.784 * pow(m * h / 1.0e8, 1. / 3.);
    rv2 = pow(OmM * Delta_C / (OmMz * 18. * pow(pi, 2.)), -1. / 3.);
    rv3 = (10. / zp) / h * 1000.;
    r_vir = rv1 * rv2 * rv3;                                          // pc
    rho_cr = rho_cr0 * (OmL + OmM * pow(zp, 3.) + OmR * pow(zp, 4.)); // msun/pc^3

    res = pow(rho_cr * delta_c, 2.0);
    res = res * pow(r_vir / c / (1. + c), 3.0);
    res = res / 3.0 * (pow(1. + c, 3.) - 1.0);
    res = res * 1.0e18; // convert to msun^2/Mpc^3

    return res;
}

double Interp_Fast(double *Tab, double xmin, double xmax, int nx, double x)
{
    // Interpolate, x axis must be linear in [xmin xmax] with size nx
    // Speed: 1E-8 seconds per call for nx=100

    double y1, y2, y, dx, x1, x2;
    int id1, id2;

    dx = (xmax - xmin) / ((double)(nx - 1));
    id1 = (int)floor((x - xmin) / dx);
    // Allow some redundencies
    if (id1 < 0)
    {
        if (id1 == -1)
        {
            id1 = 0;
        }
        else
        {
            printf("crash imminent @ Interp_Fast, overflow detected, xmin = %E, xmax = %E, x = %E, id1 = %d\n", xmin, xmax, x, id1);
            LOG_ERROR("Interpolation overflow");
            Throw(ValueError);
        }
    }
    if (id1 > nx - 2)
    {
        if (id1 <= nx)
        {
            id1 = nx - 2;
        }
        else
        {
            printf("crash imminent @ Interp_Fast, overflow detected, xmin = %E, xmax = %E, x = %E, id1 = %d\n", xmin, xmax, x, id1);
            LOG_ERROR("Interpolation overflow");
            Throw(ValueError);
        }
    }
    id2 = id1 + 1;
    x1 = xmin + ((double)id1) * dx;
    x2 = x1 + dx;
    y1 = Tab[id1];
    y2 = Tab[id2];
    y = (y2 - y1) * (x - x1) / (x2 - x1) + y1;

    return y;
}

double Conditional_HMF(double growthf, double M, double Delta, double sigma2, struct CosmoParams *cosmo_params)
{
    // Return conditional HMF dn/dm, unit: 1/(Msun Mpc^3)
    // ---- inputs ----
    // growthf: growth factor
    // M: halo mass in Msun
    // Delta: Overdensity
    // sigma2: sigma2

    double OmegaMh2 = cosmo_params->OMm * pow(cosmo_params->hlittle, 2);
    double r;
    double RhoM = OmegaMh2 * 2.775E11; // Average matter density in Msun/Mpc^3
    float dNdM_old = dNdM_conditional((float)growthf, (float)log(M), 30.0, Deltac, (float)Delta, (float)sigma2);
    double dNdM_double = (double)dNdM_old;
    double RhoM_gird = RhoM * (1 + Delta);

    // This should be positive anyway, check where the - sign comes from
    return fabs(-RhoM_gird * dNdM_double / M / sqrt(2 * PI));
    // printf("Using debug version for Conditional_HMF\n");
    // return fabs(-RhoM * dNdM_double / M / sqrt(2 * PI));

}

double BoostFactor(double z, double growthf, double delta, double sigma2, struct CosmoParams *cosmo_params, int hmf_model)
{
    // Get the boost factor all in one go, set growthf externaly to avoid repeated computations
    int idx;
    double m, dndm, lm1, lm2, dlm, lm, RhoM, RhoCr, dfcoll_dm, Halo, dHalo_dm, RhoM_grid, Ln10, HI, RhoDM, res, fcoll;

    RhoCr = 2.775e11 * pow(cosmo_params->hlittle, 2); // Critical density, msun/Mpc^3
    RhoM = cosmo_params->OMm * RhoCr;
    RhoDM = (cosmo_params->OMm - cosmo_params->OMb) * RhoCr;
    RhoM_grid = hmf_model == -1 ? RhoM * (1. + delta) : RhoM;
    Ln10 = log(10.);

    // mass range
    lm1 = -5.99;
    lm2 = 19.99;
    dlm = (lm2 - lm1) / ((double)Boost_nmh - 1.0);

    fcoll = 0.0;
    Halo = 0.0;

    for (idx = 0; idx < Boost_nmh; idx++)
    {
        lm = lm1 + ((double)idx) * dlm;
        m = pow(10.0, lm);

        // 1 - hmf
        if (hmf_model == -1)
        {
            dndm = Conditional_HMF(growthf, m, delta, sigma2, cosmo_params);
        }
        else if (hmf_model == 0)
        {
            dndm = dNdM(growthf, m);
        }
        else if (hmf_model == 1)
        {
            dndm = dNdM_st(growthf, m);
        }
        else if (hmf_model == 2)
        {
            dndm = dNdM_WatsonFOF(growthf, m);
        }
        else if (hmf_model == 3)
        {
            dndm = dNdM_WatsonFOF_z(z, growthf, m);
        }
        else
        {
            LOG_ERROR("Wrong choice of hmf_model!");
            Throw(ValueError);
        }

        // 2 - fcoll
        dfcoll_dm = m * dndm / RhoM_grid;
        fcoll += m * dfcoll_dm * dlm * Ln10;

        // 3 - Halo
        HI = HaloProfile_Integrator(z, m); // \int dr r^2 \rho^2_{dm}
        dHalo_dm = 4. * PI * dndm * HI / (pow(RhoDM, 2.0) * pow(1. + z, 3.));
        Halo += m * dHalo_dm * dlm * Ln10;
    }

    res = pow(1. - fcoll, 2.) + Halo;
    
    res = fmax(1., res);
    /*
    if (print_debug_info && hmf_model != -1)
    {
        FILE *OutputFile;
        OutputFile = fopen("Fcoll_tmp.txt", "a");
        fprintf(OutputFile, "%f   %E\n", z, fcoll);
        fclose(OutputFile);
    }
    */
    return res;
}

double Read_EFF_Tab(int idxz, int idxE, double *Tab, int UseLog)
{    /* 
    A way to read 2D EFF table, I can't get pointer to work
    ---- inputs ----
    idxz : z index
    idxE : Ek index
    Tab : EFF Table
    UseLog : return log of fc, useful if u want to do subsequent interpolation in log space
    */
    int id;
    double r;
    if ((idxE > Ek_axis_Size - 1) || (idxz > Redshift_Size - 1))
    {
        printf("Error: array index overflow.\n");
        exit(1);
    }
    id = idxE * Redshift_Size + idxz;
    r = Tab[id];
    if (r < 0)
    {
        printf("Error: fc is negative.\n");
        exit(1);
    }
    if (r < 1E-100)
    {
        r=1E-100;
    }
    if (UseLog)
    {
        r=log(r);
    }
    return r;
}

int Find_Sign(double x)
{
    // Find sign of a number
    int result;
    if (x > 0.0)
    {
        result = 1;
    }
    else if (x < 0.0)
    {
        result = -1;
    }
    else
    {
        // maybe x is too close to 0?
        result = 0;
    }
    return result;
}

int Find_Index_1D(double *Tab, double x, int nx)
{
    /* 
    From table Tab with length nx, find closest left element index for x
    return -1 if not in range
    */
    double x1, x2, x3;
    int id1, id2, id3, Stop, s1, s2, s3, idx, count;
    Stop = 0;
    id1 = 0;
    id3 = nx - 1;

    count = 0;
    while (Stop == 0)
    {
        count = count + 1;
        id2 = (int)round((((double)id1) + ((double)id3)) / 2.0);
        if (id3 == id1 + 1)
        {
            idx = id1;
            Stop = 1;
        }

        x1 = Tab[id1];
        x2 = Tab[id2];
        x3 = Tab[id3];

        s1 = Find_Sign(x - x1);
        s2 = Find_Sign(x - x2);
        s3 = Find_Sign(x - x3);
        if (s1 != s2)
        {
            id3 = id2;
        }
        else if (s2 != s3)
        {
            id1 = id2;
        }
        if ((s1 == s3) || (count > 50))
        {
            idx = -1;
            Stop = 1;
        }
    }
    return idx;
}

double Interp_2D(double Ek, double z, double *Tab)
{
    /* 
    Interpolate (in log) from a 2D table
    -- inputs --
        Ek : kinetic energy
        z : redshift
        Tab : EFF table
     */
    int eid1, eid2, zid1, zid2;
    double le, le1, le2, lz, lz1, lz2, f11, f12, f21, f22, f1, f2, F1, F2, f;
    eid1 = Find_Index_1D(Ek_GeV_axis, Ek, Ek_axis_Size);
    zid1 = Find_Index_1D(Redshift_axis, z, Redshift_Size);
    if (zid1 == -1)
    {
        return 0.0;
    }
    if (eid1 == -1)
    {
        printf("Error: Ek not in range, exitting, debugging info:.\n");
        printf("Ek1 = %E, Ek = %E, Ek2 = %E\n", Ek_GeV_axis[eid1], Ek, Ek_GeV_axis[eid1 + 1]);
        printf("min(Ek) = %E, Ek = %E, max(Ek) = %E\n", Ek_GeV_axis[0], Ek, Ek_GeV_axis[Ek_axis_Size - 1]);
        exit(1);
    }

    eid2 = eid1 + 1;
    zid2 = zid1 + 1;

    le = log10(Ek);
    le1 = log10(Ek_GeV_axis[eid1]);
    le2 = log10(Ek_GeV_axis[eid2]);
    lz = log10(1+z);
    lz1 = log10(1+Redshift_axis[zid1]);
    lz2 = log10(1+Redshift_axis[zid2]);

    f11 = Read_EFF_Tab(zid1, eid1, Tab, 1);
    f12 = Read_EFF_Tab(zid1, eid2, Tab, 1);
    f21 = Read_EFF_Tab(zid2, eid1, Tab, 1);
    f22 = Read_EFF_Tab(zid2, eid2, Tab, 1);

    // fix E1
    f1 = f11;
    f2 = f21;
    F1 = (f2 - f1) / (lz2 - lz1) * (lz - lz1) + f1;
    // fix E2
    f1 = f12;
    f2 = f22;
    F2 = (f2 - f1) / (lz2 - lz1) * (lz - lz1) + f1;

    f = (F2 - F1) / (le2 - le1) * (le - le1) + F1;
    f = exp(f);
    return f;
}

double Interp_EFF(double z, double m, int particle, int channel, int HMF, int PS, int UseBoost)
{
    /*
    Interpolate DM deposition efficiencies
    -- inputs --
    z : z
    m : DM mass in GeV
    particle : particle that DM annihilate into
        0 - gamma gamma
        1 - electron + positron
    channel : deposition channel
        0 - HIon
        1 - LyA
        2 - Heat
    HMF : HMF model, follows UserParams
    PS : POWER_SPECTRUM in UserParams
    UseBoost : Whether to use Halo Boost
        0 - no, ignore HMF and PS
        1 - yes
    */
    double Ek, r;
    double EFF_Tab[Redshift_Size*Ek_axis_Size];
    CopyEFF(EFF_Tab, particle, channel, HMF, PS, UseBoost);
    Ek = (particle==1) ? m-0.5109989E-3 : m;
    if ((HMF < 0) || (HMF > 3))
    {
        printf("Error: unsupported HMF, HMF=%d\n", HMF);
        exit(1);
    }
    if ((PS < 0) || (PS > 4))
    {
        printf("Error: unsupported PS, PS=%d\n", PS);
        exit(1);
    }
    if (Ek < 0.0)
    {
        printf("Error: Ek is negative, Ek=%E\n", Ek);
        exit(1);
    }
    r = Interp_2D(Ek, z, EFF_Tab);
    return r;
}

double PeeblesFactor(double z, double xe, double Tk, double Delta)
{
    /*
    Compute Peebles C factor following Dodelson's Modern Cosmology II (see p108)
    Analytic behavior of C:
    xe->1: C->1
    xe->0: None
    T->0: C->1
    T->inf: C->0
    */
    double E0, alpha_fsc, hbar, LightSpeed, me, np0, x, fx1, fx2, fx, beta2, pi, L2y, nH, H, La1, La2, La, result, xH;
    alpha_fsc = 1/137.03599976;
    E0 = 13.6 * 1.602176634E-19;
    hbar = 1.0545718002693302e-34;
    LightSpeed = 299792458.0;
    me = 9.109383632044565e-31;
    np0 = 0.1901567053460595;// Number density of H nuclei today
    pi = 3.14159265358979323846264338327;
    
    x = E0/(Tk * 1.38064852E-23);
    fx1 = log(x);
    fx2 = x * exp(x/4.0);
    fx = fx1/fx2;
    beta2 = 9.78 * pow(E0/(2 * pi), 1.5) * pow(alpha_fsc, 2.0) /(pow(me, 0.5) * LightSpeed * hbar) * fx;
    L2y = 8.227;
    xH = xe > 1.0 ? 0: 1-xe; // Avoid numerical overflow
    nH = np0 * pow(1.0+z, 3.0) * xH * (1+Delta);
    //H = 2.192695336552484e-18 * pow(0.69026731839 + 0.30964168161*pow(1.0+z, 3.0) + 9.1E-5 * pow(1.0+z, 4.0), 0.5);
    H = hubble(z);
    La1 = H * pow(3.0*E0, 3.0);
    La2 = nH * pow(8*pi, 2.0) * pow(hbar*LightSpeed, 3.0);
    La = La1/La2;
    // result = (La + L2y)/(La + L2y + beta2); // This can give numerical issue when xe->0
    result = 1 - beta2/(La + L2y + beta2); 
    // C->1 if xe->1, but this is kinda problemtic because if T -> inf then we instead expect C->0, however in fact T can not go to such extremes so in the end we should be safe
    if (xe > 1 - 1.0E-50)
    {
        result = 1.0;
    }
    if ((result > 1) || (result < 0))
    {
        result = 0.0/0.0;//Sanity check, if not within 0-1 then give NaN
    }
    return result;
}

double Compute_JLyA(double z, double dEdVdt_LyA)
{
    /*
    Compute J_LyA from DM following 1603.06795 and 10.1088/1475-7516/2024/01/005
    */
    double H, LightSpeed, nu_LyA, Q, ELyA, r;
    Q = 1.602176634E-19;
    LightSpeed = 299792458.0;
    nu_LyA = 2.466349E+15; // Frequency in Hz of a LyA photon with energy 10.2 eV
    H = hubble(z);
    ELyA = 10.2 * Q; // LyA energy 10.2 eV in J
    r = LightSpeed * dEdVdt_LyA / (4.0 * PI * H * nu_LyA * ELyA) * 1.0E-4; //1.0E-4 converts m^-2 to cm^-2
    if (isfinite(r) == 0)
    {
        printf("DarkSide: JLyA is infinite or nan, crash imminent, J_LyA = %E\n", r);
    }
    // Test the impact of LyA, give noisy shout out to make sure user knows this
    // r = 0.0;
    // printf("Test version with JLyA = 0, current z = %f\n", z);
    return r;
}

double EoR_Rate_DM(double z, double Boost, double Delta, struct AstroParams *astro_params, struct FlagOptions *flag_options, struct UserParams *user_params, double xe, double Tk, double dt_dzp, int GetIon)
{
    /*
    Get dxe/dz and dT/dz from DM injection
    Doing everything in SI units
    IMPORTANT NOTE: 21cmFAST defines xe as ne/(nH + nHe)
    ---- inputs ----
    GetIon : choose output
        0 - dT/dt
        1 - dxe/dt
        2 - J_LyA
    */
    double RhoC_c2_2, Pann, Q, P27, dEdVdt_Inj, nH, dxe_dt, dT_dt, fHe, dxe_dz, dT_dz, B, ntot, kB, f_HIon, f_Heat, f_LyA, Peebles, J_LyA;
    if (flag_options->USE_HALO_BOOST)
    {
        B = Boost;
    }
    else
    {
        B = 1.0;
    }

    // RhoC_c2_2 = [OmC * Rho_cr *c^2]^2
    RhoC_c2_2 = 4.0610249075234377e-20;
    Q = 1.602176634E-19;
    kB = 1.38064852E-23;
    fHe = 0.08112582781456953;

    P27 = (double)astro_params->Pann27;
    Pann = P27 * 1.0E-27 * 1.0E-6 / (1.0E9 * Q); // m^3/s/J

    dEdVdt_Inj = B * Pann * RhoC_c2_2 * pow(1.0 + z, 6.0); // Again this is in SI unit

    nH = 0.19015670534605955 * pow(1.0 + z, 3.0) * (1.0 + Delta);
    Peebles = PeeblesFactor(z, xe, Tk, Delta);
    // First do ionisation
    if (user_params->DM_Dep_Method == 0)
    {// SSCK
        f_HIon = (1.0 - xe) / 3.0;
        f_LyA = (1.0 - xe) / 3.0;
        f_Heat = (1.0 + 2 * xe) / 3.0;
    }
    else
    {// Tracy, technically in doing so we are missing the feedbacks and inhomogeneuity. To avoid unphysical xe>1 case I am adding a calibration factor of (1-xe) to HIon and LyA
        f_HIon = (1.0 - xe)*Interp_EFF(z, astro_params->mdm, user_params->DM_ANN_Channel, 0, user_params->HMF, user_params->POWER_SPECTRUM, flag_options->USE_HALO_BOOST);
        f_LyA = (1.0 - xe)*Interp_EFF(z, astro_params->mdm, user_params->DM_ANN_Channel, 1, user_params->HMF, user_params->POWER_SPECTRUM, flag_options->USE_HALO_BOOST);
        f_Heat = Interp_EFF(z, astro_params->mdm, user_params->DM_ANN_Channel, 2, user_params->HMF, user_params->POWER_SPECTRUM, flag_options->USE_HALO_BOOST);
    }
    // This is somewhat inaccurate, we are not doing Helium here but we are assuming that xe is shared by H and He. should mention this in our paper
    dxe_dt = f_HIon * dEdVdt_Inj / (13.6 * Q * nH * (1+fHe)) + (1 - Peebles) * f_LyA * dEdVdt_Inj/ (10.2 * Q * nH * (1+fHe));

    ntot = nH * (1. + xe + fHe + xe * fHe); // Neutral H and ionized H (proton), e from H, neutral and neutral He, e from He, note that He is singly ionized
    dT_dt = f_Heat * dEdVdt_Inj * 2.0 / (3.0 * kB * ntot);
    
    dxe_dz = dxe_dt * dt_dzp;
    dT_dz = dT_dt * dt_dzp;

    if (GetIon == 1)
    {
        return dxe_dz;
    }
    else if (GetIon == 0)
    {
        return dT_dz;
    }
    else if (GetIon == 2)
    {
        J_LyA = Compute_JLyA(z, f_LyA * dEdVdt_Inj);
        return J_LyA;
    }
    else
    {
        LOG_ERROR("Error: wrong choice of GetIon, must be [0, 1, 2]!");
        Throw(ValueError);
    }
}
