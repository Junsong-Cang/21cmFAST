reload_boost = 0
reload_EFF = 1
print_data = 1

'''
Compute boost factor and deposition efficiency from DM annihilation, this largely follows 10.1088/1475-7516/2022/03/012
Pipeline:
    Compute Boost factors with GetBoost
    Compute EFF with GetEFF which calls GetEFF_Kernel
    main prints out interpolation table and functions to c_products/Tables.h
    Move c_products/Tables.h to 21cmFAST src
'''
import numpy as np
import h5py
from PyLab import Hubble, TimeNow, os
import py21cmfast as p21c
from joblib import Parallel, delayed
main_path = '/Users/cangtao/FileVault/GitHub/21cmFAST/src/py21cmfast/DarkSide/'
if __name__ != '__main__' and (reload_boost or reload_EFF or print_data):
    raise Exception('You are importing this as a lib, in this case thou shalt not reload anything')

def GetBoost(HMF=0, POWER_SPECTRUM=2):
    # Compute Boost factor
    zmin = 4
    zmax = 100
    LC_Quantities = ('brightness_temp','Ts_box','xH_box','Tk_box','Boost_box')
    GLB_Quantities = ('brightness_temp','Ts_box','xH_box','Tk_box', 'Boost_box')

    user_params = p21c.UserParams(
        HII_DIM = 50,
        N_THREADS = 1,
        USE_INTERPOLATION_TABLES = True,
        HMF = HMF,
        POWER_SPECTRUM = POWER_SPECTRUM,
        BOX_LEN = 200)
    astro_params = p21c.AstroParams(
        Pann27 = 1,
        NU_X_THRESH = 500.0)
    flag_options = p21c.FlagOptions(
        USE_MASS_DEPENDENT_ZETA = True,
        INHOMO_RECO = True,
        USE_TS_FLUCT = True,
        USE_HALO_BOOST = True,
        INHOMO_HALO_BOOST = False)
    t1 = TimeNow()
    with p21c.global_params.use(Z_HEAT_MAX = zmax):
        lc = p21c.run_lightcone(
            redshift=zmin,
            max_redshift=zmax,
            astro_params=astro_params, 
            flag_options=flag_options,
            user_params = user_params,
            lightcone_quantities=LC_Quantities,
            global_quantities=GLB_Quantities)
        z = lc.node_redshifts
        b = lc.global_Boost
    t2 = TimeNow()
    print('Time used:', t2 - t1)
    return z, b

if reload_boost:
    HMF = [0, 1, 2, 3]
    PS = [0, 1, 2, 3, 4] # PS=5 doesn't work because the k range is not enough for our small halo mass down to 1E-6
    
    params = []
    for idh, hmf in enumerate(HMF):
        for idp, ps in enumerate(PS):
            params.append([hmf, ps])
    params = np.array(params)
    def kernel(idx):
        param = params[idx]
        hmf, ps = param[0], param[1]
        z, b = GetBoost(HMF=hmf, POWER_SPECTRUM=ps)
        np.savez('tmp_'+str(idx)+'.npz', z=z, b=b)
    swap = Parallel(n_jobs = 12)(delayed(kernel)(idx) for idx in np.arange(0, len(params)))
    lh, lp = len(HMF), len(PS)
    z = np.load('tmp_0.npz')['z'][::-1]
    nz = len(z)
    B = np.zeros([lh, lp, nz])
    idx = 0
    for idh in np.arange(0, lh):
        for idp in np.arange(0, lp):
            #B[idh, idp, :] = np.load('tmp_0.npz')['b'][::-1] # Previous buggy version, this is so fucked up
            B[idh, idp, :] = np.load('tmp_'+str(idx)+'.npz')['b'][::-1]
            idx = idx + 1
    np.savez(main_path + 'data/BoostTemplates.npz', z = z, HMF = HMF, PS = PS, B = B)
    os.system('rm tmp_*.npz')

def Load_Data():
    '''
    Load required datasets
    Note that axis doc in TransferFunctions.h5 are written for MatLab
    '''

    # 1. Boost factor
    BoostData = np.load(main_path + 'data/BoostTemplates.npz')
    z = BoostData['z']
    HMF = BoostData['HMF']
    PS = BoostData['PS']
    B = BoostData['B']
    
    # 2. Transfer function
    Transfer_File = main_path + 'data/TransferFunctions.h5'
    f = h5py.File(Transfer_File, 'r')
    TF = f['T_Array'][:]
    TF_z = f['Axis/z'][:]
    TF_E = f['Axis/Kinetic_Energy_GeV'][:]
    f.close()
    LnZp = np.log(1+TF_z)
    dLnZp = LnZp[0:len(LnZp)-1] - LnZp[1:len(LnZp)]
    dLnZp = np.sum(dLnZp)/len(dLnZp)
    
    # Organize data
    BoostData = {'z': z, 'B':B, 'HMF' : HMF, 'POWER_SPECTRUM':PS}
    TransferData = {'z': TF_z, 'E': TF_E, 'TF': TF, 'dLnZp': dLnZp}
    r = {'BoostData': BoostData, 'TransferData': TransferData}
    # print(TF.shape)
    # print('dLnZp=', dLnZp)
    return r

DataSet = Load_Data()

def GetEFF_Kernel(idxE=20, idxz=0, spec = 'e', channel='HIon', HMF = 0, PS=2, UseBoost = 1):
    '''
    Compute deposition efficiency for annihilating DM
    -- inputs --
    channel:
        HIon
        LyA
        Heat
    '''
    CalibrateB = 0
    # Energy is in GeV
    m = 0.5109989E-3 if spec == 'e' else 0
    Ek = DataSet['TransferData']['E'][idxE]
    z = DataSet['TransferData']['z'][idxz]
    B = 10**np.interp(x=z, xp = DataSet['BoostData']['z'], fp = np.log10(DataSet['BoostData']['B'][HMF, PS, :]), left=np.nan, right=0)
    H = Hubble(z)
    sid = 0 if spec == 'y' else 1
    if channel == 'HIon':
        cid = 0
    elif channel == 'LyA':
        cid = 2
    elif channel == 'Heat':
        cid = 3
    else:
        raise Exception('Wrong choice of channel')
    if not UseBoost: B = 1
    Prefix = H * Ek / (B * (m + Ek) * (1+z)**3)
    #Prefix = H * Ek / (B * ( Ek) * (1+z)**3) # I think this what's used in HyRec, it's the fraction injected KINETIC energy
    
    # Now the next bit
    zi = DataSet['TransferData']['z']
    
    Bi = 10**np.interp(x=zi, xp = DataSet['BoostData']['z'], fp = np.log10(DataSet['BoostData']['B'][HMF, PS, :]), left=np.nan, right=0)
    if CalibrateB == 1:
        # By logic of 0th order integration, using either the left or right points are both ok
        dx = DataSet['TransferData']['dLnZp']
        x0 = np.log(1+zi)
        x1 = x0+dx
        zi_ = np.exp(x1)
        Bi = 10**np.interp(x=zi_, xp = DataSet['BoostData']['z'], fp = np.log10(DataSet['BoostData']['B'][HMF, PS, :]), left=np.nan, right=0)
    if not UseBoost: Bi = 1
    Hi = Hubble(zi)
    T = DataSet['TransferData']['TF'][:, idxE, idxz, cid, sid]
    fc = Prefix*np.sum((1+zi)**3 * Bi * T / Hi)
    return fc

def GetEFF():
    '''
    Compute all EFF fc, format: fc[UseBoost, Species, Channel, HMF, PS, Energy, redshift]
    '''
    UseBoost = [0, 1]
    spec = ['y', 'e']
    channel = ['HIon', 'LyA', 'Heat']
    HMF = [0, 1, 2, 3]
    PS = [0, 1, 2, 3, 4]
    Ek = DataSet['TransferData']['E']
    z = DataSet['TransferData']['z']
    nb = len(UseBoost)
    ns = len(spec)
    nc = len(channel)
    nh = len(HMF)
    lp = len(PS)
    ne = len(Ek)
    nz = len(z)
    
    fc = np.zeros([nb, ns, nc, nh, lp, ne, nz])
    t1 = TimeNow()
    for idb, UseBoost_ in enumerate(UseBoost):
        for ids, spec_ in enumerate(spec):
            for idc, channel_ in enumerate(channel):
                for idh, HMF_ in enumerate(HMF):
                    for idp, PS_ in enumerate(PS):
                        for ide in np.arange(0, ne):
                            for idz in np.arange(0, nz):
                                fc[idb, ids, idc, idh, idp, ide, idz] = GetEFF_Kernel(
                                    idxE=ide, idxz=idz, spec=spec_, channel=channel_, HMF=HMF_, PS=PS_, UseBoost=UseBoost_)
    t2 = TimeNow()
    print('Time used:', t2 - t1)
    np.savez(main_path + 'data/EEE_Table.npz', fc = fc, HMF=HMF,PS=PS, Ek=Ek, z=z)

if reload_EFF:GetEFF()

def Load_EFF_Data(spec = 'e', channel='HIon', HMF = 0, PS=2, UseBoost = 0, process = 'ann'):
    '''
    Read pre-computed EFF Table
    output:
        name: name of dataset, to be used for c. Example: EFF_ann_e_HIon_HMG
        fc: EFF, format fc[Energy, redshift]
    '''
    # First get name
    s1 = 'EFF_'+process+'_'+spec+'_'+channel
    s2 = '_'+str(HMF) + str(PS) if UseBoost else '_HMG'
    name = s1 + s2
    # Now EFF
    sid = 0 if spec=='y' else 1
    if channel == 'HIon':
        cid = 0
    elif channel == 'LyA':
        cid = 1
    elif channel == 'Heat':
        cid = 2
    else:
        raise Exception('Wrong channel')
    fc = np.load(main_path + 'data/EEE_Table.npz')['fc'][UseBoost, sid, cid, HMF, PS, :, :]
    return name, fc

def PrintData_Kernel(name, fc, filename):
    '''
    Print EFF to a file, example:
    double EFF_ann_e_HIon_HMG[2520] = {
      9.7416E-04, 9.5136E-04, ...};
    '''
    fc_shape = np.shape(fc)
    fc_size = str(int(fc_shape[0] * fc_shape[1]))
    ne = fc_shape[0]
    nz = fc_shape[1]
    file = open(filename, 'a')
    head = 'double '+name+'['+fc_size+'] = {'
    print(head, file=file)
    for ide in np.arange(0, ne):
        s = '  '
        for idz in np.arange(0, nz):
            f = fc[ide, idz]
            f_str = "{0:.4E}".format(f)
            f_str = f_str + '};' if ide==ne-1 and idz==nz-1 else f_str+', '
            s+=f_str
        print(s, file=file)
    file.close()

def Print_C_funtion(filename):
    '''
    Called after printing (and defining) table data, print c function that copies an EFF table to another array, the latter is used
    for interpolation based on user input of z, E, PS, HMF, channel, etc
    '''
    cmd = 'cat c_products/swap.c >> '+filename # defines the function that copies array
    os.system(cmd)
    file = open(filename, 'a')
    specs = ['y', 'e']
    channels = ['HIon', 'LyA', 'Heat']
    print('void CopyEFF(double *Tab, int particle, int channel, int HMF, int PS, int UseBoost)', file=file)
    print('{ // Fill Tab with EFF', file=file)
    # print('  double EFF[Redshift_Size*Ek_axis_Size];', file=file)
    # raise Exception('Forgot UseBoost option in c, need to update')
    # HMG
    print('  if (UseBoost == 0)', file=file)
    print('  {', file=file)
    for spec in [0, 1]:
        for channel in [0, 1, 2]:
            print('    if ((particle == '+str(spec)+')&&'+'(channel=='+str(channel)+'))', file=file)
            print('    {', file=file)
            name, fc = Load_EFF_Data(spec=specs[spec], channel=channels[channel], UseBoost=0)
            print('      CopyArray(Tab, '+name+');', file=file)
            print('    }', file=file)
            # print(name)
    print('  }', file=file)
    print('  else', file = file)
    print('  {', file=file)
    for spec in [0, 1]:
        for channel in [0, 1, 2]:
            for HMF in [0, 1, 2, 3]:
                for PS in [0, 1, 2, 3, 4]:
                    print('    if ((particle == '+str(spec)+')&&(channel=='+str(channel)+
                          ')&&(HMF=='+str(HMF)+')&&(PS=='+str(PS)+'))',
                          file=file)
                    print('    {', file=file)
                    name, fc = Load_EFF_Data(spec=specs[spec], channel=channels[channel], UseBoost=1, HMF=HMF, PS=PS)
                    print('      CopyArray(Tab, '+name+');', file=file)
                    print('    }', file=file)
    print('  }', file=file)
    print('}', file=file)
    
    file.close()

def main(filename = main_path + 'c_products/Tables.h'):
    '''
    The main function, load and print all tables and functions
    '''
    EFF = np.load('data/EEE_Table.npz')
    Ek_GeV = EFF['Ek']
    z = EFF['z']
    # Leave some redundencies
    z[0] = z[0] * 0.999
    z[-1] = z[-1] * 1.001
    Ek_GeV[0] = Ek_GeV[0] * 0.999
    Ek_GeV[-1] = Ek_GeV[-1] * 1.001
    
    nz = len(z)
    ne = len(Ek_GeV)
    file = open(filename, 'w')
    print('// EFF data, this is not meant to be human-readable', file=file)
    print('#define Redshift_Size    '+str(nz), file=file)
    print('#define Ek_axis_Size    '+str(ne), file=file)
    # Print z
    s = 'double Redshift_axis[Redshift_Size] = {'
    for idz, z_ in enumerate(z):
        s += "{0:.4E}".format(z_)
        if idz != nz-1: s += ','
    s += '};'
    print(s, file=file)
    
    # Print Ek
    s = 'double Ek_GeV_axis[Ek_axis_Size] = {'
    for ide, Ek_ in enumerate(Ek_GeV):
        s += "{0:.4E}".format(Ek_)
        if ide != ne-1: s += ','
    s += '};'
    print(s, file=file)
    file.close()

    # Print EFF
    for spec in ['e', 'y']:
        for channel in ['HIon', 'LyA', 'Heat']:
            # BoostLess
            name, fc = Load_EFF_Data(spec=spec,channel=channel, UseBoost=0)
            PrintData_Kernel(name=name,fc=fc,filename=filename)
            for HMF in [0, 1, 2, 3]:
                for PS in [0, 1, 2, 3, 4]:
                    name, fc = Load_EFF_Data(spec=spec,channel=channel, HMF=HMF, PS=PS, UseBoost=1)
                    PrintData_Kernel(name=name,fc=fc,filename=filename)
    
    # Now that we printed table definitions, let's print functions that prepares them for interpolation
    Print_C_funtion(filename)

if print_data:main()
