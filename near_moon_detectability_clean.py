import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
import pickle


#define whether to save plots or outputs
save_plots = False
save_outs = False

#define cams
cams = ['CMOS', 'SWIR']

#define colourbar variables
cmos_vmin = 10
cmos_vmax = 18
swir_vmin =  8
swir_vmax = 15
cmos_cmap = 'inferno'
swir_cmap = 'viridis'
cmos_cmap_steps = 16
swir_cmap_steps = 14

#define arrays in which to store regime one and two data
df1s=[]
df2s=[]


#loop through cams
for cam in cams:

    #read in ZP and background mean vs rms fits
    data = np.load(cam+'_fit_params.npz')
    zp = data['zp']
    rms_popt = data['popt']
    rms_pcov = data['pcov']

    #print ZP and fit params
    print(f'ZP = {np.round(zp,3)}')
    print(f'rms_popt = {rms_popt}')
    print(f'rms_pcov = {rms_pcov}')

    #for SWIR, also read in background cap values    
    if cam=='SWIR':
        data = np.load(cam+'_param_caps.npz')
        bg_mean_cap = data['bg_mean']
        bg_rms_cap = data['bg_rms']
        
        #print background cap values
        print(f'bg_mean cap = {np.round(bg_mean_cap,3)}')
        print(f'bg_rms cap = {np.round(bg_rms_cap,3)}')
    

    #define array of moon illumination and separation values to test    
    illums_fine = np.linspace(20,100,40) #%
    fovs_fine = np.linspace(2.5,15,26) #degrees

    #create a 2D array of combinations to test    
    illums = np.repeat(illums_fine, len(fovs_fine))
    fovs = np.tile(fovs_fine, int(len(illums)/len(fovs_fine)))
    fine_vals = np.array([illums, fovs])
    dense_points = np.stack([fine_vals[0].ravel(), fine_vals[1].ravel()], -1)

    #define functions used to convert background mean to rms
    #SWIR
    def rms_func1(x, a,b,c,d,e):
        return a + b*x + c*x**2 + d*x**3 + e*x**4
    #CMOS
    def rms_func2(x, a,b,c):
            return a*x**b + c
    
    #read in interpolation fits
    #background mean
    with open(cam+'_bg_mean_rbf_spline.pkl', 'rb') as f:
        #read fit
        bg_mean_rbf = pickle.load(f)
        #generate background mean values for each illumination/separation combination
        bg_mean_rbf_data = 10**bg_mean_rbf(dense_points).reshape(fine_vals[0].shape)
        #for SWIR, cap background mean values
        if cam=='SWIR':
            bg_mean_rbf_data = [min(x, bg_mean_cap) for x in bg_mean_rbf_data]
    #background rms
    with open(cam+'_bg_rms_rbf_spline.pkl', 'rb') as f:
        #read fit
        bg_rms_rbf = pickle.load(f)
        #generate background rms values for each illumination/separation combination
        bg_rms_rbf_data = 10**bg_rms_rbf(dense_points).reshape(fine_vals[0].shape)
        #for SWIR, cap background rms values
        if cam=='SWIR':
            bg_rms_rbf_data = [min(x, bg_rms_cap) for x in bg_rms_rbf_data]
        #overwrite fitted background rms values with ones scaled from background mean
        if cam=='CMOS':
            bg_rms_rbf_fit_data = rms_func2(bg_mean_rbf_data, *rms_popt)
        elif cam=='SWIR':
            bg_rms_rbf_fit_data = [10**rms_func1(np.log10(d), *rms_popt) for d in bg_mean_rbf_data]
    
    #surface plots for background stats as functions of moon illumination and separation
    for rbf_data in [bg_mean_rbf_data, bg_rms_rbf_fit_data]:
        x = illums_fine
        y = fovs_fine
        X, Y = np.meshgrid(x, y)
        Z = np.asarray(rbf_data).reshape(len(x),len(y))
        Z = np.transpose(Z)
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(X, Y, Z, cmap=cm.viridis)
        ax.set_zscale('log')
        ax.set_xlabel('moon illumination (%)')
        ax.set_ylabel('moon separation (deg)')
        ax.set_zlabel('value')
        plt.show()
    

    #define FWHM, pixel scale (camera specific) and pixel buffer
    if cam=='CMOS':
        fwhm=6.01
        pix_scale = 0.98
    elif cam=='SWIR':
        fwhm=3.55
        pix_scale = 3.12
    buff = int(fwhm*3)
    
    
    

    #define arrays of moon illumination and separation values to test
    illums = np.linspace(20, 100, 100) #%
    seps = np.linspace(2.5, 15, 100) #degrees
    #define arrays of target rates and number of frames to test    
    rates = np.linspace(0, 50, 100) #''/s
    n_frames = np.linspace(1, 50, 100)

    #define exposure time and signal-to-noise ratio limit    
    exp_time = 1
    SNR_lim = 10

    #initialise output arrays    
    ii1=[]
    ss1=[]
    rr=[]
    mm1=[]
    ii2=[]
    ss2=[]
    nn=[]
    mm2=[]

    #loop through illuminations
    for i in illums:
        #loop through separations
        for s in seps:
            #calculate background mean
            bg_mean = 10**bg_mean_rbf([[i,s]])[0]
            #if using SWIR, cap the value
            if cam=='SWIR':
                bg_mean = min(bg_mean, bg_mean_cap)
            #calulate the background rms, by converting from the background mean
            if cam=='CMOS':
                bg_rms = rms_func2(bg_mean, *rms_popt)
            elif cam=='SWIR':
                bg_rms = 10**rms_func1(np.log10(bg_mean), *rms_popt)
            #cycle through target rates
            for r in rates:
                #convert rate from ''/s to pixels/s
                r2 = r/pix_scale
                #calculate the flux of a target, based on SNR, aperture area and target rate
                req_flux = (SNR_lim**2 + np.sqrt(SNR_lim**4 + 4*buff*bg_rms**2*SNR_lim**2*(r2+buff/2))) / 2
                #convert flux to magnitude
                m = zp - 2.5 * np.log10(req_flux/exp_time)
                #save results to arrays
                ii1.append(i)
                ss1.append(s)
                rr.append(r)
                mm1.append(m)
            #cycle through number of frames
            for n in n_frames:
                #set target rate to 0
                r = 0
                #convert SNR limit, to account for multiple frames
                SNR = SNR_lim/np.sqrt(n)
                #calculate the flux of a target, based on SNR, aperture area and number of frames
                req_flux = (SNR**2 + np.sqrt(SNR**4 + 4*buff*bg_rms**2*SNR**2*(r+buff/2))) / 2
                #convert flux to magnitude
                m = zp - 2.5 * np.log10(req_flux/exp_time)            
                #save results to arrays
                ii2.append(i)
                ss2.append(s)
                nn.append(n)
                mm2.append(m)
    
    
    #create output dataframe under regime 1 (limiting magnitude as a function of moon illumination, moon separation and target rate)
    df = pd.DataFrame({'moon_illum' : ii1, 'moon_sep' : ss1, 'target_rate' : rr, 'limiting_mag' : mm1})
    #save to file
    if save_outs:
        df.to_csv(cam+'_near_moon_detectability_out1.csv', index=False)

    #plot a data cube of the results, colour points by limiting magnitude
    #use a fixed colourbar scale for CMOS or SWIR
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    if cam=='CMOS':
        a = ax.scatter(df['moon_illum'], df['moon_sep'], df['target_rate'], c=df['limiting_mag'], cmap=plt.get_cmap(cmos_cmap, cmos_cmap_steps), vmin=cmos_vmin, vmax=cmos_vmax)
    elif cam=='SWIR':
        a = ax.scatter(df['moon_illum'], df['moon_sep'], df['target_rate'], c=df['limiting_mag'], cmap=plt.get_cmap(swir_cmap, swir_cmap_steps), vmin=swir_vmin, vmax=swir_vmax)
    ax.xaxis.set_inverted(True)
    ax.yaxis.set_inverted(True)
    ax.set_xlabel('$m_i$ (%)')
    ax.set_ylabel('$m_s$ ($^\circ$)')
    ax.zaxis.set_rotate_label(False)
    ax.set_zlabel("$r$\n($''$/s)", rotation='horizontal')
    ax.view_init(elev=30, azim=45)
    if cam=='CMOS':
        plt.colorbar(a, label='$M$ ($G$-band)')
    elif cam=='SWIR':
        plt.colorbar(a, label='$M$ ($J$-band)')
    ax.plot((20,20), (15,2.5), (50,50), c='gray', zorder=1000, lw=2, alpha=0.75)
    ax.plot((20,20), (2.5,2.5), (50,0), c='gray', zorder=1000, lw=2, alpha=0.75)
    ax.plot((20,100), (2.5,2.5), (50,50), c='gray', zorder=1000, lw=2, alpha=0.75)
    plt.tight_layout()
    #save to file
    if save_plots:
        plt.savefig(f'plots/{cam}_target_rate.png')
    plt.show()

    #print min and max limiting magnitude values
    print(f"{cam}, as function of target rate:")
    print(f"\tmax mag: {max(df['limiting_mag'])}")
    print(f"\tmin mag: {min(df['limiting_mag'])}")
    
    #create output dataframe under regime 2 (limiting magnitude as a function of moon illumination, moon separation and number of frames)
    df = pd.DataFrame({'moon_illum' : ii2, 'moon_sep' : ss2, 'n_frames' : nn, 'limiting_mag' : mm2})
    #save to file
    if save_outs:
        df.to_csv(cam+'_near_moon_detectability_out2.csv', index=False)
    
    #plot a data cube of the results, colour points by limiting magnitude
    #use a fixed colourbar scale for CMOS or SWIR
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    if cam=='CMOS':
        a = ax.scatter(df['moon_illum'], df['moon_sep'], df['n_frames'], c=df['limiting_mag'], cmap=plt.get_cmap(cmos_cmap, cmos_cmap_steps), vmin=cmos_vmin, vmax=cmos_vmax)
    elif cam=='SWIR':
        a = ax.scatter(df['moon_illum'], df['moon_sep'], df['n_frames'], c=df['limiting_mag'], cmap=plt.get_cmap(swir_cmap, swir_cmap_steps), vmin=swir_vmin, vmax=swir_vmax)
    ax.xaxis.set_inverted(True)
    ax.yaxis.set_inverted(True)
    ax.set_xlabel('$m_i$ (%)')
    ax.set_ylabel('$m_s$ ($^\circ$)')
    ax.zaxis.set_rotate_label(False)
    ax.set_zlabel('$n$', rotation='horizontal')
    ax.view_init(elev=30, azim=45)
    if cam=='CMOS':
        plt.colorbar(a, label='$M$ ($G$-band)')
    elif cam=='SWIR':
        plt.colorbar(a, label='$M$ ($J$-band)')
    ax.plot((20,20), (15,2.5), (50,50), c='gray', zorder=1000, lw=2, alpha=0.75)
    ax.plot((20,20), (2.5,2.5), (50,0), c='gray', zorder=1000, lw=2, alpha=0.75)
    ax.plot((20,100), (2.5,2.5), (50,50), c='gray', zorder=1000, lw=2, alpha=0.75)
    plt.tight_layout()
    #save to file
    if save_plots:
        plt.savefig(f'plots/{cam}_n_frames.png')
    plt.show()
    
    #print min and max limiting magnitude values
    print(f"{cam}, as function of n_frames:")
    print(f"\tmax mag: {max(df['limiting_mag'])}")
    print(f"\tmin mag: {min(df['limiting_mag'])}")




#define result dataframes for each combo of cam and regime
cmos_rate = df1s[0]
swir_rate = df1s[1]
cmos_frames = df2s[0]
swir_frames = df2s[1]

#list values of r and n to use for cross-sections
rates = [50,40,30,20,10,0]
frames = [50,40,30,20,10,1]
#find closest values in list of tested r and n values
rates = [min(np.linspace(0, 50, 100), key=lambda x:abs(x-myNumber)) for myNumber in rates]
frames = [min(np.linspace(1, 50, 100), key=lambda x:abs(x-myNumber)) for myNumber in frames]

#create 6X4 array of subfigs
fig, axes = plt.subplots(nrows=6, ncols=4, sharex=True, sharey=True, figsize=(7.75,7*4/3))

#for each subfigure select correct dataframe and colourbar variables
for i, ax in enumerate(axes.flatten()):
    if i%4==0:
        df = cmos_rate
        vmin, vmax = cmos_vmin, cmos_vmax
        cmap, cmap_steps = cmos_cmap, cmos_cmap_steps
    elif i%4==1:
        df = cmos_frames
        vmin, vmax = cmos_vmin, cmos_vmax
        cmap, cmap_steps = cmos_cmap, cmos_cmap_steps
    elif i%4==2:
        df = swir_rate
        vmin, vmax = swir_vmin, swir_vmax
        cmap, cmap_steps = swir_cmap, swir_cmap_steps
    elif i%4==3:
        df = swir_frames
        vmin, vmax = swir_vmin, swir_vmax
        cmap, cmap_steps = swir_cmap, swir_cmap_steps

    #filter data to relevant cross-section only
    try:
        rate = rates[i//4]
        df = df.loc[(df['target_rate'] == rate)].reset_index(drop=True)
    except:
        frame = frames[i//4]
        df = df.loc[(df['n_frames'] == frame)].reset_index(drop=True)

    #plot data
    if i%4<2:
        a = ax.scatter(df['moon_sep'], df['moon_illum'], c=df['limiting_mag'], s=0.1, cmap=plt.get_cmap(cmap, cmap_steps), vmin=vmin, vmax=vmax)
    else:
        b = ax.scatter(df['moon_sep'], df['moon_illum'], c=df['limiting_mag'], s=0.1, cmap=plt.get_cmap(cmap, cmap_steps), vmin=vmin, vmax=vmax)

    #set shared x/y labels
    ax.set_xlim(2.5,15)
    ax.set_ylim(20,100)
    ax.xaxis.set_inverted(True)
    if i%4==0:
        ax.set_ylabel('$m_i$ (%)')
        ax.set_xticks([2.5, 6.7, 10.8, 15.0])
    if i//4==5:
        ax.set_xlabel('$m_s$ ($^\circ$)')
        ax.set_yticks([20, 60, 100])
    ax.set_title(['CMOS  ' if i%4<2 else 'SWIR  '][0] + ['$r='+str(round(rate))+"\,''$/s" if (i%4)%2==0 else '$n='+str(round(frame))+'$'][0])

#adjust subfig spacing to fit colourbars (one for cmos, one for swir)
fig.subplots_adjust(right=0.75)
cbar_ax = fig.add_axes([0.8, 0.1, 0.03, 0.8])
fig.colorbar(a, cax=cbar_ax, label='$M$ ($G$-band)')
cbar_ax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
fig.colorbar(b, cax=cbar_ax, label='$M$ ($J$-band)')
plt.subplots_adjust(left=0.1, right=0.75, top=0.96, bottom=0.06, hspace=0.35, wspace=0.2)
#save to file
if save_plots:
    plt.savefig('plots/2d_slices.png')
plt.show()


