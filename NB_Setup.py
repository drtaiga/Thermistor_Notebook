#Import the schemdraw modules for drawing circuit schematics
import schemdraw.elements as elm
import schemdraw

import numpy as np
import matplotlib.pyplot as plt

from ipywidgets import interact, interactive, fixed, interact_manual
import ipywidgets as widgets

#import scipy.optimize as fsolve
from scipy.optimize import fsolve 
import scipy.optimize as opt

#---------------------------------
#-- Voltage divider schematics
def draw_dividers():

    v_divider = schemdraw.Drawing(inches_per_unit=.5, unit=2)   
    v_divider.add(elm.DOT,open='true',label='$V_o$')
    vdr1 = v_divider.add(elm.RES, d='down', label='$R_0$')
    v_divider.add(elm.DOT)
    vdr2 = v_divider.add(elm.RES, d='down', label='$R_1$')
    v_divider.add(elm.GND, botlabel='Fig. 1A',lblofst=1)
    v_divider.add(elm.LINE, d='right', xy=vdr1.end, l=v_divider.unit)
    v_divider.add(elm.DOT, open='true', label='$V^{(A)}$')
    v_divider.here = [v_divider.here[0]+2*v_divider.unit, v_divider.here[1]]
    vdr3=v_divider.add(elm.DOT)
    v_divider.add(elm.RES, d='up', label='$R_1$')
    v_divider.add(elm.DOT, open='true', label='$V_o$')
    v_divider.add(elm.RES,d='down',xy=vdr3.start, label='$R_0$')
    v_divider.add(elm.GND, botlabel='Fig. 1B',lblofst=1)
    v_divider.add(elm.LINE, d='right', xy=vdr3.end, l=v_divider.unit)
    v_divider.add(elm.DOT, open='true', label='$V^{(B)}$')
    return v_divider.draw()

def draw_dividers2():
    
    #SchemDraw schematic of a generic voltage divider circuit
    v2_divider = schemdraw.Drawing(inches_per_unit=.5, unit=2)
    v2_divider.add(elm.DOT,open='true',label='$V_o$')
    v2dr1 = v2_divider.add(elm.RES, d='down', label='$R_0$')
    v2_divider.add(elm.DOT)
    v2dr2 = v2_divider.add(elm.RES_VAR, reverse='true',d='down', label='$R(T)$')
    v2_divider.add(elm.GND, botlabel='Fig. 2A',lblofst=1)
    v2_divider.add(elm.LINE, d='right', xy=v2dr1.end, l=v2_divider.unit)
    v2_divider.add(elm.DOT, open='true', label='$V^{(2A)}$')
    v2_divider.here = [v2_divider.here[0]+2*v2_divider.unit, v2_divider.here[1]]
    v2dr3=v2_divider.add(elm.DOT)
    v2_divider.add(elm.RES_VAR, d='up', flip='true',label='$R(T)$')
    v2_divider.add(elm.DOT, open='true', label='$V_o$')
    v2_divider.add(elm.RES,d='down',xy=v2dr3.start, label='$R_0$')
    v2_divider.add(elm.GND, botlabel='Fig. 2B',lblofst=1)
    v2_divider.add(elm.LINE, d='right', xy=v2dr3.end, l=v2_divider.unit)
    v2_divider.add(elm.DOT, open='true', label='$V^{(2B)}$')
    return v2_divider.draw()

def draw_divamp():
    
    div_amp = schemdraw.Drawing(inches_per_unit=.5, unit=2)
    
    op = div_amp.add(elm.Opamp,flip='true')
    div_amp.add(elm.LINE, d='left', xy=op.in2, l=div_amp.unit*.75)
    p1=div_amp.add(elm.DOT)
    div_amp.add(elm.LINE,l=1.0, d='left')
    div_amp.add(elm.RES_VAR,d='left',label="$R(T)$",reverse='true',flip='true')
    div_amp.add(elm.LINE,d='left',l=0.5)
    div_amp.add(elm.DOT,open='true',label='$V_o$')
    div_amp.add(elm.LINE,d='down', l=div_amp.unit*1, xy=p1.start)
    div_amp.add(elm.GND)   
    p3=div_amp.add(elm.LINE,d='left', xy=op.in1, l=div_amp.unit/4)
    div_amp.add(elm.LINE,d='right', xy=op.out,l=1)
    div_amp.add(elm.DOT,open='true',label='$V$')
    div_amp.add(elm.LINE,xy=p3.end,d='down',l=.75)
    div_amp.add(elm.GND)
    div_amp.add(elm.LINE,d='down', xy=op.vd, l=.5)
    div_amp.add(elm.GND)
    div_amp.add(elm.LINE,d='up', xy=op.vs, l=.5)
    div_amp.add(elm.VDD,label='$V_0$')
    div_amp.add(elm.LINE,d='down', xy=op.n2, l=.5)
    div_amp.add(elm.LINE,d='right',l=.25)
    div_amp.add(elm.DOT,rgtlabel='$V_{ref}$')
    return div_amp.draw()    
#--- Voltage divider schematics    
#---------------------------------


#---------------------------------
#--- Wheatstone bridge schematics
def draw_bridge():
    
    wbridge = schemdraw.Drawing(inches_per_unit=.5, unit=3)
    br1 = wbridge.add(elm.RES,theta=45, toplabel='$R_1$')
    br_top=wbridge.add(elm.DOT)
    br2 = wbridge.add(elm.RES,theta=-45, toplabel='$R_3$')
    br_right=wbridge.add(elm.DOT)
    br3 = wbridge.add(elm.RES_VAR, theta=-135, botlabel='$R(T)$', flip='true',reverse='true')
    br_bot=wbridge.add(elm.DOT,botlabel='Fig. 3A',lblofst=1)
    br4 = wbridge.add(elm.RES,theta=135, botlabel='$R_2$')
    br_left=wbridge.add(elm.DOT)
    wbridge.add(elm.LINE,d='right',xy=br_top.start,l=wbridge.unit*1.25)
    wbridge.add(elm.DOT, open=True, label='$V_T^{(3A)}$')
    wbridge.add(elm.LINE,d='right',xy=br_bot.start,l=wbridge.unit*1.25)
    wbridge.add(elm.DOT,open=True, label='$V_B^{(3A)}$')
    wbridge.add(elm.LINE,d='left',xy=br_left.start,l=wbridge.unit/4)
    wbridge.add(elm.VDD,label='$V_0$')
    wbridge.add(elm.LINE,d='right',xy=br_right.start,l=wbridge.unit/4)
    wbridge.add(elm.GND)

    wbridge.here = [wbridge.here[0]+1.5*wbridge.unit, wbridge.here[1]]
    br2_left=wbridge.add(elm.DOT)
    br5=wbridge.add(elm.RES,theta=45, toplabel='$R_1$')
    br2_top=wbridge.add(elm.DOT)
    br6 = wbridge.add(elm.Resistor(theta=-45, toplabel='$R_3$'))
    br2_right=wbridge.add(elm.Dot())
    br7 = wbridge.add(elm.RES, theta=-135, botlabel='$R_2$', flip='true',reverse='true')
    br2_bot=wbridge.add(elm.DOT,botlabel='Fig. 3B',lblofst=1)
    br8 = wbridge.add(elm.RES_VAR,theta=135, flip='true',botlabel='$R(T)$')

    wbridge.add(elm.LINE,d='left',xy=br2_left.start,l=wbridge.unit/4)
    wbridge.add(elm.VDD,label='$V_0$')
    wbridge.add(elm.LINE,d='right',xy=br2_right.start,l=wbridge.unit/4)
    wbridge.add(elm.GND)
    wbridge.add(elm.LINE,d='right',xy=br2_top.start,l=wbridge.unit*1.25)
    wbridge.add(elm.DOT, open=True, label='$V_T^{(3B)}$')
    wbridge.add(elm.LINE,d='right',xy=br2_bot.start,l=wbridge.unit*1.25)
    wbridge.add(elm.DOT, open=True, label='$V_B^{(3B)}$')

    return wbridge.draw()


def draw_bridgeamp():
    
    wwbridge = schemdraw.Drawing(inches_per_unit=.5, unit=3)
    wbr1 = wwbridge.add(elm.RES,theta=45, toplabel='$R_1$')
    wbr_top=wwbridge.add(elm.DOT)
    #wbridge.add(elm.Vdd(label='$V_0$'))
    wbr2 = wwbridge.add(elm.RES,theta=-45, toplabel='$R_3$')
    wbr_right=wwbridge.add(elm.DOT)
    wbr3 = wwbridge.add(elm.RES_VAR,theta=-135,flip='true', botlabel='$R(T)$')
    wbr_bot=wwbridge.add(elm.DOT)
    #wbridge.add(elm.Ground())
    wbr4 = wwbridge.add(elm.RES,theta=135, botlabel='$R_2$')
    wbr_left=wwbridge.add(elm.DOT)
    wwbridge.add(elm.LINE,d='right',xy=wbr_top.start,l=wwbridge.unit*1.25)
    rn1=wwbridge.add(elm.DOT,open=True, label='$V_T$')
    wwbridge.add(elm.LINE,d='right',xy=wbr_bot.start,l=wwbridge.unit*1.25)
    rn2=wwbridge.add(elm.DOT, open=True, botlabel='$V_B$')
    wwbridge.add(elm.LINE,d='left',xy=wbr_left.start,l=wwbridge.unit/4)
    wwbridge.add(elm.VDD,label='$V_0$')
    wwbridge.add(elm.LINE,d='right',xy=wbr_right.start,l=wwbridge.unit/4)
    wwbridge.add(elm.GND)
    wwbridge.add(elm.LINE,d='down',xy=rn1.start,l=wwbridge.unit*0.5)
    wwbridge.add(elm.LINE,d='right',l=wwbridge.unit*0.5)
    O1=wwbridge.add(elm.OPAMP,anchor='in2',flip='true')
    wwbridge.add(elm.LINE,d='up',xy=rn2.start,l=wwbridge.unit*0.5)
    wwbridge.add(elm.LINE,d='left', l=wwbridge.unit*0.5, xy=O1.in1)
    wwbridge.add(elm.LINE,d='up', xy=O1.vs, l=1/2)
    wwbridge.add(elm.VDD,label='$V_0$')
    wwbridge.add(elm.LINE,d='down', xy=O1.vd, l=1/2)
    wwbridge.add(elm.GND)
    wwbridge.add(elm.LINE,d='right', xy=O1.out,l=1)
    wwbridge.add(elm.DOT,open='true',label='$V$')
    wwbridge.add(elm.LINE,d='down', xy=O1.n2, l=.5)
    wwbridge.add(elm.LINE,d='right',l=.25)
    wwbridge.add(elm.DOT,rgtlabel='$V_{ref}$')
    
    return wwbridge.draw()    
    
#--- Wheatstone bridge schematics
#---------------------------------


#---------------------------------
# Voltage divider: thermistor-to-ground configuration (Configuration 'A')
def div_tg(T, Vin, res0, B, R25, T25):
    #R25 = 1.0e4
    #B = 3977.0
    res = R25 * np.exp(B*(1.0/T-1.0/T25))
    f=res/(res0+res) * Vin
    return f

def ddt_div_tp(T, Vin, res0, B, R25, T25):
    
    res = R25 * np.exp(B*(1.0/T-1.0/T25))
    f=Vin*B*res0/T**2.0/(res0+res)**2.0 * res

    return f

def ddt_div_tg(T, Vin, res0, B, R25, T25):
    
    res = R25 * np.exp(B*(1.0/T-1.0/T25))
    f= Vin*(B*res**2.0 /(T**2.0 * (res0+res)**2.0) - \
        B*res/(T**2.0 * (res0+res)))

    return f

# Voltage divider: thermistor-to-power configuration (Configuration 'B')
def div_tp(T, Vin, res0, B, R25, T25):
    res = R25 * np.exp(B*(1.0/T-1.0/T25))
    f = res0/(res0+res) * Vin
    return f

# Bridge: thermistor-to-ground configuration
def b_tg(T, Vin, res0, B, R25, T25, rho):
    res = R25 * np.exp(B*(1.0/T-1.0/T25))
    f = -(1/(1.0 + rho)-res/(res0+res)) * Vin
    return  f

# Bridge: thermistor-to-power configuration
def b_tp(T, Vin, res0, B, R25, T25, rho):
    res = R25 * np.exp(B*(1.0/T-1.0/T25))
    f = -(rho/(1.0+rho) - res0/(res0+res)) * Vin
    return f

def div_plot(R0):
    R25 = 1.0e4
    B = 3977.0
    T25 = 273.15 + 25.0
    V_in=3.3
    # Input temperatures in Kelvin
    T25 = 273.15 + 25.0
    temp_K = np.arange(250.0, 310.0, 0.2)
    fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    out_tg=np.clip(div_tg(temp_K, V_in, R0, B, R25, T25),0,V_in)
    out_tp=np.clip(div_tp(temp_K, V_in, R0, B, R25, T25),0,V_in)
    
    t_infl = opt.fsolve(f_inflection, 280., args=R0)
    
    slope_tg=ddt_div_tg(t_infl, V_in, R0, B, R25, T25)
    slope_tp=ddt_div_tp(t_infl, V_in, R0, B, R25, T25)
    
    lin_tg = div_tg(t_infl, V_in, R0, B, R25, T25) + \
        slope_tg*(temp_K - t_infl)
    lin_tp = div_tp(t_infl, V_in, R0, B, R25, T25) + \
        slope_tp*(temp_K - t_infl)
    
    
    ax.plot(temp_K,out_tg,label='$V^{(2A)}$')
    ax.plot(temp_K,out_tp,label='$V^{(2B)}$')
    plt.plot(t_infl, np.clip(div_tg(t_infl, V_in, R0, B, R25, T25),0,V_in), \
             marker='o',color="gray" )
    plt.plot(t_infl, np.clip(div_tp(t_infl, V_in, R0, B, R25, T25),0,V_in), \
             marker='o', color="gray" )
    plt.plot(temp_K,lin_tg,':', color="gray")
    plt.plot(temp_K,lin_tp,':', color="gray")
    plt.ylim(0.,4.0)
    plt.legend()
    plt.title('Outputs from the Thermistor Voltage Divider Configurations \n (10K Thermistor: Vishay Model Booo )')
    plt.ylabel('Divider Circuit Output (Volts)')
    plt.xlabel('Temperature (Kelvin)')
    return plt.show()

def divider_plot():

    return widgets.interact(div_plot, \
        R0=widgets.IntSlider(min=5000, max=100000, step=500, value=17900.,description=r'\(R_0 (\Omega) \)'))


# Bridge: thermistor-to-ground configuration
def bridge_tg(T, Vin, res0, B, R25, T25,rho):
    res = R25 * np.exp(B*(1.0/T-1.0/T25))
    f=(1.0/(rho+1.0)-res/(res0+res)) * Vin
    return  f

# Bridge: thermistor-to-power configuration
def bridge_tp(T, Vin, res0, B, R25, T25,rho):
    res = R25 * np.exp(B*(1.0/T-1.0/T25))
    f=(1.0/(rho+1.0) - res0/(res0+res)) * Vin
    return f



def br_plot(R0, RHO):
    
    R25 = 1.0e4
    B = 3977.0
    T25 = 273.15 + 25.0
    V_in=3.3
    
    temp_K = np.arange(250.0, 310.0, 0.2)
    
    fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    
    out_tg=np.clip(bridge_tg(temp_K, V_in, R0, \
                                 B, R25, T25,RHO), 0,V_in)
    out_tp=np.clip(bridge_tp(temp_K, V_in, R0, \
                                 B, R25, T25,RHO),0,V_in)


    
    t_infl = opt.fsolve(f_inflection, 280., args=R0)
    
    slope_tg=-ddt_div_tg(t_infl, V_in, R0, B, R25, T25)
    slope_tp=-ddt_div_tp(t_infl, V_in, R0, B, R25, T25)
    
    lin_tg = bridge_tg(t_infl, V_in, R0, B, R25, T25, RHO) + \
        slope_tg*(temp_K - t_infl)
    lin_tp = bridge_tp(t_infl, V_in, R0, B, R25, T25, RHO) + \
        slope_tp*(temp_K - t_infl)

    ax.plot(temp_K,out_tg,label='$V^{(3A)}$')
    ax.plot(temp_K,out_tp,label='$V^{(3B)}$')
    
    plt.plot(t_infl, np.clip(bridge_tg(t_infl, V_in, R0, B, R25, T25,RHO),0,V_in), \
             marker='o',color="gray" )
    plt.plot(t_infl, np.clip(bridge_tp(t_infl, V_in, R0, B, R25, T25, RHO),0,V_in), \
             marker='o', color="gray" )
    plt.plot(temp_K,lin_tg,':', color="gray")
    plt.plot(temp_K,lin_tp,':', color="gray")
    #ax.plot(temp_K,out_div_tp,label='$V_{div}$')
    plt.ylim(0.,3.5)
    plt.legend()
    plt.title('Outputs from the Amplified Voltage Divider Configurations')
    plt.ylabel('Divider Circuit Output (Volts)')
    plt.xlabel('Temperature (Kelvin)')
    
    return plt.show()

def bridge_plot2():

    return widgets.interact(br_plot, \
        R0=widgets.IntSlider(min=5000, max=100000, step=100, value=17500.,description=r'\(R_0 (\Omega) \)'), \
        RHO=widgets.FloatSlider(min=0, max=0.5, step=.005, value=0,description=r'\(\rho\)'))


def amp_plot(R0,V_ref,RHO,A_G):
    R25 = 1.0e4
    B = 3977.0
    T25 = 273.15 + 25.0
    V_in=3.3
    temp_K = np.arange(250.0, 310.0, 0.2)
    #temp_K = np.arange(250.0, 310.0, 0.2)
    fig, ax = plt.subplots(1, 1, figsize=(6,4))
    out_tg=np.clip(A_G*b_tg(temp_K, V_in, R0, \
                                 B, R25, T25,RHO)+V_ref,0.0,V_in)
    out_tp=np.clip(A_G*b_tp(temp_K, V_in, R0, \
                                 B, R25, T25,RHO)+V_ref,0.0,V_in)
    out_div_tp=np.clip(A_G*div_tp(temp_K, V_in, R0, B, R25, T25)+V_ref,0,V_in)
    
    #ax.plot(temp_K,out_tg,label='$V_{TG}$')
    ax.plot(temp_K,out_tp,label='$V_{br}$')
    ax.plot(temp_K,out_div_tp,label='$V_{div}$')
    plt.ylim(0.,3.5)
    plt.legend()
    plt.title('Outputs from the Amplified Voltage Divider Configurations')
    plt.ylabel('Divider Circuit Output (Volts)')
    plt.xlabel('Temperature (Kelvin)')
    
    return plt.show()

def ampcircuits_plot():
    
    return widgets.interact(amp_plot, \
        R0=widgets.IntSlider(min=5000, max=100000, step=100,value=17500.,description=r'\(R_0 (\Omega) \)'), \
        V_ref=widgets.FloatSlider(min=0, max=3.3, step=.01,value=0,description=r'\(V_{\it ref} (V)\)'), \
        RHO=widgets.FloatSlider(min=0, max=1, step=.01,value=0,description=r'\(\rho\)'), \
        A_G=widgets.FloatSlider(min=1, max=2, step=.01,value=0,description=r'\(A_{G}\)'))


#---------------------------------

def therm_res(T, B, R25, T25):
    
    #B-Parameter equation for thermistor resistance
    
    res = R25 * np.exp(B*(1.0/T-1.0/T25))
    
    return res

def d_therm_res(T, B, R25, T25):
    
    f= - B/T**2.0 * therm_res(T,B,R25,T25)
    
    return f

def thermres_plot():
    
    #Define the temperature range in Kelvin
    temp_K = np.arange(250.0, 310.0, 0.2)

    #Define thermistor parameters

    R25 = 1.0e4
    T25 = 273.15 + 25.0
    B = 3977.0
    
    r_t = therm_res(temp_K, B, R25, T25)
    
    fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    ax.plot(temp_K,r_t/1000.)
    plt.title('Temperature Dependence of Thermistor Resistance \n (10K Thermistor: Vishay Model )');
    plt.ylabel('Thermistor Resistance ($k\Omega$)');
    plt.xlabel('Temperature (Kelvin)');   
    
    return plt.show()

def inflection():
    
    R25 = 1.0e4
    T25 = 273.15 + 25.0
    B = 3977.0
    
    #Define the temperature range in Kelvin
    temp_K = np.arange(250.0, 310.0, 0.2)

    
    res0= (B-2.0*temp_K)/(B+2.0*temp_K)*therm_res(temp_K,B,R25,T25)
    
    fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    ax.plot(temp_K,res0/1000.,temp_K,therm_res(temp_K,B,R25,T25)/1000.)
    plt.title('Temperature Dependence of Thermistor Resistance \n (10K Thermistor: Vishay Model )');
    plt.ylabel('Thermistor Resistance ($k\Omega$)');
    plt.xlabel('Temperature (Kelvin)');   
    
    return plt.show()

def f_inflection(tkel, *data):
    
    R0=data
    
    R25 = 1.0e4
    T25 = 273.15 + 25.0
    B = 3977.0
    
    f = (B-2.0*tkel)/(B+2.0*tkel)*therm_res(tkel,B,R25,T25)-R0
    
    return f

def f_prime(tkel):

    R25 = 1.0e4
    T25 = 273.15 + 25.0
    B = 3977.0
    
    fp = - (2.0 + 2.0*(B-2.0*tkel)/(B+2.0*tkel) + B*(B-2.0*tkel)/tkel**2.0) \
        * therm_res(tkel,B,R25,T25) / (B+2.0*tkel)
    
    return fp
    


