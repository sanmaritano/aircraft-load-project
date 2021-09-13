from math import sqrt
from inputs import *
from vn_diagram import *
import matplotlib.pyplot as plt
import pandas as pd
i = 0

show_empennage = str(input('Do you want to show empennage data? (Y/N): '))

'''
	Note 1: This code was made following Mario Lott's methodology.
'''

diving_rho = dynamic_pressure(air_density, diving_speed)

manuever_rho = dynamic_pressure(air_density, manuever_speed)

stall_rho = dynamic_pressure(air_density, stall_speed)

cruising_rho = dynamic_pressure(air_density, cruising_speed)


x_wf = quarter_chord - x_CG_foward                         #AC-CGfw Distance

x_wa = quarter_chord - x_CG_aft                         # AC-CGaft Distance

def aerodynamic_moment (cm_max, aerodynamic_chord, wing_surface, rho):
	cma = cm_max*aerodynamic_chord*wing_surface*rho
	return cma

cm_stall = aerodynamic_moment(cm_max, aerodynamic_chord, wing_surface, stall_rho)/GRAVITY
cm_manuever = aerodynamic_moment(cm_max, aerodynamic_chord, wing_surface, manuever_rho)/GRAVITY
cm_diving = aerodynamic_moment(cm_max, aerodynamic_chord, wing_surface, diving_rho)/GRAVITY

def balanceTailLoad (cma, maximumTakeOffWeight, n , leverArm, x):

  	variable1 = cma
  	variable2 = n*maximumTakeOffWeight*x
  	variable3 = leverArm - x

  	return (variable1 - variable2)/variable3

# Note 2: fw means FOWARD, aft means rear region, back I think.
  
fzt_fw_down = balanceTailLoad(cm_diving, mtow, cn_max_pos, lever_arm, x=x_wf) 
fzt_fw_up = balanceTailLoad(cm_stall, mtow, cn_max_neg, lever_arm, x=x_wf)

fzt_tras_baixo = balanceTailLoad(cm_diving, mtow, cn_max_pos, lever_arm, x=x_wa)
fzt_tras_cima = balanceTailLoad(cm_stall, mtow, cn_max_neg, lever_arm, x=x_wa)


# Unchecked Manuever


K = 0.75                                             # Elevator's Hinge Moment


def n_cor(n, press, K):                               # FAR 423 - also Appendix A
    w = 4.8 + (0.534*(n*press))
    W = (K*w*n)/4.4
    return W


w_distribuited = n_cor(cn_max_pos, w_distribuited, K)


# Checked Manuevers

def pitching_accelaration(n, v):
	variable = (39/v)*n*(n-1.5)
	return variable


pitching_manuever = pitching_accelaration(cn_max_pos, manuever_speed*1.94384)
pitching_diving = pitching_accelaration(cn_max_pos, diving_speed*1.94384)


# Inercial Loads, use SolidWorks to take it.

def inercial_loads(n, pitchingSpeed, leverArm, gravityAcceleration):
    nt = n - ((pitchingSpeed*leverArm)/gravityAcceleration)
    return nt


nt_manuever_min = inercial_loads(CN_STALL, pitching_manuever, lever_arm, GRAVITY)
nt_manuever_max = inercial_loads(cn_max_pos, pitching_manuever, lever_arm, GRAVITY)
nt_diving_min = inercial_loads(CN_STALL, pitching_diving, lever_arm, GRAVITY)
nt_diving_max = inercial_loads(cn_max_pos, pitching_diving, lever_arm, GRAVITY)


Fzti_manuever_min = nt_manuever_min*tail_weight
Fzti_manuever_max = nt_manuever_max*tail_weight
Fzti_diving_min = nt_diving_min*tail_weight
Fzti_diving_max = nt_diving_max*tail_weight


# Resultantes Aerodinâmicas

def manuever_moment (pitchingSpeed, inercialMomentIxx):
	M = round(pitchingSpeed*inercialMomentIxx, 5)
	return M

def manuever_load_vary(manueverMoment, leverArm):
	Fztm = - round(manueverMoment/leverArm, 5)
	return Fztm

def load_factor_vary(manueverTailLoadVary, maximumTakeOffWeight):
	ncg = manueverTailLoadVary/maximumTakeOffWeight
	return ncg

def requested_load_factor(n, loadFactorVary):
	nb = n - loadFactorVary
	return nb

def polardrag_load_coef(requestedLoadFactor, maximumTakeOffWeight, wingSurface, dynamicPressure):
	cza = (requestedLoadFactor*maximumTakeOffWeight)/(wingSurface*dynamicPressure)
	return cza

def btl_vary(BalanceLoadCoefficient, wingSurface, dynamicPressure):
	Fztb = - BalanceLoadCoefficient*wingSurface*dynamicPressure
	return Fztb

def total_empennage_load (Fztm, Fztb, Fzti):
	Fzt = Fztm + Fztb + Fzti
	return Fzt

# Positive Manuever

pos_Va_moment = manuever_moment(pitching_manuever, inercial_moment_xx)
delta_Fzt_man = manuever_load_vary(pos_Va_moment, lever_arm)
pos_deltanz_man = load_factor_vary(delta_Fzt_man, mtow)
pos_nz_man = requested_load_factor(CN_STALL, pos_deltanz_man)
pos_cza_man = polardrag_load_coef(pos_nz_man, mtow, wing_surface, manuever_rho)
pos_polar_coefficient_man = polar_drag(pos_cza_man)
pos_btl_vary_man = btl_vary(pos_polar_coefficient_man, wing_surface, manuever_rho)
pos_fzt_man = total_empennage_load(delta_Fzt_man, pos_btl_vary_man, Fzti_manuever_min)

coluna1 = [pos_Va_moment, delta_Fzt_man, pos_deltanz_man, pos_nz_man,\
	 pos_cza_man, pos_polar_coefficient_man, pos_btl_vary_man, pos_fzt_man]


# Positive Diving


pos_Vd_moment = manuever_moment(pitching_diving, inercial_moment_xx)
delta_fzt_dive = manuever_load_vary(pos_Vd_moment, lever_arm)
pos_deltanz_dive = load_factor_vary(delta_fzt_dive, mtow)
pos_nz_dive = requested_load_factor(CN_STALL, pos_deltanz_dive)
pos_cza_dive = polardrag_load_coef(pos_nz_dive, mtow, wing_surface, diving_rho)
pos_polar_coefficient_dive = polar_drag(pos_cza_dive)
pos_btl_vary_dive = btl_vary(pos_polar_coefficient_dive, wing_surface, diving_rho)
pos_fzt_dive = total_empennage_load(delta_fzt_dive, pos_btl_vary_dive, Fzti_diving_min)


coluna2 = [pos_Vd_moment, delta_fzt_dive, pos_deltanz_dive, pos_nz_dive,\
	 pos_cza_dive, pos_polar_coefficient_dive, pos_btl_vary_dive, pos_fzt_dive]


neg_Va_moment = manuever_moment(-pitching_manuever, inercial_moment_xx)
neg_delta_Fzt_man = manuever_load_vary(neg_Va_moment, lever_arm)
neg_deltanz_man = load_factor_vary(neg_delta_Fzt_man, mtow)
neg_nz_man = requested_load_factor(cn_max_pos, neg_deltanz_man)
neg_cza_man = polardrag_load_coef(neg_nz_man, mtow, wing_surface, manuever_rho)
neg_polar_coefficient_man = polar_drag(neg_cza_man)
neg_btl_vary_man = btl_vary(neg_polar_coefficient_man, wing_surface, manuever_rho)
neg_fzt_man = total_empennage_load(neg_delta_Fzt_man, neg_btl_vary_man, Fzti_manuever_max)


coluna3 = [neg_Va_moment, neg_delta_Fzt_man, neg_deltanz_man, neg_nz_man,\
	 neg_cza_man, neg_polar_coefficient_man, neg_btl_vary_man, neg_fzt_man]


neg_Vd_moment = manuever_moment(-pitching_diving, inercial_moment_xx)
neg_delta_Fzt_dive = manuever_load_vary(neg_Vd_moment, lever_arm)
neg_deltanz_dive = load_factor_vary(neg_delta_Fzt_dive, mtow)
neg_nz_dive = requested_load_factor(cn_max_pos, neg_deltanz_dive)
neg_cza_dive = polardrag_load_coef(neg_nz_dive, mtow, wing_surface, diving_rho)
neg_polar_coefficient_man = polar_drag(neg_cza_dive)
neg_btl_vary_man = btl_vary(neg_polar_coefficient_man, wing_surface, diving_rho)
neg_fzt_dive = total_empennage_load(neg_delta_Fzt_dive, neg_btl_vary_man, Fzti_diving_max)	


coluna4 = [neg_Vd_moment, neg_delta_Fzt_dive, neg_deltanz_dive, neg_nz_dive,\
	 neg_cza_dive, neg_polar_coefficient_man, neg_btl_vary_man, neg_fzt_dive]

df = pd.DataFrame(list(zip(coluna1, coluna2)), 
	columns=['Va', 'Vd'])
df.to_excel('Cargas Resultantes.xlsx')

design_load = pos_fzt_dive
Fztm = design_load*h_stab_surface

# Gust Loads - Far 23.425(d) - Apply Imperial measure system.


mass_ratio_stabilizer = vertical_mass_ratio(w_distribuited, h_stab_chord, 
										h_alpha, air_density*CONSTANT_SLUG, GRAVITY)

kg = (0.88*mass_ratio_stabilizer)/(5.3+mass_ratio_stabilizer)

def stabilizer_gust_load (factor, kg_factor, Ude, V, sh, ah, derivada_e):
	variable = factor*kg_factor*Ude*V*sh*ah*derivada_e
	return variable/498

Lh1 = stabilizer_gust_load(CONSTANT_LBF, kg, Uc_gust,
		 cruising_speed*1.94, h_stab_surface*CONSTANT_AREA, h_alpha, downwash_coefficient)

Lh2 = stabilizer_gust_load(CONSTANT_LBF, kg, Ud_gust, 
			diving_speed*1.94, h_stab_surface*CONSTANT_AREA, h_alpha, downwash_coefficient)


# Critical Conditions

# FAR 23 - Apendix B.6 - uses distribuition loads to project horizontal load

distribuited_stabilizer_load = abs(design_load/h_stab_surface)

# Vertical Stabilizer Loads


w1 = 0.534*(cn_max_pos*w_distribuited)                              # Far 23 pag 356
w2 = 3.36*(sqrt(cn_max_pos*w_distribuited))
W1 = (K*w1*cn_max_pos)/(4.4)
W2 = (K*w2*cn_max_pos)/(4.4)
W3 = W2

Lv = W2*0.116


# Vertical Stabilizer Gust Loads

mi_lat = lateral_mass_ratio(inercial_moment_yy*23.73036, air_density*CONSTANT_SLUG, v_stab_chord*CONSTANT_FT, v_alpha, v_stab_surface*CONSTANT_AREA)

kgt = (0.88*mi_lat)/(5.3+mi_lat)      

def v_gust_load(factor, kgt, Ude, Speed, alphav, vt_s):
	Lvt = (factor*kgt*Ude*Speed*alphav*vt_s)
	return Lvt/498

verticalTailGustLoad = v_gust_load(CONSTANT_LBF, kgt, Uc_gust, cruising_speed*CONSTANT_KNOTS, v_alpha, v_stab_surface*CONSTANT_AREA)


# Força de Profundor e Leme

'''	FAR 23 - Appendix A23.11
	The most severe elevator and rudder 
	loads should be further considered as being 
	distributed parabolically from three times 
	the mean loading of the surface (stabilizer 
	and elevator, or fin and rudder) at the 
	leadign edge of the elevator and rudder, respectivily,
	to zero at the trailing edge according 
'''

def firstSurfaceLoading(w, localChord, surfaceChord):

	E = (localChord/surfaceChord)*CONSTANT_FT
	d = (surfaceChord-localChord)*CONSTANT_FT

	variable1 = 2*w
	variable2 = (2-E-(3*d))
	variable3 = 1-E
	return variable1*(variable2/variable3)*CONSTANT_LBF

def secondSurfaceLoading(w, localChord, surfaceChord):

	E = (localChord/surfaceChord)*CONSTANT_FT
	d = (surfaceChord-localChord)*CONSTANT_FT

	variable1 = 2*w
	variable2 = (3*d+E-1)
	return variable1*variable2*CONSTANT_LBF


maximumElevatorLoad = firstSurfaceLoading(w_distribuited*0.65, elevator_chord, h_stab_chord)*elevator_surface
minimumElevatorLoad = secondSurfaceLoading(w_distribuited, elevator_chord, h_stab_chord)*elevator_surface

maximumRudderLoad = firstSurfaceLoading(W2*0.65, rudder_chord, v_stab_chord)*rudder_surface
minimumRudderLoad = secondSurfaceLoading(W2, rudder_chord, v_stab_chord)*rudder_surface

maximumRudderLoad1 = firstSurfaceLoading(W1*0.65, rudder_chord, v_stab_chord)*rudder_surface
minimumRudderLoad1 = secondSurfaceLoading(W1, rudder_chord, v_stab_chord)*rudder_surface


vectorElevator = [-maximumElevatorLoad, -minimumElevatorLoad, 0]
vectorEChord = [0, (h_stab_chord-elevator_chord), h_stab_chord]
vectorRudder = [maximumRudderLoad, -minimumRudderLoad, 0]
vectorRChord = [0, (v_stab_chord-rudder_chord), h_stab_chord]
vectorRudder1 = [maximumRudderLoad1, -minimumRudderLoad1, 0]
vectorRChord1 = [0, (v_stab_chord-rudder_chord), h_stab_chord]


omega = (cn_max_pos*w_distribuited)*0.466
w_distribuited = K*omega*(cn_max_pos/4.4)

La = round(w_distribuited*aileron_surface, 4)


# Flaps:

omega = ((2.3131/1.6)*w_distribuited)*0.64

Ff = omega*flap_surface


# Relatório

if show_empennage == str('Y'):

	print(f'\nMCA [N.m]: {cm_manuever} kgf.m')
	print(f'\nFzt Diving CG foward: {fzt_fw_down} kgf')
	print(f'\nFzt Diving CG aft: {fzt_fw_up} kgf')
	print(f'\nFzt Stall CG foward: {fzt_tras_baixo} kgf')
	print(f'\nFzt Stall CG aft: {fzt_tras_cima} kgf')
	print(f'\nInercial Moment Ixx = {inercial_moment_xx}')
	print(f'\nHorizontal stabilizer total weight : {tail_weight} Kg')
	print(f'\nInercial loads for nz = 1 at Va: \nnt = {Fzti_manuever_min}')
	print(f'\nInercial loads for nCA = 2.5 at Va: \nnt = {Fzti_manuever_max}')
	print(f'\nInercial loads for nCA = 1 at Vd:  \nnt = {Fzti_diving_min}')
	print(f'\nInercial loads for nCA = 2.5 at Vd: \nnt = {Fzti_diving_max}')

	print()
	print(f'\n\nHORIZONTAL STABILIZER LOADS:')
	print(f'\nXwa = {x_wa} m\nXwf = {x_wf} m')
	print(f'\nMaximum load factor: {cn_max_pos}')
	print(f'\nHorizontal tail mean force: {Fztm} Kgf')
	print(f'\nManuever pitching acceleration: \nθa = {pitching_manuever} rad/s²')
	print(f'\nDiving pitching acceleration: \nθd = {pitching_diving} rad/s²')
	print(f'\nVertical gust load at 50 ft/s:{Lh1} Kgf\n\nVertical gust load at 25ft/s {Lh2} Kgf')
	print(f'\nLateral gust load: {verticalTailGustLoad} Kgf')

	print(f'\nAerodynamics Resultant Loads Table: ')
	print()
	print(coluna1)
	print(coluna2)
	print(coluna3)
	print(coluna4)

	print()
	print(f'\nVERTICAL STABILIZER LOADS:')

	print(f'\nAbrupt deflection yaw - FAR 23.441(a)(1):\n{W2} kgf/m²')
	print(f'\nComeback deflection yaw - FAR 23.441(a)(2):\n{W3} kgf/m²')
	print(f'\n15 degrees deflection yaw - FAR 23.441(a)(3):\n{W1} kgf/m²')

	print(f'\nLateral mass ratio = {mi_lat}')
	print(f'\nLateral gust alleviation factor = {kgt}')
	print(f'\nHorizontal tail distribuited load: {distribuited_stabilizer_load} Kgf')
	print(f'\nVertical tail distribuited load: {Lv} Kgf')

	print()
	print(f'\nHorizontal tail leading edge load: {maximumElevatorLoad} Kgf')
	print(f'\nElevators hinge moment line load: {minimumElevatorLoad} Kgf')
	print(f'\nVertical tail leading edgeload: {maximumRudderLoad} Kgf')
	print(f'\nRudders hinge moment line load: {minimumRudderLoad} Kgf')

	print()
	print(f'Aileron loads: {La} Kgf')
	print(f'Flap loads: {round(2*Ff, 4)} Kgf')

	fig = plt.figure(1)


	plt.subplot(2, 1, 1)
	plt.plot(vectorRChord1, vectorRudder1, '.-')
	plt.title('Rudder load distribuition - Abrupt e 15 degrees')
	plt.ylabel('Load (kgf)')
	plt.subplot(2, 1, 2)
	plt.plot(vectorRChord, vectorRudder, '.-')
	plt.xlabel('Vertical empennage chord (m)')
	plt.ylabel('Load (kgf)')

	fig = plt.figure(2)
	plt.plot(vectorEChord, vectorElevator, '.-')
	plt.title('Elevator load distribuition')
	plt.xlabel('Horizontal empennage chord (m)')
	plt.ylabel('Load  (kgf)')

	plt.tight_layout()
	plt.show()

