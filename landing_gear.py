"""
- Conditions number module.
- Some datas used are in imperial measurement system.
- Final Datas in Metric measurement system (S.I.).
- Methology: PAZMANY, L. Landing Gear Design for Light Aircraft - Part I
- Regulation Criteria: FAR23

"""


from math import sqrt, atan
from inputs import *
from vn_diagram import *


def balance(mtow, fuel_weight):

    # FAR 23.473 - Design Landing Weight:

    W = mtow - (0.25*fuel_weight)
    used_weight = mtow * (0.95)
    return max(W, used_weight)

design_landing_weight = round(balance(mtow, fuel_weight), 3)

print(f'Weight used for landing gears project: {design_landing_weight}')

# Velocity in m/s

FT_FACTOR = 0.092903

converted_area = wing_surface/FT_FACTOR

def descentLoad (S_wing_converted, Vs, cl_stall):
    #Imperial measure system (ft², mph)
    L = S_wing_converted*(Vs**2)*cl_stall
    return L/390

K_load = descentLoad(converted_area, stallSpeed*CONSTANT_KNOTS, cl_stall)

print(K_load)

def load_ratio(K_load, designLandingWeight):
    r = K_load/designLandingWeight
    return r

ratio = round(load_ratio(K_load, design_landing_weight*CONSTANT_KGF), 3)
    
if ratio > 2/3:

    print(f'{ratio} ration doesnt respondes FAR23.473(b) criteria.')

else:

    print(f'{ratio} ratio succeeded FAR23.473(b) criteria.')


    def drop_height(designLandingWeight, S_wing):

        def equivalent_diving_speed (designLandingWeight, S_wing):
            # Relative Diving Velocity, metric measure (m/s)

            Vdm = 4.4*(((designLandingWeight)/(S_wing))**(1/4)) 
            return Vdm

        Vd = equivalent_diving_speed (designLandingWeight, S_wing)

        print('-'*80)
        print(f'Velocidade Descendente (ft/s): {Vd}')

        #FAR FAR23.473 uses imperial measurement system, ft/s.

        if Vd < 7 or Vd > 10:
            print("Não condiz com a norma.")
        else:
            print("Está dentro dos parâmetros da Norma FAR23.473(b).")


        variable1 = (designLandingWeight/S_wing)**0.5
        return variable1 * 3.6 

    equivalent_height = drop_height(design_landing_weight, wing_surface)

    #FAR 23.725 uses imperial measurement system, inches.
    

    if equivalent_height < 9.2 or equivalent_height > 18.7:

        print("Total Failure. FAR 23.725 is not being followed.")

    else:

        print(f'Process succeeded. FAR 23.725 approves. {equivalent_height}')

    metric_height = equivalent_height*CONSTANT_IN


    def drop_test_time(g, h):
        t = (2*h)/g
        return sqrt(t)

    time = drop_test_time(GRAVITY, metric_height)


    def shock_absortion_stroke(W, 
                            S, nz, 
                            tire_deflection,
                            tire_efficiency,
                            strut_efficiency):

        # Assuming gravity in ft/s²
        # ds: Strut Deflection - TPP
        # dt: Tire Deflection - TPP

        ds_1 = 0.3*((W/S)**0.5) - tire_deflection*((nz*tire_efficiency)-0.333)
        ds_2 = (nz*strut_efficiency)-0.33
        return ds_1/ds_2


    shock = shock_absortion_stroke(design_landing_weight*CONSTANT_KGF, converted_area, cn_max_pos, tire_deflection, tiree_efficiency, strut_efficiency)
    metric_shock = shock*CONSTANT_IN

    print(f'Shock defleciton (ft): {shock}')
    print(f'Shock defleciton (m): {metric_shock}')

    def effective_weight(W, h, ratio, shock):

        variable2 = 1-ratio
        variable3 = h + (variable2*shock)
        variable4 = (h+shock)
        division1 = variable3/variable4
        return division1*W

    effectWeight = effective_weight(design_landing_weight, metric_height, ratio, metric_shock)

    print(f'Effective Weight: {effectWeight}')

    
def components_height(shock, CG_Wheel_Difference):

    # to metric measure system

    h = CG_Wheel_Difference - shock

    return h

# coming soon

h1 = components_height(metric_shock, CG_Wheel_Difference=nose_strut_height)
h2 = components_height(metric_shock, CG_Wheel_Difference=landgear_Height)

K = 0.25


# According to FAR23.474(g) and 725(d) minimum load factor = 2, 
# calculate with dv/dt (acceleration of drop test) in g's,
# plus 1. If there's no dynamic reaction, assume 1g + 1;

variable = equivalent_height/shock

ng = (variable*CONSTANT_IN)+1


# According to FAR23.723(e), Limit inertia Load Factor:
# According to FAR23.723(f), the value may not be more than limit inercia load factor

n = (ng*(effectWeight/design_landing_weight)) + ratio

H = K * n * design_landing_weight
V = (n - ratio) * design_landing_weight

print(f'\nGround load factor: {ng, n}')
print(f'\nLanding gear deflection: {h2}')

print(f'\nVertical load: {V}')
print(f'\nHorizontal load: {H}')

# Imperial Measure System

cg_to_landgear = 0.048
cg_to_nosestrut = 0.272

d = cg_to_nosestrut + cg_to_landgear

semi_spanwise = wing_spanwise/2

def three_wheels_landing(a, b, d, H, V):

    #Condição
    
    #Somatório dos Momentos
    
    Vn = V*(a/d)
    VN = V*(b/d)
    HN = H*(b/d)
    Hn = H*(a/d)
    return VN,HN,(Vn+VN), (Hn+HN)

def two_wheels_landing(V, H, a, d):
    Vn = V*(a/d)
    Hn = H*(a/d)
    return Vn, Hn

def one_wheel_landing(V, H, a, d):
    Vn = V*(a/d)
    Hn = H*(a/d)
    return 2*Vn ,2*Hn

def nose_steering_angle (d, semiSpanwise):
    delta = d/semiSpanwise
    return atan(delta)


t_r=three_wheels_landing(cg_to_nosestrut, cg_to_landgear,d, H, V)
d_r=two_wheels_landing(V, H, cg_to_nosestrut, d)
u_r=one_wheel_landing(V, H, cg_to_nosestrut, d)
steering = nose_steering_angle(d, semi_spanwise)

print('-'*80)
print(f'Nose strut load (Vertical) | Nose strut load (Horizontal) | Level landing load (Vertical) | Level landing load (Horizontal)')
print(t_r)
print('-'*80)
print(f' Landing gear load (Vertical) | Landing gear load (Horizontal)')
print(d_r)
print('-'*80)
print(f'One wheel landing load (Vertical) | One wheel landing load (Horizontal) ')
print(u_r)
print('-'*80)
print(f'Steering Angle 𝛿 (rad)')
print(steering)
print('-'*80)
print(f'Drop Test Height (m)')
print(metric_height)
print('-'*80)
print(f'Drop Test Time (s)')
print(time)
print('-'*80)




