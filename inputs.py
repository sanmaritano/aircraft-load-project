from math import pow

# Constants

GRAVITY = 9.81

#Contant to metric system

CONSTANT_LBF = 0.453592

CONSTANT_KGF = 2.20462

CONSTANT_TORQUE = 725.25

CONSTANT_SLUG = 0.00194

CONSTANT_AREA = 10.7639

CONSTANT_FT = 3.28084

CONSTANT_IN = 0.0254

CONSTANT_KNOTS = 1.94384

CN_STALL = 1

CL_STALL = 1

CN_NULL = 0


air_density = 1

mtow = 1

# Wing + Aileron

# Derivate with dCl/da on wing
alpha = 4

wing_spanwise = 2

aerodynamic_chord = 0.5 

quarter_chord = 0.25*aerodynamic_chord

root_chord = 0.5

wing_surface = wing_spanwise*aerodynamic_chord

aileron_surface = 0.03

deflection_aileron_up = 25

deflection_aileron_down = 25

flap_surface = 0.03

cl_max = 2

clmax_flapped = 2.5

cd_max = 0.5

cm_max = -0.3


#Horizontal Empennage

h_stab_chord = 0.25

h_stab_length = 1

h_stab_surface = h_stab_chord*h_stab_length

elevator_chord = 0.1

elevator_surface = 0.03

h_alpha = 2

downwash_coefficient = 0.5

# Vertical Empennage

v_stab_chord = 0.25

v_stab_length = 0.5

v_stab_surface = v_stab_chord*v_stab_length

rudder_chord = 0.1

rudder_surface = 0.075

lever_arm = 1

# Derivate with dCl/da on VS

v_alpha = 1

x_CG_position = 0.5

x_CG_foward = 0.45

x_CG_aft = 0.55

# Structural

cn_max_pos = 2.5

cn_max_neg = -1

cn_ult_pos = cn_max_pos * 1.5

cn_ult_neg = cn_max_neg * 1.5

inercial_moment_xx = 0.5
inercial_moment_yy = 0.25

tail_weight = 0.25


# aluminum efficiency
tire_efficiency = 0.5
# Deflection in foot
tire_deflection = 0.5
# sandwich matrix efficiency
strut_efficiency = 0.5

nose_strut_height = 0.5  # m
landgear_Height = 0.5  # m

# Engine

fuel_weight = 0.25  # kgf

Vc_max = 25
# Velocidade máxima em voo nivelado
Vc_min = 5
# Velocidade mínima em voo nivelado

cruise = 0

dive = 0

Uc_gust = 50

Ud_gust = 25

engine_mass = 0.5

maximum_power = 2.5    #HP

maximum_rot = 20000        #RPM


def dynamic_pressure(airDensity, speed):
    q = airDensity * pow(speed, 2)
    return q/2

def polar_drag (cl):
    '''
    Enter your theorical polar drag like this one:
    '''
    cd = 0.05 + pow(0.25, 2)
    return cd

def weight_distribuition(mtow, wing_surface):

    return mtow/wing_surface


w_distribuited = weight_distribuition(mtow, wing_surface)


def vertical_mass_ratio(W, chord, a, air_density, g):
    variable1 = (2 * W)
    variable2 = (air_density * chord * a * g)
    return variable1/variable2

def gust_alleviation_factor(mass_ratio):
    return (0.88 * mass_ratio) / (5.3 + mass_ratio)

def lateral_mass_ratio(Iyy, air_density, verticalStabilizerChord, alphav, verticalStabilizerSurface):
    variable1 = 2*Iyy
    variable2 = air_density*verticalStabilizerChord*alphav*verticalStabilizerSurface
    return variable1/variable2