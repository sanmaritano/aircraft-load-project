from math import sqrt
from matplotlib import pyplot as plt
from inputs import *
import pandas as pd

show_vn = str(input('Do you want to show VN data? (Y/N): '))

    
# Speed Equations


if cruise == 0:

    def cruising_speed_equation(mtow, g, wing_area, maximum_cruise_speed):

        first_cruising_equation = (2.4 * \
            sqrt((mtow * g) / wing_area)) * 0.544

        second_cruising_equation = 0.9 * maximum_cruise_speed

        if show_vn == str('Y'):
            
            print("Cruising speed may not be more than 90% of maximum cruise speed.")

            if first_cruising_equation <= max(first_cruising_equation, second_cruising_equation) <= 0.9 * maximum_cruise_speed:
                print("\nCruising Speed is according to FAR23.335(a)")
            else:
                print("\nBack to analysis.")

        return max(first_cruising_equation, second_cruising_equation)
        

    cruising_speed = cruising_speed_equation(mtow, GRAVITY, wing_surface, Vc_max)


else:

    cruising_speed = cruise


if dive == 0:

    def diving_speed_equation(cruising_speed):

        firstDivingEquation = 1.25 * cruising_speed
        

        return max(firstDivingEquation, firstDivingEquation)

        
    diving_speed = diving_speed_equation(cruising_speed)


else:
    diving_speed = dive


def speed_curve(n, g, W, rho, cl_max):
            variable1 = (2 * abs(n) * W * g)
            variable2 = (rho * cl_max)
            return sqrt(variable1/variable2)

# Maximum Load Factor limits curve - definition of stall region

def curves():

    
    n = 0
    vetor1=[]
    vetor2=[]

    
    while n < cn_ult_pos:
        
        positive_speed = speed_curve(n, GRAVITY, w_distribuited, air_density, cl_max)

        vetor2.append(positive_speed)
        vetor1.append(n)

        n += 0.05

        plt.plot(vetor2, vetor1, color='darkorange')
    
    if show_vn == str('Y'):
        print(f'\nPositive Curve (1):\n{vetor1} \n{vetor2}') 

    df1 = pd.DataFrame(list(zip(vetor1, vetor2)), 
    columns=['negativeLoad', 'negativeSlopSpeed'])
    df1.to_excel('vn.xlsx')

    n = 0
    vetor3=[]
    vetor4=[]

    while n > (cn_ult_neg-0.05):

        negative_speed = speed_curve(n, GRAVITY, w_distribuited, air_density, cl_max)

        vetor3.append(n)
        vetor4.append(negative_speed)
        
        n -= 0.05
        
        plt.plot(vetor4, vetor3, color='darkorange')

    df2 = pd.DataFrame(list(zip(vetor3, vetor4)), 
    columns=['negativeLoad', 'negativeSlopSpeed'])
    df2.to_excel('vn.xlsx')
    
    if show_vn == str('Y'):
        print(f'\nNegative Curve (2): \n{vetor4}')
        

curves()           


stall_speed = speed_curve(CN_STALL, GRAVITY, w_distribuited, air_density, cl_max)


stall_flap_speed = speed_curve(CN_STALL, GRAVITY, w_distribuited, air_density, clmax_flapped)

manuever_speed = stall_speed * sqrt(cn_max_pos)

def flap_speed_equation(stall_flap_speed, stall_speed):
    
    # According to FAR23.345(2-b)
    # flapSpeed may be the higher value: 140% of Stall Speed or 180% Flapped Stall Speed.

    variable1 = 1.4*stall_speed
    variable2 = 1.8*stall_flap_speed

    if show_vn == str('Y'):
        print(variable1, variable2)
        
    return min(variable1, variable2)

n = CN_STALL

flap_speed = flap_speed_equation(stall_flap_speed, stall_speed)

def flap_load_factor (rho, mtow, wing_area, flap_speed, cl_max, g):
    variable1 = rho*wing_area*(flap_speed**2)*cl_max
    variable2 = 2*mtow*g
    return variable1/variable2


n_flap = flap_load_factor(air_density, mtow, wing_surface, flap_speed, clmax_flapped, GRAVITY)

n = CN_STALL

vetor3 = []
vetor4 = []

while n < n_flap:

    flapCurve = speed_curve(n, GRAVITY, w_distribuited, air_density, clmax_flapped)
    vetor3.append(flapCurve)
    vetor4.append(n)
    n += 0.05

plt.plot(vetor3, vetor4, color='magenta')


if show_vn == str('Y'):     

    def plotX():

        ult_pos_speed = speed_curve(cn_ult_pos, GRAVITY, w_distribuited, air_density, cl_max)

        vector_v_ult1 = [ult_pos_speed, diving_speed]
        vector_n_ult1 = [cn_ult_pos, cn_ult_pos]

        plt.plot(vector_v_ult1, vector_n_ult1,   '.-', color='darkorange', label='Stall Line')

        print(f"\nUltimate Positive Line (3): \n{vector_v_ult1} \n {vector_n_ult1}")


        ult_neg_speed = speed_curve(cn_ult_neg, GRAVITY, w_distribuited, air_density, cl_max)

        vector_v_ult2 = [ult_neg_speed, diving_speed]
        vector_n_ult2 = [cn_ult_neg, cn_ult_neg]
        plt.plot(vector_v_ult2,  vector_n_ult2, '.-', color='darkorange' )

        print(f'\nUltimate Negative Line (4): \n{vector_v_ult2}\n{vector_n_ult2}')

        vector_v_lim1 = [manuever_speed, diving_speed]
        vector_n_lim1 = [cn_max_pos, cn_max_pos]

        plt.plot(vector_v_lim1, vector_n_lim1, '.-', color = 'b', label='Manuever Envelope')

        print(f'\nLimit Positive Line (5): \n{vector_v_lim1}\n{vector_n_lim1}')

        vector_v_lim2 = [stall_speed, cruising_speed, diving_speed]
        vector_n_lim2 = [(-1) * CN_STALL, (-1) * CN_STALL, CN_NULL]

        plt.plot(vector_v_lim2, vector_n_lim2, '.-', color='b')

        print(f'\nLimit Negative Line (6): \n{vector_v_lim2}\n{vector_n_lim2}')

        vector_x_null = [stall_flap_speed, flap_speed]
        vector_y_null = [0,0]

        plt.plot(vector_x_null, vector_y_null, '.-', color='magenta', label='Flap Envelope')

    plotX()

    # Gusts Calculation

        
    mass_ratio = vertical_mass_ratio(mtow, aerodynamic_chord, alpha, air_density, GRAVITY)

    kg = gust_alleviation_factor(mass_ratio)


    def GustEquation(v, alpha, kg, Ude, weight_distribuition):
        
        variable1 = (v * 1.9438) * alpha * kg * Ude
        variable2 = 498 * weight_distribuition

        return variable1/variable2


    cruisingPositiveGust = 1 + GustEquation(cruising_speed, alpha, kg, Uc_gust, w_distribuited)

    cruisingNegativeGust = 1 - GustEquation(cruising_speed, alpha, kg, Uc_gust, w_distribuited)

    divingPositiveGust = 1 + GustEquation(diving_speed, alpha, kg, Ud_gust, w_distribuited)

    divingNegativeGust = 1 - GustEquation(diving_speed, alpha, kg, Ud_gust, w_distribuited)
        

    minimum_gust_speed = stall_speed * sqrt(cruisingPositiveGust)

    # Gust Vectors

    def gustPlot():

        positive_vector_cruising_gust = [1, cruisingPositiveGust]

        up_vector_cruising_speed = [0, cruising_speed]

        down_vector_cruising_speed = [diving_speed, cruising_speed]

        negative_vector_cruising_gust = [1, cruisingNegativeGust]

        down_vector_cruising_speed = [diving_speed, cruising_speed]

        vector_diving_speed = [0, diving_speed]

        n_rajada_mergulho_pos_vector = [1, divingPositiveGust]

        n_rajada_mergulho_neg_vector = [1, divingNegativeGust]

        n_rajada_combinado_cima = [divingPositiveGust, cruisingPositiveGust]

        n_rajada_combinado_baixo = [divingNegativeGust, cruisingNegativeGust]

        # Plot Rajadas

        plt.plot(up_vector_cruising_speed, positive_vector_cruising_gust, color='lime', label='Gust Envelope')
        plt.plot(up_vector_cruising_speed, negative_vector_cruising_gust, color='lime')
        plt.plot(vector_diving_speed, n_rajada_mergulho_pos_vector, color='lime')
        plt.plot(vector_diving_speed, n_rajada_mergulho_neg_vector, color='lime')
        plt.plot(down_vector_cruising_speed, n_rajada_combinado_cima, color='lime')
        plt.plot(down_vector_cruising_speed, n_rajada_combinado_baixo, color='lime')

    gustPlot()

    def plotY():

        vector_v_diving1 = [diving_speed, diving_speed]
        vector_n_diving1 = [cn_ult_neg, cn_ult_pos]

        plt.plot(vector_v_diving1, vector_n_diving1,   '.-', color='darkorange')

        print(f'\nDiving Speed Vertical Line(7): \n{vector_v_diving1}\n{vector_n_diving1}')

        vector_v_diving2 = [diving_speed, diving_speed]
        vector_n_diving2 = [0, cn_max_pos]

        plt.plot(vector_v_diving2, vector_n_diving2,   '.-', color='b')

        vector_stall = [stall_speed, stall_speed, stall_speed]
        vector_stall2 = [CN_STALL, CN_NULL, -CN_STALL]
        
        plt.plot(vector_stall, vector_stall2, '.-', color='b')

        print(f'\nStall Speed Line (8): \n{vector_stall}\n{vector_stall2}')

        vector_v_cruising = [cruising_speed, cruising_speed, cruising_speed]
        vector_n_cruising = [cn_max_pos, CN_NULL, (-1) * CN_STALL]
        plt.plot(vector_v_cruising, vector_n_cruising, '.-', color='b')
            
        print(f'\nCruising Speed Vertical Line (9): \n{vector_v_cruising}\n{vector_n_cruising}')

        vector_v_minstall = [stall_flap_speed, stall_flap_speed]
        vector_n_minstall = [0, CN_STALL]
        plt.plot(vector_v_minstall, vector_n_minstall, '.-', color='magenta')

        vector_v_flap = [flap_speed, flap_speed]
        vector_n_minstall = [0, n_flap]
        plt.plot(vector_v_flap, vector_n_minstall, '.-', color='magenta')

    plotY()

    def txt():

        print(f'\nNegative Curve (2): \n{vetor4}')

        print(f'\nMass Ratio: {mass_ratio}')

        print(f'\nGust Alleviation Factor: {kg}')

        print(f'\nUltimate Load Factor+: {cn_ult_pos}\nUltimate Load Factor-: {cn_ult_neg}')

        print(f'\nflap load factor: {n_flap}')

        print(f"\nCruising Speed: {cruising_speed} m/s")

        print(f'\nDiving Speed: {diving_speed} m/s')

        print(f'\nStall Speed: {stall_speed} m/s')

        print(f'\nFlap Speed: {flap_speed} m/s')

        print(f'\nManuever Speed: {manuever_speed} m/s')

        print(f'\nMinimum Gust Speed: {minimum_gust_speed} m/s')

        print(f'\nCruising Speed at Positive Gust:{cruisingPositiveGust}\nCruising Speed at Negative Gust:{cruisingNegativeGust}')

        print(f'\nDiving Speed at Positive Gust:{divingPositiveGust}\nDiving Speed at Negative Gust:{divingNegativeGust}')

    txt()

    plt.tight_layout()
    plt.title('V-N Envelope')
    plt.xlabel('Speed (m/s)')
    plt.ylabel('Load factor')

    plt.legend()
    plt.grid()
    plt.show()
    
