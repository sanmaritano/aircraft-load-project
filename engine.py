from inputs import *
from vn_diagram import *

print('follow FAR 23.361 for your engine load conditions.')

print()
K = int(input("Select the factor for mean torque according to your engine: "))

def torque (K, factor, maximumPower, RPM):
  variable1 = (maximumPower/RPM)
  return K*factor*variable1

f = torque(K, CONSTANT_TORQUE, maximum_power, maximum_rot)
print(f'Limit torque: {round(f, 4)} Kgf.m')

# carga limite no motor

propeller_efficiency = 0.5304 

fz = cn_max_pos*engine_mass

fxp = round(propeller_efficiency*(maximum_power/manuever_speed)*76 , 4)

#Cargas Laterais
ny = cn_max_pos/3
fy = ny * engine_mass

print(f'\nX axis load = {fxp} Kgf\nY axis load = {round(fy, 4)} Kgf\
\nZ axis load = {fz} Kgf')
