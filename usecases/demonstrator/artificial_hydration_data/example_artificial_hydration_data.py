import yaml
import fenics_concrete

data_file = 'artificial_hydration_data.yaml'
#Example 1:
# read file and access artificial data:
with open(data_file) as file:
    hydration_data = yaml.safe_load(file)

# data is given in dictionary
# data[mix ratio: 0/.2/0.5/.8/1][temperature: 20/40/60][time/heat]
# it is assumed that this is hydration data for two distinct mixes
# mix 1: mix ration = 0 and mix 2: mix ratio = 1
# there are 3 intermediate mixes with 20/80, 50/50 and 80/20 ratio between mix 1 and 2
# for each of the 5 mixes there are 3 temperature measurements, each at 20, 40 and 60 degree
# for each temperature there is a list with the time and the heat values

# loop over all data, print lists
for mix_r in hydration_data:
    for temp in hydration_data[mix_r]:
        print(mix_r,temp,'time:',hydration_data[mix_r][temp]['time'])
        print(mix_r,temp,'heat:',hydration_data[mix_r][temp]['heat'])

# Example 2:
# use the hydration model with a specific set of parameters
# these parameters are roughly correct for mix 1 (mix_r=0)
# mix 2 will have different values but I assume them to be in a similar order of magnitude
# alpha_max and T_ref should NOT be estimated or changed
# for now we assume  them to be constant, you can also try to use alpha_max = 1 to see what happens

# Some value within which I would expect the dataset to be
parameter = fenics_concrete.Parameters()  # using the current default values

parameter['B1'] = 2.916E-4  # in 1/s (le 0, < 0.1)
parameter['B2'] = 0.0024229  # - (le 0, smaller 1)
parameter['eta'] = 5.554  # something about diffusion (should be larger 0)
parameter['alpha_max'] = 0.875  # also possible to approximate based on equation with w/c (larger 0 and max 1)
parameter['E_act'] = 47002   # activation energy in Jmol^-1 (no relevant limits)
parameter['T_ref'] = 25  # reference temperature in degree celsius
parameter['igc'] = 8.3145  # ideal gas constant in [J/K/mol], CONSTANT!!!
parameter['zero_C'] = 273.15  # in Kelvin, CONSTANT!!!
parameter['Q_pot'] = 500e3 # potential heat per weight of binder in J/kg

# this is the minimal time step used in the simulation
# using a larger value will increase the speed but decrease the accuracy
dt = 300 # value in seconds

# this is the simulated temperature, needs to be adjusted depending on the temperature of the experimental data
T = 30
# this is the list of measured time data as given by the experiments
time_list = [0,5000,10000,20000,100000]

# initiate material problem, for this the "fenics_concrete" conda package needs to be installed
# use: 'mamba install -c etamsen fenics_concrete"
problem = fenics_concrete.ConcreteThermoMechanical()

# get the hydration function
# this might change in the future to make it more easily accessible but for now it should work like this
hydration_fkt = problem.get_heat_of_hydration_ftk()
# the results are a heat list and a degree of hydration list, which you can ignore for now
heat_list, doh_list= hydration_fkt(T, time_list, dt, parameter)

# the results!!!
print(heat_list)

# x = np.array([0.3])
# b = np.array([2.916,2.4229,5.554,5])
# # convert the above to run the pretraining
# x = torch.tensor(x).reshape(1,-1)
# b = torch.tensor(b).reshape(1,-1)
# # convert the above to a 1x4 array

# # reshape the above to 2d tensor
# #b = torch.tensor(b).reshape(1,-1)
# #x = torch.tensor(x).reshape(1,-1)

# nn_pretrained = pretrain_nn_mean(x, b, epochs=100, lr=1e-3, hidden_dim=10)