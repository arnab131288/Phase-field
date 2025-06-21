"""
Configuration module for Phase Field Simulation.
Defines all simulation parameters and constants derived from them.
"""
params = {
    'NX': 400,                           # Mesh size in x-direction 
    'NY': 100,                           # Mesh size in y-direction
    'DELTA_X': 0.03,                     # Grid spacing in x-direction
    'DELTA_Y': 0.03,                     # Grid spacing in y-direction
    'DELTA_T': 1e-5,                     # Temporal spacing
    'TOTAL_TIMESTEPS': 50000,            # Total simulation timesteps
    'TIME_BREAK': 1000,                  # Timesteps after which output csv files are written
    'TAU': 0.0003,                       # Interface relaxation parameter  
    'EPSILON': 0.01,                     # Interface width
    'T_E': 1.0,                          # Eqbm temp 
    'ALPHA': 0.9,                        # Coupling constant controlling the driving force
    'GAMMA': 10.0,                       # Constant controlling strength of driving force
    'A': 0.01,                           # Magnitude of noise strength
    'K': 0.9,                            # Coupling constant related to latent heat
    'DELTA': 0.0,                        # Strength of anisotropy
    'FOLD': 4,                           # Symmetry of anisotropy
    'T_FILL': 0.0,                       # Initial filling of the temp field
    'PHI_FILL': 1.0,                     # Initial filling of the phi field 
    'OUT_DIR': '/home/arnab/workspace/kobayashi/python/output'                  # Path to save output directory
}

# Derived constants
rdX = 1.0 / params['DELTA_X']
rdY = 1.0 / params['DELTA_Y']
dx2 = rdX * rdX
dy2 = rdY * rdY

# Exported constants for numba routines (avoid dict import issues with numba)
NX      = params['NX']
NY      = params['NY']
DELTA_X = params['DELTA_X']
DELTA_Y = params['DELTA_Y']
DELTA_T = params['DELTA_T']
TOTAL_TIMESTEPS = params['TOTAL_TIMESTEPS']
TIME_BREAK = params['TIME_BREAK']
TAU     = params['TAU']
EPSILON = params['EPSILON']
T_E     = params['T_E']
ALPHA   = params['ALPHA']
GAMMA   = params['GAMMA']
A       = params['A']
K       = params['K']
DELTA   = params['DELTA']
FOLD    = params['FOLD']
T_FILL  = params['T_FILL']
PHI_FILL = params['PHI_FILL']

