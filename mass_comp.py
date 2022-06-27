import sys
import Hirogen_Functions
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('masses.csv')

df = df.dropna()

df = df[df.lgm_tot_p50 > 0]
df = df[df.lgm_fib_p50 > 0]

sf_port_mass = df['SF_logMass']
pass_port_mass = df['logMass']
mpa_mod_mass = df['lgm_tot_p50']
mpa_fib_mass = df['lgm_fib_p50']



fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2, figsize = (16,14), sharex=True, sharey=True)

ax1.scatter(sf_port_mass, mpa_mod_mass)
ax2.scatter(sf_port_mass, mpa_fib_mass)
ax3.scatter(pass_port_mass, mpa_mod_mass)
ax4.scatter(pass_port_mass, mpa_fib_mass)

ax1.set_xlabel(r'Portsmouth Starforming Mass, $\log(M_{\odot})$')
ax1.set_ylabel(r'MPA-JHU Model Mass, $\log(M_{\odot})$')

ax2.set_xlabel(r'Portsmouth Starforming Mass, $\log(M_{\odot})$')
ax2.set_ylabel(r'MPA-JHU Fibre Mass, $\log(M_{\odot})$')

ax3.set_xlabel(r'Portsmouth Passive Mass, $\log(M_{\odot})$')
ax3.set_ylabel(r'MPA-JHU Model Mass, $\log(M_{\odot})$')

ax4.set_xlabel(r'Portsmouth Passive Mass, $\log(M_{\odot})$')
ax4.set_ylabel(r'MPA-JHU Fibre Mass, $\log(M_{\odot})$')

plt.show()