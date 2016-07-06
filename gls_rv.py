import numpy as np
import matplotlib.pyplot as plt

# this requires Astropy 1.2.1
from astropy.stats import LombScargle


time, rv, rv_error = np.loadtxt('EPIC-9792_SOPHIE.rdb', skiprows=2, unpack=True)


fig = plt.figure()
ax = fig.add_subplot(111)
ax.errorbar(time, rv, rv_error, fmt='o')
ax.set_xlabel('Time [days]')
ax.set_ylabel('RV [km/s]')
plt.show()



# using astropy's Lomb-Scargle periodogram is sooo easy...
frequency, power = LombScargle(time, rv, rv_error).autopower()



fig = plt.figure()
ax = fig.add_subplot(111)
ax.semilogx(1/frequency, power, 'k-', lw=2)
ax.set_xlabel('Period [days]')
ax.set_ylabel('Power')
plt.show()


# determine the period of maximum power
ind = power.argmax()
print 'period of maximum power:', (1/frequency)[ind], 'days'