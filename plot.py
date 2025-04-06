from matplotlib import pyplot as plt
import pandas as pd

df = pd.read_csv('./casova_narocnost.csv')

# Plot the data with 'n' as logarithmic scale
plt.scatter(df['n'], df['elapsedTime'], marker='+')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Počet rovnic')
plt.grid(which='both', linestyle='--', linewidth=0.5, alpha=.5)
# plt.tight_layout()
plt.ylabel('Časová náročnost (s)')
plt.savefig('casova_narocnost.png')
