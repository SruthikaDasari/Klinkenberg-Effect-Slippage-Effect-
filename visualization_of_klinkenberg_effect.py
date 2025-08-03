import numpy as np
import matplotlib.pyplot as plt

print("=== KLINKENBERG EFFECT ANALYSIS ===")

mean_pressure = float(input("Enter mean pressure (Pm): "))
gas_perm = float(input("Enter gas permeability (Kog): "))
guess_liquid_perm = float(input("Enter initial guess for liquid permeability (Kl): "))

def klinkenberg_function(Kl, Pm, Kog):
    return 6.9 * Kl**0.64 + Pm * Kl - Pm * Kog

def klinkenberg_derivative(Kl, Pm):
    return 4.416 * Kl**(-0.36) + Pm

Kl = guess_liquid_perm
tolerance = 1e-7
max_iterations = 1000
iteration = 0

while abs(klinkenberg_function(Kl, mean_pressure, gas_perm)) > tolerance and iteration < max_iterations:
    derivative = klinkenberg_derivative(Kl, mean_pressure)
    if derivative == 0:
        print("Zero derivative encountered. Exiting.")
        break
    Kl = Kl - (klinkenberg_function(Kl, mean_pressure, gas_perm) / derivative)
    if Kl <= 0:
        Kl = 1e-5
    iteration += 1

if iteration == max_iterations:
    print("\nWarning: Newton-Raphson did not converge within the limit.")
else:
    print(f"\nEstimated Liquid Permeability (Kl): {Kl:.4f} after {iteration} iterations")

lab_Kl = float(input("Enter lab-measured liquid permeability (Klo): "))
difference = Kl - lab_Kl
print(f"Difference from lab value: {difference:.4f}")

b = 6.9 * Kl**(-0.36)

inv_Pm = np.linspace(0.01, 1.0, 100)
Kog = Kl + b * Kl * inv_Pm

plt.figure(figsize=(10, 6))
plt.plot(inv_Pm, Kog, label="Kog = Kl + b*Kl*(1/Pm)", color='blue')
plt.plot(inv_Pm, np.full_like(inv_Pm, Kl), label="Klo = Kl (constant)", color='green', linestyle='--')
plt.axvline(1/mean_pressure, color='purple', linestyle='--', label=f'Current 1/Pm = {1/mean_pressure:.3f}')
plt.scatter([1/mean_pressure], [gas_perm], color='red', zorder=5, label=f"Measured Kog = {gas_perm}")
plt.title("Klinkenberg Effect: Gas Slippage Correction Curve")
plt.xlabel("1 / Mean Pressure (1/Pm)")
plt.ylabel("Permeability")
plt.grid(True)
plt.axhline(0, color='black', linewidth=0.5)
plt.axvline(0, color='black', linewidth=0.5)
plt.legend()
plt.tight_layout()
plt.show()
