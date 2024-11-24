import numpy as np
import matplotlib.pyplot as plt

# Define translation constants
c_x = 2  # Rightward translation
c_y = 10  # Upward translation

# Define x values (only negative for the third quadrant)
x = np.linspace(-100, -0.1, 500)  # Negative values only (avoid x = 0 to prevent division by zero)

# Compute the derivative of 1/x
f_prime = -1 / x**2

# Filter to retain only the desired third quadrant
# Flip x to positive (absolute value for log scale) and apply translation
x_translated = -x + c_x  # Flip negative x to positive and translate
f_prime_translated = f_prime + c_y  # Translate upward

# Plot the result
plt.figure(figsize=(10, 6))
plt.plot(x_translated, f_prime_translated, label="Translated Derivative", color="blue")

# Set log-log scales
plt.xscale("log")
plt.yscale("log")
import numpy as np
import matplotlib.pyplot as plt

# Define translation constants
c_x = 2  # Rightward translation
c_y = 10  # Upward translation

# Define x values (negative only for the third quadrant)
x = np.linspace(-100, -0.1, 500)  # Strictly negative x-values

# Compute the derivative of 1/x
f_prime = -1 / x**2

# Transform x and y for the desired translation
x_translated = -x + c_x  # Flip x to positive and translate
f_prime_translated = f_prime + c_y  # Translate upward

# Explicitly filter out any unintended positive x-values
mask = x_translated > 0  # Ensure we only include positive x after translation

x_translated = x_translated[mask]
f_prime_translated = f_prime_translated[mask]

# Plot the result
plt.figure(figsize=(10, 6))
plt.plot(x_translated, f_prime_translated, label="Translated Derivative", color="blue")

# Set log-log scales
plt.xscale("log")
plt.yscale("log")

# Add labels and title
plt.xlabel("x (log scale)")
plt.ylabel("f'(x) (log scale)")
plt.title("Translated Derivative of f(x) = 1/x (Log-Log Scale)")
plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.legend()
plt.show()

# Add labels and title
plt.xlabel("x (log scale)")
plt.ylabel("f'(x) (log scale)")
plt.title("Translated Derivative of f(x) = 1/x (Log-Log Scale)")
plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.legend()
plt.show()
