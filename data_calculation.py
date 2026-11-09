import csv
import math

n = 600
m = 1
delta_a, delta_d = 0.2, 0.1

data = [
    ["True", 12.0, 5.0],
    ["True", 12.0, 3.6],
    ["True", 12.0, 3.2],
    ["Virtual", 26.0, 14.0],
    ["Virtual", 26.0, 9.6],
    ["Virtual", 26.0, 8.1],
]

results = []
for row in data:
    spec, a, d = row
    lam = d / (m * n * math.sqrt(a**2 + d**2)) * 10**6
    delta_lam = lam * ((delta_d / d) + ((a * delta_a + d * delta_d) / (a**2 + d**2)))
    error = (delta_lam/lam) * 100
    results.append([spec, a, delta_a, d, delta_d, round(lam,2), round(delta_lam,2), round(error,2)])

for i in range(3):
    n = (3, 4, 5)
    ridberg = (1.096 * (1 / (4) - 1 / (n[i] **2)))**(-1) * 100
    results.append(["Ridberg"," ", " ", " ", " ", round(ridberg,2) , " "])

with open("results.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Spectrum","a","Delta_a","d","Delta_d","Lambda","Delta_Lambda","Error"])
    writer.writerows(results)
