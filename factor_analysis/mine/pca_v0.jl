using DelimitedFiles
using Statistics
using LinearAlgebra
using MultivariateStats
using Plots
pyplot(reuse=false, size=(600,450))
closeall()
clearconsole()

# Define working path
path = "/home/leniac/Work/2019/modos_sup/factor_analysis/"
cd(path)

## Shaggy
# Load data
X = readdlm("MatrizDatos.dat")'
times = readdlm("ti2.dat")[:]
ω = 1:size(X,1)
# Data augmentation
# X2 = repeat(X, 1, 3)

# Train a PCA model
M = fit(PCA, X; maxoutdim=5, method=:svd)

# apply FactorAnalysis model to testing set
Yte = transform(M, X)

# reconstruct testing observations (approximately)
Xr = reconstruct(M, Yte)

plot(ω, X'); gui()
plot(ω, Xr'); gui()

W = projection(M)

plot(ω, Yte')
plot(times, W)
