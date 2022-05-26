using DelimitedFiles
# using Statistics
# using LinearAlgebra
# using MultivariateStats
using Plots
pyplot(reuse=false, size=(600,450))
using LaTeXStrings
closeall()
clearconsole()

# Define working path
path = "/home/leniac/Work/2019/modos_sup/factor_analysis/"
cd(path)
include("PFAMalinowski.jl")
using Main.PFAMalinowski: PFA, PFATT

# Load data
X = readdlm("MatrizDatos.dat")
Xerr = readdlm("MatrizError.dat")
times = readdlm("ti2.dat")[:]
ω = 1:size(X,1)

## PFA preliminary: determine the number of components
# If you know them beforehand ignore this section and select k

# PFA to the error matrix
_ = PFA(Xerr)

# # PFA to the data matrix
_ = PFA(X)

## PFA to the data matrix with k factors: abstract
# Number of fators chosen
k = 2
results1 = PFA(X, k)

## Target transformation: real

# Spectra for which you have the targets
spectra_target = [1; size(X, 2)]

# Initial concentration of the selected spectra
ci = [1.0; 0.0]

# Perform transformation
results2 = PFATT(spectra_target, ci, X, results1, k)

## Plots
wnstr = L"Wavenumber [cm$^{-1}$]"
abstrs = "Absorbance [a. u.]"
timestr = "Time [min]"

plt1 = plot(ω, results1.matrixResults.Xreprod, title=("Reproduced abstract data"), yaxis=(abstrs), xaxis=(wnstr), legend=false);

plt2 = plot(ω, results2.realFac, title=("Real principal factors"), yaxis=(abstrs), xaxis=(wnstr), legend=false);

plt3 = plot(times, results2.realConc', line=(false), marker=(:auto), title=("Real principal loadings/concentrations"), yaxis=("[units]"), xaxis=(timestr), legend=false);

plt4 = plot(ω, results2.realSpectraReprod, title=(@sprintf("Reproduced real data, RMS = %.2f%s", results2.rms*100, "%")), yaxis=(abstrs), xaxis=(wnstr), legend=false);

plot(plt1, plt2, plt3, plt4, layout = (2,2)); gui()
