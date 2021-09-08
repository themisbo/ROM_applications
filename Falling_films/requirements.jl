using Pkg

dependencies = [
    "FFTW",
    "Parameters",
    "DifferentialEquations",
    "Plots",
    "Statistics",
    "PyCall"
]

Pkg.add(dependencies)