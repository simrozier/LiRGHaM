########################################
# LOADING ALL NEEDED PACKAGES
########################################

using Cuba # To compute 3D integrals.
using BenchmarkTools # To have precise benchmark measurements
using ArgParse # To be able to parse command-line arguments
# using HDF5 # To be able to dump/read HDF5 files
# using GSL # To have access to Bessel functions
using Interpolations # To make interpolations
using StaticArrays # To have access to fast static arrays
using LinearAlgebra # To have access to core algebra functions
# using Base.Threads # To parallelise over threads
using HypergeometricFunctions # To use the 2F1 hypergeometric function in the DF. From https://github.com/JuliaMath/HypergeometricFunctions.jl. Copyright (c) 2018 JuliaApproximation
using QuadGK # To compute integrals using the Gauss-Kronrod quadrature
using SpecialFunctions # Methods for Bessel functions
using FunctionZeros # Methods for the zeros of Bessel functions
using Jacobi # For Gegenbauer polynomials
# using Plots # For plotting
using DelimitedFiles # For dropping data
using Roots # Numerical root finding
# using SphericalHarmonics # Computing spherical harmonics
using Cubature # Computing 2D integrals
using AssociatedLegendrePolynomials  # Time and memory efficient spherical harmonics