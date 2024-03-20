@with_kw struct theory
    Nc::Float64
    Nf::Float64
    g2::Float64
    mQ::Float64 
    # IR power law coefficient
    κGl::Float64 = 1
    κGh::Float64 = 0
    # generic Casimirs
    CA::Float64
    CF::Float64 
    T::Float64  
    # general γ aus Muta 3.4.18
    γ::Float64 = (-13CA + 8*T*Nf)/(22CA - 8*T*Nf)
    δ::Float64 = (-1-γ)/2
    # exponents for three-gluon vertex model
    a::Float64 # scaling - for decoupling choose a = -1
    α::Float64 = 1-a/δ -2a
    β::Float64 = -1 - a
end
# Casimirs SU(N) fundamental 
#CA::Float64 = Nc
#CF::Float64 = (Nc^2-1)/(2Nc)
#T::Float64  = 1/2
# Casimirs SU(N) adjoint 
#CA::Float64 = Nc
#CF::Float64 = Nc
#T::Float64  = Nc
# Casimirs Sp(2N) = Sp(Nc) fundamental 
#CA::Float64 = Nc/2+1
#CF::Float64 = (Nc+1)/4
#T::Float64  = 1/2
# Casimirs SO(N) fundamental 
#CA::Float64 = Nc-2
#CF::Float64 = (Nc-1)/2
#T::Float64  = 1
@with_kw struct parameter
    #renormalization
    µ::Float64 = 45000
    ZUV::Float64 = 1 # value of dressing at gluon subtraction point
    GIR::Float64 = 15 # value of the ghost dressing at the lowest external momentum
    # iteration
    maxiterNR::Int = 50 # maximal number of iterations using Newton's method
    linesearch::Int = 10 # maximal number of steps in line search
    epsilonNR::Float64 = 10^(-11) #convergence criterion for Newton's method
    quarkrelax::Float64 = 1 # relaxation parameter for quark iteration
    h::Float64 = 10^(-9) #step size for forward differentiation in Jacobian
    propiter::Int = 50
    # integration and interpolation
    nr::Int = 200 # radial integration points
    na::Int = 50  # angular integration points
    nrQ::Int = 300 # radial integration points
    naQ::Int = 100  # angular integration points
    nCheby::Int = 80  # number of Chebyshev polynomials for expansion
    nQ::Int = 500 # external momenta for quark integration
    # external momentum range and integeration range
    IR::Float64 = 10^(-8) # IR cutoff
    UV::Float64 = 50000 # UV cutoff
    ϵ::Float64 = 10^(-5) # smallest external momentum
    Λ::Float64 = UV*0.95 # highest external momentum
    # quark iteration
    maxiterQ::Int = 400
    ϵQ::Float64 = 10^-14
    ϵtotal::Float64 = 10^-7
    splinesQ::Bool = true
    # 3-Gluon Vertex Model
    GluonModel::Symbol = :CSF
    hIR::Float64 = -1
    Λ3G::Float64 = 1 #quenched
end
#-------------------------------------------------------------------------------
# structures for propagator calculation
struct QuadratureProp
    p2::Array{Float64,1}
    k2::Array{Float64,2}
    q2::Array{Float64,3}
    RadWeights::Array{Float64,2} # includes jacobian of transformation to [IR,UV]
    AngWeights::Array{Float64,1}
    k2_inv::Array{Float64,2}
    q2_inv::Array{Float64,3}
    # arrays for Huber/Smekal-model
    s3::Array{Float64,3}
    s6::Array{Float64,3}
    s3_inv::Array{Float64,3}
    s6_inv::Array{Float64,3}
    damping::Array{Float64,3}
end
mutable struct RenormalizationConstants
    Z2::Float64
    Zm::Float64
    Z3::Float64
    Z3_tilde::Float64
end
mutable struct InterpolatedYM
    Gy::Array{Float64,2}
    Zy::Array{Float64,2}
    Gz::Array{Float64,3}
    Zz::Array{Float64,3}
    Gs3::Array{Float64,3}
    Zs3::Array{Float64,3}
    Gs6::Array{Float64,3}
    VertexHS::Array{Float64,3}
end
struct InterpolatedQuark
    Ay::Array{Float64,2}
    By::Array{Float64,2}
    Az::Array{Float64,3}
    Bz::Array{Float64,3}
end

mutable struct KernelsProp
    Ghost::Array{Float64,3}
    GhostLoop::Array{Float64,3}
    GluonLoop::Array{Float64,3}
    QuarkLoop::Array{Float64,3}
end
#------------------------------------------------------------------------------
# structures for vertex calculation
@with_kw struct parameterGGV
    nr::Int = 480 # needs to be a multiple of 6
    na1::Int = 25
    na2::Int = 25
    IR::Float64 = 10^(-8)
    UV::Float64 = 10^(13)
    maxiter::Int = 1
    converged::Float64 = 10^(-4)
    sysitermax::Int = 20 # maximum # of iterations for system
    eps_sys::Float64 = 10^-5 # convergence criterion for system
    rel::Float64 = 1 # relaxation parameter
    # spacing of external momenta
    Δxy::Float64   = 0.20
    Δcosθ::Float64 = 0.1111
    logϵ::Float64 = -7.01
    logΛ::Float64 = +8.00
end
struct QuadratureGGV
    a::Array{Float64,1}
    cosθ₂::Array{Float64,1}
    cosθ₁::Array{Float64,1}
    RadWeights::Array{Float64,1}
    θ₁Weights::Array{Float64,1}
    θ₂Weights::Array{Float64,1}
end
mutable struct InterpolationPointsGGV  #structur containing interpolation points
    I1::Array{Float64,1}
    I2::Array{Float64,2}
    I3::Array{Float64,3}
    I4::Array{Float64,2}
    I5::Array{Float64,3}
    I6::Array{Float64,1}
    newI1::Array{Float64,2}
    newI2::Array{Float64,3}
    newI6::Array{Float64,2}
end
mutable struct KernelsGGV
    Abelian::Array{Float64,3}
    NonAbelian::Array{Float64,3}
end
mutable struct IndexEqualGGV
    xind1::Array{Int64,3} # trilinear index
    yind1::Array{Int64,3} # trilinear index
    zind::Array{Int64,3}  # trilinear index
    indI2::Array{Int64,2} # bilinear index
    indI4::Array{Int64,2} # bilinear index
    newIndI1::Array{Int64,2} # bilinear index
    newIndI6::Array{Int64,2} # bilinear index
    indI1::Array{Int64,1} # 1dim index
    indI6::Array{Int64,1} # 1dim index
    indp2I1::Array{Int64,1} # for propagators
    indp2I2::Array{Int64,2} # for propagators
    indp2I3::Array{Int64,3} # for propagators
end
mutable struct InterpolatedGGV
    AI1::Array{Float64,2}
    AI2::Array{Float64,3}
    AI3::Array{Float64,2}
    AI4::Array{Float64,3}
    GI1::Array{Float64,1}
    GI2::Array{Float64,2}
    GI3::Array{Float64,3}
    ZI1::Array{Float64,1}
    ZI2::Array{Float64,2}
    ZI3::Array{Float64,3}
end
