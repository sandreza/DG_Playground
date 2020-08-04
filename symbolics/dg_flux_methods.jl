abstract type AbstractFluxMethod end
abstract type AbstractBoundaryCondition end
struct NeglectFlux  <: AbstractFluxMethod end
struct Central <: AbstractFluxMethod end
struct Upwind  <: AbstractFluxMethod end
struct Rusanov{𝒯} <: AbstractFluxMethod
    α::𝒯
end
struct RusanovBC{𝒯} <: AbstractFluxMethod
    α::𝒯
end
struct Slider{𝒯, 𝒮} <: AbstractFluxMethod
    α::𝒯
    v::𝒮
end
# Boundary Conditions
struct Dirichlet{𝒯} <: AbstractBoundaryCondition
    left::𝒯
    right::𝒯
end
# Boundary Conditions
struct Dirichlet2{𝒯} <: AbstractBoundaryCondition
    left::𝒯
    right::𝒯
end
struct FluxBC{𝒯} <: AbstractBoundaryCondition
    left::𝒯
    right::𝒯
end
struct Inflow{𝒯} <: AbstractBoundaryCondition
    in::𝒯
end
struct Inflow2{𝒯} <: AbstractBoundaryCondition
    in::𝒯
end
struct Outflow{𝒯} <: AbstractBoundaryCondition
    out::𝒯
end
struct Neumann{𝒯} <: AbstractBoundaryCondition
    left::𝒯
    right::𝒯
end
struct Periodic <: AbstractBoundaryCondition end
struct NoFlux   <: AbstractBoundaryCondition end
struct FreeFlux <: AbstractBoundaryCondition end