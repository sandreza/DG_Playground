"""
"""

export ProductDomain
export IntervalDomain, PeriodicIntervalDomain
export ProductDomain
export RectangularDomain
export BoxDomain
export SphereSurfaceDomain

export UnitInterval, UnitSquare, UnitCube
export SphericalAnnulus

export get_geometric_dimension, get_topological_dimension

"""
"""
abstract type AbstractDomain end

function get_geometric_dimension end
function get_topological_dimension end

"""
"""
struct IntervalDomain{FT} <: AbstractDomain
    length::FT
    function IntervalDomain(a, b)
        @assert a < b
        length = abs(b - a)
        return new{typeof(length)}(length)
    end
end

"""
"""
struct PeriodicIntervalDomain{FT} <: AbstractDomain
    length::FT
    function PeriodicIntervalDomain(length)
        @assert length > 0
        return new{typeof(length)}(length)
    end
end

get_geometric_dimension(::Union{IntervalDomain, PeriodicIntervalDomain}) = 1
get_topological_dimension(::Union{IntervalDomain, PeriodicIntervalDomain}) = 1

"""
"""
struct ProductDomain{DT} <: AbstractDomain
    domains::DT
end
×(args::AbstractDomain...) = ProductDomain{typeof(args)}(args)

function get_geometric_dimension(pd::ProductDomain)
    return sum(map(get_geometric_dimension, pd.domains))
end
function get_topological_dimension(pd::ProductDomain)
    return sum(map(get_topological_dimension, pd.domains))
end

"""
"""
struct SphereSurfaceDomain{FT} <: AbstractDomain
    r::FT
    SphereSurfaceDomain(r) = new{typeof(r)}(r)
end

get_geometric_dimension(::SphereSurfaceDomain) = 2
get_topological_dimension(::SphereSurfaceDomain) = 3


RectangularDomain(a, b, c, d) = IntervalDomain(a, b) × IntervalDomain(c, d)
BoxDomain(a, b, c, d, e, f) = IntervalDomain(a, b) × IntervalDomain(c, d) × IntervalDomain(e, f)

UnitInterval() = IntervalDomain(0, 1)
UnitSquare() = RectangularDomain(0, 1, 0, 1)
UnitCube() = BoxDomain(0, 1, 0, 1, 0, 1)

SphericalAnnulus(r, H) = SphereSurfaceDomain(r) × IntervalDomain(0, H)

# TODO: Think about metadata for labeling boundaries
