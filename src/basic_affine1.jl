"""
Affine form with center `c`, affine components `γ` and error `Δ`. This
is the AF1 form described in Messine 2002.

Variant where Δ is an interval
"""
struct BasicAffine1{N,T<:Real,S<:AbstractLinearization}
    c::T               # mid-point
    γ::SVector{N,T}    # affine terms
    Δ::Interval{T}     # error term
end

function Base.show(io::IO, C::BasicAffine1{N,T,S}) where {N,T,S}
    print(io, "⟨", C.c, "; ", C.γ, "; ", C.Δ, "⟩")
end

# ==(C::BasicAffine1ine, D::BasicAffine1ine) = C.c == D.c && C.γ == D.γ

"""
Make a `BasicAffine1` based on an interval, which is number `i` of `n` total variables
specifying a linearization mode of type `S`.
"""
function BasicAffine1(X::Interval{T}, ::Val{N}, i::Int, ::S) where {N,T,S}
    c = mid(X)
    r = radius(X)
    γ = SVector{N,T}(ntuple(j->i==j ? r : zero(r), N))
    return BasicAffine1{N,T,S}(c, γ, Interval{T}(0))
end

"""
Make a `BasicAffine1` based on an interval, which is number `i` of `n` total variables
using the default linearization mode.
"""
function BasicAffine1(X::Interval{T}, ::Val{N}, i::Int) where {N,T}
    return BasicAffine1(X, Val(N), i, DEFAULT_LIN())
end


## Promotion rules
promote_rule(::Type{BasicAffine1{N,T,S}}, ::Type{BasicAffine1{N,Q,S}}) where {N, T<:Real, Q<:Real, S} =
    BasicAffine1{N, promote_type(T,Q), S}

promote_rule(::Type{BasicAffine1{N,T,S}}, ::Type{Q}) where {N, T<:Real, Q<:Real, S} =
    BasicAffine1{N, promote_type(T,Q), S}

promote_rule(::Type{BigFloat}, ::Type{BasicAffine1{N,T,S}}) where {N, T<:Real,S} =
    BasicAffine1{N, promote_type(T,BigFloat), S}


# convert methods:
convert(::Type{BasicAffine1{N,T,S}}, x::Bool) where {N,T,S} = convert(BasicAffine1{N,T,S}, Int(x))
convert(::Type{BasicAffine1{N,T,S}}, x::Real) where {N,T,S} = BasicAffine1{N,T,S}(atomic(Interval{T}, x))
convert(::Type{BasicAffine1{N,T,S}}, x::T) where {N,T<:Real,S} = BasicAffine1{N,T,S}(x)
convert(::Type{BasicAffine1{N,T,S}}, x::Interval{T}) where {N,T,S} = x
convert(::Type{BasicAffine1{N,T,S}}, x::Interval) where {N,T,S} = BasicAffine1{N,T,S}(atomic(Interval{T}, x))

for f in (:+,:-)
    @eval function ($f)(x::BasicAffine1{N,T,S}, y::BasicAffine1{N,T,S}) where {N,T<:Real,S<:AbstractLinearization}
        BasicAffine1{N,T,S}(($f).(x.c, y.c), ($f).(x.γ,y.γ), ($f).(x.Δ,y.Δ))
    end
end


interval(C::BasicAffine1{N,T,S}) where {N,T,S} = C.c + sum(abs.(C.γ))*(-1..1) + C.Δ


function *(x::BasicAffine1{N,T,S}, y::BasicAffine1{N,T,S}) where {N,T,S}
    c = x.c * y.c

    γ = x.c .* y.γ + y.c .* x.γ

    Δ = (x.γ ⋅ y.γ) * (0..1)  # ϵ_i^2

    if N > 1
        Δ += sum(x.γ[i] * y.γ[j] for i in 1:N, j in 1:N if i ≠ j) * (-1..1)  # ϵ_i * ϵ_j
    end

    Δ += (x.c + sum(abs.(x.γ))*(-1..1)) * y.Δ
    Δ += (y.c + sum(abs.(y.γ))*(-1..1)) * x.Δ

    Δ += x.Δ * y.Δ

    return BasicAffine1(c, γ, Δ)

end

# Scalar Operations
*(x::BasicAffine1{N,T,S}, α::T) where {N,T<:Real,S<:AbstractLinearization} = BasicAffine1{N,T,S}(α*x.c, α.*x.γ, α*x.Δ)
*(α::T, x::BasicAffine1{N,T,S}) where {N,T<:Real,S<:AbstractLinearization} = x * α

+(x::BasicAffine1{N,T,S}, α::T) where {N,T<:Real,S<:AbstractLinearization} = BasicAffine1{N,T,S}(α+x.c, x.γ, x.Δ)
+(α::T, x::BasicAffine1{N,T,S}) where {N,T<:Real,S<:AbstractLinearization} = x + α

-(x::BasicAffine1{N,T,S}) where {N,T<:Real,S<:AbstractLinearization} = BasicAffine1{N,T,S}(-x.c, .-(x.γ), -x.Δ)
-(x::BasicAffine1{N,T,S}, α::T) where {N,T<:Real,S<:AbstractLinearization} = BasicAffine1{N,T,S}(x.c - α, x.γ, x.Δ)
-(α::T, x::BasicAffine1{N,T,S}) where {N,T<:Real,S<:AbstractLinearization} = α + (-x)

/(x::BasicAffine1{N,T,S}, α::T) where {N,T<:Real,S<:AbstractLinearization} = BasicAffine1{N,T,S}(x.c/α, x.γ/α, x.Δ/α)
/(α::T, x::BasicAffine1{N,T,S}) where {N,T<:Real,S<:AbstractLinearization} = α*inv(x)


function power_by_squaring(x, p::Integer)
    if p == 1
        return copy(x)
    elseif p == 0
        return one(x)
    elseif p == 2
        return x*x
    elseif p < 0
        isone(x) && return copy(x)
        isone(-x) && return iseven(p) ? one(x) : copy(x)
        Base.throw_domerr_powbysq(x, p)
    end
    t = trailing_zeros(p) + 1
    p >>= t
    x = x
    while (t -= 1) > 0
        x *= x
    end
    y = x
    while p > 0
        t = trailing_zeros(p) + 1
        p >>= t
        while (t -= 1) >= 0
            x *= x
        end
        y *= x
    end
    return y
end

# Ran into the following error so re-defined power by squaring:
#=
ERROR: LoadError: MethodError: Cannot `convert` an object of type
  BasicAffine1{2, Float64, MinRange} to an object of type
  Union{}
Closest candidates are:
  convert(::Type{Union{}}, ::Any) at essentials.jl:203
  convert(::Type{var"#s53"} where var"#s53"<:Polynomials.AbstractPolynomial, ::Any) at C:\Users\wilhe\.julia\packages\Polynomials\1aa8e\src\common.jl:309
  convert(::Type{T}, ::Any) where T<:VecElement at baseext.jl:8
  ...
Stacktrace:
  [1] convert(#unused#::Core.TypeofBottom, x::BasicAffine1{2, Float64, MinRange})
    @ Base .\essentials.jl:203
  [2] to_power_type(x::BasicAffine1{2, Float64, MinRange})
    @ Base .\intfuncs.jl:240
  [3] power_by_squaring(x_::BasicAffine1{2, Float64, MinRange}, p::Int64)
    @ Base .\intfuncs.jl:255
  [4] ^(x::BasicAffine1{2, Float64, MinRange}, n::Int64)
    @ AffineArithmetic ~\Desktop\Package Development\AffineArithmetic.jl\src\basic_affine1.jl:14
=#
function ^(x::BasicAffine1{N,T,S}, n::Integer) where {N,T<:Real,S<:AbstractLinearization}

    invert = false

    if n < 0
        invert = true
        n = -n
    end

    result = power_by_squaring(x, n)

    if invert
        result = inv(result)
    end

    return result
end

Base.literal_pow(::typeof(^), x::BasicAffine1, ::Val{p}) where {T,p} = x^p

"General formula for affine approximation of nonlinear functions"
function affine_approx(x::BasicAffine1{N,T,S}, α, ζ, δ) where {N,T<:Real,S<:AbstractLinearization}

    c = α * x.c + ζ
    γ = α .* x.γ
    δ += α * x.Δ  # interval

    return BasicAffine1{MinRange,N,T}(c, γ, δ)
end
# to_power_type(x) = convert(Base._return_type(*, Tuple{typeof(x), typeof(x)}), x)
# to_power_type(x) = convert(Base._return_type(*, Tuple{typeof(x), typeof(x)}), x)
#=
function power_by_squaring(x_, p::Integer)
    x = to_power_type(x_)
    if p == 1
        return copy(x)
    elseif p == 0
        return one(x)
    elseif p == 2
        return x*x
    elseif p < 0
        isone(x) && return copy(x)
        isone(-x) && return iseven(p) ? one(x) : copy(x)
        throw_domerr_powbysq(x, p)
    end
    t = trailing_zeros(p) + 1
    p >>= t
    while (t -= 1) > 0
        x *= x
    end
    y = x
    while p > 0
        t = trailing_zeros(p) + 1
        p >>= t
        while (t -= 1) >= 0
            x *= x
        end
        y *= x
    end
    return y
end
=#
function Base.sqrt(x::BasicAffine1{N,T,MinRange}, X=interval(x)) where {N,T<:Real,S<:AbstractLinearization}

    a, b = X.lo, X.hi

    # min-range:  de Figuereido book, pg. 64
    α = 1/(2*√b)
    ζ = (√b)/2
    δ = ((√b - √a)^2/(2*√b))*(-1..1)

    return affine_approx(x, α, ζ, δ)
end


function Base.inv(x::BasicAffine1{N,T,MinRange}, X=interval(x)) where {N,T<:Real,S<:AbstractLinearization}

    a, b = X.lo, X.hi

    # min-range:  de Figuereido book, pg. 70
    α = -1/(b^2)
    d = interval(1/b - α*b, 1/a - α*a)
    ζ = mid(d)
    δ = radius(d) * (-1..1)

    if a < 0
        ζ = -ζ
    end

    return affine_approx(x, α, ζ, δ)
end

function Base.inv(x::BasicAffine1{N,T,Chebyshev}, X=interval(x)) where {N,T<:Real,S<:AbstractLinearization}

    a, b = X.lo, X.hi

    # min-range:  de Figuereido book, pg. 70
    α = -1/(b^2)
    d = interval(1/b - α*b, 1/a - α*a)
    ζ = mid(d)
    δ = radius(d) * (-1..1)

    if a < 0
        ζ = -ζ
    end

    return affine_approx(x, α, ζ, δ)
end

# Directly applies Lemma 1 of Rump 2015 to convex and concave functions of
# fixed monotonicity
f_inverse = Dict{Symbol,Any}()

f_inverse[:exp]   = :log
f_inverse[:exp2]  = :log2
f_inverse[:exp10] = :log10
f_inverse[:expm1] = :log1p

f_inverse[:log]   = :log
f_inverse[:log2]  = :exp2
f_inverse[:log10] = :exp10
f_inverse[:log1p] = :expm1

#f_inverse[:sqrt] = # TODO
f_inverse[:acosh] = :cosh

# Define montonicity
f_curv = Dict{Symbol,Any}()

f_curv[:exp]   = :convex
f_curv[:exp2]  = :convex
f_curv[:exp10] = :convex
f_curv[:expm1] = :convex

f_curv[:log]   = :concave
f_curv[:log2]  = :concave
f_curv[:log10] = :concave
f_curv[:log1p] = :concave

#f_curv[:sqrt] = :concave
f_curv[:acosh] = :concave

for f in (:exp, :exp2, :exp10, :expm1, :log, :log2, :log10, :log1p, :acosh)
    pval = f_curv[f] == :convex ? diffrule(:Base, f, :a) : diffrule(:Base, f, :b)
    finv = f_inverse[f]
    @eval function ($f)(x::BasicAffine1{N,T,MinRange}, X=interval(x)) where {N,T<:Real,S<:AbstractLinearization}

        a, b = X.lo, X.hi

        p = $pval
        q = (($f)(a) + ($f)(b) - p*(a + b))/2
        Δ = abs(($f)(b) - ($f)(a) - p*(b - a))/2

        return affine_approx(x, α, ζ, δ)
    end
    @eval function ($f)(x::BasicAffine1{N,T,Chebyshev}, X=interval(x)) where {N,T<:Real,S<:AbstractLinearization}

        a, b = X.lo, X.hi

        p = (($f)(b) - ($f)(a))/(b - a)
        ξ = ($finv)(p)
        q = (($f)(a) + ($f)(ξ) - p*(a + ξ))/2
        Δ = abs(($f)(ξ) - ($f)(a) - p*(ξ - a))/2

        return affine_approx(x, α, ζ, δ)
    end
end
