
# Twice Differentiable, Increasing functions of fix monotonicity (Rump 2015)
f_inverse = Dict{Symbol,Expr}()

f_inverse[:exp]   = :log
f_inverse[:exp2]  = :log2
f_inverse[:exp10] = :log10
f_inverse[:expm1] = :log1p

f_inverse[:log]   = :log
f_inverse[:log2]  = :exp2
f_inverse[:log10] = :exp10
f_inverse[:log1p] = :expm1

f_inverse[:sqrt] = # TODO
f_inverse[:acosh] = :cosh

# Define montonicity
f_curv = Dict{Symbol,Expr}()

f_curv[:exp]   = :convex
f_curv[:exp2]  = :convex
f_curv[:exp10] = :convex
f_curv[:expm1] = :convex

f_curv[:log]   = :concave
f_curv[:log2]  = :concave
f_curv[:log10] = :concave
f_curv[:log1p] = :concave

f_curv[:sqrt] = :concave
f_curv[:acosh] = :concave

for f in (:exp, :exp2, :exp10, :expm1, :log, :log2, :log10, :log1, :sqrt, :acosh)
    pval = f_curv[f] == :convex ? diffrule(:Base, f, :a) : diffrule(:Base, f, :b)
    finv = f_inverse[f]
    @eval function ($f)(x::Aff{N,T,MinRange}, X=interval(x)) where {N,T<:Real}

        a, b = X.lo, X.hi

        p = $pval
        q = (($f)(a) + ($f)(b) - p*(a + b))/2
        Δ = abs(($f)(b) - ($f)(a) - p*(b - a))/2

        return affine_approx(x, α, ζ, δ)
    end
    @eval function ($f)(x::Aff{N,T,Chebyshev}, X=interval(x)) where {N,T<:Real}

        a, b = X.lo, X.hi

        p = (($f)(b) - ($f)(a))/(b - a)
        ξ = ($finv)(p)
        q = (($f)(a) + ($f)(ξ) - p*(a + ξ))/2
        Δ = abs(($f)(ξ) - ($f)(a) - p*(ξ - a))/2

        return affine_approx(x, α, ζ, δ)
    end
end
