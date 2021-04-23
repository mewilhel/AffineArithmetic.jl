module AffineArithmetic

using IntervalArithmetic, StaticArrays, LinearAlgebra, DiffRules
using DiffRules: diffrule

import IntervalArithmetic: interval
import Base: +, -, *, /, ^, ==, exp, exp2, exp10, expm1, log, log2,
             log10, log1p, acosh, zero, one, range, show

abstract type AbstractLinearization end
struct MinRange <: AbstractLinearization end
struct Chebyshev <: AbstractLinearization end
const DEFAULT_LIN = MinRange

include("basic_affine1.jl")
#include("basic_affine2.jl")
#include("mixed_affine.jl")

const Affine = BasicAffine1

export Affine, BasicAffine1, MinRange, Chebyshev
export +, -, *, /, ^, ==, exp, exp2, exp10, expm1, log, log2,
       log10, log1p, acosh, zero, one, range, show

end
