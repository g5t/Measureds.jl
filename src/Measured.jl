"""
Values which have been determined experimentally must have some uncertainty. 
The `Measured` immutable type is intended to carry measured values and their associated variances
through mathematical operations such that the uncertainty in any value derived from experimental
measurements is known.

In general, a function of any number of variables, f(xᵢ), has an uncertainty with variance
    δf²=∑ᵢ|∂f/∂xᵢ|² δx²ᵢ
where δx²ᵢ is the variance of xᵢ. It is therefore possible to calculate the variance for any 
mathematical operation that has well-defined partial derivatives. Simple mathematical operations
have already been extended to correctly include variance calculation. Algorithms that operate on 
`Real` values via simple mathematical operations will operate on `Measured` values as well, with
correct propogation of uncertainty. More complicated algorithms may need to be extended to support
use of `Measured` values with correct error propagation.

`Measured` objects have two fields, `val` and `var`, the measured value and variance, respectively.
The value, variance, and square-root of the variance (the error) of a `Measured` object `m` can be
accessed through the functions `val(m)`, `var(m)`, and `err(m)`, respectively.
"""
immutable Measured <: Real
    val::Float64
    var::Float64
end
Measured(a::Real=0)=Measured(a,0.)
@vectorize_2arg Real Measured
@vectorize_1arg Real Measured
Measured(m::Measured)=m # no need to worry about copying because Measured values are immutable
@vectorize_1arg Measured Measured

countingerror(a::Measured)=countingerror(a.val)
countingerror{T<:AbstractFloat}(a::T)=a==0?Measured(a,one(T)):Measured(a,a)
countingerror(a::AbstractArray)=map(countingerror,a)
countingerror!{T<:Measured}(a::AbstractArray{T})= (for i=1:length(a); a[i]=countingerror(a[i]); end)


zero(::Type{Measured})=Measured(0)
one(::Type{Measured})=Measured(1)

# comparisons, only including variance when comparing two Measured values
<(a::Measured,b::Measured)=isless(a,b)
>(a::Measured,b::Measured)=isless(b,a)
<(a::Measured,b::Real)=isless(a,b)
<(a::Real,b::Measured)=isless(a,b)
>(a::Measured,b::Real)=isless(b,a)
>(a::Real,b::Measured)=isless(b,a)
<=(a::Measured,b::Measured)=isless(a,b)|isapprox(a,b)
>=(a::Measured,b::Measured)=isless(b,a)|isapprox(a,b)
<=(a::Measured,b::Real)=a.val <= b
<=(a::Real,b::Measured)=a <= b.val
>=(a::Measured,b::Real)=a.val >= b
>=(a::Real,b::Measured)=a >= b.val
isless(a::Measured,b::Measured)=isless(a.val+sqrt(a.var),b.val-sqrt(b.var))
isapprox(a::Measured,b::Measured)=!(isless(a,b)|isless(b,a))
isless{T<:AbstractFloat}(a::Measured,b::T)=isless(a.val,b)
isless{T<:AbstractFloat}(a::T,b::Measured)=isless(a,b.val)
isless{T<:Real}(a::Measured,b::T)=isless(a.val,b)
isless{T<:Real}(a::T,b::Measured)=isless(a,b.val)
isapprox{T<:Real}(a::Measured,b::T)=isapprox(a.val,b)
isapprox{T<:Real}(a::T,b::Measured)=isapprox(a,b.val)

Base.typemax(::Type{Measured})=Base.typemax(Float64)
Base.typemin(::Type{Measured})=Base.typemin(Float64)


# functions with one input value
fdf1=[:cos  (x,vx)->Base.abs(Base.sin(x)).^2.*vx
      :sin  (x,vx)->Base.abs(Base.cos(x)).^2.*vx
      :tan  (x,vx)->Base.abs(Base.sec(x).^2).^2.*vx
      :sec  (x,vx)->Base.abs(Base.sec(x).*Base.tan(x)).^2.*vx
      :csc  (x,vx)->Base.abs(Base.cot(x).*Base.csc(x)).^2.*vx
      :cot  (x,vx)->Base.abs(Base.csc(x).^2).^2.*vx
      :sinc (x,vx)->Base.abs(Base.cos(x)./x-Base.sin(x)./x.^2).^2.*vx
      :acos (x,vx)->vx./(1-x.^2)
      :asin (x,vx)->vx./(1-x.^2)
      :atan (x,vx)->vx./(1+x.^2).^2
      :cosh (x,vx)->Base.abs(Base.sinh(x)).^2.*vx
      :sinh (x,vx)->Base.abs(Base.cosh(x)).^2.*vx
      :tanh (x,vx)->Base.abs(Base.sech(x).^2).^2.*vx
      :sech (x,vx)->Base.abs(Base.sech(x).*Base.tanh(x)).^2.*vx
      :csch (x,vx)->Base.abs(Base.coth(x).*Base.csch(x)).^2.*vx
      :coth (x,vx)->Base.abs(Base.csch(x).^2).^2.*vx
      :sqrt (x,vx)->vx./x/2
      :exp  (x,vx)->vx.*Base.exp(2*x)
      :log  (x,vx)->vx.*Base.abs(1./x).^2
      :-    (x,vx)->vx
      :abs  (x,vx)->vx
     ]
for (f,df) in zip(fdf1[:,1],fdf1[:,2])
    @eval Base.$f(a::Measured)=Measured(Base.$f(a.val),$df(a.val,a.var))
    eval(quote 
        function Base.$f{T<:Measured}(a::Array{T,1})
            out=Array(T,size(a))
            for i=1:length(a); out[i]=Measured(Base.$f(a[i].val),$df(a[i].val,a[i].var)); end
            return out
        end
    end)
    eval(quote 
        function Base.$f{T<:Measured}(a::Array{T,2})
            out=Array(T,size(a))
            for j=1:size(a,2), i=1:size(a,1)
                out[i,j]=Measured(Base.$f(a[i,j].val),$df(a[i,j].val,a[i,j].var)); 
            end
            return out
        end
    end)
end
# functions with two input values
fdf2=[:+    (x,vx,y,vy)->vx+vy
      :.+   (x,vx,y,vy)->vx+vy
      :-    (x,vx,y,vy)->vx+vy
      :.-   (x,vx,y,vy)->vx+vy
      :*    (x,vx,y,vy)->Base.abs(y).^2.*vx+abs(x).^2.*vy
      :.*   (x,vx,y,vy)->Base.abs(y).^2.*vx+abs(x).^2.*vy
      :/    (x,vx,y,vy)->Base.abs(1./y).^2.*vx + abs(x./y.^2).^2.*vy
      :./   (x,vx,y,vy)->Base.abs(1./y).^2.*vx + abs(x./y.^2).^2.*vy
      :atan2 (x,vx,y,vy)->y.^2*(vy.*(x./y).^2+vx)./(x.^2+y.^2).^2
      ]
for (f,df) in zip(fdf2[:,1],fdf2[:,2])
   @eval Base.$f(a::Measured,b::Measured)=Measured(Base.$f(a.val,b.val),$df(a.val,a.var,b.val,b.var))
end
@vectorize_2arg Measured Base.atan2

# Special functions which are remapped to corresponding functions
Base.log10(x::Measured)=Base.log(x)./Base.log(10.)
@vectorize_1arg Measured Base.log10

types=[Integer,Rational,Real]
for a in types
    @eval ^(x::Measured,n::$a)=Measured(x.val.^n, abs(n*x.val.^(n-1)).^2.*x.var)
end

Base.isnan(x::Measured)=Base.isnan(x.val)
Base.isfinite(x::Measured)=Base.isfinite(x.val)
Base.signbit(a::Measured)=Base.signbit(a.val)

"""
    indexin(a::Array{Measured},b::Array{Measured})
    
Returns the highest (linear) index of `b` for each value of `a` that is within uncertainty of
one or more mbembers of `b`. The output vector contains 0 wherever `a` is not within
uncertainty of a member of `b`.
"""
function Base.indexin{T<:Measured}(a::AbstractArray{T},b::AbstractArray{T})
    ind=Array{Int}(size(a)...)
    vb=vec(b)
    for i=1:length(a)
        isvec=map(x->isapprox(a[i],x),vb)
        ind[i]=any(isvec)?findlast(isvec):0
    end
    return ind
end
    

function discernable{T<:Measured}(x::AbstractArray{T,1})
    xval,xerr=val(x),err(x)
    # BitArrays are **not** initialized to zero! Use "falses" instead
    found=falses(size(x)...)
    du=falses(size(x)...)
    while sum(found)<length(x)
        i=findfirst(!found)
        found |= abs(xval-xval[i]).<=xerr[i] # points within error have been found (including the current point)
        du[i]=true  # and the this point we choose to be the unique one
    end
    return du
end
discernablyUnique{T<:Measured}(x::AbstractArray{T,1})=val(x)[discernable(x)] # kept for legacy use
Base.unique{T<:Measured}(x::AbstractArray{T,1})=x[discernable(x)]
@deprecate discernablyUnique(x) val(unique(x))

## this is dumb: (squaring a complex number causes promotion of the number 
#Base.isinteger(::Measured)=true

Base.size(x::Measured)=tuple() # like any single value, the size of a Measured is ()

# promote_array_type used by .+ and .-
#Base.promote_array_type{T<:Real}(F,::Type{T},::Type{Measured})=Measured
#Base.promote_array_type(F,::Type{Measured},::Type{Measured})=Measured


# any operation on a Measured and any other number should yield a Measured
Base.promote_op(f,::Type{Measured},::Type{Measured})=Measured
Base.promote_op(f,::Type{Measured},::Type{Union{}})=Measured
Base.promote_op(f,::Type{Union{}},::Type{Measured})=Measured
# even in Arrays
Base.promote_array_type{T<:AbstractFloat}(F,::Type{Measured},::Type{T})=Measured
Base.promote_array_type{T<:Real}(F,::Type{Measured},::Type{T})=Measured
Base.promote_array_type{S<:Measured}(F, ::Type{S}, ::Type, ::Type) = S
Base.promote_array_type{S<:Measured, A<:AbstractFloat}(F, ::Type{S}, ::Type{A}, ::Type) = S

# promote_type seems to be used by .* and ./ broadcasts
Base.promote_type(::Type{Measured},::Type{Measured})= Measured
Base.promote_type(::Type{Measured},::Type{Union{}})= Measured
Base.promote_type(::Type{Union{}},::Type{Measured})= Measured
#XXX#Base.promote_type{T<:Real}(::Type{Measured},::Type{T})= Measured
#XXX#Base.promote_type{T<:Real}(::Type{T},::Type{Measured})= Measured

# promote seems to be used by + and -
Base.promote(a::Measured,b::Measured)=(a,b)
Base.promote(a::Real,b::Measured)=(Measured(a),b)
Base.promote(a::Measured,b::Real)=(a,Measured(b))

# promote_rule appears to be used by .^
Base.promote_rule(::Type{Measured},::Type{Measured})= Measured
Base.promote_rule{R<:Real}(::Type{Measured},::Type{R})= Measured
Base.promote_rule(::Type{Bool},::Type{Measured})=Measured
Base.promote_rule(::Type{Base.MPFR.BigFloat},::Type{Measured})=Measured
Base.promote_rule(::Type{Base.Irrational},::Type{Measured})=Measured
Base.promote_rule{R<:Real}(::Type{R},::Type{Measured})= Measured
#XXX#Base.promote_rule{M<:Measured}(::Type{M},::Type{M})= M
#XXX#Base.promote_rule{R<:Real,M<:Measured}(::Type{M},::Type{R})= M
#XXX#Base.promote_rule{R<:Real,M<:Measured}(::Type{R},::Type{M})= M

#XXX#Base.convert(::Type{Bool},x::Measured)=convert(Bool,x.val)
#XXX#Base.convert(::Type{Integer},x::Measured)=convert(Integer,x.val)
#XXX#Base.convert(::Type{Rational{Base.GMP.BigInt}},x::Measured)=convert(Rational{GMP.BigInt},x.val)
#XXX#Base.convert{T<:Integer}(::Type{Rational{T}},x::Measured)=convert(Rational{T},x.val)
Base.convert(::Type{Measured},x::Measured)=x
#Base.convert{S<:Real}(::Type{S},x::Measured)=convert(S,x.val)
Base.convert(::Type{Measured},x::Real)=Measured(convert(Float64,x),0.) # allows for, e.g., Measured[1,2,3,4]
Base.convert{T<:AbstractFloat}(::Type{T},x::Measured)=convert(T,x.val) # used in creating TripleAxis object from measured angles # TODO think about moving TripleAxis to using Measured values

val(m::Measured)=m.val
var(m::Measured)=m.var
err(m::Measured)=sqrt(m.var)
@vectorize_1arg Measured val
@vectorize_1arg Measured var
@vectorize_1arg Measured err
#errforplot{T<:Measured}(m::Array{T,1})=hcat(err(m),err(m))'

pre(x)=(x<0?-1:1)*floor(Integer,abs(x))
zeropad(a::Real,b::Int)= ( n=a>0?round(Int,floor(log10(a)))+1:1; b-=n; (b>0?"$("0"^b)":"")*"$a" )
innerpost(x,y)=round(Integer,abs(x-pre(x))/10.0^y)
post(x,y)=zeropad(innerpost(x,y),abs(y))
function showMeasured(io::IO,m::Measured,compact::Bool=true)
    v,u = val(m),err(m)
    if u>0 # some uncertainty
        pv=(v!=0)?round(Integer,floor(log10(abs(v)))):0 # power of first digit in the value
        pu=round(Integer,floor(log10(u))) # power of first digit in the error
        dv=round(Integer,v/10.0^pu) # the value rounded to the precision of 1 digit in the error
        if dv==0 && v!=0 # the value isn't zero, but the rounded value is
            while dv==0; pu=pu-1; dv=round(Integer,v/10.0^pu); end # rounded value is zero but actual value is non-zero, increase the precision
        end
        du=round(Integer,u/10.0^pu) # the error rounded to the same precision
        
        rv=dv*10.0^pu # the error-precision-rounded value (now a Float, so printing is complicated)
        # now actually do the formatting
        if (pv<0&&pu<0)||(v==0&&pu<0) # we need to try it both ways
            str1="$dv($du)ᴇ$pu"
            str2=((v<0)?"-":"")*"$(pre(rv)).$(post(rv,pu))($du)"
            str=(length(str1)<(length(str2)-3))?str1:str2 # pick str2 preferentially
        elseif pu>0 # two ways to try as well
            tdv=round(v/10.0^pv)
            tv=(v-tdv*10.0^pv)/10.0.^pv
            str1=((v<0)?"-":"")*"$(pre(tdv))."*"$(post(tv,pv-pu))($du)ᴇ$pv"
            str2="$(dv*10^pu)($(du*10^pu))"
            str=(length(str1)<(length(str2)-3))?str1:str2 # pick str2 preferentially
        elseif pu<0
            str="$dv"[1:end+pu]*"."*"$dv"[1+end+pu:end]*"($du)"
        else # pu==0
            str="$dv($du)"
        end
        print(io,str)
    else
        compact ? showcompact(io,v) : show(io,v)
    end
end    
Base.show(io::IO,m::Measured)=showMeasured(io,m,false)
Base.showcompact(io::IO,m::Measured)=showMeasured(io,m,true)
function Base.alignment(x::Measured)
    m = match(r"^(.*?)((?:[(±].*)?)$", sprint(Base.showcompact_lim, x))
    m === nothing ? (length(sprint(Base.showcompact_lim, x)), 0) : (length(m.captures[1]), length(m.captures[2]))
end

approx(a::Measured)=Measured(round(val(a),-floor(Int,log10(err(a)))),var(a))

