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
correct propogation of uncertainty. More complicated algorithms support use of `Measured` values,
however they must be extended for correct error propagation.

`Measured` objects have two fields, `val` and `var`, the measured value and variance, respectively.
The value, variance, and square-root of the variance (the error, or uncertainty) of a `Measured` object `m` can be
accessed through the functions `value(m)`, `variance(m)`, and `uncertainty(m)`, respectively.
"""
abstract type Measured <: Number end #<: Real
# TODO make Measured{As,S}ymmetric{T<:AbstractFloat} <: Measured{T}; val::T ...
# with logic to ensure the internal type can represent val to the accuracy of sqrt(var,vap,vam)
# then, values like 10(2) could easily be represented by Float16 thereby saving memory
#
# This will necessitate use of eps(T) for T in [Float16,Float32,Float64,BigFloat]
# And ultimately may be more trouble than it's worth
immutable MeasuredAsymmetric <: Measured
    val::Float64
    vap::Float64
    vam::Float64
    MeasuredAsymmetric(a::Float64,b::Float64,c::Float64)=new(a,abs(b),abs(c))
end
immutable MeasuredSymmetric  <: Measured
    val::Float64
    var::Float64
    MeasuredSymmetric(a::Float64,b::Float64)=new(a,abs(b))
end
const MS = MeasuredSymmetric
const MA = MeasuredAsymmetric

MeasuredAsymmetric(a::MA)=a
MeasuredAsymmetric(a::MS)=MA(a.val,a.var,a.var)
MeasuredAsymmetric(a::Real,b=zero(Float64),c=zero(Float64))=MA(Float64(a),Float64(b),Float64(c))

MeasuredSymmetric(a::MS)=a
MeasuredSymmetric(::MA)=throw(InexactError())
MeasuredSymmetric(a::Real,b=zero(Float64))=MS(Float64(a),Float64(b))
# the above was {T<:Real}(a::T,b=zero(T)) which fails for T<:Irrational since convert(Irrational,0) fails


countingerror(a::Measured)=countingerror(a.val)
countingerror{T<:AbstractFloat}(a::T)=a==0?MA(a,one(T),zero(T)):MS(a,a)
countingerror(a::AbstractArray)=map(countingerror,a)
countingerror!{T<:Measured}(a::AbstractArray{T})= (for i=1:length(a); a[i]=countingerror(a[i]); end)

symmetricerror(a::Measured,r)=symmetricerror(a.val,r)
symmetricerror{T<:AbstractFloat}(a::T,r)=MS(a,r^2)
symmetricerror(a::AbstractArray,r)=map(x->symmetricerror(x,r),a)
symmetricerror!{T<:Measured}(a::AbstractArray{T},r)=(for i=1:length(a); a[i]=symmetricerror(a[i],r); end)

asymmetricerror(a::Measured,r,p)=asymmetricerror(a.val,r,p)
asymmetricerror{T<:AbstractFloat}(a::T,r,p)=MA(a,r^2,p^2)
asymmetricerror(a::AbstractArray,r,p)=map(x->asymmetricerror(x,r,p),a)
asymmetricerror!{T<:Measured}(a::AbstractArray{T},r,p)=(for i=1:length(a); a[i]=asymmetricerror(a[i],r,p); end)

replaceerror!(a,b)=symmetricerror!(a,b)
replaceerror!(a,b,c)=asymmetricerror!(a,b,c)


# Generic constructors take (value, uncertainty)
#Measured(a::Real,b::Real=0)      =MS(a,b.^2) # symmetric by default
#Measured(a::Real,b::Real,c::Real)=MA(a,b.^2,c.^2)
Measured(a::Real)                =MS(a)    # symmetric by default
Measured(a::Real,b::Real)        =MS(a,b)  # two values -> symmetric
Measured(a::Real,b::Real,c::Real)=MA(a,b,c)# three values -> asymmetric
# Generic array constructors take (value, variance)
# three arrays -> asymmetric
Measured{T<:Real,R<:Real,L<:Real}(a::Array{T},b::Array{R},c::Array{L})=map(MA,a,b,c)
# two arrays -> symmetric; +1 constant -> asymmetric
Measured{T<:Real,R<:Real}(a::Array{T},b::Array{R}        )=map(MS,a,b)
Measured{T<:Real,R<:Real}(a::Array{T},b::Array{R},c::Real)=map((x,y)->MA(x,y,c),a,b)
# one array +1 constant -> symmetric; +2 constants -> asymmetric
Measured{T<:Real}(a::Array{T},b::Real,c::Real)=map(x->MA(x,b,c),a)
Measured{T<:Real}(a::Array{T},b::Real        )=map(x->MS(x,b),a)
# generalized "copy" constructor
Measured(m::Measured)=m # no need to worry about copying because Measured values are immutable
# Specific array constructors take (value, variance)
# asymmetric
MA{T<:Real,R<:Real,L<:Real}(a::Array{T},b::Array{R},c::Array{L})=map(MA,a,b,c)
MA{T<:Real,R<:Real}(a::Array{T},b::Array{R},c::Real=zero(T))=map((x,y)->MA(x,y,c),a,b)
MA{T<:Real,R<:Real}(a::Array{T},b::Real,c::Array{R})=map((x,y)->MA(x,b,y),a,c)
MA{T<:Real}(a::Array{T},b::Real=zero(T),c::Real=zero(T))=map(x->MA(x,b,c),a)
# symmetric
MS{T<:Real,R<:Real}(a::Array{T},b::Array{R})=map(MS,a,b)
MS{T<:Real}(a::Array{T},b::Real=zero(T))=map(x->MS(x,b),a)
MS{T<:Real}(a::Real,b::Array{T})=map(x->MS(a,x),b)

Base.zero(::Type{MS})=MS(0,0)
Base.one(::Type{MS})=MS(1,0)
Base.zero(::Type{MA})=MA(0,0,0)
Base.one(::Type{MA})=MA(1,0,0)


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
isless(a::MS,b::MS)=isless(a.val+sqrt(a.var),b.val-sqrt(b.var))
isless(a::MA,b::MA)=isless(a.val+sqrt(a.vap),b.val-sqrt(b.vam))
isless(a::MS,b::MA)=isless(a.val+sqrt(a.var),b.val-sqrt(b.vam))
isless(a::MA,b::MS)=isless(a.val+sqrt(a.vap),b.val-sqrt(b.var))
isless(::Measured,::Measured)=throw(InexactError()) # This shouldn't be possible as the four combinations are defined
isapprox(a::Measured,b::Measured)=!(isless(a,b)|isless(b,a))
isless{T<:AbstractFloat}(a::Measured,b::T)=isless(a.val,b)
isless{T<:AbstractFloat}(a::T,b::Measured)=isless(a,b.val)
isless{T<:Real}(a::Measured,b::T)=isless(a.val,b)
isless{T<:Real}(a::T,b::Measured)=isless(a,b.val)
isapprox{T<:Real}(a::Measured,b::T)=isapprox(a.val,b)
isapprox{T<:Real}(a::T,b::Measured)=isapprox(a,b.val)

function fuzzyindexin{T<:Measured,R<:Measured}(a::AbstractArray{T},b::AbstractArray{R})
    bidx=1:length(b)
    aidx=Array{Int}(size(a)...)
    for i=1:length(a); j=findfirst(x->isapprox(x,a[i]),b); aidx[i]=j>0?bidx[j]:0; end
    return aidx
end
fuzzyindexin(a::AbstractArray, b::AbstractArray)=Base.indexin(a,b) # for regular types just use indexin

# unitary/one-argument functions
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
      :deg2rad (x,vx)-> vx/180/180*pi*pi
      :rad2deg (x,vx)-> vx/pi/pi*180*180
     ]
for (f,df) in zip(fdf1[:,1],fdf1[:,2])
   @eval Base.$f(a::MS)=MS(Base.$f(a.val),$df(a.val,a.var))
   @eval Base.$f(a::MA)=MA(Base.$f(a.val),$df(a.val,a.vap),$df(a.val,a.vam))
end
fdf1_inverting=[:-    (x,vx)->vx] # inverting function(s) flip +/- to -/+ variances
for (f,df) in zip(fdf1_inverting[:,1],fdf1_inverting[:,2])
   @eval Base.$f(a::MS)=MS(Base.$f(a.val),$df(a.val,a.var))
   @eval Base.$f(a::MA)=MA(Base.$f(a.val),$df(a.val,a.vam),$df(a.val,a.vap))
end
fdf1_abs=[:abs   (x,vx)->vx
          :abs2  (x,vx)->4x^2*vx
         ]
for (f,df) in zip(fdf1_abs[:,1],fdf1_abs[:,2]) # abs functions flip +- to -= variance if the input value is negative
  @eval Base.$f(a::MS)=MS(Base.$f(a.val),$df(a.val,a.var))
  @eval Base.$f(a::MA)=MA(Base.$f(a.val),$df(a.val,Base.signbit(a.val)?a.vam:a.vap),$df(a.val,Base.signbit(a.val)?a.vap:a.vam))
end
# functions with two input values
fdf2=[:+     (x,vx,y,vy)->vx+vy
      :-     (x,vx,y,vy)->vx+vy
      :*     (x,vx,y,vy)->Base.abs(y).^2.*vx+Base.abs(x).^2.*vy
      :/     (x,vx,y,vy)->Base.abs(1./y).^2.*vx + Base.abs(x./y.^2).^2.*vy
      :atan2 (x,vx,y,vy)->y.^2*(vy.*(x./y).^2+vx)./(x.^2+y.^2).^2
      :^     (x,vx,y,vy)->Base.abs(y*(x)^(y-1)).^2*vx+Base.abs(x^y*log(x)).^2*vy
      ]
for (f,df) in zip(fdf2[:,1],fdf2[:,2])
   @eval Base.$f(a::MS,b::MS)=MS(Base.$f(a.val,b.val),$df(a.val,a.var,b.val,b.var))
   @eval Base.$f(a::MA,b::MA)=MA(Base.$f(a.val,b.val),$df(a.val,a.vap,b.val,b.vap),$df(a.val,a.vam,b.val,b.vam))
   @eval Base.$f(a::MS,b::MA)=MA(Base.$f(a.val,b.val),$df(a.val,a.var,b.val,b.vap),$df(a.val,a.var,b.val,b.vam))
   @eval Base.$f(a::MA,b::MS)=MA(Base.$f(a.val,b.val),$df(a.val,a.vap,b.val,b.var),$df(a.val,a.vam,b.val,b.var))
end
fdf3=[:.^    (x,vx,y,vy)->Base.abs(y.*(x).^(y-1)).^2.*vx+Base.abs(x.^y.*log(x)).^2.*vy]
for (f,df) in zip(fdf3[:,1],fdf3[:,2])
   @eval Base.broadcast(::typeof(Base.$f),a::MS,b::MS)=MS(Base.$f(a.val,b.val),$df(a.val,a.var,b.val,b.var))
   @eval Base.broadcast(::typeof(Base.$f),a::MA,b::MA)=MA(Base.$f(a.val,b.val),$df(a.val,a.vap,b.val,b.vap),$df(a.val,a.vam,b.val,b.vam))
   @eval Base.broadcast(::typeof(Base.$f),a::MS,b::MA)=MA(Base.$f(a.val,b.val),$df(a.val,a.var,b.val,b.vap),$df(a.val,a.var,b.val,b.vam))
   @eval Base.broadcast(::typeof(Base.$f),a::MA,b::MS)=MA(Base.$f(a.val,b.val),$df(a.val,a.vap,b.val,b.var),$df(a.val,a.vam,b.val,b.var))
end

# Special functions which are remapped to corresponding functions
Base.broadcast(::typeof(*),x::Measured,y::Measured)=x*y
Base.broadcast(::typeof(/),x::Measured,y::Measured)=x/y
log10(x::Measured)=log(x)./log(10.)

types=[Integer,Rational,Real]
powvar(x,v,n)=abs(n*x.^(n-1)).^2.*v
for a in types
    @eval ^(x::MS,n::$a)=MS(x.val.^n, powvar(x.val,x.var,n))
    @eval ^(x::MA,n::$a)=MA(x.val.^n, powvar(x.val,x.vap,n), powvar(x.val,x.vam,n))
end




Base.isnan(x::Measured)=Base.isnan(x.val)
Base.isfinite(x::Measured)=Base.isfinite(x.val)
Base.signbit(a::Measured)=Base.signbit(a.val)

value(m::Measured)=m.val
value(m)=m
variance(m::MS)=m.var
variance(m::MA)=(m.vap+m.vam)/2
variance(m)=0*m
positivevariance(m::MS)=m.var
negativevariance(m::MS)=m.var
positivevariance(m::MA)=m.vap
negativevariance(m::MA)=m.vam
positivevariance(m)=0m
negativevariance(m)=0m
uncertainty(m::Measured)=sqrt(variance(m))
uncertainty(m)=0m
negativeuncertainty(m::Measured)=sqrt(negativevariance(m))
positiveuncertainty(m::Measured)=sqrt(positivevariance(m))
negativeuncertainty(m)=0m
positiveuncertainty(m)=0m
errforplot{T<:MS}(m::Vector{T})=uncertainty.(m)
errforplot{T<:Measured}(m::Vector{T})=sqrt.(hcat(negativevariance.(m),positivevariance.(m)))'


# vap(x) and vam(x) of MeasuredSymmetric are defined equivalent to var(x)
function discernable{T<:Measured}(x::AbstractArray{T,1})
    xval,xerp,xerm=value(x),positiveuncertainty(x),-negativeuncertainty(x)
    # BitArrays are **not** initialized to zero! Use "falses" instead
    found=falses(size(x)...)
    du=falses(size(x)...)
    we=falses(size(x)...) # within error logical vector
    while sum(found)<length(x)
        i=findfirst(!found)
        dx=xval-xval[i]
        we[dx.>0]=dx[dx.>0].<=xerp[i]
        we[dx.<0]=dx[dx.<0].>=xerm[i]
        we[dx.==0]=true
        found |= we # points within error have been found
        du[i]=true  # and the first point we choose to be the unique one
    end
    return du
end
discernablyUnique{T<:Measured}(x::AbstractArray{T,1})=value(x)[discernable(x)] # kept for legacy use
Base.unique{T<:Measured}(x::AbstractArray{T,1})=x[discernable(x)]
@deprecate discernablyUnique(x) value(unique(x))

Base.size(::Measured)=tuple() # like any single value, the size of a Measured is ()

# promote_array_type used by .+ and .-
#Base.promote_array_type{T<:Real}(F,::Type{T},::Type{Measured})=Measured
#Base.promote_array_type(F,::Type{Measured},::Type{Measured})=Measured

#Base.promote_array_type{T<:AbstractFloat}(F,::Type{Measured},::Type{T})=Measured
#Base.promote_array_type{T<:Real}(F,::Type{Measured},::Type{T})=Measured


# promote_type seems to be used by .* and ./ broadcasts
Base.promote_type(::Type{MA},::Type{MS})=MA
Base.promote_type(::Type{MS},::Type{MA})=MA
Base.promote_type(::Type{Measured},::Type{MS})=Measured
Base.promote_type(::Type{Measured},::Type{MA})=Measured
Base.promote_type{M<:Measured}(::Type{M},::Type{M})= M
Base.promote_type{M<:Measured}(::Type{M},::Type{Union{}})= M
Base.promote_type{M<:Measured}(::Type{Union{}},::Type{M})= M
#Base.promote_type{T<:Real,M<:Measured}(::Type{M},::Type{T})= M
#Base.promote_type{T<:Real,M<:Measured}(::Type{T},::Type{M})= M

# promote seems to be used by + and -
Base.promote(a::MA,b::MS)=(a,MA(b.val,b.var,b.var))
Base.promote(a::MS,b::MA)=(MA(a.val,a.var,a.var),b)
#Base.promote{M<:Measured}(a::M,b::M)=(a,b)
#Base.promote(a::Measured,b::Measured)=throw(InexactError()) # this also shouldn't be possible
Base.promote(a::MS,b::MS)=(a,b)
Base.promote(a::MA,b::MA)=(a,b)
Base.promote{T<:Measured}(a::Real,b::T)=(T(a),b)
Base.promote{T<:Measured}(a::T,b::Real)=(a,T(b))

# promote_rule appears to be used by .^
Base.promote_rule{M<:Measured}(::Type{M},::Type{M})= M
Base.promote_rule{M<:Measured,T<:Measured}(a::Type{M},b::Type{T})=Base.promote_type(a,b)
Base.promote_rule{T<:Real,M<:Measured}(::Type{M},::Type{T})= M

Base.promote_rule{M<:Measured}(::Type{Bool},::Type{M})=M
Base.promote_rule{M<:Measured}(::Type{Base.MPFR.BigFloat},::Type{M})=M
Base.promote_rule{M<:Measured}(::Type{Base.Irrational},::Type{M})=M

Base.promote_rule{T<:Real,M<:Measured}(::Type{T},::Type{M})= M

#Base.convert(::Type{Bool},x::Measured)=convert(Bool,x.val)
#Base.convert(::Type{Integer},x::Measured)=convert(Integer,x.val)
#Base.convert(::Type{Rational{Base.GMP.BigInt}},x::Measured)=convert(Rational{GMP.BigInt},x.val)
#Base.convert{T<:Integer}(::Type{Rational{T}},x::Measured)=convert(Rational{T},x.val)
Base.convert(::Type{MA},a::MS)=MA(a.val,a.var,a.var)
Base.convert(::Type{MS},::MA)=throw(InexactError())
Base.convert(::Type{MS},a::MS)=a
Base.convert(::Type{MS},x::Real)=MS(x)
Base.convert(::Type{MA},a::MA)=a
Base.convert(::Type{MA},x::Real)=MA(x)
#Base.convert{M<:Measured}(::Type{M},x::M)=x

Base.convert(::Type{Measured},a::MS)=a
Base.convert(::Type{Measured},a::MA)=a
Base.convert(::Type{Measured},x::Real)=MS(x) # allows for, e.g., Measured[1,2,3,4]

Base.convert{T<:AbstractFloat}(::Type{T},x::Measured)=convert(T,x.val)

# TODO replace val, var, vap, vam, err with names that are not defined in Base
# NOTE this is only realy a problem with var
#
# my suggestions:
#       val -> value
#       var -> variance
#       vap -> positivevariance
#       vam -> negativevariance
#       err -> uncertainty
#
# I would NOT change the internal fieldnames, just the accessor function names

#val(m::Measured)=value(m)
#var(m::Measured)=variance(m)
#vap(m::Measured)=positivevariance(m)
#vam(m::Measured)=negativevariance(m)
#err(m::Measured)=uncertainty(m)



sigpower(value)=(value!=0)?round(Integer,floor(log10(abs(value)))):0 # power of first digit in the value
roundtopower(v,pu)=round(Integer,v/10.0^pu)
pre(x)=(x<0?-1:1)*floor(Integer,abs(x))
zeropad(a::Real,b::Int)= ( n=a>0?round(Int,floor(log10(a)))+1:1; b-=n; (b>0?"$("0"^b)":"")*"$a" )
innerpost(x,y)=round(Integer,abs(x-pre(x))/10.0^y)
post(x,y)=zeropad(innerpost(x,y),abs(y))
function showMeasured(io::IO,m::MS,compact::Bool=true)
    v,u = value(m),uncertainty(m)
    if u>0&&~(isnan(v)||isinf(v)) # some uncertainty
        pv=sigpower(v); pu=sigpower(u); dv=roundtopower(v,pu) # round v to the power of the first digit of uncertainty
        0==v||(while 0==(dv=roundtopower(v,pu)); pu-=1; end) # if v!=0, make sure the rounded value is non-zero as well
        du=roundtopower(u,pu); # the error rounded to the same precision as the value

        rv=dv*10.0^pu # the error-precision-rounded value (now a Float, so printing is complicated)
        # now actually do the formatting
        if (pv<0&&pu<0)||(v==0&&pu<0) # we need to try it both ways
            str1="$dv($du)e$pu"
            str2=((v<0)?"-":"")*"$(pre(rv)).$(post(rv,pu))($du)"
            str=(length(str1)<(length(str2)-3))?str1:str2 # pick str2 preferentially
        elseif pu>0 # two ways to try as well
            tdv=round(v/10.0^pv)
            tv=(v-tdv*10.0^pv)/10.0.^pv
            #str1=((v<0)?"-":"")*"$(pre(tdv))."*"$(post(tv,pv-pu))($du)e$pv" # produces, e.g., 1.000(5)e8 when 1.234(5)e8 is expected
            str1=((v<0)?"-":"")*"$(pre(tdv))."*"$(post(tv,pu-pv))($du)e$pv" # switched pv-pu to pu-pv. does this cause another problem?
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
function showMeasured(io::IO,a::MA,compact::Bool=true)
    v,p,m = value(a),positiveuncertainty(a),negativeuncertainty(a)
    if p>0||m>0&&~(isnan(v)||isinf(v)) # some uncertainty
        pv=sigpower(v)
        pu=min(sigpower(p),sigpower(m))
        dv=roundtopower(v,pu) # round v to the power of the first digit of uncertainty
        0==v||(while 0==(dv=roundtopower(v,pu)); pu-=1; end) # if v!=0, make sure the rounded value is non-zero as well
        dp=roundtopower(p,pu); dm=roundtopower(m,pu)# the positive and negative error rounded to the same precision as the value

        rv=dv*10.0^pu # the error-precision-rounded value (now a Float, so printing is complicated)
        # now actually do the formatting
        if (pv<0&&pu<0)||(v==0&&pu<0) # we need to try it both ways
            str1="$dv($dp,$dm)e$pu"
            str2=((v<0)?"-":"")*"$(pre(rv)).$(post(rv,pu))($dp,$dm)"
            str=(length(str1)<(length(str2)-3))?str1:str2 # pick str2 preferentially
        elseif pu>0 # two ways to try as well
            tdv=round(v/10.0^pv)
            tv=(v-tdv*10.0^pv)/10.0.^pv
            #str1=((v<0)?"-":"")*"$(pre(tdv))."*"$(post(tv,pv-pu))($dp,$dm)e$pv" # produces, e.g., 1.000(5,6)e8 when 1.234(5,6)e8 is expected
            str1=((v<0)?"-":"")*"$(pre(tdv))."*"$(post(tv,pu-pv))($dp,$dm)e$pv" # switched pv-pu to pu-pv. does this cause another problem?
            str2="$(dv*10^pu)($(dp*10^pu),$(dm*10^pu))"
            str=(length(str1)<(length(str2)-3))?str1:str2 # pick str2 preferentially
        elseif pu<0
            str="$dv"[1:end+pu]*"."*"$dv"[1+end+pu:end]*"($dp,$dm)"
        else # pu==0
            str="$dv($dp,$dm)"
        end
        print(io,str)
    else
        compact ? showcompact(io,v) : show(io,v)
    end
end


Base.show(io::IO,m::Measured)=showMeasured(io,m,false)
Base.showcompact(io::IO,m::Measured)=showMeasured(io,m,true)
function Base.alignment(x::Measured)
    s = sprint(Base.showcompact_lim, x)
    m = match(r"^(.*?)((?:[\.(±].*)?)$", s)
    m === nothing ? (length(s), 0) : (length(m.captures[1]), length(m.captures[2]))
end


# extend Base.conv:
function Base.conv{T<:Measured}(u::AbstractVector{T},v::AbstractVector{T})
  u0=value.(u); up=positiveuncertainty.(u); um=negativeuncertainty.(u)
  v0=value.(v); vp=positiveuncertainty.(v); vm=negativeuncertainty.(v)
  cmm=conv(u0-um,v0-vm); cm0=conv(u0-um,v0   ); cmp=conv(u0-um,v0+vp);
  c0m=conv(u0   ,v0-vm); c00=conv(u0   ,v0   ); c0p=conv(u0   ,v0+vp);
  cpm=conv(u0+up,v0-vm); cp0=conv(u0+up,v0   ); cpp=conv(u0+up,v0+vp);
  mat=hcat(cmm,cm0,cmp,c0m,c00,c0p,cpm,cp0,cpp)
  cp=maximum(mat,2).-c00
  cm=c00.-minimum(mat,2)
  all(cp.==cm)?MeasuredSymmetric(c00,cp.^2):MeasuredAsymmetric(c00,cp.^2,cm.^2)
end
Base.conv{T<:Measured}(u::AbstractVector{T},v::AbstractVector)=Base.conv(u,Measured(v))
Base.conv{T<:Measured}(u::AbstractVector,v::AbstractVector{T})=Base.conv(Measured(u),v)

# extend Base.issymmetric to enable efficient tests of generic Measured arrays
Base.issymmetric(a::MS)=true
Base.issymmetric(a::MA)=false
# extend Base.isinf to act on the stored value
Base.isinf(a::Measured)=Base.isinf(a.val)

allSymmetric{T<:Measured}(a::Array{T})=all(issymmetric,a)
allAsymmetric{T<:Measured}(a::Array{T})=all(x->~issymmetric(x),a)
export allSymmetric,allAsymmetric

all_sym_or_asym{T<:Measured}(a::Array{T})=map(all(issymmetric,a)?MS:MA,a)
function all_sym_or_asym{T<:Measured}(a::Array{T},b::Array{T})
  if all(issymmetric,a)&&all(issymmetric,b)
    ia=map(MS,a); ib=map(MS,b)
  else
    ia=map(MA,a); ib=map(MA,b)
  end
  return (ia,ib)
end
function all_sym_or_asym{T<:Measured}(a::Array{T},b::Array{T}...)
    # we need to unpack b from a Tuple to an array:
    b=[b...]
    println(typeof(a),": ",size(a))
    println(typeof(b),": ",size(b))
    for i=1:length(b)
        println("  ",typeof(b[i]),": ",size(b[i]))
    end
  if all(issymmetric,a)&&all(x->all(issymmetric,x),b)
    ia=map(MS,a); ib=Array{Array{MS},1}(length(b));
    for i=1:length(b); ib[i]=map(MS,b[i]); end
  else
    ia=map(MA,a); ib=Array{Array{MA},1}(length(b));
    for i=1:length(b); ib[i]=map(MA,b[i]); end
  end
  return (ia,ib...)
end

function meanextrema{T<:Measured}(v::Vector{T})
    m=mean(v)
    (n,x)=extrema(value.(v))
    errp=maximum([positiveuncertainty(m),x-value(m)])
    errm=maximum([negativeuncertainty(m),value(m)-n])
    errp==errm&&typeof(m)<:MS?MS(value(m),errp^2):MA(value(m),errp^2,errm^2)
end
export meanextrema
