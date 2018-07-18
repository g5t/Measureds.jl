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
use of `Measured` values.

`Measured` objects have two fields, `val` and `var`, the measured value and variance, respectively.
The value, variance, and square-root of the variance (the error) of a `Measured` object `m` can be
accessed through the functions `val(m)`, `var(m)`, and `err(m)`, respectively.
"""
abstract Measured <: Real
immutable AsymmetricMeasured <: Measured
    val::Float64
    varp::Float64
    varm::Float64
end
immutable SymmetricMeasured <: Measured
    val::Float64
    var::Float64
end
Measured(a::Real=0)=SymmetricMeasured(a,0.) # by default create symmetric uncertainty
@vectorize_2arg Real Measured
@vectorize_1arg Real Measured
Measured(m::Measured)=m # no need to worry about copying because Measured values are immutable
@vectorize_1arg Measured Measured

zero(::Type{Measured})=Measured(0)
one(::Type{Measured})=Measured(1)

+(a::SymmetricMeasured,b::SymmetricMeasured)=SymmetricMeasured(a.val+b.val, a.var+b.var)
-(a::SymmetricMeasured,b::SymmetricMeasured)=SymmetricMeasured(a.val-b.val, a.var+b.var)
-(a::SymmetricMeasured)=SymmetricMeasured(-a.val,a.var)
/(x::SymmetricMeasured,y::SymmetricMeasured)=SymmetricMeasured(x.val./y.val, abs(1./y.val).^2.*x.var + abs(x.val./y.val.^2).^2.*y.var)
*(x::SymmetricMeasured,y::SymmetricMeasured)=SymmetricMeasured(x.val.*y.val, abs(y.val).^2.*x.var + abs(x.val).^2.*y.var)

+(a::AsymmetricMeasured,b::AsymmetricMeasured)=AsymmetricMeasured(a.val+b.val, a.varp+b.varp, a.varm+b.varm)
-(a::AsymmetricMeasured,b::AsymmetricMeasured)=AsymmetricMeasured(a.val-b.val, a.varp+b.varp, a.varm+b.varm)
-(a::AsymmetricMeasured)=AsymmetricMeasured(-a.val,a.varm,a.varp) # insofar as negating is reflecting about zero, this makes sense
/(x::AsymmetricMeasured,y::AsymmetricMeasured)=AsymmetricMeasured(x.val./y.val, abs(1./y.val).^2.*x.varp + abs(x.val./y.val.^2).^2.*y.varp, abs(1./y.val).^2.*x.varm + abs(x.val./y.val.^2).^2.*y.varm)
*(x::AsymmetricMeasured,y::AsymmetricMeasured)=AsymmetricMeasured(x.val.*y.val, abs(y.val).^2.*x.varp + abs(x.val).^2.*y.varp, abs(y.val).^2.*x.varm + abs(x.val).^2.*y.varm)

# let promotions deal with Asymmetric*Symmetric, etc.
.*(x::Measured,y::Measured)=x*y
./(x::Measured,y::Measured)=x/y

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
isless(a::SymmetricMeasured,b::SymmetricMeasured)=isless(a.val+sqrt(a.var),b.val-sqrt(b.var))
isless(a::AsymmetricMeasured,b::AsymmetricMeasured)=isless(a.val+sqrt(a.varp),b.val-sqrt(b.varm))
isless(a::SymmetricMeasured,b::AsymmetricMeasured)=isless(a.val+sqrt(a.var),b.val-sqrt(b.varm))
isless(a::AsymmetricMeasured,b::SymmetricMeasured)=isless(a.val+sqrt(a.varp),b.val-sqrt(b.var))
isapprox(a::Measured,b::Measured)=!(isless(a,b)|isless(b,a))
isless{T<:AbstractFloat}(a::Measured,b::T)=isless(a.val,b)
isless{T<:AbstractFloat}(a::T,b::Measured)=isless(a,b.val)
isless{T<:Real}(a::Measured,b::T)=isless(a.val,b)
isless{T<:Real}(a::T,b::Measured)=isless(a,b.val)
isapprox{T<:Real}(a::Measured,b::T)=isapprox(a.val,b)
isapprox{T<:Real}(a::T,b::Measured)=isapprox(a,b.val)

^(x::SymmetricMeasured,n::Integer) =SymmetricMeasured(x.val.^n, abs(n*x.val.^(n-1)).^2.*x.var)
^(x::SymmetricMeasured,n::Rational)=SymmetricMeasured(x.val.^n, abs(n*x.val.^(n-1)).^2.*x.var)
^(x::SymmetricMeasured,n::Real    )=SymmetricMeasured(x.val.^n, abs(n*x.val.^(n-1)).^2.*x.var)
cos(x::SymmetricMeasured)=SymmetricMeasured(cos(x.val), abs(sin(x.val)).^2.*x.var)
sin(x::SymmetricMeasured)=SymmetricMeasured(sin(x.val), abs(cos(x.val)).^2.*x.var)
acos(x::SymmetricMeasured)=SymmetricMeasured(acos(x.val),x.var/(1-x.val^2))
asin(x::SymmetricMeasured)=SymmetricMeasured(asin(x.val),x.var/(1-x.val^2))
atan(x::SymmetricMeasured)=SymmetricMeasured(atan(x.val),x.var/(1+x.val^2)^2)
sqrt(x::SymmetricMeasured)=SymmetricMeasured(sqrt(x.val),x.var/2/x.val)

^(x::AsymmetricMeasured,n::Integer) =AsymmetricMeasured(x.val.^n, abs(n*x.val.^(n-1)).^2.*x.varp, abs(n*x.val.^(n-1)).^2.*x.varm)
^(x::AsymmetricMeasured,n::Rational)=AsymmetricMeasured(x.val.^n, abs(n*x.val.^(n-1)).^2.*x.varp, abs(n*x.val.^(n-1)).^2.*x.varm)
^(x::AsymmetricMeasured,n::Real)    =AsymmetricMeasured(x.val.^n, abs(n*x.val.^(n-1)).^2.*x.varp, abs(n*x.val.^(n-1)).^2.*x.varm)
cos(x::AsymmetricMeasured)=AsymmetricMeasured(cos(x.val), abs(sin(x.val)).^2.*x.varp, abs(sin(x.val)).^2.*x.varm)
sin(x::AsymmetricMeasured)=AsymmetricMeasured(sin(x.val), abs(cos(x.val)).^2.*x.varp, abs(cos(x.val)).^2.*x.varm)
acos(x::AsymmetricMeasured)=AsymmetricMeasured(acos(x.val),x.varp/(1-x.val^2),  x.varm/(1-x.val^2))
asin(x::AsymmetricMeasured)=AsymmetricMeasured(asin(x.val),x.varp/(1-x.val^2),  x.varm/(1-x.val^2))
atan(x::AsymmetricMeasured)=AsymmetricMeasured(atan(x.val),x.varp/(1+x.val^2)^2,x.varm/(1+x.val^2)^2)
sqrt(x::AsymmetricMeasured)=AsymmetricMeasured(sqrt(x.val),x.varp/2/x.val,      x.varm/2/x.val)

tan(x::Measured)=sin(x)/cos(x)

atan2(y::SymmetricMeasured,x::SymmetricMeasured)=SymmetricMeasured(atan2(y.val,x.val), x.val^2*( (y.val/x.val)^2 * x.var + y.var)/(x.val^2+y.val^2)^2 )
atan2(y::AsymmetricMeasured,x::AsymmetricMeasured)=AsymmetricMeasured(atan2(y.val,x.val), x.val^2*( (y.val/x.val)^2 * x.varp + y.varp)/(x.val^2+y.val^2)^2, x.val^2*( (y.val/x.val)^2 * x.varm + y.varm)/(x.val^2+y.val^2)^2 )
atan2(y::SymmetricMeasured,x::AsymmetricMeasured)=AsymmetricMeasured(atan2(y.val,x.val), x.val^2*( (y.val/x.val)^2 * x.varp + y.var)/(x.val^2+y.val^2)^2, x.val^2*( (y.val/x.val)^2 * x.varm + y.var)/(x.val^2+y.val^2)^2 )
atan2(y::AsymmetricMeasured,x::SymmetricMeasured)=AsymmetricMeasured(atan2(y.val,x.val), x.val^2*( (y.val/x.val)^2 * x.var + y.varp)/(x.val^2+y.val^2)^2, x.val^2*( (y.val/x.val)^2 * x.var + y.varm)/(x.val^2+y.val^2)^2 )

log(x::SymmetricMeasured)=SymmetricMeasured(log(x.val), abs(1./x.val).^2 .* x.var)
log(x::AsymmetricMeasured)=AsymmetricMeasured(log(x.val), abs(1./x.val).^2 .* x.varp, abs(1./x.val).^2 .* x.varm)

log10(x::Measured)=log(x)./log(10.)
@vectorize_1arg Measured sqrt
@vectorize_1arg Measured cos
@vectorize_1arg Measured sin
@vectorize_1arg Measured tan
@vectorize_1arg Measured acos
@vectorize_1arg Measured asin
@vectorize_1arg Measured atan
@vectorize_2arg Measured atan2
@vectorize_1arg Measured log
@vectorize_1arg Measured log10

Base.abs(x::SymmetricMeasured)=SymmetricMeasured(Base.abs(x.val),x.var)
Base.abs(x::AsymmetricMeasured)=AsymmetricMeasured(Base.abs(x.val),x.varp,xvarm)
Base.isnan(x::Measured)=Base.isnan(x.val)
Base.isfinite(x::Measured)=Base.isfinite(x.val)
Base.signbit(a::Measured)=Base.signbit(a.val)

#function discernablyUnique{T<:Measured}(x::AbstractArray{T,1})
#    xval=val(x)
#    xerr=err(x)
#    # BitArrays are **not** initialized to zero! Use "falses" instead
#    found=falses(size(x)...)
#    du=falses(size(x)...)
#    while sum(found)<length(x)
#        i=findfirst(!found)
#        withinerror = abs(xval-xval[i]) .<= xerr[i] # <= to ensure at least the current point is selected (if xerr[i]==0)
#        found = found | withinerror
#        du[i]=true
#    end
#    xval[du]
#end

function discernablyUnique{T<:Measured}(x::AbstractArray{T,1})
    found=falses(size(x)...)
    du=falses(size(x)...)
    while sum(found)<length(x)
        i=findfirst(!found)
        found = found | map(isapprox,x,x[i]) # take into account (possibly) assymetric errors
        du[i]=true
    end
    xval[du]
end
        

## this is dumb: (squaring a complex number causes promotion of the number 
#Base.isinteger(::Measured)=true

Base.size(x::Measured)=tuple() # like any single value, the size of a Measured is ()

# promote_array_type used by .+ and .-
#Base.promote_array_type{T<:Real}(F,::Type{T},::Type{Measured})=Measured
#Base.promote_array_type(F,::Type{Measured},::Type{Measured})=Measured

Base.promote_array_type{T<:AbstractFloat,M<:Measured}(F,::Type{M},::Type{T})=M
Base.promote_array_type{T<:Real,M<:Measured}(F,::Type{M},::Type{T})=M


# promote_type seems to be used by .* and ./ broadcasts
Base.promote_type(::AsymmetricMeasured,::SymmetricMeasured)=AsymmetricMeasured
Base.promote_type{M<:Measured}(::Type{M},::Type{M})= M
Base.promote_type{M<:Measured}(::Type{M},::Type{Union{}})= M
Base.promote_type{M<:Measured}(::Type{Union{}},::Type{M})= M
Base.promote_type{T<:Real,M<:Measured}(::Type{M},::Type{T})= M
Base.promote_type{T<:Real,M<:Measured}(::Type{T},::Type{M})= M

# promote seems to be used by + and -
Base.promote(a::AsymmetricMeasured,b::SymmetricMeasured)=(a,AsymmetricMeasured(b.val,b.var,b.var))
Base.promote(a::SymmetricMeasured,b::AsymmetricMeasured)=(AsymmetricMeasured(a.val,a.var,a.var),b)
Base.promote{M<:Measured}(a::M,b::M)=(a,b)
Base.promote(a::Real,b::Measured)=(Measured(a),b)
Base.promote(a::Measured,b::Real)=(a,Measured(b))

# promote_rule appears to be used by .^
Base.promote_ruleM<:Measured}(::Type{M},::Type{M})= M
Base.promote_rule{T<:Real,M<:Measured}(::Type{M},::Type{T})= M
Base.promote_rule{T<:Real,M<:Measured}(::Type{T},::Type{M})= M

Base.convert(::Type{Bool},x::Measured)=convert(Bool,x.val)
Base.convert(::Type{Integer},x::Measured)=convert(Integer,x.val)
Base.convert(::Type{Rational{Base.GMP.BigInt}},x::Measured)=convert(Rational{GMP.BigInt},x.val)
Base.convert{T<:Integer}(::Type{Rational{T}},x::Measured)=convert(Rational{T},x.val)
Base.convert(::Type{AsymmetricMeasured},x::SymmetricMeasured)=AsymmetricMeasured(x.val,x.var,x.var)
Base.convert(::Type{SymmetricMeasured},x::AsymmetricMeasured)=isapprox(x.varp,x.varm)?SymmetricMeasured(x.val,x.varp):Throw(InexactError())
Base.convert{M<:Measured}(::Type{M},x::M)=x
#Base.convert{S<:Real}(::Type{S},x::Measured)=convert(S,x.val)
Base.convert(::Type{Measured},x::Real)=Measured(convert(Float64,x)) # allows for, e.g., Measured[1,2,3,4]
Base.convert{T<:AbstractFloat}(::Type{T},x::Measured)=convert(T,x.val)



val(m::Measured)=m.val
var(m::SymmetricMeasured)=m.var
var(m::AsymmetricMeasured)=(m.varp+m.varm)/2
err(m::Measured)=sqrt(var(m))
@vectorize_1arg Measured val
@vectorize_1arg Measured var
@vectorize_1arg Measured err


pre(x)=(x<0?-1:1)*floor(Integer,abs(x))
zeropad(a::Real,b::Int)= ( n=a>0?round(Int,floor(log10(a)))+1:1; b-=n; (b>0?"$("0"^b)":"")*"$a" )
post(x,y)=zeropad(round(Integer,abs(x-pre(x))/10.0^y),abs(y))
function showMeasured(io::IO,m::Measured,compact::Bool=true)
    v,u = val(m),err(m)
    if u>0 # some uncertainty
#        if compact
            pv=(v!=0)?round(Integer,floor(log10(abs(v)))):0 # power of first digit in the value
            pu=round(Integer,floor(log10(u))) # power of first digit in the error
            dv=round(Integer,v/10.0^pu) # the value rounded to the precision of 1 digit in the error
            if dv==0 && v!=0 # the value isn't zero, but the rounded value is
                while dv==0; pu=pu-1; dv=round(Integer,v/10.0^pu); end # rounded value is zero but actual value is non-zero, increase the precision
            end
            du=round(Integer,u/10.0^pu) # the error rounded to the same precision
            # now actually do the formatting
            if (pv<0&&pu<0) # we need to try it both ways
                str1="$dv($du)ᴇ$pu"
                str2=((v<0)?"-":"")*"0."*"$(post(v,pu))($du)"
                str=(length(str1)<(length(str2)-3))?str1:str2 # pick str2 preferentially
            elseif pu>0 # two ways to try as well
                tdv=round(v/10.0^pv)
                tv=(v-tdv*10.0^pv)/10.0.^pv
                str1=((v<0)?"-":"")*"$(pre(tdv))."*"$(post(tv,pv-pu))($du)ᴇ$pv"
                str2="$(dv*10^pu)($(du*10^pu))"
                str=(length(str1)<(length(str2)-3))?str1:str2 # pick str2 preferentially
            elseif pu<0
                str="$(pre(v))"*"."*"$(post(v,pu))($du)"
            else # pu==0
                str="$dv($du)"
            end
            print(io,str)
#        else
#            show(io,v)
#            print(io," ± ")
#            showcompact(io,u)
#        end
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

## Plotting of Measured arrays
#function PyPlot.plot(x::AbstractArray{Measured,1},y::AbstractArray{Measured,1};marker="o",color="r",linestyle=" ",fillstyle="white")
#    (fillstyle == "white" || fillstyle == "w") ? (facecolor="white";fillstyle="full";):(facecolor=color);
#    h=PyPlot.errorbar(val(x),val(y),xerr=err(x),yerr=err(y))
#    PyPlot.setp(h[1],color=facecolor,marker=marker,markeredgecolor=color,markerfacecoloralt=color,fillstyle=fillstyle,linestyle=linestyle)
#    PyPlot.setp(h[2],marker=" ",color=color)
#    PyPlot.setp(h[3],color=color)
#    return h
#end
#PyPlot.plot{T<:Real}(x::AbstractArray{T,1},y::AbstractArray{T,1},o...)=PyPlot.plot(Measured(x),Measured(y),o...)
##plot{T<:Real}(x::AbstractArray{T,1},y::AbstractArray{Measured,1},o...)=plot(Measured(x),y,o...)
##plot{T<:Real}(x::AbstractArray{Measured,1},y::AbstractArray{T,1},o...)=plot(x,Measured(y),o...)

# Plotting of Measured arrays
function PyPlot.errorbar(x::AbstractArray{Measured,1},y::AbstractArray{Measured,1};
          marker="o",color=[1,0,0],linestyle="None",facecolor="white",fillstyle="full",masked=false,k...)
    (color,facecolor)=decodecolor(color,facecolor)
    masked && (color=color/3+2[1,1,1]/3;facecolor=facecolor/3+2[1,1,1]/3)
    h=PyPlot.errorbar(val(x),val(y),xerr=err(x),yerr=err(y))
    PyPlot.setp(h[1],color=facecolor,marker=marker,markeredgecolor=color,markerfacecoloralt=color,fillstyle=fillstyle,linestyle=linestyle)
    PyPlot.setp(h[2],marker=" ",color=color)
    PyPlot.setp(h[3],color=color)
    return h
end
#PyPlot.plot{T<:Real}(x::AbstractArray{T,1},y::AbstractArray{T,1},o...)=PyPlot.plot(Measured(x),Measured(y),o...)

function PyPlot.plot(x::AbstractArray{Measured,1},y::AbstractArray{Measured,1};
                    marker="",color="b",linestyle="-",facecolor="white",alpha=0.5,arealinestyle="",hatch="",fillstyle="None",k...)
    (color,facecolor)=decodecolor(color,facecolor)
    hatch!=""? (fillfacecolor="None";edgecolor=color;filllinewidth=0.):(fillfacecolor=color;edgecolor="None";filllinewidth=0.)
    #xlm=transpose(val(x).+hcat(-err(x),err(x)))[:]
    xlm=transpose(hcat(val(x),val(x)))[:] # ignore x-uncertainties for now TODO find a good way to include x-uncertainties in determining shaded area
    yll=transpose(val(y).-hcat(err(y),err(y)))[:]
    ymm=transpose(val(y).+hcat(err(y),err(y)))[:]
    hf=PyPlot.fill_between(xlm,yll,ymm,facecolor=fillfacecolor,alpha=alpha,hatch=hatch,edgecolor=edgecolor,linewidth=filllinewidth)
    h=PyPlot.plot(val(x),val(y),linestyle=linestyle,color=color)
end
PyPlot.plot{T<:Real}(x::AbstractArray{T,1},y::AbstractArray{Measured,1};k...)=PyPlot.plot(Measured(collect(x)),y;k...)
PyPlot.errorbar{T<:Real}(x::AbstractArray{T,1},y::AbstractArray{Measured,1};k...)=PyPlot.errorbar(Measured(collect(x)),y;k...)

   


#end # module vv



