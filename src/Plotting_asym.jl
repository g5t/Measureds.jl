# In these functions setting markeredgecolor and markerfacecoloralt in
# the initial call to errorbar, plot, etc., is neccessary to avoid the error
#     /.../matplotlib/lines.py:1150: UnicodeWarning:
#       Unicode unequal comparison failed to convert both arguments to Unicode
#      - interpreting them as being unequal
# which seems to arrise from a comparison of a (unicode?) string to a Vector{Float}
# which could be fixed, but not by me.

# Plotting of Measured arrays
function errorbar{T<:Measured,R<:Measured}(x::AbstractVector{T},y::AbstractVector{R};
          marker="o",color=[1,0,0],linestyle="None",facecolor="white",fillstyle="full",masked=false,k...)
    @assert size(x)==size(y) "x:$(size(x)) and y:$(size(y)) do not match"
    (color,facecolor)=decodecolor(color,facecolor)
    masked && (color=color/3+2[1,1,1]/3;facecolor=facecolor/3+2[1,1,1]/3)
    h=plt.errorbar(value.(x),value.(y),xerr=errforplot(x),yerr=errforplot(y),color=facecolor,marker=marker,markeredgecolor=color,markerfacecoloralt=color,fillstyle=fillstyle,linestyle=linestyle)
    # plt.setp(h[1],color=color,marker=marker,markeredgecolor=color,markerfacecoloralt=facecolor,fillstyle=fillstyle,linestyle=linestyle)
    plt.setp(h[2],marker=" ",color=color)
    plt.setp(h[3],color=color)
    return h
end
errorbar{T<:Real,R<:Measured}(x::AbstractVector{T},y::AbstractVector{R};k...)=errorbar(Measured(collect(x),zero(T)),y;k...)

circle(o...;k...)=PyPlot2TikZ.PyPlot.matplotlib[:patches][:Circle](o...;k...)
ellipse=PyPlot2TikZ.PyPlot.matplotlib[:patches][:Ellipse] # find a way to generalize this?


function errorbar{T<:Measured,R<:Measured,S<:Measured}(x::AbstractVector{T},y::AbstractVector{R},z::AbstractVector{S};
        color="k",linestyle="None",facecolor="white",fillstyle="full",masked=false,
        edgestyle="-",vmin=nothing,vmax=nothing,norm=nothing,scalezx=1,scalezy=1,k...)
    @assert size(x)==size(y)==size(z) "x:$(size(x)), y:$(size(x)) and z:$(size(y)) do not match"
    (color,facecolor)=decodecolor(color,facecolor)
    masked && (color=color/3+2[1,1,1]/3;facecolor=facecolor/3+2[1,1,1]/3)

    vx=value.(x)
    vy=value.(y)

    h=plt.errorbar(vx,vy,xerr=errforplot(x),yerr=errforplot(y),color=color,marker="",linestyle=linestyle)
    plt.setp(h[2],marker=" ",color=color)
    plt.setp(h[3],color=color)

    # now the "fun" part
    ax=gca().obj # find a way to generalize this?
    #circle=PyPlot2TikZ.PyPlot.matplotlib[:patches][:Circle] # find a way to generalize this?
    # XXX maybe use CirclePolygon instead?

    # here we need to deal with the color map
    minz=isa(vmin,Real)? vmin : minimum(z)
    maxz=isa(vmax,Real)? vmax : maximum(z)
    plotr=isa(norm,Function)? norm(z): (z-value(minz))/value(maxz-minz)

    r0=value.(plotr)
    rn=r0-negativeuncertainty.(plotr)
    rp=r0+positiveuncertainty.(plotr)

    r0[r0.<0]=0
    rn[rn.<0]=0
    rp[rp.<0]=0
    kes=["-","--",":","-."]
    wes=findfirst(kes.==edgestyle)
    nes=circshift(kes,-1)[ length(kes)>=wes>0 ? wes : 1]

    c0=[ellipse([xi,yi],width=scalezx*ri,height=scalezy*ri,angle=0,edgecolor=color,facecolor=facecolor,linestyle=edgestyle) for (xi,yi,ri) in zip(vx,vy,r0)]
    cm=[ellipse([xi,yi],width=scalezx*ri,height=scalezy*ri,angle=0,edgecolor=color,facecolor="none",linestyle=nes) for (xi,yi,ri) in zip(vx,vy,rn)]
    cz=[ellipse([xi,yi],width=scalezx*ri,height=scalezy*ri,angle=0,edgecolor=color,facecolor="none",linestyle=nes) for (xi,yi,ri) in zip(vx,vy,rp)]

    for i=1:length(c0)
        ax[:add_artist](c0[i])
        ax[:add_artist](cm[i])
        ax[:add_artist](cz[i])
        c0[i][:set_clip_box](ax[:bbox])
        cm[i][:set_clip_box](ax[:bbox])
        cz[i][:set_clip_box](ax[:bbox])
    end
    return h
end

# Define a 3D errorbar for use with, e.g., waterfallerrorbar plotting
function errorbar3D(x::AbstractVector{T},y::AbstractVector{R},z::AbstractVector{S};
    marker="o",color="k",linestyle="None",facecolor="white",fillstyle="full",masked=false,k...) where {T<:Measured,R<:Measured,S<:Measured}
    @assert size(x)==size(y)==size(z) "x:$(size(x)), y:$(size(x)) and z:$(size(y)) do not match"
    (color,facecolor)=decodecolor(color,facecolor)
    masked && (color=color/3+2[1,1,1]/3;facecolor=facecolor/3+2[1,1,1]/3)

    vx=value.(x); vy=value.(y); vz=value.(z); nn=NaN*ones(vx)
    px=vec(hcat(vx,vx,nn,vx,vx,nn,vx-negativeuncertainty.(x),vx+positiveuncertainty.(x),nn)')
    py=vec(hcat(vy,vy,nn,vy-negativeuncertainty.(y),vy+positiveuncertainty.(y),nn,vy,vy,nn)')
    pz=vec(hcat(vz-negativeuncertainty.(z),vz+positiveuncertainty.(z),nn,vz,vz,nn,vz,vz,nn)')

    plt.plot3D(px,py,pz,color=color,linestyle="-")
    plt.plot3D(vx,vy,vz,marker=marker,color=color,linestyle=linestyle,facecolor=facecolor,fillstyle=fillstyle)
end
function errorbar3D(x::AbstractVector{T},y::AbstractVector{R},z::AbstractVector{S};
    marker="o",color="k",linestyle="None",facecolor="white",fillstyle="full",masked=false,k...) where {T<:Measured,R<:Measured,S}
    @assert size(x)==size(y)==size(z) "x:$(size(x)), y:$(size(x)) and z:$(size(y)) do not match"
    (color,facecolor)=decodecolor(color,facecolor)
    masked && (color=color/3+2[1,1,1]/3;facecolor=facecolor/3+2[1,1,1]/3)

    vx=value.(x); vy=value.(y); vz=value.(z); nn=NaN*ones(vx)
    px=vec(hcat(vx,vx,nn,vx-negativeuncertainty.(x),vx+positiveuncertainty.(x),nn)')
    py=vec(hcat(vy-negativeuncertainty.(y),vy+positiveuncertainty.(y),nn,vy,vy,nn)')
    pz=vec(hcat(vz,vz,nn,vz,vz,nn)')

    plt.plot3D(px,py,pz,color=color,linestyle="-")
    plt.plot3D(vx,vy,vz,marker=marker,color=color,linestyle=linestyle,facecolor=facecolor,fillstyle=fillstyle)
end
function errorbar3D(x::AbstractVector{T},y::AbstractVector{R},z::AbstractVector{S};
    marker="o",color="k",linestyle="None",facecolor="white",fillstyle="full",masked=false,k...) where {T<:Measured,R,S<:Measured}
    @assert size(x)==size(y)==size(z) "x:$(size(x)), y:$(size(x)) and z:$(size(y)) do not match"
    (color,facecolor)=decodecolor(color,facecolor)
    masked && (color=color/3+2[1,1,1]/3;facecolor=facecolor/3+2[1,1,1]/3)

    vx=value.(x); vy=value.(y); vz=value.(z); nn=NaN*ones(vx)
    px=vec(hcat(vx,vx,nn,vx-negativeuncertainty.(x),vx+positiveuncertainty.(x),nn)')
    py=vec(hcat(vy,vy,nn,vy,vy,nn)')
    pz=vec(hcat(vz-negativeuncertainty.(z),vz+positiveuncertainty.(z),nn,vz,vz,nn)')

    plt.plot3D(px,py,pz,color=color,linestyle="-")
    plt.plot3D(vx,vy,vz,marker=marker,color=color,linestyle=linestyle,facecolor=facecolor,fillstyle=fillstyle)
end
function errorbar3D(x::AbstractVector{T},y::AbstractVector{R},z::AbstractVector{S};
    marker="o",color="k",linestyle="None",facecolor="white",fillstyle="full",masked=false,k...) where {T,R<:Measured,S<:Measured}
    @assert size(x)==size(y)==size(z) "x:$(size(x)), y:$(size(x)) and z:$(size(y)) do not match"
    (color,facecolor)=decodecolor(color,facecolor)
    masked && (color=color/3+2[1,1,1]/3;facecolor=facecolor/3+2[1,1,1]/3)

    vx=value.(x); vy=value.(y); vz=value.(z); nn=NaN*ones(vx)
    px=vec(hcat(vx,vx,nn,vx,vx,nn)')
    py=vec(hcat(vy,vy,nn,vy-negativeuncertainty.(y),vy+positiveuncertainty.(y),nn)')
    pz=vec(hcat(vz-negativeuncertainty.(z),vz+positiveuncertainty.(z),nn,vz,vz,nn)')

    plt.plot3D(px,py,pz,color=color,linestyle="-")
    plt.plot3D(vx,vy,vz,marker=marker,color=color,linestyle=linestyle,facecolor=facecolor,fillstyle=fillstyle)
end
function errorbar3D(x::AbstractVector{T},y::AbstractVector{R},z::AbstractVector{S};
    marker="o",color="k",linestyle="None",facecolor="white",fillstyle="full",masked=false,k...) where {T<:Measured,R,S}
    @assert size(x)==size(y)==size(z) "x:$(size(x)), y:$(size(x)) and z:$(size(y)) do not match"
    (color,facecolor)=decodecolor(color,facecolor)
    masked && (color=color/3+2[1,1,1]/3;facecolor=facecolor/3+2[1,1,1]/3)

    vx=value.(x); vy=value.(y); vz=value.(z); nn=NaN*ones(vx)
    px=vec(hcat(vx-negativeuncertainty.(x),vx+positiveuncertainty.(x),nn)')
    py=vec(hcat(vy,vy,nn)')
    pz=vec(hcat(vz,vz,nn)')

    plt.plot3D(px,py,pz,color=color,linestyle="-")
    plt.plot3D(vx,vy,vz,marker=marker,color=color,linestyle=linestyle,facecolor=facecolor,fillstyle=fillstyle)
end
function errorbar3D(x::AbstractVector{T},y::AbstractVector{R},z::AbstractVector{S};
    marker="o",color="k",linestyle="None",facecolor="white",fillstyle="full",masked=false,k...) where {T,R,S<:Measured}
    @assert size(x)==size(y)==size(z) "x:$(size(x)), y:$(size(x)) and z:$(size(y)) do not match"
    (color,facecolor)=decodecolor(color,facecolor)
    masked && (color=color/3+2[1,1,1]/3;facecolor=facecolor/3+2[1,1,1]/3)

    vx=value.(x); vy=value.(y); vz=value.(z); nn=NaN*ones(vx)
    px=vec(hcat(vx,vx,nn)')
    py=vec(hcat(vy,vy,nn)')
    pz=vec(hcat(vz-negativeuncertainty.(z),vz+positiveuncertainty.(z),nn)')

    plt.plot3D(px,py,pz,color=color,linestyle="-")
    plt.plot3D(vx,vy,vz,marker=marker,color=color,linestyle=linestyle,facecolor=facecolor,fillstyle=fillstyle)
end
function errorbar3D(x::AbstractVector{T},y::AbstractVector{R},z::AbstractVector{S};
    marker="o",color="k",linestyle="None",facecolor="white",fillstyle="full",masked=false,k...) where {T,R<:Measured,S}
    @assert size(x)==size(y)==size(z) "x:$(size(x)), y:$(size(x)) and z:$(size(y)) do not match"
    (color,facecolor)=decodecolor(color,facecolor)
    masked && (color=color/3+2[1,1,1]/3;facecolor=facecolor/3+2[1,1,1]/3)

    vx=value.(x); vy=value.(y); vz=value.(z); nn=NaN*ones(vx)
    px=vec(hcat(vx,vx,nn)')
    py=vec(hcat(vy-negativeuncertainty.(y),vy+positiveuncertainty.(y),nn)')
    pz=vec(hcat(nn,vz,vz,nn)')

    plt.plot3D(px,py,pz,color=color,linestyle="-")
    plt.plot3D(vx,vy,vz,marker=marker,color=color,linestyle=linestyle,facecolor=facecolor,fillstyle=fillstyle)
end

"""
    plot(x::Vector,y::Vector{Measured})
Overlay a `fill_between` and a `plot` of a `x::Vector{Union{Real,Measured}}` and a vector of
`Measured` `y` values. Uses keyword inputs: `marker`, `color`, `linestyle`, `facecolor`, `alpha`,
`arealinestyle`, `hatch`, `fillstyle`; ignoring any additional keyword arguments.
As currently implemented, there is no way to produce different `fill_between` fill- and `plot`
line-colors since both represent the uncertainty and value of `y`. One could rewite this function
if they desire such functionality.
"""
# for now strip the uncertainty from x if it is Measured as there is no (easy) way to include it in plotting
plot{T<:Measured}(x::AbstractVector{T},y::AbstractVector{T};k...)=plot(value.(x),y;k...)
plot{T<:Measured,R<:Measured}(x::AbstractVector{T},y::AbstractVector{R};k...)=plot(value.(x),y;k...)
function plot{T<:Real,R<:Measured}(x::AbstractVector{T},y::AbstractVector{R};
                                   marker="",
                                   color="b",
                                   linestyle="-",
                                   facecolor="white",
                                   alpha=0.5,
                                   arealinestyle="",
                                   hatch="",
                                   fillstyle="None",
                                   k...)
    @assert size(x)==size(y) "x:$(size(x)) and y:$(size(y)) do not match"
    (color,facecolor)=decodecolor(color,facecolor)
    hatch!=""?(fillfacecolor="None";edgecolor=color;filllinewidth=0.):(fillfacecolor=color;edgecolor="None";filllinewidth=0.)
    if any(uncertainty(y).>0) # only plot the uncertainty if any of it is non-zero
    hf=plt.fill_between(x,value.(y).-negativeuncertainty.(y),value.(y).+positiveuncertainty.(y),
                        facecolor=fillfacecolor,alpha=alpha,hatch=hatch,edgecolor=edgecolor,linewidth=filllinewidth)
    end
    h=plt.plot(x,value.(y),linestyle=linestyle,color=color)
end

# possibly useless scatter plotting of two Measured arrays and one (possibly Measured) Real array
errorscatter{T<:Measured,R<:Measured,S<:Measured}(x::AbstractVector{T},y::AbstractVector{R},z::AbstractVector{S};k...)=errorscatter(x,y,value.(z);k...)
function errorscatter{T<:Measured,R<:Measured,S<:Real}(x::AbstractVector{T},y::AbstractVector{R},z::AbstractVector{S};
                                                  edgecolors="face",cmap=nothing,norm=nothing,vmin=nothing,vmax=nothing,k...)
    @assert length(x)==length(y)==length(z)
    # we need to know which axis to insert our ellipses; gca creates an axis if none exist
    ax=gca().obj # find a way to generalize this?

    # here we need to deal with the color map
    minz=isa(vmin,Real)? vmin : minimum(z)
    maxz=isa(vmax,Real)? vmax : maximum(z)
    plotz=isa(norm,Function)? norm(z): (z-minz)/(maxz-minz)
    plotc=isa(cmap,PyPlot2TikZ.PyPlot.ColorMap)? cmap : PyPlot2TikZ.PyPlot.cm_get_cmap("inferno")

    fc=map(plotc,plotz)
    typeof(edgecolors)<:AbstractString && "face"==edgecolors && (edgecolors=copy(fc))

    drawedges=false; singleedgecolor=false
    if typeof(edgecolors)<:Vector
        # first check if it's a color spec (3 or 4 values)
        if eltype(edgecolors)<:Real && 2<length(edgecolors)<5
            drawedges=true; singleedgecolor=true;
        elseif length(edgecolors)==length(fc) && (eltype(edgecolors)<:AbstractString || (eltype(edgecolors)<:Vector && all(q->(eltype(q)<:Real)&(2<length(q)<5), edgecolors)))
            drawedges=true
        end
    elseif typeof(edgecolors)<:AbstractString
        drawedges=true; singleedgecolor=true;
    end
    drawedges||(edgecolors="none"; singleedgecolor=true)

    xm=value.(x)-negativeuncertainty.(x)
    xp=value.(x)+positiveuncertainty.(x)
    ym=value.(y)-negativeuncertainty.(y)
    yp=value.(y)+positiveuncertainty.(y)

    xc=(xp+xm)/2
    xw= xp-xm
    yc=(yp+ym)/2
    yw= yp-ym

    ells=[ ellipse([xi,yi],width=wi,height=hi,angle=0) for (xi,wi,yi,hi) in zip(xc,xw,yc,yw)]

    for i=sortperm(plotz) # always draw from lowest to highest intensity
        ax[:add_artist](ells[i])
        ells[i][:set_clip_box](ax[:bbox]) # this seems unneccessary
        ells[i][:set_facecolor](fc[i]) # alpha is included in the colormap, no need for [:set_alpha]()
    end

    if singleedgecolor
        foreach(q->q[:set_edgecolor](edgecolors),ells)
    else
        for i=1:length(ells); ells[i][:set_edgecolor](edgecolors[i]); end
    end

    # finally, expand the axis limits if necessary so that we can see the added points
    (xmin,xmax)=getp(gca(),"xlim")
    (ymin,ymax)=getp(gca(),"ylim")
    xmin=minimum(vcat(xmin,xm))
    xmax=maximum(vcat(xmax,xp))
    ymin=minimum(vcat(ymin,ym))
    ymax=maximum(vcat(ymax,yp))
    setp(gca(),"xlim",[xmin,xmax],"ylim",[ymin,ymax])
    return nothing
end

function scatter{T<:Measured,R<:Measured}(x::AbstractVector{T},y::AbstractVector{R};s=50,c="r",k...)
    @assert length(x)==length(y)
#    sok=(typeof(s)<:AbstractArray&eltype(s)<:Number&length(x)==length(s))|typeof(s)<:Number
#    @assert sok "keyword `s` must be either a number or an array of numbers the same length as `x` and `y`"
#    cok=(typeof(c)<:AbstractArray&(eltype(c)<:Number|eltype(c)<:AbstractString)&length(x)==length(c))|typeof(c)<:Number|typeof(c)<:AbstractString
#    @assert cok "keyword `c` must be either a number, string, or array of either the same length as `x` and `y`"

    # we're dealing with s and c here just in case either is an Array{Measured} or a Mesured
    # everything else can be handled by the PyPlot2TikZ.PyPlot.scatter (FIXME we want just PyPlot2TikZ.scatter to work)
    PyPlot2TikZ.PyPlot.scatter(value.(x),value.(y);s=value.(s),c=value.(c),k...)
end
