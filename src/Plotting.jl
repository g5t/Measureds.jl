
# Plotting of Measured arrays
function errorbar(x::AbstractArray{Measured,1},y::AbstractArray{Measured,1};
          marker="o",color=[1,0,0],linestyle="None",facecolor="white",fillstyle="full",masked=false,k...)
    (color,facecolor)=decodecolor(color,facecolor)
    masked && (color=color/3+2[1,1,1]/3;facecolor=facecolor/3+2[1,1,1]/3)
    h=errorbar(val(x),val(y),xerr=err(x),yerr=err(y))
    setp(h[1],color=facecolor,marker=marker,markeredgecolor=color,markerfacecoloralt=color,fillstyle=fillstyle,linestyle=linestyle)
    setp(h[2],marker=" ",color=color)
    setp(h[3],color=color)
    return h
end
errorbar{T<:Real}(x::AbstractArray{T,1},y::AbstractArray{Measured,1};k...)=errorbar(Measured(collect(x)),y;k...)

#function plot(x::AbstractArray{Measured,1},y::AbstractArray{Measured,1};
#                    marker="",color="b",linestyle="-",facecolor="white",alpha=0.5,arealinestyle="",hatch="",fillstyle="None",k...)
#    (color,facecolor)=decodecolor(color,facecolor)
#    hatch!=""? (fillfacecolor="None";edgecolor=color;filllinewidth=0.):(fillfacecolor=color;edgecolor="None";filllinewidth=0.)
#    #xlm=transpose(val(x).+hcat(-err(x),err(x)))[:]
#    # ignore x-uncertainties for now TODO find a good way to include x-uncertainties in determining shaded area
#    hf=fill_between(val(x),val(y).-err(y),val(y).+err(y),facecolor=fillfacecolor,alpha=alpha,hatch=hatch,edgecolor=edgecolor,linewidth=filllinewidth)
#    h=plot(val(x),val(y),linestyle=linestyle,color=color)
#end
#plot{T<:Real}(x::AbstractArray{T,1},y::AbstractArray{Measured,1};k...)=plot(Measured(collect(x)),y;k...)
#

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
plot(x::AbstractArray{Measured,1},y::AbstractArray{Measured,1};k...)=plot(val(x),y;k...)
function plot{T<:Real}(x::AbstractArray{T,1},y::AbstractArray{Measured,1};
                       marker="",color="b",linestyle="-",facecolor="white",alpha=0.5,arealinestyle="",hatch="",fillstyle="None",k...)
    (color,facecolor)=decodecolor(color,facecolor)
    hatch!=""?(fillfacecolor="None";edgecolor=color;filllinewidth=0.):(fillfacecolor=color;edgecolor="None";filllinewidth=0.)
    hf=fill_between(x,val(y).-err(y),val(y).+err(y),facecolor=fillfacecolor,alpha=alpha,hatch=hatch,edgecolor=edgecolor,linewidth=filllinewidth)
    h=plot(x,val(y),linestyle=linestyle,color=color)
end

# possibly useless scatter plotting of two Measured arrays and one (possibly Measured) Real array
function scatter{T<:Real}(x::AbstractArray{Measured,1},y::AbstractArray{Measured,1},z::AbstractArray{T,1};
                          edgecolors="face",cmap=nothing,norm=nothing,vmin=nothing,vmax=nothing,k...
    @assert length(x)==length(y)==length(z)
    # we need to know which axis to insert our ellipses; gca creates an axis if none exist
    ax=gca().obj # find a way to generalize this?
    ellipse=PyPlot2TikZ.PyPlot.matplotlib[:patches][:Ellipse] # find a way to generalize this?
    
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
        elseif length(edgecolors)==length(fc) && (eltype(edgecolors)<:AbstractString || (eltype(edgecolors)<:Vector && all(x->(eltype(x)<:Real)&(2<length(x)<5), edgecolors)))
            drawedges=true
        end
    elseif typeof(edgecolors)<:AbstractString
        drawedges=true; singleedgecolor=true;
    end
    drawedges||(edgecolors="none"; singleedgecolor=true)
    
    xm=value(x)-uncertainty(x)
    xp=value(x)+uncertainty(x)
    ym=value(y)-uncertainty(y)
    yp=value(y)+uncertainty(y)
    
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
        foreach(x->x[:set_edgecolor](edgecolors),ells)
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
