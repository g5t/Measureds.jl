__precompile__(true)
module Measureds
    import Base: (==),(+),(-),(/),(*),(^),(./),(.*),(%),(<),(>),(<=),(>=),
    sqrt,cos,sin,tan,acos,asin,atan,atan2,log,log10,zero,one,prod,convert,
    sort,abs,isless,isapprox,isnan,sum,getindex,setindex!,cross,dot,norm,
    round,ceil,floor,size,ndims,convert,promote_rule,sum,bin,var,get,find,unique

    using PyPlot2TikZ # to allow for dispatch overloading of PyPlot plotting of Measureds
    import PyPlot2TikZ: plot,errorbar #,scatter FIXME once scatter is implemented properly in PyPlot2TikZ
    const plt=PyPlot2TikZ

    export Measured,MeasuredSymmetric,MeasuredAsymmetric
    export discernablyUnique,discernable,isless,isapprox,fuzzyindexin
    export val,var,vap,vam,err # these should all be deprecated
    export value,variance,positivevariance,negativevariance,uncertainty,positiveuncertainty,negativeuncertainty
    export decodecolor
    export countingerror,countingerror!,symmetricerror,asymmetricerror
    export symmetricerror!,asymmetricerror!,replaceerror!
    export @m_str

    include("DecodeColor.jl")
#    include("Measured.jl") # datatype definition, display, and mathmatical function overloading       ### Measured.jl
#    include("Plotting.jl") # overloading of PyPlot functions for Measured values                      ### Measured.jl
    include("Measured_asym.jl") # datatype definition, display, and mathmatical function overloading ### Measured_asym.jl
    include("Plotting_asym.jl") # overloading of PyPlot functions for Measured values                ### Measured_asym.jl

    include("Parse.jl") # ultimately overloading of Base.parse(::Measured)


    @deprecate(val,value)
    @deprecate(var,variance)
    @deprecate(vap,positivevariance)
    @deprecate(vam,negativevariance)
    @deprecate(err,uncertainty)
end # module
