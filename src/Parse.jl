"""
Array{,1}
`Measureds._syntax` is a `Vector` of `Tuple{Regex,Vector{Int},Vararg{AbstractArray{Int,1},N}}` that contains
regular expressions which match valid string representations of values and their symmetric or
asymmetric errors. In such a string, a relative uncertainty is indicated by one or two numbers in
parentheses following a value, e.g., `1.234(5)` indicates an uncertainty of plus or minus `5` in the
last digit of `1.234`; alternatively the same uncertainty can be specified as an absolute value as
`1.234±0.005`.

There are ten valid syntax included in `Measureds._syntax`, in order they are:

| syntax | sym | rel |
|:---:|:---:|:---:|
| `val`(`err`)                                   | symmetric  | relative |
| `val`(`err`)×10`pow`                           | symmetric  | relative |
| `val`±`err`                                    | symmetric  | absolute |
| (`val`±`err`)×10`pow`                          | symmetric  | absolute |
| `val`×10`vpow`±`err`×10`epow`                  | symmetric  | absolute |
| `val`(`pos`,`neg`)                             | asymmetric | relative |
| `val`(`pos`,`neg`)×10`pow`                     | asymmetric | relative |
| `val`±(`pos`,`neg`)                            | asymmetric | absolute |
| (`val`±(`pos`,`neg`))×10`pow`                  | asymmetric | absolute |
| `val`×10`vpow`±(`pos`×10`ppow`,`neg`×10`npow`) | asymmetric | absolute |

where `val` represents its value, `err` symmetric uncertainty, `pow` overall exponent,
`vpow` value-specific exponent, `epow` uncertainty-exponent, `pos` positive uncertainty,
`neg` negative uncertainty, `ppow` postive-uncertainty-exponent, and `npow` negative-uncertainty
exponent. In the above table, each `±` can alternatively be `+/-` or `+-`, and each `×10` can
alternatively be `x10`, `e`, `E`, or `ᴇ`. The exponents `pow`, `vpow`, `epow`, `ppow`, and `npow`
must be integers (unsigned integers are taken to represent positive exponents) and can optionally
be raised characters, e.g., `×10⁻³`.

The first integer vector of each `Tuple` in `Measured._syntax` provides the access order of the
matched (and stored) substrings for each part of the `Measured` value and its uncertainty.
If present, the second integer vector indicates which of the accessed substrings represent
exponents.
"""
#const _syntax=[(r"^([0-9,\.]+)\(([0-9]+)\)$",[1,2]),
#               (r"^([0-9,\.]+)\(([0-9]+)\)(?:(?:x|×)10|(?:e|E|ᴇ))((?:\+|\-)?[0-9]+)$",[1,2,3,3],[3,4]),
#               (r"^([0-9,\.]+)(?:\+/?\-|±)([0-9,\.]+)$",[1,2]),
#               (r"^\(([0-9,\.]+)(?:\+/?\-|±)([0-9,\.]+)\)(?:(?:x|×)10|(?:e|E|ᴇ))((?:\+|\-)?[0-9]+)$",[1,2,3,3],[3,4]),
#               (r"^([0-9,\.]+)(?:e|E|ᴇ|[x,×]10)((?:\+|\-)?[0-9]+)(?:\+/?\-|±)([0-9,\.]+)(?:e|E|ᴇ|[x,×]10)((?:\+|\-)?[0-9]+)$",[1,3,2,4],[3,4]),
#               (r"^([0-9,\.]+)\(\+?([0-9]+),\-?([0-9]+)\)$",[1,2,3]),
#               (r"^([0-9,\.]+)\(\+?([0-9]+),\-?([0-9]+)\)(?:e|E|ᴇ|[x,×]10)((?:\+|\-)?[0-9]+)$",[1,2,3,4,4,4],4:6),
#               (r"^([0-9,\.]+)(?:\+/?\-|±)\(\+?([0-9,\.]+),\-?([0-9,\.]+)\)$",[1,2,3]),
#               (r"^\(([0-9,\.]+)(?:\+/?\-|±)\(\+?([0-9,\.]+),\-?([0-9,\.]+)\)\)(?:e|E|ᴇ|[x,×]10)((?:\+|\-)?[0-9]+)$",[1,2,3,4,4,4],4:6),
#               (r"^([0-9,\.]+)(?:e|E|ᴇ|[x,×]10)((?:\+|\-)?[0-9]+)(?:\+/?\-|±)\(\+?([0-9,\.]+)(?:e|E|ᴇ|[x,×]10)((?:\+|\-)?[0-9]+),\-?([0-9,\.]+)(?:e|E|ᴇ|[x,×]10)((?:\+|\-)?[0-9]+)\)$",[1,3,5,2,4,6],4:6)]
#((?:\+|\-|⁻|⁺)?[0-9,⁰,¹,²,³,⁴,⁵,⁶,⁷,⁸,⁹]+)
const _syntax=[(r"^([0-9,\.]+)\(([0-9]+)\)$",[1,2]),
               (r"^([0-9,\.]+)\(([0-9]+)\)(?:(?:x|×)10|(?:e|E|ᴇ))((?:\+|\-|⁻|⁺)?[0-9,⁰,¹,²,³,⁴,⁵,⁶,⁷,⁸,⁹]+)$",[1,2,3,3],[3,4]),
               (r"^([0-9,\.]+)(?:\+/?\-|±)([0-9,\.]+)$",[1,2]),
               (r"^\(([0-9,\.]+)(?:\+/?\-|±)([0-9,\.]+)\)(?:(?:x|×)10|(?:e|E|ᴇ))((?:\+|\-|⁻|⁺)?[0-9,⁰,¹,²,³,⁴,⁵,⁶,⁷,⁸,⁹]+)$",[1,2,3,3],[3,4]),
               (r"^([0-9,\.]+)(?:e|E|ᴇ|[x,×]10)((?:\+|\-|⁻|⁺)?[0-9,⁰,¹,²,³,⁴,⁵,⁶,⁷,⁸,⁹]+)(?:\+/?\-|±)([0-9,\.]+)(?:e|E|ᴇ|[x,×]10)((?:\+|\-|⁻|⁺)?[0-9,⁰,¹,²,³,⁴,⁵,⁶,⁷,⁸,⁹]+)$",[1,3,2,4],[3,4]),
               (r"^([0-9,\.]+)\(\+?([0-9]+),\-?([0-9]+)\)$",[1,2,3]),
               (r"^([0-9,\.]+)\(\+?([0-9]+),\-?([0-9]+)\)(?:e|E|ᴇ|[x,×]10)((?:\+|\-|⁻|⁺)?[0-9,⁰,¹,²,³,⁴,⁵,⁶,⁷,⁸,⁹]+)$",[1,2,3,4,4,4],4:6),
               (r"^([0-9,\.]+)(?:\+/?\-|±)\(\+?([0-9,\.]+),\-?([0-9,\.]+)\)$",[1,2,3]),
               (r"^\(([0-9,\.]+)(?:\+/?\-|±)\(\+?([0-9,\.]+),\-?([0-9,\.]+)\)\)(?:e|E|ᴇ|[x,×]10)((?:\+|\-|⁻|⁺)?[0-9,⁰,¹,²,³,⁴,⁵,⁶,⁷,⁸,⁹]+)$",[1,2,3,4,4,4],4:6),
               (r"^([0-9,\.]+)(?:e|E|ᴇ|[x,×]10)((?:\+|\-|⁻|⁺)?[0-9,⁰,¹,²,³,⁴,⁵,⁶,⁷,⁸,⁹]+)(?:\+/?\-|±)\(\+?([0-9,\.]+)(?:e|E|ᴇ|[x,×]10)((?:\+|\-|⁻|⁺)?[0-9,⁰,¹,²,³,⁴,⁵,⁶,⁷,⁸,⁹]+),\-?([0-9,\.]+)(?:e|E|ᴇ|[x,×]10)((?:\+|\-|⁻|⁺)?[0-9,⁰,¹,²,³,⁴,⁵,⁶,⁷,⁸,⁹]+)\)$",[1,3,5,2,4,6],4:6)]
const _syntax_has_relerr=BitArray([1,1,0,0,0,1,1,0,0,0])
function lowerexponents(a::AbstractString)
    s=split(a,"") # we need an array of characters
    raised=Dict("⁰"=>"0","¹"=>"1","²"=>"2","³"=>"3","⁴"=>"4","⁵"=>"5","⁶"=>"6","⁷"=>"7","⁸"=>"8","⁹"=>"9","⁻"=>"-","⁺"=>"+")
    for k in keys(raised)
        any(s.==k) && (s[s.==k]=raised[k])
    end
    join(s)
end
splitmeasuredstring{T<:Integer}(str::AbstractString,typ::Regex,pos::AbstractVector{T})=map(string,match(typ,str).captures[pos]) # map(string ensures the eltype is <:AbstractString and not a Union
function splitmeasuredstring{T<:Integer}(str::AbstractString,typ::Regex,pos::AbstractVector{T},epos::AbstractVector{T})
    cap=splitmeasuredstring(str,typ,pos)
    cap[epos]="E".*cap[epos]
    return cap
end
function relative_to_absolute_uncertainty!(sstring::Vector)
    ne=fld(length(sstring),2) # sstring is [#1,#2,#3...,exp1,exp2,exp3,...] we just care about #s here as all exp *should* be the same
    @assert all(sstring[ne+2:end].==sstring[ne+1]) "All relative errors should have the same exponent as the value!"
    cv=split(sstring[1],"") # split the first string, the value, into characters
    # number of digits beyond the decimal place, or 0 if 0 or multiple "." exist
    p=sum(cv.==".")==1?length(cv)-findfirst(cv,"."):0
    0==p && (return sstring)
    for i=2:ne # now loop over the relative uncertainties
        ssi=sstring[i]
        n=length(ssi)
        sstring[i]= n<p? "0."*"0"^(p-n)*ssi : ssi[1:(n-p)]*"."*ssi[(n-p+1):end]
    end
    return sstring
end
function decodemeasuredstring(str::AbstractString)
    matched=0; for i=1:length(_syntax); ismatch(_syntax[i][1],str) && (matched=i; break); end
    if 0==matched
#        try tmp=parse(str); catch problem # eval(parse(str)) is bad since it could interpret arbitrary code
#            typeof(problem)<:TypeError?error("decodemeasuredstring can't deal with $str, please expand."):rethrow(problem)
#        end
        tmp=parse(str)
        isa(tmp,Number)||error("decodemeasuredstring can't deal with $str, please expand.")
        return Measured(tmp) # tmp is a regular number, with no uncertainty
    end
    splitstr=splitmeasuredstring(lowerexponents(str),_syntax[matched]...)
    # by here we know that str represents some sort of measured value.
    # if l=length(splitstr) <= 3 then there is no exponent information, fake it now
    length(splitstr)<4 && ( splitstr=vcat(splitstr,repeat([""],outer=size(splitstr))) )
    # now splitstr is either ["val","err","val_exp","err_exp"] or ["val","pos","neg","val_exp","pos_exp","err_exp"]
    _syntax_has_relerr[matched] && relative_to_absolute_uncertainty!(splitstr) # convert err|pos,neg from relative to absolute if necessary
    n=fld(length(splitstr),2)
    vals=map(parse,splitstr[1:n].*splitstr[n+1:2n]) # ["1.23","0.01","E5","E5"] becomes ["1.23E5","0.01E5"], which is then parsed
    vals[2:n].^=2 # was [value,err|pos,neg] but needs to be [value,variance(s)...]
    return Measured(vals...) # let the default constructor decide if its symmetric or antisymmetric
end

Base.parse(::Type{Measured},s::AbstractString)=decodemeasuredstring(s)
# also allow for creation of Measured values by, e.g., m"5(1)"
macro m_str(p)
    decodmeasuredstring(p)
end
