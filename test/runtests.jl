using Measureds
using Base.Test

# test parsing of measured strings
# special characters for Measured values as strings: ± ᴇ ×
ms1="1.234(5)"
ms2="1.234±0.005"
m1=parse(Measured,ms1)
m2=parse(Measured,ms2)
@test m1==m2
ms1="1.234(5)ᴇ-3"
ms2="1.234e-3+-5x10-6"
ms3="(1.234±0.005)×10-3"
m1=parse(Measured,ms1)
m2=parse(Measured,ms2)
m3=parse(Measured,ms3)
@test m1==m2==m3

# now test asymmetric uncertainties
ms1="1.234(5,67)"
ms2="1.234±(0.005,0.067)"
m1=parse(Measured,ms1)
m2=parse(Measured,ms2)
@test m1==m2
ms1="1.234(5,67)E-3"
ms2="1.234x10-3+/-(5e-6,6.7E-5)"
ms3="(1.234+-(0.005,0.067))×10-3"
m1=parse(Measured,ms1)
m2=parse(Measured,ms2)
m3=parse(Measured,ms3)
@test m1==m2==m3
