import sys
import os
sys.path.append( os.environ['AFQMCLAB_DIR']+"/scripts/supercubic" )
from setHoping import *

#Method Parameter
dt                                   = 0.005
thermalSize                          = 1000
writeNumber                          = 600
measureNumberPerWrite                = 2
measureSkipStep                      = 5
walkerSizePerThread                  = 70
decompTypeU                          = "densitySpin"  # "densityCharge", "densitySpin"
decompTypeU1                         = "densitySpin"  # "densityCharge", "densitySpin" 
decompTypeU2J                        = "densitySpin"  # "densityCharge", "densitySpin"
decompTypeJ                          = "spin"         # "charge", "spin"
forceType                            = "dynamicForce" # "dynamicForce", "constForce"
forceCap                             = 1.5
initialPhiTFlag                      = "setFromModel" #"setFromModel", "setRandomly", "readFromFile"
initialWalkerFlag                    = "sampleFromPhiT"  #"setFromModel", "setRandomly", "sampleFromPhiT","readFromFile","readAllWalkers"
mgsStep                              = 10
popControlStep                       = 10
initPopControlMaxSize                = 0
logEnergyCap                         = 1.62635
ET                                   = -86.6066
backGroundETInit                       = "EstimateFromPhiTWalker" #"EstimateFromPhiTWalker", "EstimateFromPhiT", "readFromFile"
ETAdjustStep                         = 10
ETAdjustMaxSize                      = 200
ETAndBackGroundGrowthEstimateStep    = 5
ETAndBackGroundGrowthEstimateMaxSize = 1000
seed                                 = 985456376

#write method_param
f = open('afqmc_param', 'w')
f.write(" {:<36} {:<26.18e} \n".format("dt", dt ) )
f.write(" {:<36} {:<26} \n".format("thermalSize", thermalSize) )
f.write(" {:<36} {:<26} \n".format("writeNumber", writeNumber) )
f.write(" {:<36} {:<26} \n".format("measureNumberPerWrite", measureNumberPerWrite) )
f.write(" {:<36} {:<26} \n".format("measureSkipStep", measureSkipStep) )
f.write(" {:<36} {:<26} \n".format("walkerSizePerThread",walkerSizePerThread) )
f.write(" {:<36} {:<26} \n".format("decompTypeU",decompTypeU) )
f.write(" {:<36} {:<26} \n".format("decompTypeU1",decompTypeU1) )
f.write(" {:<36} {:<26} \n".format("decompTypeU2J", decompTypeU2J))
f.write(" {:<36} {:<26} \n".format("forceType", forceType) )
f.write(" {:<36} {:<26.18e} \n".format("forceCap", forceCap) )
f.write(" {:<36} {:<26} \n".format("initialPhiTFlag",initialPhiTFlag) )
f.write(" {:<36} {:<26} \n".format("initialWalkerFlag",initialWalkerFlag) )
f.write(" {:<36} {:<26} \n".format("mgsStep", mgsStep) )
f.write(" {:<36} {:<26} \n".format("popControlStep", popControlStep) )
f.write(" {:<36} {:<26} \n".format("initPopControlMaxSize", initPopControlMaxSize) )
f.write(" {:<36} {:<26.18e} \n".format("logEnergyCap", logEnergyCap) )
f.write(" {:<36} {:<26.18e} \n".format("ET", ET) )
f.write(" {:<36} {:<26} \n".format("backGroundETInit", backGroundETInit) )
f.write(" {:<36} {:<26} \n".format("ETAdjustStep", ETAdjustStep) )
f.write(" {:<36} {:<26} \n".format("ETAdjustMaxSize", ETAdjustMaxSize) )
f.write(" {:<36} {:<26} \n".format("ETAndBackGroundGrowthEstimateStep", ETAndBackGroundGrowthEstimateStep) )
f.write(" {:<36} {:<26} \n".format("ETAndBackGroundGrowthEstimateMaxSize", ETAndBackGroundGrowthEstimateMaxSize) )
f.write(" {:<36} {:<26} \n".format("seed", seed) )
f.close()

#Model Parameter
latt_n     = [4,3]
ktwist     = [0.012,0.034]
t1         = 1.0
U          = 4.0
mu         = 0.0
Ntot       = 10
UpDnFlag   = 0    # 0 up=dn, 1 up=conj(dn) ==> different twist

numberkana = 2
site_i     = [0, 1]
site_j     = [1, 2]
U1         = [0, 0]
U2         = [0, 0]
J          = [0, 0]

#Set lattice information
latt = Latt_class( latt_n )
up_i, up_j, up_K = HubbardNearestNeighborHopping(latt, ktwist, t1)
if UpDnFlag == 0:
    dn_i = up_i; dn_j = up_j; dn_K = up_K
elif UpDnFlag == 1:
    dn_i, dn_j, dn_K = HubbardNearestNeighborHopping(latt, -np.array(ktwist), t1)
else:
    print( "WRONG!!! Do not know UpDnFlag!!!" )
    sys.exit(1)

Tmatrix = np.zeros((2 * latt.L, 2 * latt.L), dtype='complex', order='F')
for i in range( len(up_K) ):
    Tmatrix[up_i[i], up_j[i]] += up_K[i]
for i in range( len(dn_K) ):
    Tmatrix[dn_i[i] + latt.L, dn_j[i] + latt.L] += dn_K[i]

f = open("model_param", 'w')
f.write( '{:16d} \n'.format(latt.L) )
f.write( '{:16d} \n'.format(Ntot) )
f.write( '{:16d} \n'.format(numberkana) )
f.write( '{:>16} \n'.format(decompTypeJ) )
for i in range( 2*latt.L ):
    for j in range( 2*latt.L ):
        f.write( '{:26.18e} {:26.18e} \n'.format(Tmatrix[j, i].real, Tmatrix[j, i].imag))
for i in range( latt.L ):
    f.write( '{:26.18e} \n'.format( U ) )
for i in range( numberkana ):
    f.write( '{:16d} \n'.format(site_i[i]) )
for i in range( numberkana ):
    f.write( '{:16d} \n'.format(site_j[i]) )
for i in range( numberkana ):
    f.write( '{:26.18e} \n'.format(U1[i]) )
for i in range( numberkana ):
    f.write( '{:26.18e} \n'.format(U2[i]) )
for i in range( numberkana ):
    f.write( '{:26.18e} \n'.format(J[i]) )
for i in range( numberkana ):
    f.write( '{:26.18e} \n'.format(0.0) )
f.close()

#Check Kmatrix
# Kmatrix  = Tmatrix.copy()
# for i in range(numberkana):
#     ki=site_i[i]
#     Kmatrix[ki, ki] = Kmatrix[ki, ki] - J[i]*0.5
#     Kmatrix[ki+latt.L, ki+latt.L] = Kmatrix[ki+latt.L, ki+latt.L] - J[i]*0.5
#
#     kj=site_j[i]
#     Kmatrix[kj, kj] = Kmatrix[kj, kj] - J[i]*0.5
#     Kmatrix[kj+latt.L, kj+latt.L] = Kmatrix[kj+latt.L, kj+latt.L] - J[i]*0.5
#
# w, v = np.linalg.eigh(Kmatrix)
# print( "{:<26.18f}".format( np.sum( w[0:Ntot] ) ) )
