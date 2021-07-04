import sys
import os
import subprocess
sys.path.append( os.environ['AFQMCLAB_DIR']+"/scripts/supercubic" )
from lattClass import *
from rham295K import *

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
backGroundETInit                     = "EstimateFromPhiTWalker" #"EstimateFromPhiTWalker", "EstimateFromPhiT", "readFromFile"
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
latt_n  = [1,1,1]
latt    = Latt_class( latt_n )
U_one   = 2.3
J_one   = 0.35
Ntot    = latt.L*2*4   # 2 Ru per unit cell, 4 electrons per Ru.

numberkana = latt.L*2*3 # 2 Ru per unit cell, 3 kanamori per Ru
site_i = []; site_j = []
for i in range(latt.L*2):
    site_i.append( 0+i*3 ); site_i.append( 0+i*3 ); site_i.append( 1+i*3 )
    site_j.append( 1+i*3 ); site_j.append( 2+i*3 ); site_j.append( 2+i*3 )

U1  = [U_one-2*J_one]*latt.L*6
U2  = [U_one-3*J_one]*latt.L*6
J   = [J_one]*latt.L*6

#Total state is 2*3*latt.L, additional 2 for soc
Tmatrix = np.zeros(( 12 * latt.L, 12 * latt.L), dtype='complex', order='F')
jumpx = [0] if latt.n[0]==1 else [-1,0,1]
jumpy = [0] if latt.n[1]==1 else [-1,0,1]
jumpz = [0] if latt.n[2]==1 else [-1,0,1]

if latt.n[0]==2: jumpx = [0, 1]
if latt.n[1]==2: jumpy = [0, 1]
if latt.n[2]==2: jumpz = [0, 1]

for i in range(latt.L):
    coor_i = latt.coor(i)
    coor_j = [0,0,0]
    for dx in jumpx:
        coor_j[0] = latt.bound( coor_i[0]+dx, latt.n[0]  )
        for dy in jumpy:
            coor_j[1] = latt.bound( coor_i[1]+dy, latt.n[1]  )
            for dz in jumpz:
                coor_j[2] = latt.bound( coor_i[2]+dz, latt.n[2]  )
                j = latt.index( coor_j )
                m = np.array( Hopping[(dx, dy, dz)] )
                Tmatrix[ 6*i:6*(i+1), 6*j:6*(j+1) ] += m[0:6, 0:6]

Tmatrix[ 6*latt.L:12*latt.L, 6*latt.L:12*latt.L ] = Tmatrix[ 0:6*latt.L, 0:6*latt.L ]

f = open("model_param", 'w')
f.write( '{:16d} \n'.format(6*latt.L) )
f.write( '{:16d} \n'.format(Ntot) )
f.write( '{:16d} \n'.format(numberkana) )
f.write( '{:>16} \n'.format(decompTypeJ) )
for i in range( 12*latt.L ):
    for j in range( 12*latt.L ):
        f.write( '{:26.18e} {:26.18e} \n'.format(Tmatrix[j, i].real, Tmatrix[j, i].imag))
for i in range( 6*latt.L ):
    f.write( '{:26.18e} \n'.format( U_one ) )
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
#     Kmatrix[ki+12*latt.L, ki+12*latt.L] = Kmatrix[ki+12*latt.L, ki+12*latt.L] - J[i]*0.5
#
#     kj=site_j[i]
#     Kmatrix[kj, kj] = Kmatrix[kj, kj] - J[i]*0.5
#     Kmatrix[kj+12*latt.L, kj+12*latt.L] = Kmatrix[kj+12*latt.L, kj+12*latt.L] - J[i]*0.5
#
# print "Difference between Kmatrix and Tmatrix: "
# print (Kmatrix-Tmatrix).diagonal()
#
# w, v = np.linalg.eigh(Kmatrix)
# print "Eigen of Kmatrix"
# print(w)
#
# print "FE total energy"
# print( "{:<26.18f}".format( np.sum( w[0:Ntot] ) ) )
