import subprocess
import sys
import os

if len(sys.argv) < 2:
    sys.exit("Missing arguments!!! Example: python dataAnalyze.py blockSize")
blockSize  = int( sys.argv[1] )

#Read Lattice Size
f = open("model_param", 'r')
firstLine =  f.readline()
L =  int(firstLine)
secondLine =  f.readline()
thirdLine =  f.readline()
numberOfKana = int( thirdLine )
f.close()

if os.path.isfile('./TNum.dat'):
    print( "\033[1m" "Calculate TAverage." "\033[0m" )
    subprocess.call( os.environ['AFQMCLAB_DIR']+"/bin/NumDenErrorAnalysis TNum.dat den.dat TAverage.dat", shell=True)

if os.path.isfile('./UNum.dat'):
    print( "\033[1m" "Calculate UAverage." "\033[0m" )
    subprocess.call( os.environ['AFQMCLAB_DIR']+"/bin/NumDenErrorAnalysis UNum.dat den.dat UAverage.dat", shell=True)

if os.path.isfile('./U1Num.dat'):
    print( "\033[1m" "Calculate U1Average." "\033[0m" )
    subprocess.call( os.environ['AFQMCLAB_DIR']+"/bin/NumDenErrorAnalysis U1Num.dat den.dat U1Average.dat", shell=True)

if os.path.isfile('./U2Num.dat'):
    print( "\033[1m" "Calculate U2Average." "\033[0m" )
    subprocess.call( os.environ['AFQMCLAB_DIR']+"/bin/NumDenErrorAnalysis U2Num.dat den.dat U2Average.dat", shell=True)

if os.path.isfile('./JNum.dat'):
    print( "\033[1m" "Calculate JAverage." "\033[0m" )
    subprocess.call( os.environ['AFQMCLAB_DIR']+"/bin/NumDenErrorAnalysis JNum.dat den.dat JAverage.dat", shell=True)

if os.path.isfile('./HNum.dat'):
    print( "\033[1m" "Calculate HAverage." "\033[0m" )
    subprocess.call( os.environ['AFQMCLAB_DIR']+"/bin/NumDenErrorAnalysis HNum.dat den.dat HAverage.dat", shell=True)

if os.path.isfile('./NupTotNum.dat'):
    print( "\033[1m" "Calculate NupTotAverage." "\033[0m" )
    subprocess.call( os.environ['AFQMCLAB_DIR']+"/bin/NumDenErrorAnalysis NupTotNum.dat den.dat NupTotAverage.dat", shell=True)

if os.path.isfile('./NdnTotNum.dat'):
    print( "\033[1m" "Calculate NdnTotAverage." "\033[0m" )
    subprocess.call( os.environ['AFQMCLAB_DIR']+"/bin/NumDenErrorAnalysis NdnTotNum.dat den.dat NdnTotAverage.dat",shell=True)

if os.path.isfile('./SplusTotNum.dat'):
    print( "\033[1m" "Calculate SplusTotAverage." "\033[0m" )
    subprocess.call( os.environ['AFQMCLAB_DIR']+"/bin/NumDenErrorAnalysis SplusTotNum.dat den.dat SplusTotAverage.dat", shell=True)

if os.path.isfile('./SminusTotNum.dat'):
    print( "\033[1m" "Calculate SminusTotAverage." "\033[0m" )
    subprocess.call( os.environ['AFQMCLAB_DIR']+"/bin/NumDenErrorAnalysis SminusTotNum.dat den.dat SminusTotAverage.dat", shell=True)

if os.path.isfile('./NupNum.dat'):
    print( "\033[1m" "Calculate NupAverage." "\033[0m" )
    subprocess.call( os.environ['AFQMCLAB_DIR']+"/bin/NumArrayDenErrorAnalysis NupNum.dat den.dat Nup_{1:d}_Average.dat {0:d} {1:d}".format(L, blockSize),shell=True)

if os.path.isfile('./NdnNum.dat'):
    print( "\033[1m" "Calculate NdnAverage." "\033[0m" )
    subprocess.call( os.environ['AFQMCLAB_DIR']+"/bin/NumArrayDenErrorAnalysis NdnNum.dat den.dat Ndn_{1:d}_Average.dat {0:d} {1:d}".format(L, blockSize),shell=True)

if os.path.isfile('./SplusNum.dat'):
    print( "\033[1m" "Calculate SplusAverage." "\033[0m" )
    subprocess.call( os.environ['AFQMCLAB_DIR']+"/bin/NumArrayDenErrorAnalysis SplusNum.dat den.dat Splus_{1:d}_Average.dat {0:d} {1:d}".format(L, blockSize),shell=True)

if os.path.isfile('./SminusNum.dat'):
    print( "\033[1m" "Calculate SminusAverage." "\033[0m" )
    subprocess.call( os.environ['AFQMCLAB_DIR']+"/bin/NumArrayDenErrorAnalysis SminusNum.dat den.dat Sminus_{1:d}_Average.dat {0:d} {1:d}".format(L, blockSize),shell=True)

if os.path.isfile('./kanamoriBgNum.dat'):
    print( "\033[1m" "Calculate kanamoriBgAverage." "\033[0m" )
    subprocess.call( os.environ['AFQMCLAB_DIR']+"/bin/NumArrayDenErrorAnalysis kanamoriBgNum.dat den.dat kanamoriBg_{1:d}_Average.dat {0:d} {1:d}".format(numberOfKana, blockSize),shell=True)