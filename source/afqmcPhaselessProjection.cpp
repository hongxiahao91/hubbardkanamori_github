//
// Created by boruoshihao on 7/8/17.
//

#include "../include/afqmcPhaseless.h"
#include "../include/afqmcWalkerPop.h"

using namespace std;
using namespace tensor_hao;

void AfqmcPhaseless::projectExpHalfDtK()
{
    WalkerRight walkerTemp;
    for(int i = 0; i < method.walkerSizePerThread; ++i)
    {
        if( walkerIsAlive[i] )
        {
            oneBodyWalkerRightOperation.applyToRight(expHalfDtK, walker[i], walkerTemp);
            walker[i] = move(walkerTemp);
        }
    }
}

void AfqmcPhaseless::projectExpMinusHalfDtK()
{
    WalkerRight walkerTemp;
    for(int i = 0; i < method.walkerSizePerThread; ++i)
    {
        if( walkerIsAlive[i] )
        {
            oneBodyWalkerRightOperation.applyToRight(expMinusHalfDtK, walker[i], walkerTemp);
            walker[i] = move(walkerTemp);
        }
    }
}

void AfqmcPhaseless::projectExpMinusDtKExpMinusDtV()
{
    complex<double> logOverlap;
    double logDiff;
    double phase;
    WalkerRight walkerTemp;
    complex<double> auxForce;
    for(int i = 0; i < method.walkerSizePerThread; ++i)
    {
        if( walkerIsAlive[i] )
        {
            WalkerWalkerOperation walkerWalkerOperation(phiT, walker[i]);
            logOverlap = walkerWalkerOperation.returnLogOverlap();

            //Cap weight real
            logDiff = ( logOverlap - logOverlapBackup[i] ).real();
            if( (logDiff-method.logEnergyCap)>1e-12  )  //如果<phiT|phi_new>/<phiT|phi_old> > EnergyCap
            {
                cout<<"Cap weight in projection: "
                    <<setw(25)<<i+MPIRank()*method.walkerSizePerThread
                    <<setw(25)<<logDiff<<setw(25)<<method.logEnergyCap<<endl;
                walker[i].addLogw( method.logEnergyCap  - logDiff ); //log<phiT|phi_new> = log<phiT|phi_old> + logEngCap, in other words, <phiT|phi_new>/<phiT|phi_old> = EnergyCap
                logOverlap += ( method.logEnergyCap - logDiff );
            }

            //Force weight to be real
            phase = logOverlap.imag();
            //if (phase > 1e-8) cout << "weight: "<< phase <<endl;
            if (cos(phase) <= 0.0) { walkerIsAlive[i] = false; continue; }
            walker[i].addLogw( log(cos(phase)) - complex<double>(0, phase) );

            logOverlapBackup[i] = logOverlap + log(cos(phase)) - complex<double>(0, phase);

            if (method.forceType == "constForce")
            {
                twoBodyAux = expMinusDtV.sampleAuxFromForce(constForce);
                twoBodySample = expMinusDtV.getTwoBodySampleFromAuxForce(twoBodyAux, constForce);
                auxForce = expMinusDtV.calculateAuxForce(twoBodyAux, constForce);
            }
            else if (method.forceType == "dynamicForce")
            {
                dynamicForce = mixedMeasure.getForce(expMinusDtV, walkerWalkerOperation, method.forceCap);
                twoBodyAux = expMinusDtV.sampleAuxFromForce(dynamicForce);
                twoBodySample = expMinusDtV.getTwoBodySampleFromAuxForce(twoBodyAux, dynamicForce);
                auxForce = expMinusDtV.calculateAuxForce(twoBodyAux, dynamicForce);
            }
            else
            {
                cout << "Error!!! Do not know method.forceType " << method.forceType << endl;
                exit(1);
            }

            //Cos projection
            phase = auxForce.imag();
            //if (phase > 1e-8) cout << "Phase: "<< phase <<endl;
            if (cos(phase) <= 0.0) { walkerIsAlive[i] = false; continue; }
            walker[i].addLogw( log( cos(phase) ) );

            twoBodySampleWalkerRightOperation.applyToRight(twoBodySample, walker[i], walkerTemp);

            walkerTemp.addLogw(method.dt * method.ET);

            oneBodyWalkerRightOperation.applyToRight(expMinusDtK, walkerTemp, walker[i]);
        }
    }
}

void AfqmcPhaseless::projectOneStep(size_t &mgsIndex, size_t &popControlIndex)
{
    projectExpMinusDtKExpMinusDtV();

    mgsIndex++;
    if (mgsIndex == method.mgsStep)
    {
        modifyGM();
        mgsIndex = 0;
    }

    popControlIndex++;
    if (popControlIndex == method.popControlStep)
    {
        popControl();
        popControlIndex = 0;
    }
}

void AfqmcPhaseless::modifyGM()
{
    for(int i = 0; i < method.walkerSizePerThread; ++i)
    {
        if( walkerIsAlive[i] ) walker[i].stabilize();
    }
}

void AfqmcPhaseless::popControl(double popWeightCap)
{
    vector<double> weightPerThread( method.walkerSizePerThread );
    complex<double> logOverlap; double overlapReal;
    double logDiff;
    for(int i = 0; i < method.walkerSizePerThread; ++i)
    {
        if( walkerIsAlive[i] )
        {
            walkerWalkerOperation.set(phiT, walker[i]);

            logOverlap = walkerWalkerOperation.returnLogOverlap();

            //Cap weight real
            if( logOverlapBackup.size() > 0 )
            {
                logDiff = ( logOverlap - logOverlapBackup[i] ).real();
                if( (logDiff-method.logEnergyCap)>1e-12  )
                {
                    cout<<"Cap weight in projection before popControl: "
                        <<setw(25)<<i+MPIRank()*method.walkerSizePerThread
                        <<setw(25)<<logDiff<<setw(25)<<method.logEnergyCap<<endl;
                    walker[i].addLogw( method.logEnergyCap  - logDiff );
                    logOverlap += ( method.logEnergyCap - logDiff );
                }
            }

            overlapReal = (exp(logOverlap)).real();
            if (overlapReal <= 0.0) { weightPerThread[i] = 0.0; walkerIsAlive[i] = false; }
            else { weightPerThread[i] = overlapReal; }

            walker[i].addLogw( -logOverlap );
        }
        else
        {
            weightPerThread[i] = 0.0;
        }
    }

    checkAndResetWalkerIsAlive();
    resetLogOverlapBackup();

    vector<double> weight;
#ifdef MPI_HAO
    if( MPIRank()==0 ) weight.resize( method.walkerSize );
    MPI_Gather( weightPerThread.data(), method.walkerSizePerThread, MPI_DOUBLE_PRECISION,
                weight.data(), method.walkerSizePerThread, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD );
#else
    weight = move( weightPerThread );
#endif


    if( MPIRank() == 0 ) popCheck(weight);
    if( MPIRank() == 0 ) capWeight(weight, popWeightCap);

    vector<int> table;
    if( MPIRank()==0 ) table=popConfiguration( MPISize(), weight );

    vector<AfqmcWalkerPop> walkerPop;
    walkerPop.reserve( method.walkerSizePerThread );
    for(int i=0; i<method.walkerSizePerThread; i++) walkerPop.push_back( AfqmcWalkerPop(walker[i]) );
    populationControl(walkerPop, table);

    if( MPIRank()==0 ) cout<<endl;
}

void AfqmcPhaseless::checkAndResetWalkerIsAlive()
{
    size_t aliveWalkerPerThread(0);
    for (int i = 0; i < method.walkerSizePerThread; ++i)
    {
        if( walkerIsAlive[i] ) aliveWalkerPerThread++;
    }

    size_t aliveWalker = MPISum(aliveWalkerPerThread);

    if( MPIRank()==0 )
    {
        cout<<"Total number of walker is "<<method.walkerSize<<"."<<endl;
        cout<<"Currently "<<aliveWalker<<" walkers are still alive."<<endl;
        cout<<"Currently "<<method.walkerSize-aliveWalker<<" walkers are killed."<<endl;

        if( aliveWalker == 0) { cout<<"Error!!! All walkers are killed!"<<endl; exit(1); }
    }

    for (int i = 0; i < method.walkerSizePerThread ; ++i)
    {
        walkerIsAlive[i] = true;
    }
}

void AfqmcPhaseless::resetLogOverlapBackup()
{
    logOverlapBackup.resize(method.walkerSizePerThread);
    for (int i = 0; i < method.walkerSizePerThread ; ++i)
    {
        logOverlapBackup[i] = 0.0;
    }
}
