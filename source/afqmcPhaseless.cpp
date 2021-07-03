//
// Created by boruoshihao on 7/8/17.
//

#include "../include/afqmcPhaseless.h"

using namespace std;
using namespace tensor_hao;

AfqmcPhaseless::AfqmcPhaseless() { }

AfqmcPhaseless::~AfqmcPhaseless() { }

void AfqmcPhaseless::run()
{
    initialParameters();

    initialPhiT();
    
    initialWalker();

    initialMeasure();

    if( std::abs(method.dt) < 1e-12 )
    {
        estimateMemory();

        measureWithoutProjection();
    }
    else
    {
        initialExpOneBody();

        estimateMemory();

        measureWithProjection();
    }

    prepareStop();
}

void AfqmcPhaseless::initialParameters()
{
    if( MPIRank()==0 ) method.read("afqmc_param");
    if( MPIRank()==0 ) method.print();
    MPIBcast(method);

    randomHaoInit(method.seed, 1);
    if( method.seed != -1 ) randomHaoSave();

    if( MPIRank()==0 ) model.read("model_param");
    MPIBcast(model);

    expMinusDtV = model.returnExpMinusAlphaV( method.dt, method.decompTypeU, method.decompTypeU1, method.decompTypeU2J );
    constForce = expMinusDtV.readForce("constForce_param");
}

void AfqmcPhaseless::initialMeasure()
{
    mixedMeasure.setModelWalker(model, phiT);
}

void AfqmcPhaseless::initialExpOneBody()
{
    //Init backGround to init exp( OneBoy )
    if( method.backGroundETInit == "readFromFile" )
    {
        //Already read from the model.
    }
    else if( method.backGroundETInit == "EstimateFromPhiTWalker")
    {
        addPureMeasure();
        adjustETAndBackGroundThenResetMeasurement();
    }
    else if( method.backGroundETInit == "EstimateFromPhiT")
    {
        addPhiTMeasure();
        adjustETAndBackGroundThenResetMeasurement();
    }
    else
    {
        cout<<"Error!!! do not know the type of backGroundETInit!"<<endl;
        exit(1);
    }

    expMinusDtK     = model.returnExpMinusAlphaK(  method.dt     );
    expMinusHalfDtK = model.returnExpMinusAlphaK(  method.dt*0.5 );
    expHalfDtK      = model.returnExpMinusAlphaK( -method.dt*0.5 );
}

void AfqmcPhaseless::estimateMemory()
{
    double mem(0.0);
    mem += model.getMemory();
    mem += expMinusHalfDtK.getMemory()+expHalfDtK.getMemory()+expMinusDtK.getMemory();
    mem += expMinusDtV.getMemory();
    mem += constForce.getMemory()*2.0;

    twoBodyAux = expMinusDtV.sampleAuxFromForce(constForce);
    mem += twoBodyAux.getMemory();

    twoBodySample = expMinusDtV.getTwoBodySampleFromAuxForce(twoBodyAux, constForce);
    mem += twoBodySample.getMemory();

    mem += phiT.getMemory();
    mem += ( walker[0].getMemory()+1.0 ) * method.walkerSizePerThread;

    walkerWalkerOperation.set(phiT, walker[0]);
    mixedMeasure.addMeasurement(walkerWalkerOperation, 1.0);
    mem += walkerWalkerOperation.getMemory();
    mem += mixedMeasure.getMemory();
    mixedMeasure.reSet();

    //Make a slightly big estimation for uncounted memory.
    mem*=1.2;
    if(MPIRank()==0)
    {
        cout<<"Memory need for this program is roughly: "<<mem/1e9<<"G per process."<<endl;
        cout<<"Please make sure available memory is larger than this.\n"<<endl;
    }
}

void AfqmcPhaseless::measureWithoutProjection()
{
    if( MPIRank() == 0 ) cout<<"Measure without projection."<<endl;

    addPureMeasure();
    writeAndResetMeasurement();
}

void AfqmcPhaseless::measureWithProjection()
{
    if( MPIRank() == 0 ) cout<<"Start the projection..."<<endl;

    double beta;

    beta = (method.thermalSize+method.writeNumber*method.measureNumberPerWrite*method.measureSkipStep)*method.dt;
    if( MPIRank() == 0 ) cout<<"Total beta will be "<<beta<<endl;

    projectExpMinusHalfDtK();

    size_t mgsIndex(0), popControlIndex(0);

    if( MPIRank() == 0 ) cout<<"\nThermalize..."<<endl;

    size_t numberOfGrowthMeasure = (method.ETAndBackGroundGrowthEstimateMaxSize-method.ETAdjustMaxSize)/method.ETAndBackGroundGrowthEstimateStep;

    for (size_t i = 0; i < method.thermalSize; ++i)
    {
        if ( i<method.ETAdjustMaxSize )
        {
            if ( (i+1)%method.ETAdjustStep == 0 )
            {
                addMixedMeasurement();
                adjustETThenResetMeasurement();
            }
        }

        if (i < method.ETAndBackGroundGrowthEstimateMaxSize && i >= method.ETAdjustMaxSize)
        {
            if ( (i + 1 - method.ETAdjustMaxSize) % method.ETAndBackGroundGrowthEstimateStep == 0 )
            {
                addMixedMeasurement();
            }

            if ( i == (method.ETAndBackGroundGrowthEstimateMaxSize - 1) )
            {
                if(numberOfGrowthMeasure>0)
                {
                    adjustETAndBackGroundThenResetMeasurement();

                    projectExpHalfDtK();
                    expMinusDtK     = model.returnExpMinusAlphaK(  method.dt     );
                    expMinusHalfDtK = model.returnExpMinusAlphaK(  method.dt*0.5 );
                    expHalfDtK      = model.returnExpMinusAlphaK( -method.dt*0.5 );
                    projectExpMinusHalfDtK();
                }
            }
        }

        if( MPIRank() == 0 ) cout<<i*method.dt<<endl;

        if( i<method.initPopControlMaxSize )
        {
            if( method.popControlStep>0 ) popControlIndex=method.popControlStep-1;
            else popControlIndex=0;
        }

        projectOneStep(mgsIndex, popControlIndex);
    }

    if( MPIRank() == 0 ) cout<<"\nMeasure..."<<endl;

    for (size_t i = 0; i < method.writeNumber; ++i)
    {
        for (size_t j = 0; j < method.measureNumberPerWrite; ++j)
        {
            addMixedMeasurement();

            for (size_t k = 0; k < method.measureSkipStep; ++k)
            {
                beta = ( method.thermalSize+k+j*method.measureSkipStep
                         +i*method.measureSkipStep*method.measureNumberPerWrite)*method.dt;
                if (MPIRank() == 0) cout << beta << endl;
                projectOneStep(mgsIndex, popControlIndex);
            }
        }

        beta = ( method.thermalSize+(i+0.5)*method.measureNumberPerWrite*method.measureSkipStep-0.5 )*method.dt;
        if (MPIRank() == 0) writeFile( beta, "beta.dat", ios::app);
        writeAndResetMeasurement();
    }

    projectExpHalfDtK();
}

void AfqmcPhaseless::prepareStop()
{
    if( MPIRank()==0 ) method.write("afqmc_param");
    if( MPIRank()==0 ) model.writeBackGround("model_param");
    writeWalkers();
    randomHaoSave();

    if( std::abs(method.dt) >= 1e-12 )
    {
        if( MPIRank()==0 ) cout<<"twoBodySampleWalkerRightOperation information: "<<endl;
        twoBodySampleWalkerRightOperation.print();
    }
}
