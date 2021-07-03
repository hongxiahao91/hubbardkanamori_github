//
// Created by boruoshihao on 6/23/17.
//

#include "../include/afqmcPhaselessMethod.h"

using namespace std;

AfqmcPhaselessMethod::AfqmcPhaselessMethod()
{
    setDefault();
}

AfqmcPhaselessMethod::~AfqmcPhaselessMethod()
{

}


void AfqmcPhaselessMethod::read(const std::string &filename)
{
    readBySearchString(dt, "dt", filename);
    readBySearchString(thermalSize, "thermalSize", filename);
    readBySearchString(writeNumber, "writeNumber", filename);
    readBySearchString(measureNumberPerWrite, "measureNumberPerWrite", filename);
    readBySearchString(measureSkipStep, "measureSkipStep", filename);

    readBySearchString(walkerSizePerThread, "walkerSizePerThread", filename);
    walkerSize = walkerSizePerThread*MPISize();

    readBySearchString(decompTypeU, "decompTypeU", filename);
    readBySearchString(decompTypeU1, "decompTypeU1", filename);
    readBySearchString(decompTypeU2J, "decompTypeU2J", filename);


    readBySearchString(forceType, "forceType", filename);
    readBySearchString(forceCap, "forceCap", filename);
    readBySearchString(initialPhiTFlag, "initialPhiTFlag", filename);
    readBySearchString(initialWalkerFlag, "initialWalkerFlag", filename);

    readBySearchString(mgsStep, "mgsStep", filename);
    readBySearchString(popControlStep, "popControlStep", filename);
    readBySearchString(initPopControlMaxSize, "initPopControlMaxSize", filename);

    readBySearchString(logEnergyCap, "logEnergyCap", filename);

    readBySearchString(ET, "ET", filename);
    readBySearchString(backGroundETInit, "backGroundETInit", filename);
    readBySearchString(ETAdjustStep, "ETAdjustStep", filename);
    readBySearchString(ETAdjustMaxSize, "ETAdjustMaxSize", filename);
    readBySearchString(ETAndBackGroundGrowthEstimateStep, "ETAndBackGroundGrowthEstimateStep", filename);
    readBySearchString(ETAndBackGroundGrowthEstimateMaxSize, "ETAndBackGroundGrowthEstimateMaxSize", filename);

    readBySearchString(seed, "seed", filename);

    analysis();
}

void AfqmcPhaselessMethod::write(const std::string &filename)
{
    ofstream file;
    file.open(filename, ios::out|ios::trunc);
    if ( ! file.is_open() ) {cout << "Error opening file in File!!! "<<filename<<endl; exit(1);}

    file<<left<<endl;

    file<<setw(36)<<"dt "<<setw(26)<<dt<<endl;
    file<<setw(36)<<"thermalSize "<<setw(26)<<thermalSize<<endl;
    file<<setw(36)<<"writeNumber "<<setw(26)<<writeNumber<<endl;
    file<<setw(36)<<"measureNumberPerWrite "<<setw(26)<<measureNumberPerWrite<<endl;
    file<<setw(36)<<"measureSkipStep "<<setw(26)<<measureSkipStep<<endl;
    file<<endl;

    file<<setw(36)<<"walkerSizePerThread "<<setw(26)<<walkerSizePerThread<<endl;
    file<<endl;

    file<<setw(36)<<"decompTypeU "<<setw(26)<<decompTypeU<<endl;
    file<<setw(36)<<"decompTypeU1 "<<setw(26)<<decompTypeU1<<endl;
    file<<setw(36)<<"decompTypeU2J "<<setw(26)<<decompTypeU2J<<endl;

    file<<setw(36)<<"forceType "<<setw(26)<<forceType<<endl;
    file<<setw(36)<<"forceCap "<<setw(26)<<forceCap<<endl;
    file<<setw(36)<<"initialPhiTFlag "<<setw(26)<<initialPhiTFlag<<endl;
    file<<setw(36)<<"initialWalkerFlag "<<setw(26)<<initialWalkerFlag<<endl;
    file<<endl;

    file<<setw(36)<<"mgsStep "<<setw(26)<<mgsStep<<endl;
    file<<setw(36)<<"popControlStep "<<setw(26)<<popControlStep<<endl;
    file<<setw(36)<<"initPopControlMaxSize "<<setw(26)<<initPopControlMaxSize<<endl;
    file<<endl;

    file<<setw(36)<<"logEnergyCap "<<setw(26)<<logEnergyCap<<endl;
    file<<endl;

    file<<setw(36)<<"ET "<<setw(26)<<ET<<endl;
    file<<setw(36)<<"backGroundETInit "<<setw(26)<<backGroundETInit<<endl;
    file<<setw(36)<<"ETAdjustStep "<<setw(26)<<ETAdjustStep<<endl;
    file<<setw(36)<<"ETAdjustMaxSize "<<setw(26)<<ETAdjustMaxSize<<endl;
    file<<setw(36)<<"ETAndBackGroundGrowthEstimateStep "<<setw(26)<<ETAndBackGroundGrowthEstimateStep<<endl;
    file<<setw(36)<<"ETAndBackGroundGrowthEstimateMaxSize "<<setw(26)<<ETAndBackGroundGrowthEstimateMaxSize<<endl;
    file<<endl;

    file<<setw(36)<<"seed "<<setw(26)<<seed<<endl;
    file<<endl;

    file.close();
}

void AfqmcPhaselessMethod::print()
{
    cout<<left<<endl;

    cout<<setw(36)<<"AFQMC parameters: \n"<<endl;

    cout<<setw(36)<<"dt "<<setw(26)<<dt<<endl;
    cout<<setw(36)<<"thermalSize "<<setw(26)<<thermalSize<<endl;
    cout<<setw(36)<<"writeNumber "<<setw(26)<<writeNumber<<endl;
    cout<<setw(36)<<"measureNumberPerWrite "<<setw(26)<<measureNumberPerWrite<<endl;
    cout<<setw(36)<<"measureSkipStep "<<setw(26)<<measureSkipStep<<endl;
    cout<<endl;

    cout<<setw(36)<<"walkerSizePerThread "<<setw(26)<<walkerSizePerThread<<endl;
    cout<<setw(36)<<"walkerSize "<<setw(26)<<walkerSize<<endl;
    cout<<endl;

    cout<<setw(36)<<"decompTypeU "<<setw(26)<<decompTypeU<<endl;
    cout<<setw(36)<<"decompTypeU1 "<<setw(26)<<decompTypeU1<<endl;
    cout<<setw(36)<<"decompTypeU2J "<<setw(26)<<decompTypeU2J<<endl;

    cout<<setw(36)<<"forceType "<<setw(26)<<forceType<<endl;
    cout<<setw(36)<<"forceCap "<<setw(26)<<forceCap<<endl;
    cout<<setw(36)<<"initialPhiTFlag "<<setw(26)<<initialPhiTFlag<<endl;
    cout<<setw(36)<<"initialWalkerFlag "<<setw(26)<<initialWalkerFlag<<endl;
    cout<<endl;

    cout<<setw(36)<<"mgsStep "<<setw(26)<<mgsStep<<endl;
    cout<<setw(36)<<"popControlStep "<<setw(26)<<popControlStep<<endl;
    cout<<setw(36)<<"initPopControlMaxSize "<<setw(26)<<initPopControlMaxSize<<endl;
    cout<<endl;

    cout<<setw(36)<<"logEnergyCap "<<setw(26)<<logEnergyCap<<endl;
    cout<<endl;

    cout<<setw(36)<<"ET "<<setw(26)<<ET<<endl;
    cout<<setw(36)<<"backGroundETInit "<<setw(26)<<backGroundETInit<<endl;
    cout<<setw(36)<<"ETAdjustStep "<<setw(26)<<ETAdjustStep<<endl;
    cout<<setw(36)<<"ETAdjustMaxSize "<<setw(26)<<ETAdjustMaxSize<<endl;
    cout<<setw(36)<<"ETAndBackGroundGrowthEstimateStep "<<setw(26)<<ETAndBackGroundGrowthEstimateStep<<endl;
    cout<<setw(36)<<"ETAndBackGroundGrowthEstimateMaxSize "<<setw(26)<<ETAndBackGroundGrowthEstimateMaxSize<<endl;
    cout<<endl;

    cout<<setw(36)<<"seed "<<setw(26)<<seed<<endl;
    cout<<endl;
}

#ifdef MPI_HAO
void MPIBcast(AfqmcPhaselessMethod &buffer, int root, MPI_Comm const &comm)
{
    MPIBcast(buffer.dt, root, comm);
    MPIBcast(buffer.thermalSize, root, comm);
    MPIBcast(buffer.writeNumber, root, comm);
    MPIBcast(buffer.measureNumberPerWrite, root, comm);
    MPIBcast(buffer.measureSkipStep, root, comm);

    MPIBcast(buffer.walkerSizePerThread, root, comm);
    MPIBcast(buffer.walkerSize, root, comm);

    MPIBcast(buffer.decompTypeU, root, comm);
    MPIBcast(buffer.decompTypeU1, root, comm);
    MPIBcast(buffer.decompTypeU2J, root, comm);

    MPIBcast(buffer.forceType, root, comm);
    MPIBcast(buffer.forceCap, root, comm);
    MPIBcast(buffer.initialPhiTFlag, root, comm);
    MPIBcast(buffer.initialWalkerFlag, root, comm);

    MPIBcast(buffer.mgsStep, root, comm);
    MPIBcast(buffer.popControlStep, root, comm);
    MPIBcast(buffer.initPopControlMaxSize, root, comm);

    MPIBcast(buffer.logEnergyCap, root, comm);

    MPIBcast(buffer.ET, root, comm);
    MPIBcast(buffer.backGroundETInit, root, comm);
    MPIBcast(buffer.ETAdjustStep, root, comm);
    MPIBcast(buffer.ETAdjustMaxSize, root, comm);
    MPIBcast(buffer.ETAndBackGroundGrowthEstimateStep, root, comm);
    MPIBcast(buffer.ETAndBackGroundGrowthEstimateMaxSize, root, comm);

    MPIBcast(buffer.seed, root, comm);
}
#endif

void AfqmcPhaselessMethod::setDefault()
{
    dt=0.01;
    thermalSize = 200;
    writeNumber = 60;
    measureNumberPerWrite = 2;
    measureSkipStep = 5;

    walkerSizePerThread = 300 / MPISize();
    walkerSize = walkerSizePerThread * MPISize();

    decompTypeU = "None";
    decompTypeU1 = "None";
    decompTypeU2J = "None";

    forceType="dynamicForce";
    forceCap = 1.5;
    initialPhiTFlag = "readFromFile";
    initialWalkerFlag = "readFromFile";

    mgsStep = 10;
    popControlStep = 10;
    initPopControlMaxSize = 0;

    logEnergyCap = 23.0*sqrt(dt);

    ET = -10;
    backGroundETInit="EstimateFromPhiTWalker";
    ETAdjustStep = 10;
    ETAdjustMaxSize = 100;
    ETAndBackGroundGrowthEstimateStep = 5;
    ETAndBackGroundGrowthEstimateMaxSize = 200;

    seed = 985456376;
}

void AfqmcPhaselessMethod::analysis()
{
    if( std::abs(dt) < 1e-12 )
    {
        cout<<"The code will measure everything without projection!"<<endl;
        return;
    }

    if( dt < 0.0 )
    {
        cout<<"Error!!! dt must be postive!"<<endl;
        exit(1);
    }

    if( writeNumber==0 )
    {
        cout<<"Warning!!! writeNumber = 0, code will not write any measurement to disk!"<<endl;
    }

    if( measureNumberPerWrite == 0 )
    {
        cout<<"Error!!! measureNumberPerWrite = 0, code will not measure anything!"<<endl;
        exit(1);
    }

    if( measureSkipStep ==0 )
    {
        cout<<"Error!!! measureSkipStep ==0, code will not measure anything!"<<endl;
        exit(1);
    }

    if( walkerSizePerThread <= 0 )
    {
        cout<<"Error!!! walkerSizePerThread is not a positive number!"<<endl;
        exit(1);
    }

    if( forceCap < 0.0 )
    {
        cout<<"Error!!! forceCap must be positive or zero!"<<endl;
        exit(1);
    }

    if( mgsStep==0 )
    {
        cout<<"Warning!!! The code will not do mgs!"<<endl;
    }

    if( popControlStep==0 )
    {
        cout<<"Warning!!! The code will not do population control!"<<endl;
    }

    if( initPopControlMaxSize>thermalSize )
    {
        cout<<"Warning!!! initPopControlMaxSize is larger than thermalSize, only works up to thermalSize!"<<endl;
    }

    if( ETAdjustMaxSize > thermalSize )
    {
        cout << "Error!!! We should not adjust ET and backGround after thermalizing!" << endl;
        exit(1);
    }

    if( ETAndBackGroundGrowthEstimateMaxSize > thermalSize )
    {
        cout<<"Error!!! We should not use growth estimator to adjust ET and backGround after thermalizing!"<<endl;
        exit(1);
    }

    if( ETAdjustMaxSize >= ETAndBackGroundGrowthEstimateMaxSize )
    {
        cout<<"ETAndBackGroundGrowthEstimator is off!"<<endl;
    }
    else
    {
        size_t numberOfMeasure=(ETAndBackGroundGrowthEstimateMaxSize-ETAdjustMaxSize)/ETAndBackGroundGrowthEstimateStep;
        cout<<"Number Of ETAndBackGroundGrowthEstimation is "<<numberOfMeasure<<endl;

        if( numberOfMeasure == 0 ) cout<<"ETAndBackGroundGrowthEstimator is off!"<<endl;
        else cout<<"ETAndBackGroundGrowthEstimator is on!"<<endl;
    }

    if( logEnergyCap<0.0 )
    {
        cout<<"Error!!! logEnergyCap can not be negative!"<<endl;
        exit(1);
    }
}