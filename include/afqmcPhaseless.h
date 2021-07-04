//
// Created by boruoshihao on 7/8/17.
//

#ifndef AFQMCLAB_AFQMCPHASELESS_H
#define AFQMCLAB_AFQMCPHASELESS_H

#include "afqmcPhaselessDefine.h"
#include "afqmcPhaselessMethod.h"

class AfqmcPhaseless
{
 private:
    AfqmcPhaselessMethod method;
    Model model;
    OneBody expMinusDtK, expMinusHalfDtK, expHalfDtK;
    TwoBody expMinusDtV;
    TwoBodyForce dynamicForce, constForce;

    OneBodyWalkerRightOperation oneBodyWalkerRightOperation;
    TwoBodySampleWalkerRightOperation twoBodySampleWalkerRightOperation;

    TwoBodyAux twoBodyAux;
    TwoBodySample twoBodySample;

    WalkerLeft phiT;
    std::vector<WalkerRight> walker;
    std::vector<bool> walkerIsAlive;
    std::vector<std::complex<double>> logOverlapBackup;

    WalkerWalkerOperation walkerWalkerOperation;
    ModelMeasureMixed mixedMeasure;

 public:
    AfqmcPhaseless();
    ~AfqmcPhaseless();

    void run();
    void initialParameters();
    void initialMeasure();
    void initialExpOneBody();
    void estimateMemory();
    void measureWithoutProjection();
    void measureWithProjection();
    void prepareStop();

 private:
    void initialPhiT();
    void initialWalker();
    void writeWalkers();
    void checkOverlap(WalkerRight &oneWalker);
    void initialMgsAndPopControl();

    void projectExpHalfDtK();
    void projectExpMinusHalfDtK();
    void projectExpMinusDtKExpMinusDtV();
    void projectOneStep(size_t &mgsIndex, size_t &popControlIndex);
    void modifyGM();
    void popControl(double popWeightCap=0.2);
    void checkAndResetWalkerIsAlive();
    void resetLogOverlapBackup();

    void addPhiTMeasure();
    void addPureMeasure();
    void addMixedMeasurement();
    void writeAndResetMeasurement();
    void adjustETThenResetMeasurement();
    void adjustETAndBackGroundThenResetMeasurement();
};

#endif //AFQMCLAB_AFQMCPHASELESS_H