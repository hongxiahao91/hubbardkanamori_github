//
// Created by boruoshihao on 6/23/17.
//

#ifndef AFQMCLAB_AFQMCPHASELESSDEFINE_H
#define AFQMCLAB_AFQMCPHASELESSDEFINE_H

#include "afqmclab.h"

typedef Hop OneBody;

//typedef KanamoriInteractReal       TwoBody;
//typedef KanamoriInteractRealAux    TwoBodyAux;
//typedef KanamoriInteractRealForce  TwoBodyForce;
//typedef KanamoriInteractRealSample TwoBodySample;
//
//typedef HubbardKanamoriReal Model;
//typedef HubbardKanamoriRealMeasureFixSDSD ModelMeasureMixed;

typedef KanamoriInteract       TwoBody;
typedef KanamoriInteractAux    TwoBodyAux;
typedef KanamoriInteractForce  TwoBodyForce;
typedef KanamoriInteractSample TwoBodySample;

typedef HubbardKanamori Model;
typedef HubbardKanamoriMeasureFixSDSD ModelMeasureMixed;

//typedef SD2s   WalkerLeft;
//typedef SD2is  WalkerRight;
//typedef Hop2isSD2isOperation OneBodyWalkerRightOperation;
//typedef CholeskyRealSampleSD2isOperation TwoBodySampleWalkerRightOperation;
//typedef SD2sSD2isOperation WalkerWalkerOperation;
//typedef RealMaterialMoleculeMeasureFixedSD2sSD2is ModelMeasureMixed;

//typedef MDCas2s WalkerLeft;
//typedef SD2is  WalkerRight;
//typedef Hop2isSD2isOperation OneBodyWalkerRightOperation;
//typedef CholeskyRealSampleSD2isOperation TwoBodySampleWalkerRightOperation;
//typedef MDCas2sSD2isOperation WalkerWalkerOperation;
//typedef RealMaterialMoleculeMeasureFixedMDCas2sSD2is ModelMeasureMixed;

typedef SD WalkerLeft;
typedef SD  WalkerRight;
typedef HopSDOperation OneBodyWalkerRightOperation;
typedef LogHopSDOperation TwoBodySampleWalkerRightOperation;
typedef SDSDOperation WalkerWalkerOperation;


#endif //AFQMCLAB_AFQMCPHASELESSDEFINE_H
