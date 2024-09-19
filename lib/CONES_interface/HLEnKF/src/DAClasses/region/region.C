/*--------------------------------*- C++ -*----------------------------------*\
|                         __________  _   _____________                       |
|                        / ____/ __ \/ | / / ____/ ___/                       |
|                       / /   / / / /  |/ / __/  \__ \                        |
|                      / /___/ /_/ / /|  / /___ ___/ /                        |
|                      \____/\____/_/ |_/_____//____/                         |
|                                                                             |
| C.  : Coupling               |    C.O.N.Es. version : 2405                  |
| O.  : OpenFOAM (with)        |    OpenFOAM version  : OpenFOAM-org v9       |
| N.  : Numerical              |    Lucas Villanueva   /   Miguel M.Valero    |
| Es. : EnvironmentS           |    Tom Moussie        /   Sarp Er            |
|      ===============         |    Paolo Errante      /   Marcello Meldi     |
\*---------------------------------------------------------------------------*/
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
/**
 * @dir ./lib/CONES_interface/HLEnKF/src/DAClasses/region
 * @brief **Region class declarations and definitions.**
 */

/**
 * @file region.C
 * @authors Lucas Villanueva (lucas.villanueva@ensma.fr)
 * @authors Miguel M. Valero (miguel.martinez_valero@ensam.eu)
 * @authors Tom Moussie (tom.moussie@ensam.eu)
 * @authors Sarp Er (sarp.er@ensam.eu)
 * @authors Paolo Errante (paolo.errante@ensam.eu)
 * @authors Marcello Meldi (marcello.meldi@ensam.eu)
 * @brief Region class definition.
 * @version 0.1
 * @date 2024
 * 
 * @copyright Copyright (c) 2024
 */

#include <set>
#include "region.h"

/*** Region Class ***/
/** Member functions **/

/* Constructors */
region::region(/*Default*/)
{}

region::region(Eigen::ArrayXi clipIDs, Eigen::ArrayXi cellIDs): 
nClips_(clipIDs.rows()),
clipIDs_(clipIDs),
nCells_(cellIDs.rows()),
cellIDs_(cellIDs),
nProbes_(clipIDs.rows()),
probeIDs_(clipIDs)
{}

region::region(Eigen::ArrayXi clipIDs, Eigen::ArrayXi cellIDs, Eigen::ArrayXi probeIDs): 
nClips_(clipIDs.rows()),
clipIDs_(clipIDs),
nCells_(cellIDs.rows()),
cellIDs_(cellIDs),
nProbes_(probeIDs.rows()),
probeIDs_(probeIDs)
{}

/* Setters */

void region::set_nObsProbes(observations obsData){
    /* Count the number of probes of each type and the total number of observations in the region */
    if (obsData.get_obsType() == "U"){
		switch (obsData.get_velocityCaseSwitch()){
			case 1 :
			case 2 :
			case 3 :
                nObsProbesVel_ = nClips_;
                nObsProbesPres_ = 0;
                nObsProbesCf_ = 0;
                nObs_ = nObsProbesVel_;
                obsCoords_ = Eigen::MatrixXd::Zero(nObs_, 3);
                for (int i = 0 ; i < nObsProbesVel_; i++)
                {
                    obsCoords_(i,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                }
        		break;
			case 4 :
			case 5 :
			case 6 :
                nObsProbesVel_ = nClips_;
                nObsProbesPres_ = 0;
                nObsProbesCf_ = 0;
                nObs_ = 2*nObsProbesVel_;
                obsCoords_ = Eigen::MatrixXd::Zero(nObs_, 3);
                for (int i = 0 ; i < nObsProbesVel_; i++)
                {
                    obsCoords_(i,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                    obsCoords_(i+nObsProbesVel_,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i+nObsProbesVel_,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i+nObsProbesVel_,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                }
        		break;
			case 7 :
                nObsProbesVel_ = nClips_;
                nObsProbesPres_ = 0;
                nObsProbesCf_ = 0;
                nObs_ = 3*nObsProbesVel_;
                obsCoords_ = Eigen::MatrixXd::Zero(nObs_, 3);
                for (int i = 0 ; i < nObsProbesVel_; i++)
                {
                    obsCoords_(i,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                    obsCoords_(i+nObsProbesVel_,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i+nObsProbesVel_,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i+nObsProbesVel_,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                    obsCoords_(i+2*nObsProbesVel_,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i+2*nObsProbesVel_,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i+2*nObsProbesVel_,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                }
        		break;
		}
	}
	else if (obsData.get_obsType() == "p"){
        nObsProbesVel_ = 0;
        nObsProbesPres_ = nClips_;
        nObsProbesCf_ = 0;
        nObs_ = nObsProbesPres_;
        obsCoords_ = Eigen::MatrixXd::Zero(nObs_, 3);
        for (int i = 0 ; i < nObsProbesPres_; i++)
        {
            obsCoords_(i,0) = obsData.get_probeCoord(probeIDs_(i))[0];
            obsCoords_(i,1) = obsData.get_probeCoord(probeIDs_(i))[1];
            obsCoords_(i,2) = obsData.get_probeCoord(probeIDs_(i))[2];
        }
	}
	else if (obsData.get_obsType() == "Up"){
		switch (obsData.get_velocityCaseSwitch()){
			case 1 :
			case 2 :
			case 3 :
                for (int i = 0 ; i < nClips_; i++)
				{
                    if (obsData.get_probeVarLabels(clipIDs_(i))!="p") {nObsProbesVel_ = nObsProbesVel_ + 1;}
					else {nObsProbesPres_ = nObsProbesPres_ + 1;}
				}
                nObsProbesCf_ = 0;
                nObs_ = nObsProbesVel_+nObsProbesPres_;
                obsCoords_ = Eigen::MatrixXd::Zero(nObs_, 3);
                for (int i = 0; i < nObsProbesVel_; ++i){
                    obsCoords_(i,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                }
                for (int i = nObsProbesVel_; i < (nObs_); ++i){
                    obsCoords_(i,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                }
				break;
			case 4 :
			case 5 :
			case 6 :
                for (int i = 0 ; i < nClips_; i++)
				{
                    if (obsData.get_probeVarLabels(clipIDs_(i))!="p") {nObsProbesVel_ = nObsProbesVel_ + 1;}
					else {nObsProbesPres_ = nObsProbesPres_ + 1;}
				}
                nObsProbesCf_ = 0;
                nObs_ = 2*nObsProbesVel_+nObsProbesPres_;
                obsCoords_ = Eigen::MatrixXd::Zero(nObs_, 3);
                for (int i = 0; i < nObsProbesVel_; ++i){
                    obsCoords_(i+nObsProbesVel_,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i+nObsProbesVel_,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i+nObsProbesVel_,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                    obsCoords_(i+2*nObsProbesVel_,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i+2*nObsProbesVel_,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i+2*nObsProbesVel_,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                }
                for (int i = 2*nObsProbesVel_; i < (nObs_); ++i){
                    obsCoords_(i,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                }
        		break;
			case 7 :
                for (int i = 0 ; i < nClips_; i++)
				{
                    if (obsData.get_probeVarLabels(clipIDs_(i))!="p") {nObsProbesVel_ = nObsProbesVel_ + 1;}
					else {nObsProbesPres_ = nObsProbesPres_ + 1;}
				}
                nObsProbesCf_ = 0;
                nObs_ = 3*nObsProbesVel_+nObsProbesPres_;
                obsCoords_ = Eigen::MatrixXd::Zero(nObs_, 3);
                for (int i = 0; i < nObsProbesVel_; ++i){
                    obsCoords_(i,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                    obsCoords_(i+nObsProbesVel_,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i+nObsProbesVel_,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i+nObsProbesVel_,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                    obsCoords_(i+2*nObsProbesVel_,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i+2*nObsProbesVel_,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i+2*nObsProbesVel_,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                }
                for (int i = 3*nObsProbesVel_; i < (nObs_); ++i){
                    obsCoords_(i,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                }
        		break;
		}
	}
	else if (obsData.get_obsType() == "UCf"){
		switch (obsData.get_velocityCaseSwitch()){
			case 1 :
			case 2 :
			case 3 :
                for (int i = 0 ; i < nClips_; i++)
				{
                    if (obsData.get_probeVarLabels(clipIDs_(i))!="Cf") {nObsProbesVel_ = nObsProbesVel_ + 1;}
					else {nObsProbesCf_ = nObsProbesCf_ + 1;}
				}
                nObsProbesPres_ = 0;
                nObs_ = nObsProbesVel_+nObsProbesCf_;
                obsCoords_ = Eigen::MatrixXd::Zero(nObs_, 3);
                for (int i = 0; i < nObsProbesVel_; ++i){
                    obsCoords_(i,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                }
                for (int i = nObsProbesVel_; i < (nObs_); ++i){
                    obsCoords_(i,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                }
				break;
			case 4 :
			case 5 :
			case 6 :
                for (int i = 0 ; i < nClips_; i++)
				{
                    if (obsData.get_probeVarLabels(clipIDs_(i))!="Cf") {nObsProbesVel_ = nObsProbesVel_ + 1;}
					else {nObsProbesCf_ = nObsProbesCf_ + 1;}
				}
                nObsProbesPres_ = 0;
                nObs_ = 2*nObsProbesVel_+nObsProbesCf_;
                obsCoords_ = Eigen::MatrixXd::Zero(nObs_, 3);
                for (int i = 0; i < nObsProbesVel_; ++i){
                    obsCoords_(i+nObsProbesVel_,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i+nObsProbesVel_,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i+nObsProbesVel_,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                    obsCoords_(i+2*nObsProbesVel_,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i+2*nObsProbesVel_,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i+2*nObsProbesVel_,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                }
                for (int i = 2*nObsProbesVel_; i < (nObs_); ++i){
                    obsCoords_(i,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                }
        		break;
			case 7 :
                for (int i = 0 ; i < nClips_; i++)
				{
                    if (obsData.get_probeVarLabels(clipIDs_(i))!="Cf") {nObsProbesVel_ = nObsProbesVel_ + 1;}
					else {nObsProbesCf_ = nObsProbesCf_ + 1;}
				}
                nObsProbesPres_ = 0;
                nObs_ = 3*nObsProbesVel_+nObsProbesCf_;
                obsCoords_ = Eigen::MatrixXd::Zero(nObs_, 3);
                for (int i = 0; i < nObsProbesVel_; ++i){
                    obsCoords_(i,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                    obsCoords_(i+nObsProbesVel_,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i+nObsProbesVel_,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i+nObsProbesVel_,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                    obsCoords_(i+2*nObsProbesVel_,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i+2*nObsProbesVel_,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i+2*nObsProbesVel_,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                }
                for (int i = 3*nObsProbesVel_; i < (nObs_); ++i){
                    obsCoords_(i,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                }
        		break;
		}
    }
}

void region::set_nObsProbesNoClip(observations obsData){
    /* Count the number of probes of each type and the total number of observations in the region */
    
    nObsProbesVel_ = obsData.get_nObsProbesVel();
    nObsProbesPres_ = obsData.get_nObsProbesPres();
    nObsProbesCf_ = obsData.get_nObsProbesCf();
    nObs_ = obsData.get_nObs();
    obsCoords_ = Eigen::MatrixXd::Zero(nObs_, 3);

    if (obsData.get_obsType() == "U"){
		switch (obsData.get_velocityCaseSwitch()){
			case 1 :
			case 2 :
			case 3 :
                for (int i = 0 ; i < nObsProbesVel_; i++)
                {
                    obsCoords_(i,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                }
        		break;
			case 4 :
			case 5 :
			case 6 :
                for (int i = 0 ; i < nObsProbesVel_; i++)
                {
                    obsCoords_(i,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                    obsCoords_(i+nObsProbesVel_,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i+nObsProbesVel_,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i+nObsProbesVel_,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                }
        		break;
			case 7 :
                for (int i = 0 ; i < nObsProbesVel_; i++)
                {
                    obsCoords_(i,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                    obsCoords_(i+nObsProbesVel_,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i+nObsProbesVel_,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i+nObsProbesVel_,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                    obsCoords_(i+2*nObsProbesVel_,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i+2*nObsProbesVel_,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i+2*nObsProbesVel_,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                }
        		break;
		}
	}
	else if (obsData.get_obsType() == "p"){
        for (int i = 0 ; i < nObsProbesPres_; i++)
        {
            obsCoords_(i,0) = obsData.get_probeCoord(probeIDs_(i))[0];
            obsCoords_(i,1) = obsData.get_probeCoord(probeIDs_(i))[1];
            obsCoords_(i,2) = obsData.get_probeCoord(probeIDs_(i))[2];
        }
	}
	else if (obsData.get_obsType() == "Up"){
		switch (obsData.get_velocityCaseSwitch()){
			case 1 :
			case 2 :
			case 3 :
                for (int i = 0; i < nObsProbesVel_; ++i){
                    obsCoords_(i,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                }
                for (int i = nObsProbesVel_; i < (nObs_); ++i){
                    obsCoords_(i,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                }
				break;
			case 4 :
			case 5 :
			case 6 :
                for (int i = 0; i < nObsProbesVel_; ++i){
                    obsCoords_(i+nObsProbesVel_,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i+nObsProbesVel_,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i+nObsProbesVel_,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                    obsCoords_(i+2*nObsProbesVel_,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i+2*nObsProbesVel_,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i+2*nObsProbesVel_,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                }
                for (int i = 2*nObsProbesVel_; i < (nObs_); ++i){
                    obsCoords_(i,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                }
        		break;
			case 7 :
                for (int i = 0; i < nObsProbesVel_; ++i){
                    obsCoords_(i,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                    obsCoords_(i+nObsProbesVel_,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i+nObsProbesVel_,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i+nObsProbesVel_,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                    obsCoords_(i+2*nObsProbesVel_,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i+2*nObsProbesVel_,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i+2*nObsProbesVel_,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                }
                for (int i = 3*nObsProbesVel_; i < (nObs_); ++i){
                    obsCoords_(i,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                }
        		break;
		}
	}
	else if (obsData.get_obsType() == "UCf"){
		switch (obsData.get_velocityCaseSwitch()){
			case 1 :
			case 2 :
			case 3 :
                for (int i = 0; i < nObsProbesVel_; ++i){
                    obsCoords_(i,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                }
                for (int i = nObsProbesVel_; i < (nObs_); ++i){
                    obsCoords_(i,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                }
				break;
			case 4 :
			case 5 :
			case 6 :
                for (int i = 0; i < nObsProbesVel_; ++i){
                    obsCoords_(i+nObsProbesVel_,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i+nObsProbesVel_,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i+nObsProbesVel_,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                    obsCoords_(i+2*nObsProbesVel_,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i+2*nObsProbesVel_,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i+2*nObsProbesVel_,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                }
                for (int i = 2*nObsProbesVel_; i < (nObs_); ++i){
                    obsCoords_(i,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                }
        		break;
			case 7 :

                for (int i = 0; i < nObsProbesVel_; ++i){
                    obsCoords_(i,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                    obsCoords_(i+nObsProbesVel_,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i+nObsProbesVel_,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i+nObsProbesVel_,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                    obsCoords_(i+2*nObsProbesVel_,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i+2*nObsProbesVel_,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i+2*nObsProbesVel_,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                }
                for (int i = 3*nObsProbesVel_; i < (nObs_); ++i){
                    obsCoords_(i,0) = obsData.get_probeCoord(probeIDs_(i))[0];
                    obsCoords_(i,1) = obsData.get_probeCoord(probeIDs_(i))[1];
                    obsCoords_(i,2) = obsData.get_probeCoord(probeIDs_(i))[2];
                }
        		break;
		}
    }
    
}

/* Getters */

/* Destructor */
region::~region()
{}

/* Outside class but related functions */
Eigen::ArrayXi removeDuplicatesSet(const Eigen::Ref<const Eigen::ArrayXi>& arr) 
{
    std::set<int> uniqueValues(arr.data(), arr.data() + arr.size());
    Eigen::ArrayXi result(uniqueValues.size());
    std::copy(uniqueValues.begin(), uniqueValues.end(), result.data());
    return result;
}

Eigen::ArrayXi read_clipFile_info(const std::string stringRootPath, const std::string filename, float cwipiVerbose)
{
    /* 
    Input : clippingCells.txt file
    Output : Array containing the values of the clipping file
    */

    std::string clipping_file = stringRootPath + "/" + filename;
    std::cout << clipping_file << std::endl;
 
    int count_clippingCells = 0;

    std::vector<std::vector<std::string>> clips_content;
    std::vector<std::string> clips_row;
    std::string clips_line, clips_word;

    std::fstream clips_stream (clipping_file, std::ios::in);
    if(clips_stream.is_open())
    {
        while(getline(clips_stream, clips_line))
        {
            clips_row.clear();
 
            std::stringstream clips_str(clips_line);
 
            while(getline(clips_str, clips_word, ','))
                clips_row.push_back(clips_word);
            clips_content.push_back(clips_row);
            count_clippingCells = count_clippingCells + 1; 
        }
    }
    else {
        std::cerr << "Couldn't open clipping file.\n";
    }

    if (cwipiVerbose) std::cout << "clippingCells.txt file succesfully read" << std::endl;
    if (cwipiVerbose) std::cout << "count_clippingCells = " << count_clippingCells << std::endl;


    //======== Create an eigen array containing the ID of the regions ========
    Eigen::ArrayXi clippingFile(count_clippingCells); 
    for(int i=0; i<count_clippingCells; i++)
    {
        clippingFile(i) = std::stod(clips_content[i][0]);
    }

    if (cwipiVerbose) std::cout << "Array of the clippingCells created" << std::endl;
    // if (cwipiVerbose) std::cout << "clippingCells = " << clippingCells << std::endl;
    return clippingFile;
}

std::vector<int> calculate_nCells_clips(const Eigen::Ref<const Eigen::ArrayXi>& clippingFile, float cwipiVerbose)
{
    /* 
    Input : Array of the values in the clipping file
    Output : Array containing the number of cells for each clipping
    */

    if (cwipiVerbose) std::cout << "Calculating the number of cells per clip from the clip file ..." << std::endl;

    int count_cellsClipsDone = 0;   // Count the number of cells of all the clipping already done in the hyperlocalization

    int nb_clipCells = clippingFile(0);

    int clippingFileSize = clippingFile.rows();

    std::vector<int> nCells_clips;
    
    while (count_cellsClipsDone < clippingFileSize)
    {
        nb_clipCells = clippingFile(count_cellsClipsDone);
        nCells_clips.push_back(nb_clipCells);
        // if (cwipiVerbose) std::cout << "calculate_nCells_clips : pushing " << nb_clipCells << std::endl;
        count_cellsClipsDone = count_cellsClipsDone + nb_clipCells + 1;
    }

    // if (cwipiVerbose) std::cout << "calculate_nCells_clips : Size of nCells_clips = " << nCells_clips.size() << std::endl;

    return nCells_clips;
}

std::vector<Eigen::ArrayXi> calculate_cellIDs_clips(const Eigen::Ref<const Eigen::ArrayXi>& clippingFile, float cwipiVerbose)
{
    /* 
    Input : Array of the values in the clipping file
    Output : List of arrays containing the IDs of the cells for each clipping
    */
   
    if (cwipiVerbose) std::cout << "Retreiving the cells IDs for each clip from the clip file ..." << std::endl;

    std::vector<Eigen::ArrayXi> cellIDs_clips;

    int count_cellsClipsDone = 0;   // Count the number of cells of all the clipping already done in the hyperlocalization
    int count_nbClips = 0;          // Count the number of regions 
    int IDclippingCell = 0;         // Current cell index of the clipping file

    int nb_clipCells = clippingFile(0);
    Eigen::ArrayXi oneClipCells = Eigen::ArrayXi::Zero(nb_clipCells);

    int clippingFileSize = clippingFile.rows();

    for (int i=0; i<clippingFileSize; i++)
    {
        // if (cwipiVerbose) std::cout << "cellIDs : inside count_clippingCells loop i = " << i << std::endl;
        if (i==(nb_clipCells + count_cellsClipsDone)){
            count_nbClips = count_nbClips +1;
            cellIDs_clips.push_back(oneClipCells);

            if (i<(clippingFileSize-1)){                                             // If regions are not done
                count_cellsClipsDone = count_cellsClipsDone + nb_clipCells +1;       // Update the count of the clipping already performed
                nb_clipCells = clippingFile(i+1);                                    // Update the number of cells in next clipping
                oneClipCells.resize(nb_clipCells);        // Resize temp_state over the new number of cells in the clipping
            }

        }
        else{
            IDclippingCell = clippingFile(i+1);
            oneClipCells(i - count_cellsClipsDone)=IDclippingCell;
            }
    }

    return cellIDs_clips;
}

std::vector<region> clips_to_DA_regions(std::vector<int> nCells_clips, std::vector<Eigen::ArrayXi> cellIDs_clips, float cwipiVerbose)
{
    /*
        Input : Raw clipping regions
        Output : merges the intersecting ones and constructs the data assimilation regions by making 
        the connection between the observations and regions.
    */

    if (cwipiVerbose) std::cout << "Preparing the data asssimilation regions.. \n" << std::endl;

    int nClips          = nCells_clips.size();   // number of clipping regions 
    int untreatedClips  = nClips;
    int count_region    = 0;                     // variable used as a counter for the DA regions (not the eventual number of regions)
    int nRegions        = 0;                     // Number of regions in the end of the preprocessing
    bool isRegionUnique = true;                  // To check if the region is unique or if a merge will happen

    std::vector<Eigen::ArrayXi> cellIDs_regions = cellIDs_clips; // List of DA regions (start by copying clipping cellIDs)
    std::vector<Eigen::ArrayXi> region_to_clip(nClips);   // Initialize the list of regions considering 1 clip = 1 region before merging them
    for (int i=0; i<nClips; i++)
    {
        region_to_clip[i] = Eigen::ArrayXi::Ones(1);
        region_to_clip[i] = region_to_clip[i]*i;
    }

    //** Loop over the clipping regions and merge the intersecting ones together **
    while (untreatedClips > 0)
    {
        // if (cwipiVerbose) std::cout << "\nREGION - "<< count_region << "-----------------------------------------------" << std::endl;
        // if (cwipiVerbose) std::cout << "untreated clips : " << untreatedClips <<", count_region = " << count_region << std::endl;

        Eigen::ArrayXi cellsRegion = cellIDs_regions[count_region];     //  Cell IDs of the current region

        if (region_to_clip[count_region].rows() == 1)
        {
            untreatedClips -= 1;
            // if (cwipiVerbose) std::cout << "One clip in the region, so untreatedClips reduced to : " << untreatedClips << std::endl;
        }

        // if (cwipiVerbose) std::cout << "- Checking overlapping with other clips ..." << std::endl;
        isRegionUnique = true;                                          // To check if the region is unique or if a merge will happen
        for (unsigned int j=count_region+1; j<cellIDs_regions.size(); j++)
        {
            // if (cwipiVerbose) std::cout << "j = " << j << "+++++++" << std::endl; // debugg
            // if (cwipiVerbose) std::cout << "cellsRegion = " << cellsRegion << std::endl;
            // if (cwipiVerbose) std::cout << "cellIDs_regions[" << j << "] = " <<cellIDs_regions[j] << std::endl;

            Eigen::ArrayXi mergedCellIDs(cellsRegion.rows() + cellIDs_regions[j].rows());     //Merge the current clip with the others
            mergedCellIDs << cellsRegion, 
                             cellIDs_regions[j];

            Eigen::ArrayXi uniqueMergedCellIDs = removeDuplicatesSet(mergedCellIDs);          //Create an array with non-repeated cells

            if (uniqueMergedCellIDs.rows() != mergedCellIDs.rows())                           //If the length of the arrays are different some cells are repeated, overlapping => merge
            {
            // In case of overlapping!
                // if (cwipiVerbose) std::cout << "===>> THERE IS AN OVERLAPPING! MERGING CLIP-" << region_to_clip[j][0] << " with REGION-" << count_region << ".." << std::endl;
                cellIDs_regions[count_region] = uniqueMergedCellIDs;                      // update the CellIDs list with the merged one
                cellIDs_regions.erase(cellIDs_regions.begin()+j);                         // and also delete the row of the clipping being merged
                Eigen::ArrayXi temp_region(region_to_clip[count_region].rows() + region_to_clip[j].rows());
                temp_region << region_to_clip[count_region],
                               region_to_clip[j];
                region_to_clip[count_region] = temp_region;                               // add the ID of the merged clipping to the region
                region_to_clip.erase(region_to_clip.begin()+j);                           // and also delete the region of the clipping being merged
                untreatedClips -= 1;                                                      // One more Clip gone..
                isRegionUnique = false;                                                   // Boolean to mark if the region is unique or not 
                break;
            }
            // if (cwipiVerbose) std::cout << "NO OVERLAPPING FOUND" << std::endl;
        }
    
        if (isRegionUnique == true)
        {
            // if (cwipiVerbose) std::cout << "REGION IS NOW UNIQUE SKIPPING! ==========" << std::endl;
            count_region += 1;                           // if there is no overlapping, increase the region count by one
        }
    }

    nRegions = region_to_clip.size();                    // number of regions after merging the overlapping clipping regions

    Eigen::ArrayXi nCells_region = Eigen::ArrayXi::Zero(nRegions);
    Eigen::ArrayXi nClips_region = Eigen::ArrayXi::Zero(nRegions);

    std::vector<region> regionList(nRegions);            // List of DA_Regions e.g. 4th DA region can be reached by regionList[3]

    for (int i=0; i<nRegions; i++)
    {
        // nCells_region[i] = cellIDs_regions[i].rows();                                                                    // debugg
        // nClips_region[i] = region_to_clip[i].rows();                                                                     // debugg
        regionList[i] = region(region_to_clip[i], cellIDs_regions[i]);
    }

    // if (cwipiVerbose) std::cout << "\n Number of Regions : " << nRegions << std::endl;                                   // debugg
    // if (cwipiVerbose) std::cout << "\n Number of clips in each regionList : " << nClips_region << std::endl;             // debugg
    // if (cwipiVerbose) std::cout << "\n Number of Cells in each regionList : " << nCells_region << std::endl << "\n";     // debugg
    // if (cwipiVerbose) std::cout << "\n IDs of the clips for each DA regionList : " << region_to_clip << std::endl;       // debugg

    // return the list of DA_Regions with all the relevant information
    return regionList;
}

std::vector<region> initialize_DA_regions(const std::string stringRootPath, const std::string filename, float cwipiVerbose, bool clippingSwitch, int nb_cells, const observations obsData)
{
    std::vector<region> regionList;
    if (clippingSwitch)
    {
        //** Build the regions from the clipping file **
        Eigen::ArrayXi clippingFile = read_clipFile_info(stringRootPath, filename, cwipiVerbose);
        std::vector<int> nCells_clips = calculate_nCells_clips(clippingFile, cwipiVerbose);
        std::vector<Eigen::ArrayXi> cellIDs_clips = calculate_cellIDs_clips(clippingFile, cwipiVerbose);
        regionList = clips_to_DA_regions(nCells_clips, cellIDs_clips, cwipiVerbose);

        /** Initialize the count of observations probes for each region **/
        for (unsigned int i = 0; i < regionList.size(); i++)
        {
            regionList[i].set_nObsProbes(obsData);
        }
    }
    else
    {
        //** Build the only region with the total number of cells of the simulation **
        Eigen::ArrayXi clipIDs(1);
        clipIDs(0)=0;
        Eigen::ArrayXi cellIDs(nb_cells);
        for (int i=0; i<nb_cells; i++){cellIDs(i)=i;}
        Eigen::ArrayXi probesIDs(obsData.get_nObsProbes());
        for (int i=0; i<obsData.get_nObsProbes(); i++){probesIDs(i)=i;}
        regionList.push_back(region(clipIDs, cellIDs, probesIDs));
        regionList[0].set_nObsProbesNoClip(obsData);
    }

    if (cwipiVerbose) std::cout << "=>> REPORT - Data Assimilation Regions ===================" << std::endl;
    if (cwipiVerbose) std::cout << " Size of the regionList = " << regionList.size() << std::endl;
    if (cwipiVerbose){
        for (unsigned int i=0; i<regionList.size(); i++)
        {
            std::cout << "\t Region - " << i <<" :  Number of Clips = " << regionList[i].get_nClips() << "   |   Number of Cells = " << regionList[i].get_nCells();
            std::cout << "   |   clip IDs = " << regionList[i].get_clipIDs().transpose() << "   |   Number of Vel Probes = " << regionList[i].get_nObsProbesVel();
            std::cout << "   |   Number of Pres Probes = " << regionList[i].get_nObsProbesPres() << "   |   Number of Cf Probes = " << regionList[i].get_nObsProbesCf();
            std::cout << "   |   Total Number of observations = " << regionList[i].get_nObs();
            // std::cout << "   |   Cells IDs = " << regionList[i].get_cellIDs().transpose();
            std::cout << std::endl;
        }
        std::cout << std::endl << "\n";
    }

    return regionList;
}