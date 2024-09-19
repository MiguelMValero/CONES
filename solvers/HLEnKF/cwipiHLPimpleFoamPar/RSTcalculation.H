scalar deltaT = runTime.deltaTValue();
scalar currentTime = runTime.value();
scalar startTime = runTime.startTime().value();
scalar beta = deltaT/(currentTime-startTime);

Info << "startTime = " << startTime << nl << endl;
Info << "currentTime = " << currentTime << nl << endl;
Info << "deltaT = " << deltaT << nl << endl;
Info << "beta = " << beta << nl << endl;

Info << "Premier terme de UMean = " << UMean[0].component(0) << endl;
Info << "Premier terme de U2Mean = " << U2Mean[0].component((symmTensor::XX)) << endl;

U2Mean = (1-beta)*U2Mean + beta*(sqr(U)); 
UMean = (1-beta)*UMean + beta*U;

Info << "Premier terme de UMean_upt = " << UMean[0].component(0) << endl;
Info << "Premier terme de U2Mean_upt = " << U2Mean[0].component((symmTensor::XX)) << endl;