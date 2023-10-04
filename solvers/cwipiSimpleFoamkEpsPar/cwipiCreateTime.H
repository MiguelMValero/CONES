Foam::Info<< "cwipi Create time \n" << Foam::endl;

int myGlobalRank;
MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);
char rankChar[50];
sprintf(rankChar, "%i", myGlobalRank-1);
std::cout << "rankChar" << rankChar << std::endl;

std::string globalCaseName = args.globalCaseName();
std::string globalRootPath = args.rootPath() + "/" + globalCaseName; //For Uinterpolation path

char RunCase[250];
strcpy(RunCase,globalCaseName.c_str());
strcat(RunCase,"/processor");
strcat(RunCase,rankChar);

Foam::Time runTime(Foam::Time::controlDictName, args.rootPath(), RunCase);

