/*
  This file is part of the CWIPI library. 

  Copyright (C) 2011  ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mpi.h>

#include <cstring>
#include <ctime>
#include <cstdlib>

#include <bftc_mem.h>

#include <fvmc_parall.h>

#include "cwipi_config.h"
#include "cwipi_priv.h"
#include "applicationPropertiesDataBase.hxx"
#include "applicationPropertiesDataBase_i.hxx"

namespace cwipi {

  ApplicationPropertiesDataBase::ApplicationPropertiesDataBase()
    : _distantApplicationPropertiesDataBase(*(new std::map <std::string, ApplicationProperties * > ())),
      _localApplicationProperties(NULL)
  {
    //
    // Time protect : stop application if time > 2010
#if defined(PROTECT)
    char buf[5];
    time_t t = time(0);
    strftime(buf, 5, "%Y\n", localtime(&t));
    int year = atoi(buf);
    if (year > 2010) {
      std::cout << "Use of CWIPI is limited to December 31, 2010" << std::endl;
      exit(1);
    }
#endif

  }

  ApplicationPropertiesDataBase::~ApplicationPropertiesDataBase()
  {
#if defined(DEBUG) && 0
    std::cout << "destroying ApplicationPropertiesDataBase." << std::endl;
#endif
    if (_localApplicationProperties != NULL)
      delete _localApplicationProperties;

    if (!_distantApplicationPropertiesDataBase.empty()) {
      typedef std::map <std::string, ApplicationProperties * >::iterator CI;
      for (CI p = _distantApplicationPropertiesDataBase.begin();
           p != _distantApplicationPropertiesDataBase.end(); p++) {
        if (p->second != NULL)
          delete p->second;
      }
      _distantApplicationPropertiesDataBase.clear();

    }
    delete &_distantApplicationPropertiesDataBase;
  }



  MPI_Comm  ApplicationPropertiesDataBase::init(const char* applicationName,
                                                 const MPI_Comm globalComm)
  {

    // Initialize MPI
    // --------------

    int currentRank;
    int globalCommSize;
    int color = 0;

    {
      int flag;

      MPI_Initialized(&flag);
      if (!flag)
        bftc_error(__FILE__, __LINE__, 0, "MPI is not initialized\n");

      MPI_Comm_rank(globalComm, &currentRank);
      MPI_Comm_size(globalComm, &globalCommSize);

    }


    // Search applications
    // -------------------

    {
      int j = 0;
      int index = 0;
      int totalLength = 0;
      int nameLength = 0;

      std::string currentString = "";

      nameLength = strlen(applicationName) + 1;

      MPI_Allreduce (&nameLength, &totalLength, 1, MPI_INT, MPI_SUM,
                     globalComm);

//       BFT::BFTC_MALLOC(mergeNames, totalLength, char) ;
      char *mergeNames =  (char *) malloc (sizeof(char) * (totalLength));

      int *namesLength =  (int *) malloc (sizeof(int) * (globalCommSize));
      int *iproc =  (int *) malloc (sizeof(int) * (globalCommSize));

      MPI_Allgather(&nameLength,
                    1,
                    MPI_INT,
                    namesLength,
                    1,
                    MPI_INT,
                    globalComm);

      iproc[0] = 0;
      for(int i = 1; i < globalCommSize; i++)
        iproc[i] = namesLength[i-1] + iproc[i-1];

      MPI_Allgatherv((void*) const_cast <char*> (applicationName),
                     nameLength,
                     MPI_CHAR,
                     mergeNames,
                     namesLength,
                     iproc,
                     MPI_CHAR,
                     globalComm);

      free ( iproc);
      free ( namesLength);

      for (int irank = 0; irank < globalCommSize; irank++) {

        assert(index <= totalLength);

        const char *ptCurrentName = mergeNames + index;
        std::string currentName(ptCurrentName);

        if (currentString != currentName) {

          ApplicationProperties *currentApplicationProperties =
            new ApplicationProperties(currentName, globalComm);

          if (!strcmp(currentName.c_str(), applicationName)) {
            _localApplicationProperties = currentApplicationProperties;
            color = j+1;
          }
          else {
            std::pair<std::string, ApplicationProperties *>
              newPair(std::string(currentName), currentApplicationProperties);

            std::pair<std::map<std::string, ApplicationProperties *>::iterator, bool>
              p = _distantApplicationPropertiesDataBase.insert(newPair);

            if (!p.second)
              bftc_error(__FILE__, __LINE__, 0,
                        "The MPI ranks range is not continous or\n"
                        "There are two applications with the same name '%s'  \n", currentName.c_str());
          }

          currentApplicationProperties->setBeginningRank(irank);
          if (currentString != "") {
            if (!strcmp(currentString.c_str(), applicationName))
              _localApplicationProperties->setEndRank(irank-1);
            else
              _distantApplicationPropertiesDataBase[currentString]->setEndRank(irank-1);
          }
          currentString = currentName;

          j += 1;
        }

        if (currentString != "") {
          if (!strcmp(currentString.c_str(), applicationName))
            _localApplicationProperties->setEndRank(globalCommSize-1);
          else
            _distantApplicationPropertiesDataBase[currentString]->setEndRank(globalCommSize-1);
        }

        index += currentName.size() + 1;
        assert(index <= totalLength);
      }

      free ( mergeNames);

      // Create current application communicator
      // ---------------------------------------

      MPI_Comm localComm = MPI_COMM_NULL ;
      MPI_Comm_split(globalComm, color, currentRank, &localComm);
      _localApplicationProperties->setLocalComm(localComm);

      return localComm;
    }
  }


  void ApplicationPropertiesDataBase::mergeParameters(const std::string &applicationName)
  {
    _mergeIntParameters(applicationName);
    _mergeDoubleParameters(applicationName);
    _mergeStringParameters(applicationName);
  }


  void ApplicationPropertiesDataBase::_mergeIntParameters(const std::string &applicationName)
  {

    typedef std::map <std::string, int>::iterator Iterator;

    MPI_Status status;

    const std::map <std::string, ApplicationProperties * >::iterator p = _distantApplicationPropertiesDataBase.find(applicationName);
    if (p == _distantApplicationPropertiesDataBase.end())
      bftc_error(__FILE__, __LINE__, 0,
                "'%s' application not found \n", applicationName.c_str());

    std::map <std::string, int> &localControlParameters   = _localApplicationProperties->_intControlParameters;
    std::map <std::string, int> &distantControlParameters = p->second->_intControlParameters;

    int NLocalControlParameters = localControlParameters.size();
    int NDistantControlParameters = distantControlParameters.size();

    const MPI_Comm& localComm    = _localApplicationProperties->getLocalComm();
    const MPI_Comm& globalComm   = _localApplicationProperties->getGlobalComm();

    int distantBeginningRank = p->second->_beginningRank;

    //
    // Clear distant parameters copy

    distantControlParameters.clear();

    //
    // Check that all local pocesses are the same parameter values

    int localCommSize = -1;
    int currentRank = -1;

    MPI_Comm_rank(localComm, &currentRank);
    MPI_Comm_size(localComm, &localCommSize);

    if (localCommSize > 1) {

      for (Iterator p1 = localControlParameters.begin(); p1 != localControlParameters.end(); p1++) {
        int value             = p1->second;
        const char *paramName = p1->first.c_str();
        int nameSize    = p1->first.size();

        if (currentRank == 0) {
          for (int irank = 1; irank < localCommSize; irank++){

            int   distantNameSize = 0;
            char *distantParamName = NULL;
            int   distantValue;

            MPI_Recv(&distantNameSize, 1, MPI_INT, irank, 0, localComm, &status);
            if (distantNameSize != nameSize)
              bftc_error(__FILE__, __LINE__, 0,
                        "Inconsistency between local parameters\n");

            distantParamName =  (char *) malloc (sizeof(char) * (distantNameSize+1));
            MPI_Recv(distantParamName, distantNameSize+1, MPI_CHAR, irank, 0, localComm, &status);
            if (strcmp(distantParamName, paramName))
              bftc_error(__FILE__, __LINE__, 0,
                        "Inconsistency between local parameters\n");
            free ( distantParamName);

            MPI_Recv(&distantValue, 1, MPI_INT, irank, 0, localComm, &status);
            if (distantValue != value)
              bftc_error(__FILE__, __LINE__, 0,
                        "Different values for '%s' parameter for the rank 0 and the rank '%i'\n", paramName, irank);
          }
        }

        else {
          MPI_Send(&nameSize, 1, MPI_INT, 0, 0, localComm);
          MPI_Send(const_cast <char* > (paramName), nameSize+1, MPI_CHAR, 0, 0, localComm);
          MPI_Send(&value, 1, MPI_INT, 0, 0, localComm);
        }
      }
    }

    //
    // parameters exchange between beginning ranks

    if (currentRank == 0 ) {

      MPI_Sendrecv(&NLocalControlParameters,   1, MPI_INT, distantBeginningRank, 100,
                   &NDistantControlParameters, 1, MPI_INT, distantBeginningRank, 100,
                   globalComm, &status);

      int NParameterMax = MAX(NLocalControlParameters, NDistantControlParameters);

      Iterator p1 = localControlParameters.begin();

      for (int i = 0; i < NParameterMax; i++ ) {

        if (i >= NDistantControlParameters) {

          int value             = p1->second;
          const char *paramName = p1->first.c_str();
          int nameSize          = p1->first.size();

          MPI_Send(&nameSize, 1, MPI_INT, distantBeginningRank, 0, globalComm);
          MPI_Send(const_cast <char *> (paramName), nameSize+1, MPI_CHAR, distantBeginningRank, 0, globalComm);
          MPI_Send(&value, 1, MPI_INT, distantBeginningRank, 0, globalComm);

        }

        else if (p1 == localControlParameters.end()){

          int   distantNameSize = 0;
          char *distantParamName = NULL;
          int   distantValue;

          MPI_Recv(&distantNameSize, 1, MPI_INT, distantBeginningRank, 0, globalComm, &status);

          distantParamName =  (char *) malloc (sizeof(char) * (distantNameSize+1));
          MPI_Recv(distantParamName, distantNameSize+1, MPI_CHAR, distantBeginningRank, 0, globalComm, &status);

          MPI_Recv(&distantValue, 1, MPI_INT, distantBeginningRank, 0, globalComm, &status);

          distantControlParameters[std::string(distantParamName)] = distantValue;

          free ( distantParamName);

        }

        else {

          int value             = p1->second;
          const char *paramName = p1->first.c_str();
          int nameSize          = p1->first.size();

          int   distantNameSize = 0;
          char *distantParamName = NULL;
          int   distantValue;

          MPI_Sendrecv(&nameSize       , 1, MPI_INT, distantBeginningRank, 0,
                       &distantNameSize, 1, MPI_INT, distantBeginningRank, 0,
                       globalComm, &status);

          distantParamName =  (char *) malloc (sizeof(char) * (distantNameSize+1));

          MPI_Sendrecv(const_cast <char*> (paramName), nameSize+1,        MPI_CHAR, distantBeginningRank, 0,
                       distantParamName              , distantNameSize+1, MPI_CHAR, distantBeginningRank, 0,
                       globalComm, &status);

          MPI_Sendrecv(&value       , 1, MPI_INT, distantBeginningRank, 0,
                       &distantValue, 1, MPI_INT, distantBeginningRank, 0,
                       globalComm, &status);

          distantControlParameters[std::string(distantParamName)] = distantValue;

          free ( distantParamName);
        }

        if (p1 != localControlParameters.end())
          p1++;
      }
    }

    //
    // Local beginning rank send parameters to other local ranks

    if (localCommSize > 1) {

      int size = 0;
      if (currentRank == 0)
        size =  distantControlParameters.size();

      MPI_Bcast(&size, 1, MPI_INT, 0, localComm);

      Iterator p1;
      if (currentRank == 0)
        p1 = distantControlParameters.begin();

      for (int i = 0 ; i < size; i++) {

        int value;
        char *paramName = NULL;
        int nameSize;

        if (currentRank == 0) {
          value     = p1->second;
          paramName = const_cast <char *> (p1->first.c_str());
          nameSize  = p1->first.size();
          p1++;
        }

        MPI_Bcast(&nameSize, 1, MPI_INT, 0, localComm);


        if (currentRank != 0)
          paramName =  (char *) malloc (sizeof(char) * (nameSize+1));

        MPI_Bcast(paramName, nameSize+1, MPI_CHAR, 0, localComm);

        MPI_Bcast(&value, 1, MPI_INT, 0, localComm);

        if (currentRank != 0) {
          distantControlParameters[std::string(paramName)] = value;
          free ( paramName);
        }
      }
    }
  }


  void ApplicationPropertiesDataBase::_mergeDoubleParameters(const std::string &applicationName)
  {

    typedef std::map <std::string, double>::iterator Iterator;

    MPI_Status status;
    //std::string applicationNameStr = applicationName;

    const std::map <std::string, ApplicationProperties * >::iterator p = _distantApplicationPropertiesDataBase.find(applicationName);
    if (p == _distantApplicationPropertiesDataBase.end())
      bftc_error(__FILE__, __LINE__, 0,
                "'%s' application not found \n", applicationName.c_str());

    std::map <std::string, double> &localControlParameters   = _localApplicationProperties->_doubleControlParameters;
    std::map <std::string, double> &distantControlParameters = p->second->_doubleControlParameters;

    int NLocalControlParameters = localControlParameters.size();
    int NDistantControlParameters = distantControlParameters.size();

    const MPI_Comm& localComm    = _localApplicationProperties->getLocalComm();
    const MPI_Comm& globalComm   = _localApplicationProperties->getGlobalComm();

    int distantBeginningRank = p->second->_beginningRank;

    //
    // Clear distant parameters copy

    distantControlParameters.clear();
    distantControlParameters.clear();

    //
    // Check that all local pocesses are the same parameter values

    int localCommSize = -1;
    int currentRank = -1;

    MPI_Comm_rank(localComm, &currentRank);
    MPI_Comm_size(localComm, &localCommSize);

    if (localCommSize > 1) {

      for (Iterator p1 = localControlParameters.begin(); p1 != localControlParameters.end(); p1++) {
        double value             = p1->second;
        const char *paramName = p1->first.c_str();
        int nameSize    = p1->first.size();

        if (currentRank == 0) {
          for (int irank = 1; irank < localCommSize; irank++){

            int   distantNameSize = 0;
            char *distantParamName = NULL;
            double   distantValue;

            MPI_Recv(&distantNameSize, 1, MPI_INT, irank, 0, localComm, &status);
            if (distantNameSize != nameSize)
              bftc_error(__FILE__, __LINE__, 0,
                        "Inconsistency between local parameters\n");

            distantParamName =  (char *) malloc (sizeof(char) * (distantNameSize+1));
            MPI_Recv(distantParamName, distantNameSize+1, MPI_CHAR, irank, 0, localComm, &status);
            if (strcmp(distantParamName, paramName))
              bftc_error(__FILE__, __LINE__, 0,
                        "Inconsistency between local parameters\n");
            free ( distantParamName);

            MPI_Recv(&distantValue, 1, MPI_DOUBLE, irank, 0, localComm, &status);
            
            CWIPI_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")
            if (distantValue != value)
              bftc_error(__FILE__, __LINE__, 0,
                        "Different values for '%s' parameter for the rank 0 and the rank '%i'\n", paramName, irank);
            CWIPI_GCC_SUPPRESS_WARNING_POP
          }
        }

        else {
          MPI_Send(&nameSize, 1, MPI_INT, 0, 0, localComm);
          MPI_Send(const_cast <char* > (paramName), nameSize+1, MPI_CHAR, 0, 0, localComm);
          MPI_Send(&value, 1, MPI_DOUBLE, 0, 0, localComm);
        }
      }
    }

    //
    // parameters exchange between beginning ranks

    if (currentRank == 0 ) {

      MPI_Sendrecv(&NLocalControlParameters,   1, MPI_INT, distantBeginningRank, 0,
                   &NDistantControlParameters, 1, MPI_INT, distantBeginningRank, 0,
                   globalComm, &status);

      int NParameterMax = MAX(NLocalControlParameters, NDistantControlParameters);

      Iterator p1 = localControlParameters.begin();

      for (int i = 0; i < NParameterMax; i++ ) {

        if (i >= NDistantControlParameters) {

          double value             = p1->second;
          const char *paramName = p1->first.c_str();
          int nameSize          = p1->first.size();

          MPI_Send(&nameSize, 1, MPI_INT, distantBeginningRank, 0, globalComm);
          MPI_Send(const_cast <char *> (paramName), nameSize+1, MPI_CHAR, distantBeginningRank, 0, globalComm);
          MPI_Send(&value, 1, MPI_DOUBLE, distantBeginningRank, 0, globalComm);

        }

        else if (p1 == localControlParameters.end()){

          int   distantNameSize = 0;
          char *distantParamName = NULL;
          double   distantValue;

          MPI_Recv(&distantNameSize, 1, MPI_INT, distantBeginningRank, 0, globalComm, &status);

          distantParamName =  (char *) malloc (sizeof(char) * (distantNameSize+1));
          MPI_Recv(distantParamName, distantNameSize+1, MPI_CHAR, distantBeginningRank, 0, globalComm, &status);

          MPI_Recv(&distantValue, 1, MPI_DOUBLE, distantBeginningRank, 0, globalComm, &status);

          distantControlParameters[std::string(distantParamName)] = distantValue;
          free ( distantParamName);

        }

        else {

          double value             = p1->second;
          const char *paramName = p1->first.c_str();
          int nameSize          = p1->first.size();

          int   distantNameSize = 0;
          char *distantParamName = NULL;
          double   distantValue;

          MPI_Sendrecv(&nameSize       , 1, MPI_INT, distantBeginningRank, 0,
                       &distantNameSize, 1, MPI_INT, distantBeginningRank, 0,
                       globalComm, &status);

          distantParamName =  (char *) malloc (sizeof(char) * (distantNameSize+1));
          MPI_Sendrecv(const_cast <char*> (paramName), nameSize+1,        MPI_CHAR, distantBeginningRank, 0,
                       distantParamName              , distantNameSize+1, MPI_CHAR, distantBeginningRank, 0,
                       globalComm, &status);

          MPI_Sendrecv(&value       , 1, MPI_DOUBLE, distantBeginningRank, 0,
                       &distantValue, 1, MPI_DOUBLE, distantBeginningRank, 0,
                       globalComm, &status);

          distantControlParameters[std::string(distantParamName)] = distantValue;
          free ( distantParamName);

        }

        if (p1 != localControlParameters.end())
          p1++;
      }
    }

    //
    // Local beginning rank send parameters to other local ranks

    if (localCommSize > 1) {

      int size = 0;
      if (currentRank == 0)
        size =  distantControlParameters.size();

      MPI_Bcast(&size, 1, MPI_INT, 0, localComm);

      Iterator p1;
      if (currentRank == 0)
        p1 = distantControlParameters.begin();

      for (int i = 0 ; i < size; i++) {
        double value;
        char *paramName = NULL;
        int nameSize;

        if (currentRank == 0) {
          value     = p1->second;
          paramName = const_cast <char *> (p1->first.c_str());
          nameSize  = p1->first.size();
          p1++;
        }

        MPI_Bcast(&nameSize, 1, MPI_INT, 0, localComm);

        if (currentRank != 0)
          paramName =  (char *) malloc (sizeof(char) * (nameSize+1));

        MPI_Bcast(paramName, nameSize+1, MPI_CHAR, 0, localComm);

        MPI_Bcast(&value, 1, MPI_DOUBLE, 0, localComm);

        if (currentRank != 0) {
          distantControlParameters[std::string(paramName)] = value;
          free ( paramName);
        }
      }
    }
  }

  void ApplicationPropertiesDataBase::_mergeStringParameters(const std::string &applicationName)
  {

    typedef std::map <std::string, std::string>::iterator Iterator;

    MPI_Status status;

    const std::map <std::string, ApplicationProperties * >::iterator p = _distantApplicationPropertiesDataBase.find(applicationName);
    if (p == _distantApplicationPropertiesDataBase.end())
      bftc_error(__FILE__, __LINE__, 0,
                "'%s' application not found \n", applicationName.c_str());

    std::map <std::string, std::string> &localControlParameters   = _localApplicationProperties->_stringControlParameters;
    std::map <std::string, std::string> &distantControlParameters = p->second->_stringControlParameters;

    int NLocalControlParameters = localControlParameters.size();
    int NDistantControlParameters = distantControlParameters.size();

    const MPI_Comm& localComm    = _localApplicationProperties->getLocalComm();
    const MPI_Comm& globalComm   = _localApplicationProperties->getGlobalComm();

    int distantBeginningRank = p->second->_beginningRank;

    //
    // Clear distant parameters copy

    distantControlParameters.clear();

    //
    // Check that all local pocesses are the same parameter values


    int localCommSize = -1;
    int currentRank = -1;

    MPI_Comm_rank(localComm, &currentRank);
    MPI_Comm_size(localComm, &localCommSize);

    if (localCommSize > 1) {

      for (Iterator p1 = localControlParameters.begin(); p1 != localControlParameters.end(); p1++) {
        std::string value     = p1->second;
        const char *valueCStr  = value.c_str();
        int valueSize = value.size();

        const char *paramName = p1->first.c_str();
        int nameSize          = p1->first.size();

        if (currentRank == 0) {
          for (int irank = 1; irank < localCommSize; irank++){

            int   distantNameSize = 0;
            char *distantParamName = NULL;
            int   distantValueSize;
            char *distantValueCStr = NULL;

            MPI_Recv(&distantNameSize, 1, MPI_INT, irank, 0, localComm, &status);
            if (distantNameSize != nameSize)
              bftc_error(__FILE__, __LINE__, 0,
                        "Inconsistency between local parameters\n");

            distantParamName =  (char *) malloc (sizeof(char) * (distantNameSize+1));
            MPI_Recv(distantParamName, distantNameSize+1, MPI_CHAR, irank, 0, localComm, &status);
            if (strcmp(distantParamName, paramName))
              bftc_error(__FILE__, __LINE__, 0,
                        "Inconsistency between local parameters\n");
            free ( distantParamName);

            MPI_Recv(&distantValueSize, 1, MPI_INT, irank, 0, localComm, &status);
            if (distantValueSize != valueSize)
              bftc_error(__FILE__, __LINE__, 0,
                        "Inconsistency between local parameters\n");

            distantValueCStr =  (char *) malloc (sizeof(char) * (distantValueSize+1));
            MPI_Recv(distantValueCStr, distantValueSize+1, MPI_CHAR, irank, 0, localComm, &status);
            if (strcmp(distantValueCStr, valueCStr))
              bftc_error(__FILE__, __LINE__, 0,
                        "Inconsistency between local parameters\n");
            free ( distantValueCStr);
          }
        }

        else {
          MPI_Send(&nameSize, 1, MPI_INT, 0, 0, localComm);
          MPI_Send(const_cast <char* > (paramName), nameSize+1, MPI_CHAR, 0, 0, localComm);

          MPI_Send(&valueSize, 1, MPI_INT, 0, 0, localComm);
          MPI_Send(const_cast <char* > (valueCStr), valueSize+1, MPI_CHAR, 0, 0, localComm);
        }
      }
    }

    //
    // parameters exchange between beginning ranks


    if (currentRank == 0 ) {

      MPI_Sendrecv(&NLocalControlParameters,   1, MPI_INT, distantBeginningRank, 0,
                   &NDistantControlParameters, 1, MPI_INT, distantBeginningRank, 0,
                   globalComm, &status);

      int NParameterMax = MAX(NLocalControlParameters, NDistantControlParameters);

      Iterator p1 = localControlParameters.begin();

      for (int i = 0; i < NParameterMax; i++ ) {


        if (i >= NDistantControlParameters) {

          const char *valueCStr  = p1->second.c_str();
          int valueSize = p1->second.size();

          const char *paramName = p1->first.c_str();
          int nameSize          = p1->first.size();

          MPI_Send(&nameSize, 1, MPI_INT, distantBeginningRank, 0, globalComm);
          MPI_Send(const_cast <char *> (paramName), nameSize+1, MPI_CHAR, distantBeginningRank, 0, globalComm);

          MPI_Send(&valueSize, 1, MPI_INT, distantBeginningRank, 0, globalComm);
          MPI_Send(const_cast <char *> (valueCStr), valueSize+1, MPI_CHAR, distantBeginningRank, 0, globalComm);
        }

        else if (p1 == localControlParameters.end()){

          int   distantNameSize = 0;
          char *distantParamName = NULL;

          int   distantValueSize;
          char *distantValueCStr = NULL;

          MPI_Recv(&distantNameSize, 1, MPI_INT, distantBeginningRank, 0, globalComm, &status);
          distantParamName =  (char *) malloc (sizeof(char) * (distantNameSize+1));
          MPI_Recv(distantParamName, distantNameSize+1, MPI_CHAR, distantBeginningRank, 0, globalComm, &status);

          MPI_Recv(&distantValueSize, 1, MPI_INT, distantBeginningRank, 0, globalComm, &status);
          distantValueCStr =  (char *) malloc (sizeof(char) * (distantValueSize+1));
          MPI_Recv(distantValueCStr, distantValueSize+1, MPI_CHAR, distantBeginningRank, 0, globalComm, &status);

          distantControlParameters[std::string(distantParamName)] = std::string(distantValueCStr);

          free ( distantParamName);
          free ( distantValueCStr);

        }

        else {

          const char *paramName = p1->first.c_str();
          int nameSize          = p1->first.size();
          const char *valueCStr  = p1->second.c_str();
          int valueSize = p1->second.size();

          int   distantNameSize = 0;
          char *distantParamName = NULL;
          int   distantValueSize;
          char *distantValueCStr = NULL;

          MPI_Sendrecv(&nameSize       , 1, MPI_INT, distantBeginningRank, 0,
                       &distantNameSize, 1, MPI_INT, distantBeginningRank, 0,
                       globalComm, &status);

          distantParamName =  (char *) malloc (sizeof(char) * (distantNameSize+1));
          MPI_Sendrecv(const_cast <char*> (paramName), nameSize+1,        MPI_CHAR, distantBeginningRank, 0,
                       distantParamName              , distantNameSize+1, MPI_CHAR, distantBeginningRank, 0,
                       globalComm, &status);

          MPI_Sendrecv(&valueSize       , 1, MPI_INT, distantBeginningRank, 0,
                       &distantValueSize, 1, MPI_INT, distantBeginningRank, 0,
                       globalComm, &status);

          distantValueCStr =  (char *) malloc (sizeof(char) * (distantValueSize+1));
          MPI_Sendrecv(const_cast <char*> (valueCStr), valueSize+1,        MPI_CHAR, distantBeginningRank, 0,
                       distantValueCStr              , distantValueSize+1, MPI_CHAR, distantBeginningRank, 0,
                       globalComm, &status);

          distantControlParameters[std::string(distantParamName)] = std::string(distantValueCStr);
          free ( distantValueCStr);
          free ( distantParamName);

        }


        if (p1 != localControlParameters.end())
          p1++;
      }
    }

    //
    // Local beginning rank send parameters to other local ranks

    if (localCommSize > 1) {

      int size = 0;
      if (currentRank == 0)
        size =  distantControlParameters.size();

      MPI_Bcast(&size, 1, MPI_INT, 0, localComm);

      Iterator p1;
      if (currentRank == 0)
        p1 = distantControlParameters.begin();

      for (int i = 0 ; i < size; i++) {
        char *paramName = NULL;
        int nameSize;
        char *valueCStr = NULL;
        int valueSize;

        if (currentRank == 0) {
          valueCStr = const_cast <char *> (p1->second.c_str());
          valueSize = p1->second.size() ;

          paramName = const_cast <char *> (p1->first.c_str());
          nameSize  = p1->first.size();
          p1++;
        }

        MPI_Bcast(&nameSize, 1, MPI_INT, 0, localComm);

        if (currentRank != 0)
          paramName =  (char *) malloc (sizeof(char) * (nameSize+1));

        MPI_Bcast(paramName, nameSize+1, MPI_CHAR, 0, localComm);

        MPI_Bcast(&valueSize, 1, MPI_INT, 0, localComm);

        if (currentRank != 0)
          valueCStr =  (char *) malloc (sizeof(char) * (valueSize+1));

        MPI_Bcast(valueCStr, valueSize+1, MPI_CHAR, 0, localComm);

        if (currentRank != 0) {
          distantControlParameters[std::string(paramName)] = std::string(valueCStr);
          free ( paramName);
          free ( valueCStr);
        }
      }
    }
  }

  void ApplicationPropertiesDataBase::dump()
  {
    bftc_printf("\nLocal application properties\n\n");
    _localApplicationProperties->dump();

    typedef std::map <std::string, ApplicationProperties *>::iterator Iterator;

    bftc_printf("\nDistant application properties\n\n");
    for (Iterator p = _distantApplicationPropertiesDataBase.begin(); p != _distantApplicationPropertiesDataBase.end(); p++){
      p->second->dump();
      bftc_printf("\n");
    }
    bftc_printf_flush();
  }

}
