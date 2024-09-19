/*
  This file is part of the CWIPI library. 

  Copyright (C) 2021-2023  ONERA

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
#include <list>
#include <map>
#include <algorithm>    
#include <vector>       

#include "pdm_error.h"

#include <fvmc_parall.h>

#include "cwipi_config.h"
#include "codePropertiesDB.hxx"
#include "codePropertiesDB_i.hxx"

using namespace std;

namespace cwipi {

  const int CodePropertiesDB::_nIssend = 2;

  /**
   * \brief Default constructor.
   *
   */

  CodePropertiesDB::CodePropertiesDB()
    : _codePropertiesDB(*(new map <string, CodeProperties * > ())),
      _locCodePropertiesDB(*(new map <string, CodeProperties * > ()))
  {
  }

  /**
   * \brief Destructor.
   *
   */

  CodePropertiesDB::~CodePropertiesDB()
  {
#if defined(DEBUG) && 0
    cout << "destroying CodePropertiesDB." << endl;
#endif

    if (!_codePropertiesDB.empty()) {
      typedef map <string, CodeProperties * >::iterator CI;
      for (CI p = _codePropertiesDB.begin();
           p != _codePropertiesDB.end(); p++) {
        if (p->second != NULL)
          delete p->second;
      }
      _codePropertiesDB.clear();

    }
    delete &_codePropertiesDB;
    delete &_locCodePropertiesDB;
    
  }

  
  /**
   * \brief MPI Communicator Initialization.
   *
   * This function builds the current code intra-communicator from
   * the current name and the MPI communicator containing all processes of
   * all codes.
   *
   * \param [in]  globalComm      MPI communicator containing all processes 
   *                              of all codes
   * \param [in]  n_codes         Number of codes on the current rank
   * \param [in]  code_names      Codes names on the current rank  (n_codes)
   * \param [in]  is_active_rank  Current rank is active
   * \param [in]  n_param_max     Maximum number of parameters
   * \param [in]  str_size_max    Maximum size for a string
   * \param [out] intra_comms      Current codes intra-communicators  (n_codes)
   *
   */

  void 
  CodePropertiesDB::init
  (
  const MPI_Comm     globalComm,
  const int          n_codes,
  const char**       code_names, 
  const CWP_Status_t is_active_rank,
  const int          n_param_max,
  const int          str_size_max,      
  MPI_Comm           *intra_comms       
 )
  {

    // Initialize MPI
    // --------------
    
    int currentRank;
    int globalCommSize;
    
    _globalComm = globalComm;
    
    _n_param_max = n_param_max;
    _str_size_max = str_size_max;

    int flag;

    MPI_Initialized(&flag);
    if (!flag) {
      PDM_error(__FILE__, __LINE__, 0, "MPI is not initialized\n");
    }
      
    MPI_Comm_rank(globalComm, &currentRank);
    MPI_Comm_size(globalComm, &globalCommSize);

    // Search codes
    // ------------

    int index = 0;
    int properties[2]; /* The first property is the number of rank codes
                        * The second is the length of concatenated code names 
                        * The third is coupled rank state */

    properties[0] = n_codes;
    properties[1] = 0;
    for (int i = 0; i < n_codes; i++) {
      properties[1] += strlen(code_names[i]) + 1;
    }
    //properties[2] = is_active_rank;

    int *allProperties =  (int *) malloc (sizeof(int) * (2*globalCommSize));

    MPI_Allgather(&properties,
                  2,
                  MPI_INT,
                  allProperties,
                  2,
                  MPI_INT,
                  globalComm);

    char *concatenateNames =  (char *) malloc (sizeof(char) * (properties[1]));
    char *_concatenateNames = concatenateNames;
    for (int i = 0; i < n_codes; i++) {
      strcpy(_concatenateNames, code_names[i]);
      _concatenateNames += strlen(code_names[i]) + 1;
    }

    int *iproc =  (int *) malloc (sizeof(int) * (globalCommSize + 1));
    int *codesLengthName =  (int *) malloc (sizeof(int) * (globalCommSize));

    int *iproc2 =  (int *) malloc (sizeof(int) * (globalCommSize + 1));
    int *n_codes_rank =  (int *) malloc (sizeof(int) * (globalCommSize));
    
    iproc[0] = 0;
    iproc2[0] = 0;
    for(int i = 0; i < globalCommSize; i++) {
      codesLengthName[i] = allProperties[2*i+1];
      n_codes_rank[i] = allProperties[2*i];
      iproc[i+1] = allProperties[2*i+1] + iproc[i];
      iproc2[i+1] = allProperties[2*i] + iproc2[i];
    }

    int totalLength = iproc[globalCommSize];

    char *mergeNames     =  (char *) malloc (sizeof(char) * totalLength);
    int  *mergeIsCoupled =  (int  *) malloc (sizeof(int ) * globalCommSize);
    
    MPI_Allgatherv((void*) concatenateNames,
                   properties[1],
                   MPI_CHAR,
                   mergeNames,
                   codesLengthName,
                   iproc,
                   MPI_CHAR,
                   globalComm);

    MPI_Allgather((void*) &is_active_rank,
                   1,
                   MPI_INT,
                   mergeIsCoupled,
                   1,
                   MPI_INT,
                   globalComm);

    free ( concatenateNames);
    free ( codesLengthName);



    map < string, vector < int > * >  coupledRankCode;
    map < string, vector < int > * >  rankCode;
    
    int id = 1;
     
    for (int irank = 0; irank < globalCommSize; irank++) {   
      
      for (int k = 0; k < allProperties[2*irank]; k++) {
        assert(index <= totalLength);

        const char *ptCurrentName = mergeNames + index;

        string currentName = string(ptCurrentName);

        typedef map <string, CodeProperties *>::iterator Iterator;

        Iterator p = _codePropertiesDB.find(currentName);

        if (p == _codePropertiesDB.end()) {

          CodeProperties *currentCodeProperties =
            new CodeProperties(currentName, id++, irank, 
                               currentRank == irank, globalComm,
                               _n_param_max, _str_size_max);

          coupledRankCode[currentName] = new vector <int> ();
          coupledRankCode[currentName]->reserve(globalCommSize);
          rankCode[currentName] = new vector <int> ();
          rankCode[currentName]->reserve(globalCommSize);

          pair<string, CodeProperties *>
            newPair(currentName, currentCodeProperties);

          _codePropertiesDB.insert(newPair);

          if (currentRank == irank) {
            _locCodePropertiesDB.insert(newPair);
          }
        }

        else {
          
          if (currentRank == irank) {
            if (_locCodePropertiesDB.find(currentName) == _locCodePropertiesDB.end()) {
               _locCodePropertiesDB[currentName] = p->second;
               p->second->isLocalSet(true);
            }
          }
        
        }
        
        if (mergeIsCoupled[irank]) {
          coupledRankCode[currentName]->push_back(irank);
        }
      
        rankCode[currentName]->push_back(irank);

        if (currentRank == irank) {
          _codePropertiesDB[currentName]->isActiveRankset((bool) is_active_rank);
        }

        index += currentName.size() + 1;
        assert(index <= totalLength);
      }
    }

    // Create intra code communicators
    // -------------------------------

    MPI_Group globalGroup;
    MPI_Comm_group (globalComm, &globalGroup);
    
    typedef map <string, vector < int > * >::iterator Iterator;

    for (Iterator p = rankCode.begin(); 
                  p != rankCode.end(); p++){
      p->second->resize(p->second->size());
      const int *_ranks = &((*p->second)[0]);
      _codePropertiesDB[p->first]->_intraRanks = p->second;
      int _n_ranks = p->second->size();

      MPI_Group_incl (globalGroup, _n_ranks, _ranks, 
                     &(_codePropertiesDB[p->first]->_intraGroup));
      MPI_Comm_create (globalComm, 
                       _codePropertiesDB[p->first]->_intraGroup,
                       &(_codePropertiesDB[p->first]->_intraComm));

      int rootIdx = 0;
      std::vector<int> active_ranks = (*coupledRankCode[p->first]);
      /*std::vector<int> cplRankCode = (*coupledRankCode[p->first]);
      vector<int>::iterator it = std::find( cplRankCode.begin(), cplRankCode.end(), _ranks[rootIdx] );
      
      while((n_codes_rank[ _ranks[rootIdx] ] == 2 
            || it == cplRankCode.end())
            && rootIdx <_n_ranks){
        rootIdx++;
        it =  std::find( cplRankCode.begin(), cplRankCode.end(), _ranks[rootIdx] );
      }
      if(rootIdx == _n_ranks)
      PDM_error(__FILE__, __LINE__, 0, "At least one MPI process per code must be monocode.\n");*/
      _codePropertiesDB[p->first]->_rootRankInGlobalComm = active_ranks[rootIdx];
    }

    rankCode.clear();



    typedef map <string, CodeProperties *>::iterator IteratorCP;
      
    for (IteratorCP p2  = _codePropertiesDB.begin(); 
                    p2 != _codePropertiesDB.end(); 
                    p2++) {
        
      if (p2->second->_rootRankInGlobalComm == currentRank) {
        
        MPI_Win_create(p2->second->_winGlobData, 
                       4 * sizeof(int),
                       sizeof(int), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winGlob);
        
        p2->second->_winIntParamIdxNameData = 
            (int *) malloc (sizeof(int) * (n_param_max + 1));
        p2->second->_winIntParamIdxNameData[0] = 0;
        MPI_Win_create(p2->second->_winIntParamIdxNameData, 
                       (n_param_max + 1) * sizeof(int),
                       sizeof(int), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winIntParamIdxName);
        
        p2->second->_winIntParamNameData = 
            (char *) malloc (sizeof(char) * n_param_max * str_size_max);
        MPI_Win_create(p2->second->_winIntParamNameData, 
                       n_param_max * str_size_max * sizeof(char),
                       sizeof(char), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winIntParamName);
        
        p2->second->_winIntParamValueData = 
            (int *) malloc (sizeof(int) * n_param_max);
        MPI_Win_create(p2->second->_winIntParamValueData, 
                       n_param_max * sizeof(int),
                       sizeof(int), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winIntParamValue);
        
        p2->second->_winDoubleParamIdxNameData = 
            (int *) malloc (sizeof(int) * (n_param_max + 1));
        p2->second->_winDoubleParamIdxNameData[0] = 0;
        MPI_Win_create(p2->second->_winDoubleParamIdxNameData, 
                       (n_param_max + 1) * sizeof(int),
                       sizeof(int), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winDoubleParamIdxName);
        
        p2->second->_winDoubleParamNameData = 
            (char *) malloc (sizeof(char) * n_param_max * str_size_max);
        MPI_Win_create(p2->second->_winDoubleParamNameData, 
                       n_param_max * str_size_max * sizeof(char),
                       sizeof(char), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winDoubleParamName);
        
        p2->second->_winDoubleParamValueData = 
            (double *) malloc (sizeof(double) * n_param_max);
        MPI_Win_create(p2->second->_winDoubleParamValueData, 
                       n_param_max * sizeof(double),
                       sizeof(double), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winDoubleParamValue);
        
        p2->second->_winStrParamIdxNameData = 
            (int *) malloc (sizeof(int) * (n_param_max + 1));
        p2->second->_winStrParamIdxNameData[0] = 0;
        MPI_Win_create(p2->second->_winStrParamIdxNameData, 
                       (n_param_max + 1) * sizeof(int),
                       sizeof(int), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winStrParamIdxName);
        
        p2->second->_winStrParamNameData = 
            (char *) malloc (sizeof(char) * n_param_max * str_size_max);
        MPI_Win_create(p2->second->_winStrParamNameData, 
                       n_param_max * str_size_max * sizeof(char),
                       sizeof(char), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winStrParamName);
        
        p2->second->_winStrParamIdxValueData = 
            (int *) malloc (sizeof(int) * (n_param_max + 1));
        p2->second->_winStrParamIdxValueData[0] = 0;
        MPI_Win_create(p2->second->_winStrParamIdxValueData, 
                       (n_param_max + 1) * sizeof(int),
                       sizeof(int), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winStrParamIdxValue);
        
        p2->second->_winStrParamValueData = 
            (char *) malloc (sizeof(char) * n_param_max * str_size_max);
        MPI_Win_create(p2->second->_winStrParamValueData, 
                       n_param_max * str_size_max * sizeof(char),
                       sizeof(char), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winStrParamValue);
      }
      else {
        MPI_Win_create(NULL, 
                       0, 
                       sizeof(int), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winGlob);

        p2->second->_winIntParamIdxNameData = 
            (int *) malloc (sizeof(int) * (n_param_max + 1));
        MPI_Win_create(NULL, 
                       0, 
                       sizeof(int), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winIntParamIdxName);

        p2->second->_winIntParamNameData = 
            (char *) malloc (sizeof(char) * n_param_max *str_size_max);
        MPI_Win_create(NULL, 
                       0, 
                       sizeof(char), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winIntParamName);

        p2->second->_winIntParamValueData = 
            (int *) malloc (sizeof(int) * n_param_max);
        MPI_Win_create(NULL, 
                       0, 
                       sizeof(int), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winIntParamValue);

        p2->second->_winDoubleParamIdxNameData = 
            (int *) malloc (sizeof(int) * (n_param_max + 1));
        MPI_Win_create(NULL, 
                       0, 
                       sizeof(int), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winDoubleParamIdxName);

        p2->second->_winDoubleParamNameData = 
            (char *) malloc (sizeof(char) * n_param_max *str_size_max);
        MPI_Win_create(NULL, 
                       0, 
                       sizeof(char), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winDoubleParamName);

        p2->second->_winDoubleParamValueData = 
            (double *) malloc (sizeof(double) * n_param_max);
        MPI_Win_create(NULL, 
                       0, 
                       sizeof(double), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winDoubleParamValue);

        p2->second->_winStrParamIdxNameData = 
            (int *) malloc (sizeof(int) * (n_param_max + 1));
        MPI_Win_create(NULL, 
                       0, 
                       sizeof(int), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winStrParamIdxName);
        
        p2->second->_winStrParamNameData = 
            (char *) malloc (sizeof(char) * n_param_max * str_size_max);
        MPI_Win_create(NULL, 
                       0, 
                       sizeof(char), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winStrParamName);
        
        p2->second->_winStrParamIdxValueData = 
            (int *) malloc (sizeof(int) * (n_param_max + 1));
        MPI_Win_create(NULL, 
                       0, 
                       sizeof(int), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winStrParamIdxValue);
        
        p2->second->_winStrParamValueData = 
            (char *) malloc (sizeof(char) * n_param_max * str_size_max);
        MPI_Win_create(NULL, 
                       0, 
                       sizeof(char), 
                       MPI_INFO_NULL, 
                       p2->second->_globalComm, 
                       &p2->second->_winStrParamValue);
      }      

    }
         
    free ( allProperties);
    free ( iproc);
    free ( iproc2);
    free ( mergeNames);
    free ( mergeIsCoupled);
    free ( n_codes_rank);
    
    
    // Create intra coupled code group
    // -------------------------------

    for (Iterator p = coupledRankCode.begin(); 
                  p != coupledRankCode.end(); p++){
      p->second->resize(p->second->size());
      const int *_ranks = &((*p->second)[0]);
      _codePropertiesDB[p->first]->_connectableRanks = p->second;
      int _n_ranks = p->second->size();
      MPI_Group_incl (globalGroup, _n_ranks, _ranks, 
                     &(_codePropertiesDB[p->first]->_intraConnectableGroup));
      if (MPI_GROUP_EMPTY == _codePropertiesDB[p->first]->_intraConnectableGroup) {
        printf ("[%d] Group empty : %s\n", currentRank, p->first.c_str()); 
      }
      MPI_Comm_create (globalComm, 
                       _codePropertiesDB[p->first]->_intraConnectableGroup,
                       &(_codePropertiesDB[p->first]->_intraConnectableComm));
    }

    coupledRankCode.clear();
    
    for (int i = 0; i < n_codes; i++) {
      const string &nameStr = code_names[i];
      intra_comms[i] = _codePropertiesDB[nameStr]->_intraComm;
    }

  }

  /**
   * \brief Dump properties.  
   *
   */

  void 
  CodePropertiesDB::dump()
  {
    typedef map <string, CodeProperties * >::iterator CI;
    for (CI p = _codePropertiesDB.begin();
         p != _codePropertiesDB.end(); p++) {
      if (p->second != NULL)
        p->second->dump();
    }
  }

  /**
    * \brief Dump string of properties
    *
    */

  string
  CodePropertiesDB::str_dump()
  {
    string out;
    typedef map <string, CodeProperties * >::iterator CI;
    for (CI p = _codePropertiesDB.begin();
         p != _codePropertiesDB.end(); p++) {
      if (p->second != NULL) {
        out.append(p->second->str_dump());
      }
    }
    return out;
  }

}
