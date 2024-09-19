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

#include "commWithoutPart.hxx"
#include "pdm_logging.h"

using namespace std;

namespace cwipi {

  /**
   *
   * \brief Constructor.
   *
   */
  
  CommWithoutPart::CommWithoutPart()
    : Communication::Communication(), 
      _commType(CWP_COMM_PAR_WITHOUT_PART)
  {
  }

  /**
   * \brief Destructor.
   *
   */

  CommWithoutPart::~CommWithoutPart()
  {
  }

  /**
   *
   * \brief Building coupling communicator.
   *
   */

  void 
  CommWithoutPart::_cplCommCreate
  (
   CWP_Comm_t cplCodeCommType
  )
  {
    const int localRootRank = _localCodeProperties->rootRankGet();
    const int cplRootRank   = _cplCodeProperties->rootRankGet();
    
    const MPI_Comm& globalComm = _localCodeProperties->globalCommGet();

    int currentRank;
    MPI_Comm_rank(globalComm, &currentRank);

    //_isCplRank = localRootRank == currentRank;

    vector <int> cplRanks = *(_cplCodeProperties->connectableRanksGet());
    vector <int> locRanks = *(_localCodeProperties->connectableRanksGet());

    _unionCommCplRanks = new std::vector<int>(cplRanks.size());
    _unionCommLocRanks = new std::vector<int>(locRanks.size());


    //Traduction de ces rangs dans _unionCommCplRanks - Contient les rangs (unionComm) correspondants

    MPI_Group globalGroup;
    MPI_Comm_group(_localCodeProperties->globalCommGet(), &globalGroup);
      
    MPI_Group unionGroup;
    MPI_Comm_group(_unionComm, &unionGroup);      

    MPI_Group_translate_ranks(globalGroup, cplRanks.size(), &(cplRanks[0]),
                              unionGroup, &((*_unionCommCplRanks)[0]));

    MPI_Group_translate_ranks(globalGroup, locRanks.size(), &(locRanks[0]),
                              unionGroup , &((*_unionCommLocRanks)[0]));

    _cplCommCplRanks = new std::vector<int>(*_unionCommCplRanks);
    _cplCommLocRanks = new std::vector<int>(*_unionCommLocRanks);
      
    if (cplCodeCommType != CWP_COMM_PAR_WITH_PART) {

      int cplRanks2[2];
      int gap1 = 0;
      int gap2 = 1;
      
      if (_localCodeProperties->idGet() < _cplCodeProperties->idGet()) {
        gap1 = 1;
        gap2 = 0;
      }
      
      MPI_Group_translate_ranks (globalGroup, 1, &localRootRank, 
                                 unionGroup, cplRanks2 + gap1);

      MPI_Group_translate_ranks (globalGroup, 1, &cplRootRank, 
                                 unionGroup, cplRanks2 + gap2);
      
      MPI_Group_incl(unionGroup, 2, cplRanks2, &_cplGroup);
      
      MPI_Comm_create (_unionComm, _cplGroup, &_cplComm);
      
    }
  
    else {

      const int locRootRankInGlobalComm = _localCodeProperties->rootRankGet();
      int locRootRankInUnionComm;
      MPI_Group_translate_ranks(globalGroup, 1, &locRootRankInGlobalComm,
                                unionGroup, &locRootRankInUnionComm);
      int *excludeRanks = (int *) malloc(sizeof(int) * (locRanks.size() - 1));
      int j = 0;
      for (size_t i = 0; i < locRanks.size(); i++) {
        if ((*_unionCommLocRanks)[i] != locRootRankInUnionComm) {
          excludeRanks[j++] = (*_unionCommLocRanks)[i];
        }
      }

      MPI_Group_excl(unionGroup, locRanks.size()-1, excludeRanks, &_cplGroup);
      free(excludeRanks);

      MPI_Comm_create(_unionComm, _cplGroup, &_cplComm);
      
    }
  }

  /**
   *
   * \brief Synchronise
   *
   */
/*
  void
  CommWithoutPart::sync
  (
   void *tab, 
   MPI_Datatype mpiType, 
   int tabSize
  )
  {
    MPI_Bcast(tab,
              tabSize,
              mpiType,
              0,
              _localCodeProperties->intraCommGet());

  }*/

}
