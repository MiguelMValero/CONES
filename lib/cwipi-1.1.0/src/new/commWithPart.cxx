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

#include "commWithPart.hxx"

using namespace std;

namespace cwipi {

  /**
   *
   * \brief Constructor.
   *
   */

  CommWithPart::CommWithPart()
    : Communication::Communication(),
      _commType(CWP_COMM_PAR_WITH_PART)
  {
  }

  /**
   * \brief Destructor.
   *
   */

  CommWithPart::~CommWithPart()
  {
  }

  /**
   *
   * \brief Building coupling communicator.
   *
   */

  void
  CommWithPart::_cplCommCreate
  (
   CWP_Comm_t cplCodeCommType
  )
  {

    _isCplRank = true;

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

    /*
    Use the unionComm which the union of coupled ranks to compute spatial interpolation weights for partitioned and unpartitioned case.
    Use the intraComm to broadcast computation results in withOutPart case.

    _cplComm = unionComm in withPart case and contained the two code root ranks in withOutPart case.
    the _cplComm is


    */

    if (cplCodeCommType != CWP_COMM_PAR_WITH_PART) {

      //Exclusion des rangs connectable (unionComm) du code couplé pour obtenir le communicateur cplComm
      // contenant uniquement les rangs connectable du code couplé

      const int cplRootRankInGlobalComm = _cplCodeProperties->rootRankGet();
      int cplRootRankInUnionComm;
      MPI_Group_translate_ranks(globalGroup, 1, &cplRootRankInGlobalComm,
                                unionGroup, &cplRootRankInUnionComm);
      int *excludeRanks = (int *) malloc(sizeof(int) * (cplRanks.size() - 1));
      int j = 0;
      for (size_t i = 0; i < cplRanks.size(); i++) {
        if ((*_unionCommCplRanks)[i] != cplRootRankInUnionComm) {
          excludeRanks[j++] = (*_unionCommCplRanks)[i];
        }
      }

      MPI_Group_excl(unionGroup, cplRanks.size()-1, excludeRanks, &_cplGroup);
      free(excludeRanks);

      MPI_Comm_create(_unionComm, _cplGroup, &_cplComm);

    }
    else {

      _cplComm = _unionComm;
      MPI_Comm_group(_cplComm, &_cplGroup);

    }
  }

  /**
   *
   * \brief Synchronise
   *
   */

#if defined(__INTEL_COMPILER)
#pragma warning(push)
#pragma warning(disable:869)
#elif defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-value"
#endif
/*  void
  CommWithPart::sync
  (
   void *tab,
   MPI_Datatype mpiType,
   int tabSize
  )
  {
  }*/
#if defined(__INTEL_COMPILER)
#pragma warning(pop)
#elif defined(__clang__)
#pragma clang diagnostic pop
#endif
}
