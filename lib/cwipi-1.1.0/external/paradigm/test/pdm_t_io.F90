program pdm_io_fortran
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! Test du bon fonctionnement des procedures :
  ! PDM_io_open
  ! PDM_io_global_write
  ! PDM_io_par_block_write
  ! PDM_io_par_interlaced_write
  ! PDM_io_global_read
  ! PDM_io_par_block_read
  ! PDM_io_par_interlaced_read
  ! PDM_io_close
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  use iso_c_binding
  use mpi
  use pdm
  use pdm_io
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  implicit none
  integer(4)                           :: rank,size,iErr
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call MPI_INIT(iErr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, size, iErr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, iErr)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call LectureEcriture(ENDEAN=PDM_IO_LITTLEENDIAN)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call LectureEcriture(ENDEAN=PDM_IO_BIGENDIAN   )
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  call MPI_FINALIZE(iErr)
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
contains
  
  subroutine LectureEcriture(ENDEAN)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    integer(kind=4          )            :: ENDEAN
    !>
    character(80)                        :: name
    type(c_ptr)                          :: id
    integer                              :: nombeEntier
    real(8)                              :: nombreReel
    complex(8)                           :: nombreComplex
    character(len=80)                    :: message
    integer(4), pointer                  :: ptr_int(:)
    real(8), pointer                     :: ptr_r8 (:)
    complex(8), pointer                  :: ptr_c8 (:)
    character(len=:), pointer            :: buffer=>null()
    ! type(c_ptr)                          :: cptr
    integer(kind=4)                      :: s_data
    integer(kind=pdm_g_num_s)            :: n_data
    character(80), pointer               :: lignesBlock(:)
    character(80), pointer               :: lignesEntre(:)
    logical                              :: consistance
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    consistance=.true.
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  
    if( rank==0 )then
      if( rank==0 )print '(/)'
      select case(ENDEAN)
      case(PDM_IO_BIGENDIAN)    ; print '(//"PARAGIDM IO : Test BIGENDEAN")'
      case(PDM_IO_LITTLEENDIAN) ; print '(//"PARAGIDM IO : Test LITTLEENDEAN")'
      end select
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    nombeEntier=1
    nombreReel=dacos(-1d0)
    nombreComplex=cmplx(nombreReel,nombreReel,kind=8)
    write(message,'("Ecriture Globale depuis le rank: ",i3)')rank
    
    allocate(character(len=80) :: buffer)
    allocate(ptr_int(1:1))
    allocate(ptr_r8 (1:1))
    allocate(ptr_c8 (1:1))
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Ecriture
    if( rank==0 )print '(/"Tests d''ecriture")'
    
    name="toto.dat"
    call PDM_io_open(                    &
    &    nom=trim(name)                 ,& !> nom       <-- Nom du fichier
    &    fmt=PDM_IO_FMT_BIN             ,& !> fmt       <-- Fichier text ou binaire
    &    suff_t=PDM_IO_SUFF_MAN         ,& !> suff_t    <-- Type de suffixe (manuel ou automatique)
    &    suff_u=""                      ,& !> suff_u    <-- Suffixe (si suffixe manuel)
    &    s_backup=PDM_IO_BACKUP_OFF     ,& !> s_backup  <-- Active le backup d'un fichier preexistant en mode ecriture
    &    acces=PDM_IO_KIND_MPIIO_EO     ,& !> accesio   <-- Type (parallele avec mpiio, parallele sans mpiio, sequentiel)
    &    mode=PDM_IO_MOD_WRITE          ,& !> mode      <-- Mode d'acces (lecture, ecriture, lecture/ecriture)
    &    endian=ENDEAN                  ,& !> PDM_IO_LITTLEENDIAN/PDM_IO_BIGENDIAN           
    &    comm=MPI_COMM_WORLD            ,& !> msg_comm  <-- Communicateur lie au fichier
    &    prop_noeuds_actifs=-1d0        ,& !> - 1: tous les rangs participent  +1: seuls les rangs maitres de chaque noeud participent
    &    unite=id                       ,& !> unite     --> Unite du fichier
    &    ierr=iErr                       ) !> ierr      --> Indique si le fichier est de type PDM_io ou non Utiliser uniquement pour une ouverture en lecture
    
    !>>> Ecriture Globale
    if( rank==0 )print '(3x,"Ecriture Globale")'
    
    ptr_int(1)=nombeEntier
    n_data=1
    s_data=4
    call PDM_io_global_write(id,s_data,n_data,c_loc(ptr_int))
    if( rank==0 )print '(6x,"ptr_int(1)=",i0)',ptr_int(1)
    
    ptr_r8(1)=nombreReel
    n_data=1
    s_data=8
    call PDM_io_global_write(id,s_data,n_data,c_loc(ptr_r8))
    if( rank==0 )print '(6x,"ptr_r8 (1)=",e22.15)',ptr_r8(1)
    
    ptr_c8(1)=nombreComplex
    n_data= 1
    s_data=16
    call PDM_io_global_write(id,s_data,n_data,c_loc(ptr_c8))
    if( rank==0 )print '(6x,"ptr_c8 (1)=",e22.15,1x,e22.15)',ptr_c8(1)
    
    buffer=message
    n_data= 1
    s_data=80
    call PDM_io_global_write(id,s_data,n_data,c_loc(buffer))
    if( rank==0 )print '(6x,"buffer=""",a,"""")',buffer
    !<<< Ecriture Globale
    
    !>>> Ecriture Blocs
    block
      integer                              :: iRank
      integer                              :: iLine,nLines
      integer, pointer                     :: nLinesRank(:)
      character(80)                        :: ligne
      !character(80), pointer               :: lignes(:)
      integer(4), pointer                  :: iTab(:)
      integer(kind=pdm_g_num_s)            :: shift
      
      if( rank==0 )print '(3x,"Ecriture Blocs")'
      
      nLines=4+rank !> blocs de taille variable
      allocate(lignesBlock(1:nLines))
      do iLine=1,nLines
        write(ligne,'("Ecriture Bloc Rank: ",i3,2x,"ligne: ",i3)')rank,iLine
        !ligne(80:80)=C_NEW_LINE
        lignesBlock(iLine)=ligne
      enddo
      
      allocate(nLinesRank(0:size-1))
      call mpi_allgather(                   &
      &    nLines     , 1,mpi_integer      ,&
      &    nLinesRank , 1,mpi_integer      ,&
      &    MPI_COMM_WORLD                  ,&
      &    iErr                             )
      shift=sum([(nLinesRank(iRank),iRank=0,rank-1)])+1    !> debut_bloc
      
      allocate(iTab(1:1)) ; iTab(1)=1                    !> <=
      call PDM_io_par_block_write(                    &
      &    fichier=id                                ,&  !> unite
      &    t_n_composantes=PDM_STRIDE_CST_INTERLACED ,&  !> t_n_composantes
      &    n_composantes=iTab                        ,&  !> n_composantes
      &    taille_donnee=80                          ,&  !> taille_donnee
      &    n_donnees=nLines                          ,&  !> n_donnees
      &    debut_bloc=shift                          ,&  !> debut_bloc
      &    donnees=c_loc(lignesBlock)                 )  !> donnees
      deallocate(iTab)
      
      deallocate(nLinesRank)
      !deallocate(lignes)
    end block
    !<<< Ecriture Blocs
    
    !>>> Ecriture Entrelacee
    block
      ! integer                              :: iRank
      integer                              :: iLine,nLines
      character(80)                        :: ligne
      integer(kind=pdm_g_num_s), pointer   :: indirection(:)
      integer, pointer                     :: iTab(:)
      
      if( rank==0 )print '(3x,"Ecriture Entrelacee")'
      
      nLines=4
      allocate(lignesEntre(1:nLines))
      do iLine=1,nLines
        write(ligne,'("Ecriture Entrelacee Rank: ",i3,2x,"ligne: ",i3)')rank,iLine
        lignesEntre(iLine)=ligne
      enddo
      
      !> indirection
      allocate(indirection(1:nLines))
      do iLine=1,nLines
        indirection(iLine)=rank+(iLine-1)*size +1 !> Rangement Rank/iLine      
      enddo
      
      !do iRank=0,size-1
      !  if( rank==iRank )then
      !    print '(3x,"rank= ",i0)',rank
      !    do iLine=1,nLines
      !      print '(6x,"Indirection sur rank",i0,": ",i3," -> ",i3)',rank,iLine,indirection(iLine)
      !    enddo
      !  endif
      !  call mpi_barrier(MPI_COMM_WORLD,iErr)
      !enddo
      
      allocate( iTab(1:1) ) ; iTab(1)=1                 !> <=
      call PDM_io_par_interlaced_write(               &
      &    fichier        =id                        ,& !> unite
      &    t_n_composantes=PDM_STRIDE_CST_INTERLACED ,& !> t_n_composantes
      &    n_composantes  =iTab                      ,& !> *n_composantes
      &    taille_donnee  =80                        ,& !> taille_donnee
      &    n_donnees      =nLines                    ,& !> n_donnees
      &    indirection    =indirection               ,& !> *indirection
      &    donnees        =c_loc(lignesEntre)         ) !> *donnees
      deallocate(iTab)
      
      deallocate(indirection)
      
      !do iRank=0,size-1
      !  if( rank==iRank )then
      !    print '("rank= ",i0)',rank
      !    do iLine=1,nLines
      !      print '(3x,"Ecrit sur rank",i0,": """,a,"""")',rank,lignesEntre(iLine)
      !    enddo
      !  endif
      !  call mpi_barrier(MPI_COMM_WORLD,iErr)
      !enddo
    end block
    !<<< Ecriture Entrelacee
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
  
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Lecture
    if( rank==0 )print '(/"Tests de lecture")'
    
    call PDM_io_open(                    &
    &    nom=trim(name)                 ,& !> nom       <-- Nom du fichier
    &    fmt=PDM_IO_FMT_BIN             ,& !> fmt       <-- Fichier text ou binaire
    &    suff_t=PDM_IO_SUFF_MAN         ,& !> suff_t    <-- Type de suffixe (manuel ou automatique)
    &    suff_u=""                      ,& !> suff_u    <-- Suffixe (si suffixe manuel)
    &    s_backup=PDM_IO_BACKUP_OFF     ,& !> s_backup  <-- Active le backup d'un fichier preexistant en mode ecriture
    &    acces=PDM_IO_KIND_MPIIO_EO     ,& !> accesio   <-- Type (parallele avec mpiio, parallele sans mpiio, sequentiel)
    &    mode=PDM_IO_MOD_READ           ,& !> mode      <-- Mode d'acces (lecture, ecriture, lecture/ecriture)
    &    endian=ENDEAN                  ,& !> PDM_IO_LITTLEENDIAN/PDM_IO_BIGENDIAN           
    &    comm=MPI_COMM_WORLD            ,& !> msg_comm  <-- Communicateur lie au fichier
    &    prop_noeuds_actifs=-1d0        ,& !> - 1: tous les rangs participent  +1: seuls les rangs maitres de chaque noeud participent
    &    unite=id                       ,& !> unite     --> Unite du fichier
    &    ierr=iErr                       ) !> ierr      --> Indique si le fichier est de type PDM_io ou non Utiliser uniquement pour une ouverture en lecture
    
    !>>> Lecture Globale
    if( rank==0 )print '(3x,"Lecture Globale")'
    
    n_data=1
    s_data=4
    call PDM_io_global_read(id,s_data,n_data,c_loc(ptr_int))
    
    if( .not. ptr_int(1)==nombeEntier )then
      consistance=.false.
      if( rank==0 )print '(/6x,"nombreEntier=",i0)',nombeEntier
      if( rank==0 )print '( 6x,"ptr_int(1)  =",i0)',ptr_int(1)
    endif
    
    n_data=1
    s_data=8
    call PDM_io_global_read(id,s_data,n_data,c_loc(ptr_r8))
    if( .not. ptr_r8(1)==nombreReel )then
      consistance=.false.
      if( rank==0 )print '(/6x,"nombreReel=",e22.15)',nombreReel
      if( rank==0 )print '( 6x,"ptr_r8 (1)=",e22.15)',ptr_r8(1)
    endif
    
    n_data= 1
    s_data=16
    call PDM_io_global_read(id,s_data,n_data,c_loc(ptr_c8))
    if( .not. ptr_c8(1)==nombreComplex )then
      consistance=.false.
      if( rank==0 )print '(/6x,"nombreComplex=",e22.15,1x,e22.15)',nombreComplex
      if( rank==0 )print '( 6x,"ptr_c8 (1)   =",e22.15,1x,e22.15)',ptr_c8(1)
    endif
    
    n_data= 1 
    s_data=80
    call PDM_io_global_read(id, s_data, n_data, c_loc(buffer))
    if( .not. buffer==message )then
      consistance=.false.
      if( rank==0 )print '(/6x,"message: """,a,"""")',message
      if( rank==0 )print '( 6x,"buffer : """,a,"""")',buffer
    endif
    !<<< Lecture Globale
    
    !>>> Lecture Blocs
    block
      integer                              :: iRank
      integer                              :: iLine,nLines
      integer, pointer                     :: nLinesRank(:)
      ! character(80)                        :: ligne
      character(80), pointer               :: lignes(:)
      integer(kind=pdm_g_num_s)            :: shift
      integer, pointer                     :: iTab(:)
      
      if( rank==0 )print '(/3x,"Lecture Blocs")'
      
      nLines=4+rank !> blocs de taille variable
      
      allocate(nLinesRank(0:size-1))
      call mpi_allgather(                   &
      &    nLines     , 1,mpi_integer      ,&
      &    nLinesRank , 1,mpi_integer      ,&
      &    MPI_COMM_WORLD                  ,&
      &    iErr                             )
      shift=sum([(nLinesRank(iRank),iRank=0,rank-1)])+1 !> debut_bloc
      
      allocate(lignes(1:nLines))
      
      allocate(iTab(1:1)) ; iTab(1)=1                    !> <=
      call PDM_io_par_block_read(                     &
      &    fichier=id                                ,&  !> unite
      &    t_n_composantes=PDM_STRIDE_CST_INTERLACED ,&  !> t_n_composantes
      &    n_composantes=iTab                        ,&  !> n_composantes
      &    taille_donnee=80                          ,&  !> taille_donnee
      &    n_donnees=nLines                          ,&  !> n_donnees
      &    debut_bloc=shift                          ,&  !> debut_bloc
      &    donnees=c_loc(lignes)                      )  !> donnees
      deallocate(iTab)
      
      do iRank=0,size-1
        if( rank==iRank )then
          print '(6x,"rank= ",i0)',rank
          do iLine=1,nLines
            if( .not.lignesBlock(iLine)==lignes(iLine) )then
              consistance=.false.
              print '(/9x,"Ecrit sur rank",i0,": """,a,"""")',rank,lignesBlock(iLine)
              print '( 9x,"Lu    sur rank",i0,": """,a,"""")',rank,lignes     (iLine)
            endif
          enddo
        endif
        call mpi_barrier(MPI_COMM_WORLD,iErr)
      enddo
      
      deallocate(lignes)
      lignes=>null()
      
    end block
    !<<< Lecture par bloc
    
    !>>> Lecture Entrelacee
    block
      integer                              :: iRank
      integer                              :: iLine,nLines
      ! character(80)                        :: ligne
      character(80), pointer               :: lignes(:)
      integer (kind = pdm_g_num_s), pointer:: indirection(:)
      integer, pointer                     :: iTab(:)
      ! type(c_ptr)                          :: cptr
      
      if( rank==0 )print '(/3x,"Lecture Entrelacee")'
      
      nLines=4
      
      !> indirection
      allocate(indirection(1:nLines))
      do iLine=1,nLines
        indirection(iLine)=rank+(iLine-1)*size +1 !> Rangement Rank/iLine      
      enddo
      
      allocate(lignes(1:nLines))
      
      allocate( iTab(1:1) ) ; iTab(1)=1                 !> <=
      call PDM_io_par_interlaced_read(                &
      &    fichier        =id                        ,& !> unite
      &    t_n_composantes=PDM_STRIDE_CST_INTERLACED ,& !> t_n_composantes
      &    n_composantes  =iTab                      ,& !> *n_composantes
      &    taille_donnee  =80                        ,& !> taille_donnee
      &    n_donnees      =nLines                    ,& !> n_donnees
      &    indirection    =indirection               ,& !> *indirection
      &    donnees        =c_loc(lignes)              ) !> *donnees
      deallocate(iTab)
      
      deallocate(indirection)
      
      do iRank=0,size-1
        if( rank==iRank )then
          print '(6x,"rank= ",i0)',rank
          do iLine=1,nLines
            if( .not.lignesEntre(iLine)==lignes(iLine) )then
              consistance=.false.
              print '(/9x,"Ecrit sur rank",i0,": """,a,"""")',rank,lignesBlock(iLine)
              print '( 9x,"Lu    sur rank",i0,": """,a,"""")',rank,lignes     (iLine)
            endif
          enddo
        endif
        call mpi_barrier(MPI_COMM_WORLD,iErr)
      enddo
      
      deallocate(lignes)
      lignes=>null()
    end block
    !<<< Lecture Entrelacee
    
    call PDM_io_close(id)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( consistance )then
      select case(ENDEAN)
      case(PDM_IO_BIGENDIAN)    ; if( rank==0 )print '(/"BIGENDEAN Lecture/Ecriture paralleles CONSISTANTES")'
      case(PDM_IO_LITTLEENDIAN) ; if( rank==0 )print '(/"LITTLEENDIAN Lecture/Ecriture paralleles CONSISTANTES")'
      end select
    else
      select case(ENDEAN)
      case(PDM_IO_BIGENDIAN)    ; if( rank==0 )print '(/"PROBLEME :  BIGENDEAN Lecture/Ecriture paralleles INCONSISTANTES")'
      case(PDM_IO_LITTLEENDIAN) ; if( rank==0 )print '(/"PROBLEME :  LITTLEENDIAN Lecture/Ecriture paralleles INCONSISTANTES")'
      end select
    endif
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !> Lecture Fortran sur proc0
  
    !block
    !  integer       :: i,iUnit
    !  character(1)  :: char
    !  integer       :: cpt,iChar
    !  character(80) :: ligne
    !  integer       :: Reason
    !  if( rank==0 )then
    !    print '(/"Lecture Fortran sur rank 0")'
    !    
    !    select case(ENDEAN)
    !    case(PDM_IO_BIGENDIAN)
    !      open(newunit=iUnit                        ,&
    !      &    file=trim(name)                      ,&
    !      &    access='stream'                      ,&
    !      &    form='unformatted'                   ,&
    !      &    convert="BIG_ENDIAN"                 ,&
    !      &    action='read'                         )
    !    case(PDM_IO_LITTLEENDIAN)
    !      open(newunit=iUnit                        ,&
    !      &    file=trim(name)                      ,&
    !      &    access='stream'                      ,&
    !      &    form='unformatted'                   ,&
    !      &    action='read'                         )
    !    end select
    !
    !    read(iUnit)ptr_int(1) ; print '(3x,"ptr_int(1)=",i0          )',ptr_int(1)
    !    read(iUnit)ptr_r8 (1) ; print '(3x,"ptr_r8 (1)=",e22.15      )',ptr_r8 (1)
    !    read(iUnit)ptr_c8 (1) ; print '(3x,"ptr_c8 (1)=",2(e22.15,1x))',ptr_c8 (1)
    !    read(iUnit)buffer     ; print '(3x,"buffer: """,a,""""       )',buffer
    !    
    !    cpt=0
    !    iChar=0
    !    lecture: do
    !      read(iUnit,iostat=Reason)char
    !      if( Reason>0 )then
    !       stop "probleme de lecture"
    !      elseif( Reason<0 )then
    !       exit lecture
    !      endif
    !      cpt=cpt+1
    !      
    !      iChar=iChar+1
    !      !if( char==C_NEW_LINE )then
    !      if( iChar==80 )then
    !        ligne(iChar:iChar)=char
    !        print '(3x,"lecture fortran  iChar=",i0," ligne: """,a,"""")',iChar,ligne
    !        ligne=""
    !        iChar=0
    !      else
    !        ligne(iChar:iChar)=char
    !      endif
    !    enddo lecture
    !    print '("cpt=",i0)',cpt
    !    
    !    close(iUnit)
    !  endif
    !end block
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    deallocate(ptr_int)
    deallocate(ptr_r8 )
    deallocate(ptr_c8 )
    deallocate(lignesBlock)
    deallocate(lignesEntre)
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  end  subroutine LectureEcriture




end program pdm_io_fortran
