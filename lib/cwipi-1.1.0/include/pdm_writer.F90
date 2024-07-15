module pdm_writer

  use pdm

  implicit none

  !
  ! Statut

  integer, parameter :: PDM_WRITER_OFF  = 0
  integer, parameter :: PDM_WRITER_ON   = 1

  !
  ! Type de topologie

  integer, parameter :: PDM_WRITER_TOPO_CST  = 0
  integer, parameter :: PDM_WRITER_TOPO_DEFORMABLE = 1
  integer, parameter :: PDM_WRITER_TOPO_VARIABLE   = 2

  !
  ! Type d'elements géometriques

  integer, parameter :: PDM_WRITER_POINT    = 0
  integer, parameter :: PDM_WRITER_BAR2     = 1
  integer, parameter :: PDM_WRITER_TRIA3    = 2
  integer, parameter :: PDM_WRITER_QUAD4    = 3
  integer, parameter :: PDM_WRITER_POLY_2D  = 4
  integer, parameter :: PDM_WRITER_TETRA4   = 5
  integer, parameter :: PDM_WRITER_PYRAMID5 = 6
  integer, parameter :: PDM_WRITER_PRISM6   = 7
  integer, parameter :: PDM_WRITER_HEXA8    = 8
  integer, parameter :: PDM_WRITER_POLY_3D  = 9

  !
  ! Format de sortie

  integer, parameter ::  PDM_WRITER_FMT_ENSIGHT = 0

  !
  ! Format du fichier

  integer, parameter ::  PDM_WRITER_FMT_BIN   = 0
  integer, parameter ::  PDM_WRITER_FMT_ASCII = 1

  !
  ! Dimension géométrique de la sortie

  integer, parameter ::  PDM_WRITER_DIM_2 = 0
  integer, parameter ::  PDM_WRITER_DIM_3 = 1

  !
  ! Dim des variables

  integer, parameter ::  PDM_WRITER_VAR_CSTE         = 0
  integer, parameter ::  PDM_WRITER_VAR_SCALAIRE     = 1
  integer, parameter ::  PDM_WRITER_VAR_VECTOR      = 3
  integer, parameter ::  PDM_WRITER_VAR_TENSEUR_SYM  = 6
  integer, parameter ::  PDM_WRITER_VAR_TENSEUR_ASYM = 9

  !
  ! Localisation des variables

  integer, parameter ::  PDM_WRITER_VAR_VERTICES      = 0
  integer, parameter ::  PDM_WRITER_VAR_ELEMENTS     = 1
  integer, parameter ::  PDM_WRITER_VAR_PARTICULES   = 2



  interface

  !>
  !! \brief Libere un objet CS (Cedre Sortie) et retourne un pointeur NULL si pas d'erreur
  !!
  !! \param [in] cs    Pointer to \ref PDM_writer object
  !!
  !!

  subroutine PDM_writer_free (cs) &
  bind (c, name='PDM_writer_free')

    use iso_c_binding
    implicit none

    type(c_ptr), value :: cs

  end subroutine PDM_writer_free


  !>
  !! \brief Fin d'increment
  !!
  !! \param [in] cs             Pointer to \ref PDM_writer object
  !!
  !!

  subroutine PDM_writer_step_end (cs) &
  bind (c, name="PDM_writer_step_end")

    use iso_c_binding
    implicit none

    type(c_ptr), value :: cs

  end subroutine PDM_writer_step_end


  !>
  !! \brief Ecriture du maillage courant
  !!
  !! \param [in] cs              Pointer to \ref PDM_writer object
  !! \param [in] id_geom         Identificateur de l'objet geometrique
  !!
  !!

  subroutine PDM_writer_geom_write (cs,      &
                                    id_geom) &
  bind (c, name='PDM_writer_geom_write')

    use iso_c_binding
    implicit none

    type(c_ptr),    value :: cs
    integer(c_int), value :: id_geom

  end subroutine PDM_writer_geom_write


  !>
  !! \brief Liberation des donnees decrivant le maillage courant
  !! les indirections sur les numérotation absolues sont conservées
  !!
  !! \param [in] cs              Pointer to \ref PDM_writer object
  !! \param [in] id_geom         Identificateur de l'objet geometrique
  !!
  !!

  subroutine PDM_writer_geom_data_free (cs,      &
                                        id_geom) &
  bind (c, name='PDM_writer_geom_data_free')

    use iso_c_binding
    implicit none

    type(c_ptr),    value :: cs
    integer(c_int), value :: id_geom

  end subroutine PDM_writer_geom_data_free


  !>
  !! \brief Liberation des donnees decrivant le maillage courant
  !!
  !! \param [in] cs              Pointer to \ref PDM_writer object
  !! \param [in] id_geom         Identificateur de l'objet geometrique
  !!
  !!

  subroutine PDM_writer_geom_free (cs,      &
                                   id_geom) &
  bind (c, name='PDM_writer_geom_free')

    use iso_c_binding
    implicit none

    type(c_ptr),    value :: cs
    integer(c_int), value :: id_geom
  end subroutine PDM_writer_geom_free


  !>
  !! \brief Ecriture des valeurs de la variable
  !!
  !! \param [in] cs              Pointer to \ref PDM_writer object
  !! \param [in] id_var          Identificateur de la variable a ecrire
  !!
  !!

  subroutine PDM_writer_var_write (cs,     &
                                   id_var) &
  bind (c, name='PDM_writer_var_write')

    use iso_c_binding
    implicit none

    type(c_ptr),    value :: cs
    integer(c_int), value :: id_var

  end subroutine PDM_writer_var_write


  !>
  !! \brief Liberation du tableau de donnees des variables
  !!
  !! \param [in] cs              Pointer to \ref PDM_writer object
  !! \param [in] id_var          Identificateur de la variable
  !!
  !!

  subroutine PDM_writer_var_data_free (cs,     &
                                       id_var) &
  bind (c, name='PDM_writer_var_data_free')

    use iso_c_binding
    implicit none

    type(c_ptr),    value :: cs
    integer(c_int), value :: id_var

  end subroutine PDM_writer_var_data_free


  !>
  !! \brief Liberation d'une variable
  !!
  !! \param [in] cs              Pointer to \ref PDM_writer object
  !! \param [in] id_var          Identificateur de la variable
  !!
  !!

  subroutine PDM_writer_var_free (cs,     &
                                  id_var) &
  bind (c, name='PDM_writer_var_free')

    use iso_c_binding
    implicit none

    type(c_ptr),    value :: cs
    integer(c_int), value :: id_var

  end subroutine PDM_writer_var_free


  !>
  !! \brief Free format
  !!
  !!

  subroutine PDM_writer_fmt_free () &
  bind (c, name='PDM_writer_fmt_free')
  end subroutine PDM_writer_fmt_free



  !>
  !! \brief Réinitialisation des donnees decrivant le maillage courant
  !!
  !! \param [in] cs              Pointer to \ref PDM_writer object
  !! \param [in] id_geom         Identificateur de l'objet geometrique
  !!
  !!

  subroutine PDM_writer_geom_data_reset (cs,      &
                                         id_geom) &
  bind (c, name='PDM_writer_geom_data_reset')

    use iso_c_binding
    implicit none

    type(c_ptr),    value :: cs
    integer(c_int), value :: id_geom

  end subroutine PDM_writer_geom_data_reset


  end interface

  contains

  !>
  !!
  !! \brief Cree un objet CS (Cedre Sortie) et retoure un pointeur sur cet objet
  !!
  !! \param [out] cs                   Pointer to \ref PDM_writer object
  !! \param [in]  fmt                  Format de sortie
  !! \param [in]  fmt_fic              Binary or ASCII
  !! \param [in]  topologie            Indique le maillage est mobile ou non
  !! \param [in]  st_reprise           Complete les sorties des calculs precedents en reprise
  !! \param [in]  rep_sortie           Repertoire de sortie
  !! \param [in]  nom_sortie           Nom de la sortie
  !! \param [in]  pdm_mpi_com          Communicateur MSG
  !! \param [in]  acces                Type d'acces
  !! \param [in]  prop_noeuds_actifs   Proportion des noeuds actifs dans les acces au fichier
  !!                                     *  -1 : tous les processus actifs
  !!                                     *   1 : un processus par noeud
  !!                                     * 0 < val < 1 : un processus par noeud actif
  !! \param [in]  options              Options complementaires propres au format sous
  !!                                 la forme ("nom_1 = val_1 : ... : nom_n = val_n")
  !!
  !!

  subroutine PDM_writer_create (cs,                 &
                                fmt,                &
                                fmt_fic,            &
                                topologie,          &
                                st_reprise,         &
                                rep_sortie,         &
                                nom_sortie,         &
                                f_comm,             &
                                acces,              &
                                prop_noeuds_actifs, &
                                options)
    use iso_c_binding
    implicit none

    type(c_ptr)                  :: cs
    character(len = *)           :: fmt
    integer,          intent(in) :: fmt_fic
    integer,          intent(in) :: topologie
    integer,          intent(in) :: st_reprise
    character(len = *)           :: rep_sortie
    character(len = *)           :: nom_sortie
    integer,          intent(in) :: f_comm
    integer,          intent(in) :: acces
    double precision, intent(in) :: prop_noeuds_actifs
    character(len = *)           :: options

    integer(c_int)               :: c_fmt_fic
    integer(c_int)               :: c_topologie
    integer(c_int)               :: c_st_reprise
    integer(c_int)               :: c_comm
    integer(c_int)               :: c_acces
    real(c_double)               :: c_prop_noeuds_actifs

    interface
      function PDM_writer_create_c (fmt,                &
                                    fmt_fic,            &
                                    topologie,          &
                                    st_reprise,         &
                                    rep_sortie,         &
                                    nom_sortie,         &
                                    pdm_mpi_comm,       &
                                    acces,              &
                                    prop_noeuds_actifs, &
                                    options)            &
      result (cs)                                       &
      bind (c, name='PDM_writer_create')
        use iso_c_binding
        implicit none

        type(c_ptr)           :: cs
        character(c_char)     :: fmt
        integer(c_int), value :: fmt_fic
        integer(c_int), value :: topologie
        integer(c_int), value :: st_reprise
        character(c_char)     :: rep_sortie
        character(c_char)     :: nom_sortie
        integer(c_int), value :: pdm_mpi_comm
        integer(c_int), value :: acces
        real(c_double), value :: prop_noeuds_actifs
        character(c_char)     :: options

      end function PDM_writer_create_c
    end interface

    c_comm               = PDM_MPI_Comm_f2c(f_comm)
    c_fmt_fic            = fmt_fic
    c_topologie          = topologie
    c_st_reprise         = st_reprise
    c_acces              = acces
    c_prop_noeuds_actifs = prop_noeuds_actifs

    cs = PDM_writer_create_c (trim(fmt)//C_NULL_CHAR,        &
                              c_fmt_fic,               &
                              c_topologie,             &
                              c_st_reprise,            &
                              trim(rep_sortie)//C_NULL_CHAR, &
                              trim(nom_sortie)//C_NULL_CHAR, &
                              c_comm,                  &
                              c_acces,                 &
                              c_prop_noeuds_actifs,    &
                              trim(options)//C_NULL_CHAR)

  end subroutine PDM_writer_create



  !>
  !! \brief Debut d'increment
  !!
  !! \param [in] cs             Pointer to \ref PDM_writer object
  !! \param [in] physical_time  Temps
  !!

  subroutine PDM_writer_step_beg (cs,            &
                                  physical_time)
    use iso_c_binding
    implicit none

    type(c_ptr), value           :: cs
    double precision, intent(in) :: physical_time

    real(c_double)               :: c_physical_time

    interface
      subroutine PDM_writer_step_beg_c (cs,            &
                                        physical_time) &
      bind (c, name="PDM_writer_step_beg")

        use iso_c_binding
        implicit none

        type(c_ptr),    value :: cs
        real(c_double), value :: physical_time

      end subroutine PDM_writer_step_beg_c
    end interface

    c_physical_time = physical_time

    call PDM_writer_step_beg_c (cs,            &
                                c_physical_time)

  end subroutine PDM_writer_step_beg



  !>
  !! \brief Cree une nouvelle geometrie dans l'objet CS (Cedre Sortie)
  !!
  !! \param [in]  cs                Pointer to \ref PDM_writer object
  !! \param [out] id_geom           Identificateur de l'objet geom dans cs
  !! \param [in]  nom_geom          Nom de l'objet geometrique
  !!
  !!

  subroutine PDM_writer_geom_create (cs,               &
                                     id_geom,          &
                                     nom_geom,         &
                                     n_part)
    use iso_c_binding
    implicit none

    type(c_ptr), value           :: cs
    integer, intent(out)         :: id_geom
    character(len=*), intent(in) :: nom_geom
    integer, intent(in)          :: n_part

    integer(c_int)               :: c_id_geom
    integer(c_int)               :: c_n_part

    interface
      function PDM_writer_geom_create_c (cs,               &
                                         nom_geom,         &
                                         n_part)           &
      result (id_geom)                                     &
      bind (c, name='PDM_writer_geom_create')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: cs
        integer(c_int)        :: id_geom
        character(c_char)     :: nom_geom(*)
        integer(c_int), value :: n_part

      end function PDM_writer_geom_create_c
    end interface

    c_n_part           = n_part

    c_id_geom = PDM_writer_geom_create_c (cs,                    &
                                          trim(nom_geom)//C_NULL_CHAR, &
                                          c_n_part)

    id_geom = c_id_geom

  end subroutine PDM_writer_geom_create



  !>
  !! \brief Definition des coordonnees de la partition courante
  !!
  !! \param [in] cs        Pointer to \ref PDM_writer object
  !! \param [in] id_geom   Identificateur de l'objet geometrique
  !! \param [in] id_part   Indice de partition
  !! \param [in] n_som     Nombre de sommets de la partition
  !! \param [in] coords    Coordonnes des sommets
  !! \param [in] numabs    Numerotation absolue des sommets
  !!
  !!

  subroutine PDM_writer_geom_coord_set (cs,      &
                                        id_geom, &
                                        id_part, &
                                        n_som,   &
                                        coords,  &
                                        numabs,  &
                                        owner)
    use iso_c_binding
    implicit none

    type(c_ptr), value            :: cs
    integer, intent(in)           :: id_geom
    integer, intent(in)           :: id_part
    integer, intent(in)           :: n_som
    integer, intent(in)           :: owner
    double precision,     pointer :: coords(:,:)
    integer(pdm_g_num_s), pointer :: numabs(:)

    integer(c_int)                :: c_owner

    integer(c_int)                :: c_id_geom
    integer(c_int)                :: c_id_part
    integer(c_int)                :: c_n_som
    type(c_ptr)                   :: c_coords
    type(c_ptr)                   :: c_numabs

    interface
      subroutine PDM_writer_geom_coord_set_c (cs,      &
                                              id_geom, &
                                              id_part, &
                                              n_som,   &
                                              coords,  &
                                              numabs,  &
                                              owner)   &
      bind (c, name='PDM_writer_geom_coord_set')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: cs
        integer(c_int), value :: id_geom
        integer(c_int), value :: id_part
        integer(c_int), value :: n_som
        integer(c_int), value :: owner
        type(c_ptr),    value :: coords
        type(c_ptr),    value :: numabs

      end subroutine PDM_writer_geom_coord_set_c
    end interface

    c_id_geom = id_geom
    c_id_part = id_part
    c_n_som   = n_som
    c_owner   = owner

    c_coords = C_NULL_PTR
    if (associated(coords)) then
      c_coords = c_loc(coords)
    endif 
      
    c_numabs = C_NULL_PTR
    if (associated(numabs)) then
      c_numabs = c_loc(numabs)
    endif 
      

    call PDM_writer_geom_coord_set_c (cs,        &
                                      c_id_geom, &
                                      c_id_part, &
                                      c_n_som,   &
                                      c_coords,  &
                                      c_numabs,  &
                                      c_owner)

  end subroutine PDM_writer_geom_coord_set



  !>
  !! \brief Definition des coordonnees de la partition courante
  !! a partir d'un ensemble parent
  !!
  !! \param [in] cs               Pointer to \ref PDM_writer object
  !! \param [in] id_geom          Identificateur de l'objet geometrique
  !! \param [in] id_part          Indice de partition
  !! \param [in] n_som            Nombre de sommets de la partition
  !! \param [in] n_som_parent     Nombre de sommets parent
  !! \param [in] numabs           Numerotation absolue des sommets (size = n_som)
  !! \param [in] num_parent       Numerotation des sommets dans la numerotation parente (size = n_som)
  !! \param [in] coords_parent    Coordonnes des sommets parents (size = 3 * n_som_parent)
  !! \param [in] numabs_parent    Numerotation absolue des sommets parents (size = n_som_parent)
  !!
  !!

  subroutine PDM_writer_geom_coord_from_parent_set (cs,            &
                                                    id_geom,       &
                                                    id_part,       &
                                                    n_som,         &
                                                    n_som_parent,  &
                                                    numabs,        &
                                                    num_parent,    &
                                                    coords_parent, &
                                                    numabs_parent, &
                                                    owner)

    use iso_c_binding
    implicit none

    type(c_ptr), value            :: cs
    integer, intent(in)           :: id_geom
    integer, intent(in)           :: id_part
    integer, intent(in)           :: n_som
    integer, intent(in)           :: n_som_parent
    integer, intent(in)           :: owner
    integer(pdm_g_num_s), pointer :: numabs(:)
    integer(pdm_l_num_s), pointer :: num_parent(:)
    double precision,     pointer :: coords_parent(:,:)
    integer(pdm_g_num_s), pointer :: numabs_parent(:)

    integer(c_int)                :: c_id_geom
    integer(c_int)                :: c_id_part
    integer(c_int)                :: c_n_som
    integer(c_int)                :: c_n_som_parent
    integer(c_int)                :: c_owner
    type(c_ptr)                   :: c_numabs
    type(c_ptr)                   :: c_num_parent
    type(c_ptr)                   :: c_coords_parent
    type(c_ptr)                   :: c_numabs_parent

    interface
      subroutine PDM_writer_geom_coord_from_parent_set_c (cs,            &
                                                          id_geom,       &
                                                          id_part,       &
                                                          n_som,         &
                                                          n_som_parent,  &
                                                          numabs,        &
                                                          num_parent,    &
                                                          coords_parent, &
                                                          numabs_parent, &
                                                          owner) &
      bind (c, name='PDM_writer_geom_coord_from_parent_set')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: cs
        integer(c_int), value :: id_geom
        integer(c_int), value :: id_part
        integer(c_int), value :: n_som
        integer(c_int), value :: n_som_parent
        integer(c_int), value :: owner
        type(c_ptr),    value :: numabs
        type(c_ptr),    value :: num_parent
        type(c_ptr),    value :: coords_parent
        type(c_ptr),    value :: numabs_parent

      end subroutine PDM_writer_geom_coord_from_parent_set_c
    end interface

    c_id_geom      = id_geom
    c_id_part      = id_part
    c_n_som        = n_som
    c_n_som_parent = n_som_parent
    c_owner        = owner

    c_numabs = C_NULL_PTR
    if (associated(numabs)) then
      c_numabs        = c_loc(numabs)
    endif
      
    c_num_parent = C_NULL_PTR
    if (associated(num_parent)) then
      c_num_parent    = c_loc(num_parent)
    endif
      
    c_coords_parent = C_NULL_PTR
    if (associated(coords_parent)) then
      c_coords_parent = c_loc(coords_parent)
    endif
      
    c_numabs_parent = C_NULL_PTR
    if (associated(numabs_parent)) then
      c_numabs_parent = c_loc(numabs_parent)
    endif
      

    call PDM_writer_geom_coord_from_parent_set_c (cs,              &
                                                  c_id_geom,       &
                                                  c_id_part,       &
                                                  c_n_som,         &
                                                  c_n_som_parent,  &
                                                  c_numabs,        &
                                                  c_num_parent,    &
                                                  c_coords_parent, &
                                                  c_numabs_parent, &
                                                  c_owner)

  end subroutine PDM_writer_geom_coord_from_parent_set



  !>
  !! \brief Ajout d'un bloc d'elements d'un type donne
  !!
  !! \param [in] cs             Pointer to \ref PDM_writer object
  !! \param [in] id_geom        Identificateur de l'objet geometrique
  !! \param [in] t_elt          Type d'element
  !!
  !! \return   Identificateur du bloc
  !!
  !!

  subroutine PDM_writer_geom_bloc_add (cs,           &
                                       id_geom,      &
                                       t_elt,        &
                                       owner,        &
                                       id_bloc)
    use iso_c_binding
    implicit none

    type(c_ptr), value            :: cs
    integer, intent(in)           :: id_geom
    integer, intent(in)           :: t_elt
    integer, intent(in)           :: owner
    integer, intent(out)          :: id_bloc

    integer(c_int)                :: c_id_geom
    integer(c_int)                :: c_t_elt
    integer(c_int)                :: c_owner

    interface
      function PDM_writer_geom_bloc_add_c (cs,             &
                                             id_geom,      &
                                             t_elt,        & 
                                             owner)        &
      result (id_bloc)                                     &
      bind (c, name='PDM_writer_geom_bloc_add')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: cs
        integer(c_int), value :: id_geom
        integer(c_int), value :: t_elt
        integer(c_int), value :: owner
        integer(c_int)        :: id_bloc

      end function PDM_writer_geom_bloc_add_c
    end interface

    c_id_geom      = id_geom
    c_t_elt        = t_elt
    c_owner        = owner

    id_bloc = PDM_writer_geom_bloc_add_c (cs,        &
                                     c_id_geom,      &
                                     c_t_elt,        &
                                     c_owner)

  end subroutine PDM_writer_geom_bloc_add



  !>
  !! \brief Ajout d'un bloc d'elements d'un type donne dans la partition courante
  !!
  !!  - PDM_writer_POINT :
  !!
  !!   1 x
  !!
  !!  - PDM_writer_BAR2 :
  !!
  !!   1 x-------x 2
  !!
  !!  - PDM_writer_TRIA3 :
  !!
  !!   1 x-------x 3
  !!      \     /
  !!       \   /
  !!        \ /
  !!         x 2
  !!
  !!  - PDM_writer_QUAD4 :
  !!
  !!      4 x-------x 3
  !!       /       /
  !!      /       /
  !!   1 x-------x2
  !!
  !!   - PDM_writer_TETRA4 :
  !!
  !!         x 4
  !!        /|\
  !!       / | \
  !!      /  |  \
  !!   1 x- -|- -x 3
  !!      \  |  /
  !!       \ | /
  !!        \|/
  !!         x 2
  !!
  !!   - PDM_writer_PYRAMID5 :
  !!
  !!          5 x
  !!           /|\
  !!          //| \
  !!         // |  \
  !!      4 x/--|---x 3
  !!       //   |  /
  !!      //    | /
  !!   1 x-------x 2
  !!
  !!  - PDM_writer_PRSIM6 :
  !!
  !!   4 x-------x 6
  !!     |\     /|
  !!     | \   / |
  !!   1 x- \-/ -x 3
  !!      \ 5x  /
  !!       \ | /
  !!        \|/
  !!         x 2
  !!
  !!  - PDM_writer_HEXA8 :
  !!
  !!      8 x-------x 7
  !!       /|      /|
  !!      / |     / |
  !!   5 x-------x6 |
  !!     | 4x----|--x 3
  !!     | /     | /
  !!     |/      |/
  !!   1 x-------x 2
  !!
  !! \param [in] cs                  Pointer to \ref PDM_writer object
  !! \param [in] id_geom             Identificateur de l'objet geometrique
  !! \param [in] id_bloc             Identificateur du bloc
  !! \param [in] id_part             Indice de partition
  !! \param [in] t_elt               Type d'element
  !! \param [in] n_elt               Nombre d'elements dans le bloc
  !! \param [in] connec              Table de connectivite des elements
  !! \param [in] num_part            Numerotation dans la partition
  !!
  !!

  subroutine PDM_writer_geom_bloc_std_set (cs,      &
                                           id_geom, &
                                           id_bloc, &
                                           id_part, &
                                           n_elt,   &
                                           connec,  &
                                           numabs)
    use iso_c_binding
    implicit none

    type(c_ptr), value            :: cs
    integer, intent(in)           :: id_geom
    integer, intent(in)           :: id_bloc
    integer, intent(in)           :: id_part
    integer, intent(in)           :: n_elt
    integer(pdm_l_num_s), pointer :: connec(:)
    integer(pdm_g_num_s), pointer :: numabs(:)

    integer(c_int)                :: c_id_geom
    integer(c_int)                :: c_id_bloc
    integer(c_int)                :: c_id_part
    integer(c_int)                :: c_n_elt
    type(c_ptr)                   :: c_connec
    type(c_ptr)                   :: c_numabs

    interface
      subroutine PDM_writer_geom_bloc_std_set_c (cs,      &
                                                 id_geom, &
                                                 id_bloc, &
                                                 id_part, &
                                                 n_elt,   &
                                                 connec,  &
                                                 numabs)  &
      bind (c, name='PDM_writer_geom_bloc_std_set')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: cs
        integer(c_int), value :: id_geom
        integer(c_int), value :: id_bloc
        integer(c_int), value :: id_part
        integer(c_int), value :: n_elt
        type(c_ptr),    value :: connec
        type(c_ptr),    value :: numabs

      end subroutine PDM_writer_geom_bloc_std_set_c
    end interface

    c_id_geom = id_geom
    c_id_bloc = id_bloc
    c_id_part = id_part
    c_n_elt   = n_elt

    c_connec = C_NULL_PTR
    if (associated(connec)) then
      c_connec = c_loc(connec)
    endif
      
    c_numabs = C_NULL_PTR
    if (associated(numabs)) then
      c_numabs = c_loc(numabs)
    endif
      

    call PDM_writer_geom_bloc_std_set_c (cs,        &
                                         c_id_geom, &
                                         c_id_bloc, &
                                         c_id_part, &
                                         c_n_elt,   &
                                         c_connec,  &
                                         c_numabs)

  end subroutine PDM_writer_geom_bloc_std_set



  !>
  !! \brief Ajout d'un bloc de polygones dans la partition courante
  !!
  !! \param [in] cs              Pointer to \ref PDM_writer object
  !! \param [in] id_geom         Identificateur de l'objet geometrique
  !! \param [in] id_part         Indice de partition
  !! \param [in] n_elt           Nombre d'elements dans le bloc
  !! \param [in] connec_idx      Index dans la table de connectivite (dim = n_elt+1)
  !! \param [in] connec          Table de connectivite des elements (dim = connec_idx[n_elt])
  !! \param [in] numabs          Numerotation absolue des elements
  !!
  !!

  subroutine PDM_writer_geom_bloc_poly2d_set (cs,         &
                                              id_geom,    &
                                              id_bloc,    &
                                              id_part,    &
                                              n_elt,      &
                                              connec_idx, &
                                              connec,     &
                                              numabs)
    use iso_c_binding
    implicit none

    type(c_ptr), value            :: cs
    integer, intent(in)           :: id_geom
    integer, intent(in)           :: id_bloc
    integer, intent(in)           :: id_part
    integer, intent(in)           :: n_elt
    integer(pdm_l_num_s), pointer :: connec_idx(:)
    integer(pdm_l_num_s), pointer :: connec(:)
    integer(pdm_g_num_s), pointer :: numabs(:)

    integer(c_int)                :: c_id_geom
    integer(c_int)                :: c_id_bloc
    integer(c_int)                :: c_id_part
    integer(c_int)                :: c_n_elt
    type(c_ptr)                   :: c_connec_idx
    type(c_ptr)                   :: c_connec
    type(c_ptr)                   :: c_numabs

    interface
      subroutine PDM_writer_geom_bloc_poly2d_set_c (cs,         &
                                                    id_geom,    &
                                                    id_bloc,    &
                                                    id_part,    &
                                                    n_elt,      &
                                                    connec_idx, &
                                                    connec,     &
                                                    numabs)     &
      bind (c, name='PDM_writer_geom_bloc_poly2d_set')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: cs
        integer(c_int), value :: id_geom
        integer(c_int), value :: id_bloc
        integer(c_int), value :: id_part
        integer(c_int), value :: n_elt
        type(c_ptr),    value :: connec_idx
        type(c_ptr),    value :: connec
        type(c_ptr),    value :: numabs

      end subroutine PDM_writer_geom_bloc_poly2d_set_c
    end interface

    c_id_geom = id_geom
    c_id_bloc = id_bloc
    c_id_part = id_part
    c_n_elt   = n_elt

    c_connec_idx = C_NULL_PTR
    if (associated(connec_idx)) then
      c_connec_idx = c_loc(connec_idx)
    endif 
      
    c_connec = C_NULL_PTR
    if (associated(connec)) then
      c_connec     = c_loc(connec)
    endif 
      
    c_numabs = C_NULL_PTR
    if (associated(numabs)) then
      c_numabs     = c_loc(numabs)
    endif 
      

    call PDM_writer_geom_bloc_poly2d_set_c (cs,           &
                                            c_id_geom,    &
                                            c_id_bloc,    &
                                            c_id_part,    &
                                            c_n_elt,      &
                                            c_connec_idx, &
                                            c_connec,     &
                                            c_numabs)

  end subroutine PDM_writer_geom_bloc_poly2d_set



  !>
  !! \brief Ajout d'un bloc de polyedres dans la partition courante
  !!
  !! \param [in] cs              Pointer to \ref PDM_writer object
  !! \param [in] id_geom         Identificateur de l'objet geometrique
  !! \param [in] id_part         Indice de partition
  !! \param [in] n_elt           Nombre d'elements dans le bloc
  !! \param [in] n_face          Nombre de faces de chaque element (dim = n_elt)
  !! \param [in] facsom_idx      Index dans la table de connectivite des faces (dim = n_face_total+1)
  !! \param [in] facsom          Table de connectivite des faces (dim = facsom_idx[n_face_total}
  !! \param [in] cellfac_idx     Index dans la table de connectivite des cellules (dim = n_elt+1)
  !! \param [in] cellfac         Table de connectivite des elements (dim = cellfac_idx[n_elt])
  !! \param [in] numabs          Numerotation absolue des elements
  !!
  !!

  subroutine PDM_writer_geom_bloc_poly3d_set (cs,          &
                                              id_geom,     &
                                              id_bloc,     &
                                              id_part,     &
                                              n_elt,       &
                                              n_face,      &
                                              facsom_idx,  &
                                              facsom,      &
                                              cellfac_idx, &
                                              cellfac,     &
                                              numabs)
    use iso_c_binding
    implicit none

    type(c_ptr), value            :: cs
    integer, intent(in)           :: id_geom
    integer, intent(in)           :: id_bloc
    integer, intent(in)           :: id_part
    integer, intent(in)           :: n_elt
    integer, intent(in)           :: n_face
    integer(pdm_l_num_s), pointer :: facsom_idx(:)
    integer(pdm_l_num_s), pointer :: facsom(:)
    integer(pdm_l_num_s), pointer :: cellfac_idx(:)
    integer(pdm_l_num_s), pointer :: cellfac(:)
    integer(pdm_g_num_s), pointer :: numabs(:)

    integer(c_int)                :: c_id_geom
    integer(c_int)                :: c_id_bloc
    integer(c_int)                :: c_id_part
    integer(c_int)                :: c_n_elt
    integer(c_int)                :: c_n_face
    type(c_ptr)                   :: c_facsom_idx
    type(c_ptr)                   :: c_facsom
    type(c_ptr)                   :: c_cellfac_idx
    type(c_ptr)                   :: c_cellfac
    type(c_ptr)                   :: c_numabs

    interface
      subroutine PDM_writer_geom_bloc_poly3d_set_c (cs,          &
                                                    id_geom,     &
                                                    id_bloc,     &
                                                    id_part,     &
                                                    n_elt,       &
                                                    n_face,      &
                                                    facsom_idx,  &
                                                    facsom,      &
                                                    cellfac_idx, &
                                                    cellfac,     &
                                                    numabs)      &
      bind (c, name='PDM_writer_geom_bloc_poly3d_set')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: cs
        integer(c_int), value :: id_geom
        integer(c_int), value :: id_bloc
        integer(c_int), value :: id_part
        integer(c_int), value :: n_elt
        integer(c_int), value :: n_face
        type(c_ptr),    value :: facsom_idx
        type(c_ptr),    value :: facsom
        type(c_ptr),    value :: cellfac_idx
        type(c_ptr),    value :: cellfac
        type(c_ptr),    value :: numabs

      end subroutine PDM_writer_geom_bloc_poly3d_set_c
    end interface

    c_id_geom = id_geom
    c_id_bloc = id_bloc
    c_id_part = id_part
    c_n_elt   = n_elt
    c_n_face  = n_face

    c_facsom_idx = C_NULL_PTR
    if (associated(facsom_idx)) then 
      c_facsom_idx  = c_loc(facsom_idx)
    endif 
      
    c_facsom = C_NULL_PTR
    if (associated(facsom)) then 
      c_facsom      = c_loc(facsom)
    endif 
      
    c_cellfac_idx = C_NULL_PTR
    if (associated(cellfac_idx)) then 
      c_cellfac_idx = c_loc(cellfac_idx)
    endif 
      
    c_cellfac = C_NULL_PTR
    if (associated(cellfac)) then 
      c_cellfac     = c_loc(cellfac)
    endif 
      
    c_numabs = C_NULL_PTR
    if (associated(numabs)) then 
      c_numabs      = c_loc(numabs)
    endif 
      

    call PDM_writer_geom_bloc_poly3d_set_c (cs,            &
                                            c_id_geom,     &
                                            c_id_bloc,     &
                                            c_id_part,     &
                                            c_n_elt,       &
                                            c_n_face,      &
                                            c_facsom_idx,  &
                                            c_facsom,      &
                                            c_cellfac_idx, &
                                            c_cellfac,     &
                                            c_numabs)

  end subroutine PDM_writer_geom_bloc_poly3d_set




  !>
  !!
  !! \brief Ajout de cellules 3D decrites en fonctions des faces.
  !!
  !! Cette fonction détermine les types des éléments et crée des blocs regrouppant les éléments
  !! de même type. Elle retourne l'indirection vers le nouvel ordre de rangement
  !! des cellules.
  !!
  !! \param [in]  cs              Pointer to \ref PDM_writer object
  !! \param [in]  id_geom         Identificateur de l'objet geometrique
  !! \param [in]  id_part         Identificateur de partition
  !! \param [in]  n_cell          Nombre de cellules 3D ajoutées
  !! \param [in]  n_face          Nombre de faces décrites
  !! \param [in]  face_som_idx    Index de connectivite faces -> sommets
  !! \param [in]  face_som        Connectivite faces -> sommets
  !! \param [in]  cell_face_idx   Index de connectivite cellules -> faces
  !! \param [in]  cell_face       Connectivite cellules -> faces
  !! \param [in]  numabs          Numerotation absolue des cellules
  !!
  !!

  subroutine PDM_writer_geom_cell3d_cellface_add (cs,            &
                                                  id_geom,       &
                                                  id_part,       &
                                                  n_cell,        &
                                                  n_face,        &
                                                  face_som_idx,  &
                                                  face_som_nb,   &
                                                  face_som,      &
                                                  cell_face_idx, &
                                                  cell_face_nb,  &
                                                  cell_face,     &
                                                  numabs)
    use iso_c_binding
    implicit none

    type(c_ptr), value            :: cs
    integer, intent(in)           :: id_geom
    integer, intent(in)           :: id_part
    integer, intent(in)           :: n_cell
    integer, intent(in)           :: n_face
    integer(pdm_l_num_s), pointer :: face_som_idx(:)
    integer(pdm_l_num_s), pointer :: face_som_nb(:)
    integer(pdm_l_num_s), pointer :: face_som(:)
    integer(pdm_l_num_s), pointer :: cell_face_idx(:)
    integer(pdm_l_num_s), pointer :: cell_face_nb(:)
    integer(pdm_l_num_s), pointer :: cell_face(:)
    integer(pdm_g_num_s), pointer :: numabs(:)

    integer(c_int)                :: c_id_geom
    integer(c_int)                :: c_id_part
    integer(c_int)                :: c_n_cell
    integer(c_int)                :: c_n_face
    type(c_ptr)                   :: c_face_som_idx
    type(c_ptr)                   :: c_face_som_nb
    type(c_ptr)                   :: c_face_som
    type(c_ptr)                   :: c_cell_face_idx
    type(c_ptr)                   :: c_cell_face_nb
    type(c_ptr)                   :: c_cell_face
    type(c_ptr)                   :: c_numabs

    interface
      subroutine PDM_writer_geom_cell3d_cellface_add_c (cs,            &
                                                        id_geom,       &
                                                        id_part,       &
                                                        n_cell,        &
                                                        n_face,        &
                                                        face_som_idx,  &
                                                        face_som_nb,   &
                                                        face_som,      &
                                                        cell_face_idx, &
                                                        cell_face_nb,  &
                                                        cell_face,     &
                                                        numabs)        &
      bind (c, name='PDM_writer_geom_cell3d_cellface_add')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: cs
        integer(c_int), value :: id_geom
        integer(c_int), value :: id_part
        integer(c_int), value :: n_cell
        integer(c_int), value :: n_face
        type(c_ptr),    value :: face_som_idx
        type(c_ptr),    value :: face_som_nb
        type(c_ptr),    value :: face_som
        type(c_ptr),    value :: cell_face_idx
        type(c_ptr),    value :: cell_face_nb
        type(c_ptr),    value :: cell_face
        type(c_ptr),    value :: numabs

      end subroutine PDM_writer_geom_cell3d_cellface_add_c
    end interface

    c_id_geom = id_geom
    c_id_part = id_part
    c_n_cell  = n_cell
    c_n_face  = n_face

    c_face_som_idx = C_NULL_PTR
    if (associated(face_som_idx)) then
        c_face_som_idx  = c_loc(face_som_idx)
    endif
        
    c_face_som = C_NULL_PTR
    if (associated(face_som)) then
        c_face_som      = c_loc(face_som)
    endif
        
    c_cell_face_idx = C_NULL_PTR
    if (associated(cell_face_idx)) then
        c_cell_face_idx = c_loc(cell_face_idx)
    endif
        
    c_numabs = C_NULL_PTR
    if (associated(numabs)) then
        c_numabs        = c_loc(numabs)
    endif
        
    c_cell_face = C_NULL_PTR
    if (associated(cell_face)) then
        c_cell_face     = c_loc(cell_face)
    endif      

    c_face_som_nb   = C_NULL_PTR
    if (associated(face_som_nb)) then
      c_face_som_nb   = c_loc(face_som_nb)
    endif     

    c_cell_face_nb  = C_NULL_PTR
    if (associated(cell_face_nb)) then
      c_cell_face_nb  = c_loc(cell_face_nb)
    endif

    call PDM_writer_geom_cell3d_cellface_add_c (cs,              &
                                                c_id_geom,       &
                                                c_id_part,       &
                                                c_n_cell,        &
                                                c_n_face,        &
                                                c_face_som_idx,  &
                                                c_face_som_nb,   &
                                                c_face_som,      &
                                                c_cell_face_idx, &
                                                c_cell_face_nb,  &
                                                c_cell_face,     &
                                                c_numabs)

  end subroutine PDM_writer_geom_cell3d_cellface_add



  !>
  !!
  !! \brief Ajout de cellules 2D decrites en fonctions des faces.
  !!
  !! Cette fonction détermine les types des éléments et crée des blocs regrouppant les éléments
  !! de même type. Elle retourne l'indirection vers le nouvel ordre de rangement
  !! des cellules.
  !!
  !! \param [in] cs              Pointer to \ref PDM_writer object
  !! \param [in] id_geom         Identificateur de l'objet geometrique
  !! \param [in] n_cell          Nombre de cellules 3D ajoutées
  !! \param [in] n_face          Nombre de faces décrites
  !! \param [in] face_som_idx    Index de connectivite faces -> sommets
  !! \param [in] face_som        Connectivite faces -> sommets
  !! \param [in] cell_face_idx   Index de connectivite cellules -> faces
  !! \param [in] cell_face       Connectivite cellules -> faces
  !! \param [in] numabs          Numerotation absolue des cellules
  !!
  !!

  subroutine PDM_writer_geom_cell2d_cellface_add (cs,            &
                                                  id_geom,       &
                                                  id_part,       &
                                                  n_cell,        &
                                                  n_face,        &
                                                  face_som_idx,  &
                                                  face_som_nb,   &
                                                  face_som,      &
                                                  cell_face_idx, &
                                                  cell_face_nb,  &
                                                  cell_face,     &
                                                  numabs)
    use iso_c_binding
    implicit none

    type(c_ptr), value            :: cs
    integer, intent(in)           :: id_geom
    integer, intent(in)           :: id_part
    integer, intent(in)           :: n_cell
    integer, intent(in)           :: n_face
    integer(pdm_l_num_s), pointer :: face_som_idx(:)
    integer(pdm_l_num_s), pointer :: face_som_nb(:)
    integer(pdm_l_num_s), pointer :: face_som(:)
    integer(pdm_l_num_s), pointer :: cell_face_idx(:)
    integer(pdm_l_num_s), pointer :: cell_face_nb(:)
    integer(pdm_l_num_s), pointer :: cell_face(:)
    integer(pdm_g_num_s), pointer :: numabs(:)

    integer(c_int)                :: c_id_geom
    integer(c_int)                :: c_id_part
    integer(c_int)                :: c_n_cell
    integer(c_int)                :: c_n_face
    type(c_ptr)                   :: c_face_som_idx
    type(c_ptr)                   :: c_face_som_nb
    type(c_ptr)                   :: c_face_som
    type(c_ptr)                   :: c_cell_face_idx
    type(c_ptr)                   :: c_cell_face_nb
    type(c_ptr)                   :: c_cell_face
    type(c_ptr)                   :: c_numabs

    interface
      subroutine PDM_writer_geom_cell2d_cellface_add_c (cs,            &
                                                        id_geom,       &
                                                        id_part,       &
                                                        n_cell,        &
                                                        n_face,        &
                                                        face_som_idx,  &
                                                        face_som_nb,   &
                                                        face_som,      &
                                                        cell_face_idx, &
                                                        cell_face_nb,  &
                                                        cell_face,     &
                                                        numabs)        &
      bind (c, name='PDM_writer_geom_cell2d_cellface_add')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: cs
        integer(c_int), value :: id_geom
        integer(c_int), value :: id_part
        integer(c_int), value :: n_cell
        integer(c_int), value :: n_face
        type(c_ptr),    value :: face_som_idx
        type(c_ptr),    value :: face_som_nb
        type(c_ptr),    value :: face_som
        type(c_ptr),    value :: cell_face_idx
        type(c_ptr),    value :: cell_face_nb
        type(c_ptr),    value :: cell_face
        type(c_ptr),    value :: numabs

      end subroutine PDM_writer_geom_cell2d_cellface_add_c
    end interface

    c_id_geom = id_geom
    c_id_part = id_part
    c_n_cell  = n_cell
    c_n_face  = n_face

    c_face_som_idx = C_NULL_PTR
    if (associated(face_som_idx)) then
      c_face_som_idx  = c_loc(face_som_idx)
    endif
      
    c_face_som = C_NULL_PTR
    if (associated(face_som)) then
      c_face_som      = c_loc(face_som)
    endif
      
    c_cell_face_idx = C_NULL_PTR
    if (associated(cell_face_idx)) then
      c_cell_face_idx = c_loc(cell_face_idx)
    endif
      
    c_face_som_nb = C_NULL_PTR
    if (associated(face_som_nb)) then
      c_face_som_nb   = C_NULL_PTR
    endif
      
    c_cell_face = C_NULL_PTR
    if (associated(cell_face)) then
      c_cell_face     = c_loc(cell_face)
    endif
      
    c_numabs = C_NULL_PTR
    if (associated(numabs)) then
      c_numabs        = c_loc(numabs)
    endif

    c_cell_face_nb  = C_NULL_PTR
    if (associated(cell_face_nb)) then
      c_cell_face_nb  = c_loc(cell_face_nb)
    endif

    call PDM_writer_geom_cell2d_cellface_add_c (cs,              &
                                                c_id_geom,       &
                                                c_id_part,       &
                                                c_n_cell,        &
                                                c_n_face,        &
                                                c_face_som_idx,  &
                                                c_face_som_nb,   &
                                                c_face_som,      &
                                                c_cell_face_idx, &
                                                c_cell_face_nb,  &
                                                c_cell_face,     &
                                                c_numabs)

  end subroutine PDM_writer_geom_cell2d_cellface_add



  !>
  !!
  !! \brief Ajout de faces decrites en fonctions des sommets.
  !!
  !! Cette fonction détermine les types des éléments et crée des blocs regrouppant les éléments
  !! de même type. Elle retourne l'indirection vers le nouvel ordre de rangement
  !! des cellules.
  !!
  !! \param [in] cs              Pointer to \ref PDM_writer object
  !! \param [in] id_geom         Identificateur de l'objet geometrique
  !! \param [in] n_elt           Nombre de cellules 3D ajoutées
  !! \param [in] n_face          Nombre de faces décrites
  !! \param [in] face_som_idx    Index de connectivite faces -> sommets
  !! \param [in] face_som        Connectivite faces -> sommets
  !! \param [in] numabs          Numerotation absolue des faces
  !!
  !!

  subroutine PDM_writer_geom_faces_facesom_add (cs,           &
                                                id_geom,      &
                                                id_part,      &
                                                n_face,       &
                                                face_som_idx, &
                                                face_som_nb,  &
                                                face_som,     &
                                                numabs)
    use iso_c_binding
    implicit none

    type(c_ptr), value            :: cs
    integer, intent(in)           :: id_geom
    integer, intent(in)           :: id_part
    integer, intent(in)           :: n_face
    integer(pdm_l_num_s), pointer :: face_som_idx(:)
    integer(pdm_l_num_s), pointer :: face_som_nb(:)
    integer(pdm_l_num_s), pointer :: face_som(:)
    integer(pdm_g_num_s), pointer :: numabs(:)

    integer(c_int)                :: c_id_geom
    integer(c_int)                :: c_id_part
    integer(c_int)                :: c_n_face
    type(c_ptr)                   :: c_face_som_idx
    type(c_ptr)                   :: c_face_som_nb
    type(c_ptr)                   :: c_face_som
    type(c_ptr)                   :: c_numabs

    interface
      subroutine PDM_writer_geom_faces_facesom_add_c (cs,            &
                                                      id_geom,       &
                                                      id_part,       &
                                                      n_face,        &
                                                      face_som_idx,  &
                                                      face_som_nb,   &
                                                      face_som,      &
                                                      numabs)        &
      bind (c, name='PDM_writer_geom_faces_facesom_add')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: cs
        integer(c_int), value :: id_geom
        integer(c_int), value :: id_part
        integer(c_int), value :: n_face
        type(c_ptr),    value :: face_som_idx
        type(c_ptr),    value :: face_som_nb
        type(c_ptr),    value :: face_som
        type(c_ptr),    value :: numabs

      end subroutine PDM_writer_geom_faces_facesom_add_c
    end interface

    c_id_geom = id_geom
    c_id_part = id_part
    c_n_face  = n_face

    c_face_som_idx = C_NULL_PTR
    if (associated(face_som_idx)) then
      c_face_som_idx  = c_loc(face_som_idx)
    endif 
      
    c_face_som_nb = C_NULL_PTR
    if (associated(face_som_nb)) then
      c_face_som_nb   = c_loc(face_som_nb)
    endif 
      
    c_face_som = C_NULL_PTR
    if (associated(face_som)) then
      c_face_som      = c_loc(face_som)
    endif 
      
    c_numabs = C_NULL_PTR
    if (associated(numabs)) then
      c_numabs        = c_loc(numabs)
    endif 
      

    call PDM_writer_geom_faces_facesom_add_c (cs,              &
                                              c_id_geom,       &
                                              c_id_part,       &
                                              c_n_face,        &
                                              c_face_som_idx,  &
                                              c_face_som_nb,   &
                                              c_face_som,      &
                                              c_numabs)

  end subroutine PDM_writer_geom_faces_facesom_add



  !>
  !! \brief Creation d'une variable
  !!
  !! \param [in]  cs              Pointer to \ref PDM_writer object
  !! \param [out] id_var          Identificateur de l'objet variable
  !! \param [in]  st_dep_temps    Indique si la variable est dependante du temps
  !! \param [in]  id_geom         Identificateur de l'objet geometrique
  !! \param [in]  dim             Dimension de la variable
  !! \param [in]  loc             Localisation de la variable
  !! \param [in]  nom_var         Nom de la variable
  !!
  !!

  subroutine PDM_writer_var_create (cs,         &
                                    id_var,     &
                                    st_dep_tps, &
                                    dim,        &
                                    loc,        &
                                    nom_var)
    use iso_c_binding
    implicit none

    type(c_ptr), value   :: cs
    integer, intent(out) :: id_var
    integer, intent(in)  :: st_dep_tps
    integer, intent(in)  :: dim
    integer, intent(in)  :: loc
    character (len=*)    :: nom_var

    integer(c_int)       :: c_st_dep_tps
    integer(c_int)       :: c_dim
    integer(c_int)       :: c_dof_loc

    interface
      function PDM_writer_var_create_c (cs,         &
                                        st_dep_tps, &
                                        dim,        &
                                        loc,        &
                                        nom_var)    &
      result (id_var)                               &
      bind (c, name='PDM_writer_var_create')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: cs
        integer(c_int)        :: id_var
        integer(c_int), value :: st_dep_tps
        integer(c_int), value :: dim
        integer(c_int), value :: loc
        character(c_char)     :: nom_var(*)

      end function PDM_writer_var_create_c
    end interface

    c_st_dep_tps = st_dep_tps
    c_dim        = dim
    c_dof_loc    = loc

    id_var = PDM_writer_var_create_c (cs,                   &
                                      c_st_dep_tps,         &
                                      c_dim,                &
                                      c_dof_loc,            &
                                      trim(nom_var)//C_NULL_CHAR)

  end subroutine PDM_writer_var_create


  !>
  !! \brief Creation d'une variable globale constante
  !!
  !! \param [in]  cs              Pointer to \ref PDM_writer object
  !! \param [out] id_var          Identificateur de l'objet variable
  !! \param [in]  nom_var         Nom de la variable
  !! \param [in]  val_var         Valeur
  !!
  !!

  subroutine PDM_writer_cst_global_var_create (cs,      &
                                               id_var,  &
                                               nom_var, &
                                               val_var)
    use iso_c_binding
    implicit none

    type(c_ptr), value   :: cs
    integer, intent(out) :: id_var
    character (len=*)    :: nom_var
    double precision     :: val_var

    real(c_double)       :: c_var_val

    interface
      function PDM_writer_cst_global_var_create_c (cs,         &
                                                   nom_var,    &
                                                   var_val)    &
      result (id_var)                                          &
      bind (c, name='PDM_writer_cst_global_var_create')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: cs
        integer(c_int)        :: id_var
        character(c_char)     :: nom_var(*)
        real(c_double)        :: var_val

      end function PDM_writer_cst_global_var_create_c
    end interface

    c_var_val    = val_var

    id_var = PDM_writer_cst_global_var_create_c (cs,                   &
                                                 nom_var//C_NULL_CHAR, &
                                                 c_var_val)

  end subroutine PDM_writer_cst_global_var_create


  !>
  !! \brief Creation d'une variable globale constante
  !!
  !! \param [in]  cs              Pointer to \ref PDM_writer object
  !! \param [out] id_var          Identificateur de l'objet variable
  !! \param [in]  nom_var         Nom de la variable
  !! \param [in]  val_var         Valeur
  !!
  !!

  subroutine PDM_writer_cst_global_var_set (cs,      &
                                               id_var,  &
                                               val_var)
    use iso_c_binding
    implicit none

    type(c_ptr), value          :: cs
    integer                     :: id_var
    double precision            :: val_var

    real(c_double)       :: c_var_val
    integer(c_int)       :: c_id_var

    interface
      subroutine PDM_writer_cst_global_var_set_c (cs,         &
                                                   id_var,    &
                                                   var_val)    &
      bind (c, name='PDM_writer_cst_global_var_set')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: cs
        integer(c_int)        :: id_var
        real(c_double)        :: var_val

      end subroutine PDM_writer_cst_global_var_set_c
    end interface

    c_var_val    = val_var
    c_id_var     = id_var

    call PDM_writer_cst_global_var_set_c (cs,                   &
                                          c_id_var,             &
                                          c_var_val)

  end subroutine PDM_writer_cst_global_var_set


  !>
  !! \brief Mapping des noms de variable
  !!
  !! \param [in] cs              Pointer to \ref PDM_writer object
  !! \param [in] public_name     Nom Public de la variable
  !! \param [in] pivate_name     Nom privé de la variable
  !!
  !!

  subroutine PDM_writer_name_map_add (cs,           &
                                      public_name,  &
                                      private_name)
    use iso_c_binding
    implicit none

    type(c_ptr), value :: cs
    character (len=*)  :: public_name
    character (len=*)  :: private_name

    interface
      subroutine PDM_writer_name_map_add_c (cs,           &
                                            public_name,  &
                                            private_name) &
      bind (c, name='PDM_writer_name_map_add')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: cs
        character(c_char)     :: public_name(*)
        character(c_char)     :: private_name(*)

      end subroutine PDM_writer_name_map_add_c
    end interface

    call PDM_writer_name_map_add_c (cs,                        &
                                    trim(public_name)//C_NULL_CHAR,  &
                                    trim(private_name)//C_NULL_CHAR)

  end subroutine PDM_writer_name_map_add



  !>
  !! \brief Mise a jour des valeurs de la variable.
  !!
  !! Attention, les valeurs définies aux elements doivent être définies suivant l'ordre de définition des blocs !
  !!
  !! \param [in] cs              Pointer to \ref PDM_writer object
  !! \param [in] id_geom         Identificateur de l'objet geometrique
  !! \param [in] id_part         Identificateur de la partition dans l'objet geometrique
  !! \param [in] id_var          Identificateur de la variable mise à jour
  !! \param [in] val             Valeurs
  !!
  !!

  subroutine PDM_writer_var_set (cs,      &
                                 id_var,  &
                                 id_geom, &
                                 id_part, &
                                 val)
    use iso_c_binding
    implicit none

    type(c_ptr), value        :: cs
    integer, intent(in)       :: id_var
    integer, intent(in)       :: id_geom
    integer, intent(in)       :: id_part
    double precision, pointer :: val(:)

    integer(c_int)            :: c_id_var
    integer(c_int)            :: c_id_geom
    integer(c_int)            :: c_id_part
    type(c_ptr)               :: c_val

    interface
      subroutine PDM_writer_var_set_c (cs,      &
                                       id_var,  &
                                       id_geom, &
                                       id_part, &
                                       val)     &
      bind (c, name='PDM_writer_var_set')
        use iso_c_binding
        implicit none

        type(c_ptr),    value :: cs
        integer(c_int), value :: id_var
        integer(c_int), value :: id_geom
        integer(c_int), value :: id_part
        type(c_ptr),    value :: val

      end subroutine PDM_writer_var_set_c
    end interface

    c_id_var  = id_var
    c_id_geom = id_geom
    c_id_part = id_part

    c_val = C_NULL_PTR
    if (associated(val)) then
      c_val = c_loc(val)
    endif  

    call PDM_writer_var_set_c (cs,        &
                               c_id_var,  &
                               c_id_geom, &
                               c_id_part, &
                               c_val)

  end subroutine PDM_writer_var_set



  !>
  !! \brief Add a writer format
  !!
  !! Define a new format writer
  !! WARNING: has not been tested, not sure about procedure pointer interoperability
  !!
  !! \param [in] name            Name
  !! \param [in] create_fct      Customize \ref PDM_writer_create function for the new format  (or NULL)
  !! \param [in] free_fct        Customize \ref PDM_writer_free function for the new format (or NULL)
  !! \param [in] beg_step_fct    Customize \ref PDM_writer_step_beg function for the new format (or NULL)
  !! \param [in] end_step_fct    Customize \ref PDM_writer_step_end function for the new format (or NULL)
  !! \param [in] geom_create_fct Customize \ref PDM_writer_geom_create function for the new format (or NULL)
  !! \param [in] geom_write_fct  Customize \ref PDM_writer_geom_write function for the new format
  !! \param [in] geom_free_fct   Customize \ref PDM_writer_geom_free function for the new format (or NULL)
  !! \param [in] var_create_fct  Customize \ref PDM_writer_var_create function for the new format (or NULL)
  !! \param [in] var_write_fct   Customize \ref PDM_writer_var_write function for the new format
  !! \param [in] var_free_fct    Customize \ref PDM_writer_var_free function for the new format (or NULL)
  !!
  !!

  subroutine PDM_writer_fmt_add (name,            &
                                 create_fct,      &
                                 free_fct,        &
                                 beg_step_fct,    &
                                 end_step_fct,    &
                                 geom_create_fct, &
                                 geom_write_fct,  &
                                 geom_free_fct,   &
                                 var_create_fct,  &
                                 var_write_fct,   &
                                 var_free_fct)
    use iso_c_binding
    implicit none

    character (len=*)    :: name
    procedure(), pointer :: create_fct
    procedure(), pointer :: free_fct
    procedure(), pointer :: beg_step_fct
    procedure(), pointer :: end_step_fct
    procedure(), pointer :: geom_create_fct
    procedure(), pointer :: geom_write_fct
    procedure(), pointer :: geom_free_fct
    procedure(), pointer :: var_create_fct
    procedure(), pointer :: var_write_fct
    procedure(), pointer :: var_free_fct

    type(c_funptr)       :: c_create_fct
    type(c_funptr)       :: c_free_fct
    type(c_funptr)       :: c_beg_step_fct
    type(c_funptr)       :: c_end_step_fct
    type(c_funptr)       :: c_geom_create_fct
    type(c_funptr)       :: c_geom_write_fct
    type(c_funptr)       :: c_geom_free_fct
    type(c_funptr)       :: c_var_create_fct
    type(c_funptr)       :: c_var_write_fct
    type(c_funptr)       :: c_var_free_fct

    interface
      subroutine PDM_writer_fmt_add_c (name,            &
                                       create_fct,      &
                                       free_fct,        &
                                       beg_step_fct,    &
                                       end_step_fct,    &
                                       geom_create_fct, &
                                       geom_write_fct,  &
                                       geom_free_fct,   &
                                       var_create_fct,  &
                                       var_write_fct,   &
                                       var_free_fct)    &
      bind (c, name='PDM_writer_fmt_add')
        use iso_c_binding
        implicit none

        character(c_char) :: name(*)
        type(c_funptr)    :: create_fct
        type(c_funptr)    :: free_fct
        type(c_funptr)    :: beg_step_fct
        type(c_funptr)    :: end_step_fct
        type(c_funptr)    :: geom_create_fct
        type(c_funptr)    :: geom_write_fct
        type(c_funptr)    :: geom_free_fct
        type(c_funptr)    :: var_create_fct
        type(c_funptr)    :: var_write_fct
        type(c_funptr)    :: var_free_fct

      end subroutine PDM_writer_fmt_add_c
    end interface

    c_create_fct      = c_funloc(create_fct)
    c_free_fct        = c_funloc(free_fct)
    c_beg_step_fct    = c_funloc(beg_step_fct)
    c_end_step_fct    = c_funloc(end_step_fct)
    c_geom_create_fct = c_funloc(geom_create_fct)
    c_geom_write_fct  = c_funloc(geom_write_fct)
    c_geom_free_fct   = c_funloc(geom_free_fct)
    c_var_create_fct  = c_funloc(var_create_fct)
    c_var_write_fct   = c_funloc(var_write_fct)
    c_var_free_fct    = c_funloc(var_free_fct)

    call PDM_writer_fmt_add_c (trim(name)//C_NULL_CHAR, &
                               c_create_fct,      &
                               c_free_fct,        &
                               c_beg_step_fct,    &
                               c_end_step_fct,    &
                               c_geom_create_fct, &
                               c_geom_write_fct,  &
                               c_geom_free_fct,   &
                               c_var_create_fct,  &
                               c_var_write_fct,   &
                               c_var_free_fct)

  end subroutine PDM_writer_fmt_add

end module pdm_writer
