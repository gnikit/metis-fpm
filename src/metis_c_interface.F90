module metis_c_interface
   use iso_c_binding
   implicit none
   private

#if IDXTYPEWIDTH == 32
   integer, parameter, public :: idx_t = c_int32_t
#elif IDXTYPEWIDTH == 64
   integer, parameter, public :: idx_t = c_int64_t
#else
#error "Incorrect user-supplied value fo IDXTYPEWIDTH"
#endif


#if REALTYPEWIDTH == 32
   integer, parameter, public :: real_t = c_float
#elif REALTYPEWIDTH == 64
   integer, parameter, public :: real_t = c_double
#else
#error "Incorrect user-supplied value for REALTYPEWIDTH"
#endif

   integer(idx_t), parameter, public :: METIS_VER_MAJOR = 5
   integer(idx_t), parameter, public :: METIS_VER_MINOR = 1
   integer(idx_t), parameter, public :: METIS_VER_SUBMINOR = 0

   ! The maximum length of the options[] array
   integer(idx_t), parameter, public :: METIS_NOPTIONS = 40

   ! Return codes ! rstatus_et
   integer(idx_t), public :: METIS_OK              = 1
   integer(idx_t), public :: METIS_ERROR_INPUT     = -2
   integer(idx_t), public :: METIS_ERROR_MEMORY    = -3
   integer(idx_t), public :: METIS_ERROR           = -4

   ! Operation type codes ! moptype_et
   integer(idx_t), public :: METIS_OP_PMETIS = 0
   integer(idx_t), public :: METIS_OP_KMETIS = 1
   integer(idx_t), public :: METIS_OP_OMETIS = 2

   ! Options codes (i.e., options[]) ! moptions_et
   integer(idx_t), public :: METIS_OPTION_PTYPE = 0
   integer(idx_t), public :: METIS_OPTION_OBJTYPE = 1
   integer(idx_t), public :: METIS_OPTION_CTYPE = 2
   integer(idx_t), public :: METIS_OPTION_IPTYPE = 3
   integer(idx_t), public :: METIS_OPTION_RTYPE = 4
   integer(idx_t), public :: METIS_OPTION_DBGLVL = 5
   integer(idx_t), public :: METIS_OPTION_NIPARTS = 6
   integer(idx_t), public :: METIS_OPTION_NITER = 7
   integer(idx_t), public :: METIS_OPTION_NCUTS = 8
   integer(idx_t), public :: METIS_OPTION_SEED = 9
   integer(idx_t), public :: METIS_OPTION_NO2HOP = 10
   integer(idx_t), public :: METIS_OPTION_ONDISK = 11
   integer(idx_t), public :: METIS_OPTION_MINCONN = 12
   integer(idx_t), public :: METIS_OPTION_CONTIG = 13
   integer(idx_t), public :: METIS_OPTION_COMPRESS = 14
   integer(idx_t), public :: METIS_OPTION_CCORDER = 15
   integer(idx_t), public :: METIS_OPTION_PFACTOR = 16
   integer(idx_t), public :: METIS_OPTION_NSEPS = 17
   integer(idx_t), public :: METIS_OPTION_UFACTOR = 18
   integer(idx_t), public :: METIS_OPTION_NUMBERING = 19
   integer(idx_t), public :: METIS_OPTION_DROPEDGES = 20

   ! Used for command-line parameter purposes
   integer(idx_t), public :: METIS_OPTION_HELP = 21
   integer(idx_t), public :: METIS_OPTION_TPWGTS = 22
   integer(idx_t), public :: METIS_OPTION_NCOMMON = 23
   integer(idx_t), public :: METIS_OPTION_NOOUTPUT = 24
   integer(idx_t), public :: METIS_OPTION_BALANCE = 25
   integer(idx_t), public :: METIS_OPTION_GTYPE = 26
   integer(idx_t), public :: METIS_OPTION_UBVE = 27

   ! Partitioning Schemes ! mptype_et
   integer(idx_t), public :: METIS_PTYPE_RB = 0
   integer(idx_t), public :: METIS_PTYPE_KWAY = 1

   ! Graph types for meshes ! mgtype_et
   integer(idx_t), public :: METIS_GTYPE_DUAL = 0
   integer(idx_t), public :: METIS_GTYPE_NODAL = 1

   ! Coarsening Schemes ! mctype_et
   integer(idx_t), public :: METIS_CTYPE_RM = 0
   integer(idx_t), public :: METIS_CTYPE_SHEM = 1

   ! Initial partitioning schemes ! miptype_et
   integer(idx_t), public :: METIS_IPTYPE_GROW = 0
   integer(idx_t), public :: METIS_IPTYPE_RANDOM = 1
   integer(idx_t), public :: METIS_IPTYPE_EDGE = 2
   integer(idx_t), public :: METIS_IPTYPE_NODE = 3
   integer(idx_t), public :: METIS_IPTYPE_METISRB = 4

   ! Refinement schemes ! mrtype_et
   integer(idx_t), public :: METIS_RTYPE_FM = 0
   integer(idx_t), public :: METIS_RTYPE_GREEDY = 1
   integer(idx_t), public :: METIS_RTYPE_SEP2SIDED = 2
   integer(idx_t), public :: METIS_RTYPE_SEP1SIDED = 3

   ! Debug Levels ! mdbglvl_et
   integer(idx_t), public :: METIS_DBG_INFO       = 1    !< Shows various diagnostic messages
   integer(idx_t), public :: METIS_DBG_TIME       = 2    !< Perform timing analysis
   integer(idx_t), public :: METIS_DBG_COARSEN    = 4	 !< Show the coarsening progress
   integer(idx_t), public :: METIS_DBG_REFINE     = 8	 !< Show the refinement progress
   integer(idx_t), public :: METIS_DBG_IPART      = 16 	 !< Show info on initial partitioning
   integer(idx_t), public :: METIS_DBG_MOVEINFO   = 32 	 !< Show info on vertex moves during refinement
   integer(idx_t), public :: METIS_DBG_SEPINFO    = 64 	 !< Show info on vertex moves during sep refinement
   integer(idx_t), public :: METIS_DBG_CONNINFO   = 128  !< Show info on minimization of subdomain connectivity
   integer(idx_t), public :: METIS_DBG_CONTIGINFO = 256  !< Show info on elimination of connected components
   integer(idx_t), public :: METIS_DBG_MEMORY     = 2048 !< Show info related to wspace allocation

   ! Types of objectives ! mobjtype_et
   integer(idx_t), public :: METIS_OBJTYPE_CUT = 0
   integer(idx_t), public :: METIS_OBJTYPE_VOL = 1
   integer(idx_t), public :: METIS_OBJTYPE_NODE = 2

   ! C API
   !=======================================================================
   ! Contains all C functions, with exact same names as in METIS to provide
   ! an explicit interface to the C library.

   ! Graph partitioning routines
   public :: METIS_PartGraphRecursive
   public :: METIS_PartGraphKway

   ! Mesh-to-graph conversion routines
   public :: METIS_MeshToDual
   public :: METIS_MeshToNodal

   ! Mesh partitioning routines
   public :: METIS_PartMeshDual
   public :: METIS_PartMeshNodal

   ! Sparse matrix reordering routines
   public :: METIS_NodeND

   ! Utility routines
   public :: METIS_Free
   public :: METIS_SetDefaultOptions

   interface

   !> Recursive partitioning routine.
   !>
   !> This function computes a partitioning of a graph based on multilevel
   !> recursive bisection. It can be used to partition a graph into \e k
   !> parts. The objective of the partitioning is to minimize the edgecut
   !> subject to one or more balancing constraints.
   function METIS_PartGraphRecursive(nvtxs, ncon, xadj, adjncy, vwgt, vsize, adjwgt, &
                                     nparts, tpwgts, ubvec, options, edgecut, part) &
                                     bind(C, name="METIS_PartGraphRecursive") result(stat)
      import :: idx_t, real_t, METIS_NOPTIONS
      !> The number of vertices in the graph.
      integer(idx_t), intent(in) :: nvtxs
      !> The number of balancing constraints. For the standard
      !> partitioning problem in which each vertex is either unweighted
      !> or has a single weight, ncon should be 1.
      integer(idx_t), intent(in) :: ncon
      !> An array of size `nvtxs+1` used to specify the starting
      !> positions of the adjacency structure of the vertices in the
      !> `adjncy` array.
      integer(idx_t), intent(in) :: xadj(nvtxs+1)
      !> An array of size to the sum of the degrees of the
      !> graph that stores for each vertex the set of vertices that
      !> is adjancent to.
      integer(idx_t), intent(in) :: adjncy(*)
      !> An array of size `dim(nvtxs,ncon)` that stores the weights
      !> of the vertices for each constraint. The ncon weights for the
      !> ith vertex are stored in the ncon consecutive locations starting
      !> at `vwgt[i*ncon]`. When `ncon==1`, a `NULL` value can be passed indicating
      !> that all the vertices in the graph have the same weight.
      integer(idx_t), intent(in), optional :: vwgt(nvtxs,ncon)
      !> The size of the vertices for computing the total communication
      !> volume as described in Section 5.7
      integer(idx_t), intent(in), optional :: vsize(nvtxs)
      !> An array of size equal to `adjncy`, specifying the weight
      !> for each edge (i.e., `adjwgt[j]` corresponds to the weight of the
      !> edge stored in `adjncy[j]`).
      !> A `NULL` value can be passed indicating that all the edges in the
      !> graph have the same weight.
      integer(idx_t), intent(in) :: adjwgt(*)
      !> The number of desired partitions.
      integer(idx_t), intent(in) :: nparts
      !> An array of size `dim(nparts,ncon)` that specifies the
      !> desired weight for each part and constraint. The \e{target partition
      !> weight} for the ith part and jth constraint is specified
      !> at tpwgts[i*ncon+j] (the numbering of `i` and `j` starts from 0).
      !> For each constraint, the sum of the tpwgts[] entries must be
      !> 1.0 (i.e., \f$ \sum_i tpwgts[i*ncon+j] = 1.0 \f$).
      !> A `NULL` value can be passed indicating that the graph should
      !> be equally divided among the parts.
      real(real_t), intent(in), optional :: tpwgts(nparts,ncon)
      !> An array of size ncon that specifies the allowed
      !> load imbalance tolerance for each constraint.
      !> For the ith part and jth constraint the allowed weight is the
      !> `ubvec[j]*tpwgts[i*ncon+j]` fraction of the jth's constraint total
      !> weight. The load imbalances must be greater than 1.0.
      !> A `NULL` value can be passed indicating that the load imbalance
      !> tolerance for each constraint should be 1.001 (for ncon==1)
      !> or 1.01 (for ncon>1).
      real(real_t), intent(in), optional :: ubvec(ncon)
      !> The array for passing additional parameters
      !> in order to customize the behaviour of the partitioning
      !> algorithm.
      integer(idx_t), intent(in), optional :: options(METIS_NOPTIONS)
      !> Stores the cut of the partitioning.
      integer(idx_t), intent(out) :: edgecut
      !> An array of size nvtxs used to store the
      !> computed partitioning. The partition number for the ith
      !> vertex is stored in `part[i]`. Based on the numflag parameter,
      !> the numbering of the parts starts from either 0 or 1.
      integer(idx_t), intent(out) :: part(nvtxs)
      !> `METIS_OK`  indicates that the function returned normally.
      !> `METIS_ERROR_INPUT` indicates an input error.
      !> `METIS_ERROR_MEMORY` indicates that it could not allocate the required memory.
      integer(idx_t) :: stat
   end function METIS_PartGraphRecursive

   !> Multilevel k-way partitioning routine to partition a graph into k parts.
   function METIS_PartGraphKway(nvtxs, ncon, xadj, adjncy, vwgt, vsize, adjwgt, nparts, &
                                tpwgts, ubvec, options, edgecut, part) &
                                bind(C, name="METIS_PartGraphKway") result(stat)
      import :: idx_t, real_t, METIS_NOPTIONS
      !> The number of vertices in the graph.
      integer(idx_t), intent(in) :: nvtxs
      !> The number of balancing constraints. For the standard
      !> partitioning problem in which each vertex is either unweighted
      !> or has a single weight, ncon should be 1.
      integer(idx_t), intent(in) :: ncon
      !> An array of size nvtxs+1 used to specify the starting
      !> positions of the adjacency structure of the vertices in the
      !> adjncy array.
      integer(idx_t), intent(in) :: xadj(nvtxs+1)
      !> Array of size to the sum of the degrees of the
      !> graph that stores for each vertex the set of vertices that
      !> is adjancent to.
      integer(idx_t), intent(in) :: adjncy(*)
      !> Array of size `(nvtxs,ncon)` that stores the weights
      !> of the vertices for each constraint. The ncon weights for the
      !> ith vertex are stored in the ncon consecutive locations starting
      !> at vwgt[i*ncon]. When ncon==1, a NULL value can be passed indicating
      !> that all the vertices in the graph have the same weight.
      integer(idx_t), intent(in), optional :: vwgt(nvtxs,ncon)
      !> The size of the vertices for computing the total communication
      !> volume as described in Section 5.7
      integer(idx_t), intent(in), optional :: vsize(nvtxs)
      !> An array of size equal to adjncy, specifying the weight
      !> for each edge (i.e., `adjwgt[j]` corresponds to the weight of the
      !> edge stored in `adjncy[j]`).
      !> A `NULL` value can be passed indicating that all the edges in the
      !> graph have the same weight.
      integer(idx_t), intent(in), optional :: adjwgt(*)
      !> The number of desired partitions.
      integer(idx_t), intent(in) :: nparts
      !> An array of size (nparts,ncon) that specifies the
      !> desired weight for each part and constraint. The \e{target partition
      !> weight} for the ith part and jth constraint is specified
      !> at `tpwgts[i*ncon+j]` (the numbering of i and j starts from 0).
      !> For each constraint, the sum of the tpwgts[] entries must be
      !> 1.0 (i.e., \f$ \sum_i tpwgts[i*ncon+j] = 1.0 \f$).
      !> A `NULL` value can be passed indicating that the graph should
      !> be equally divided among the parts.
      real(real_t), intent(in), optional :: tpwgts(nparts,ncon)
      !> An array of size ncon that specifies the allowed
      !> load imbalance tolerance for each constraint.
      !> For the ith part and jth constraint the allowed weight is the
      !> `ubvec[j]*tpwgts[i*ncon+j]` fraction of the jth's constraint total
      !> weight. The load imbalances must be greater than 1.0.
      !> A `NULL` value can be passed indicating that the load imbalance
      !> tolerance for each constraint should be 1.001 (for ncon==1)
      !> or 1.01 (for ncon>1).
      real(real_t), intent(in), optional :: ubvec(ncon)
      !> The array for passing additional parameters
      !> in order to customize the behaviour of the partitioning
      !> algorithm.
      integer(idx_t), intent(in), optional :: options(METIS_NOPTIONS)
      !> Stores the cut of the partitioning.
      integer(idx_t), intent(out) :: edgecut
      !> An array of size nvtxs used to store the
      !> computed partitioning. The partition number for the ith
      !> vertex is stored in part[i]. Based on the numflag parameter,
      !> the numbering of the parts starts from either 0 or 1.
      integer(idx_t), intent(out) :: part(nvtxs)
      !> `METIS_OK`  indicates that the function returned normally.
      !> `METIS_ERROR_INPUT` indicates an input error.
      !> `METIS_ERROR_MEMORY` indicates that it could not allocate the required memory.
      integer(idx_t) :: stat
   end function METIS_PartGraphKway

   !> This function creates a graph corresponding to the dual of a finite element
   !> mesh.
   function METIS_MeshToDual(ne, nn, eptr, eind, ncommon, numflag, r_xadj, r_adjncy) &
                             bind(C, name="METIS_MeshToDual") result(stat)
      import :: idx_t
      !> ne is the number of elements in the mesh.
      integer(idx_t), intent(in) :: ne
      !> nn is the number of nodes in the mesh.
      integer(idx_t), intent(in) :: nn
      !> eptr is an array of size ne+1 used to mark the start and end
      !> locations in the nind array.
      integer(idx_t), intent(in) :: eptr(ne+1)
      !> eind is an array that stores for each element the set of node IDs
      !> (indices) that it is made off. The length of this array is equal
      !> to the total number of nodes over all the mesh elements.
      integer(idx_t), intent(in) :: eind(*)
      !> ncommon is the minimum number of nodes that two elements must share
      !> in order to be connected via an edge in the dual graph.
      integer(idx_t), intent(in) :: ncommon
      !> numflag is either 0 or 1 indicating if the numbering of the nodes
      !> starts from 0 or 1, respectively. The same numbering is used for the
      !> returned graph as well.
      integer(idx_t), intent(in) :: numflag
      !> r_xadj indicates where the adjacency list of each vertex is stored
      !> in r_adjncy. The memory for this array is allocated by this routine.
      !> It can be freed by calling Free().
      integer(idx_t), intent(out) :: r_xadj(nn+1)
      !> r_adjncy stores the adjacency list of each vertex in the generated
      !> dual graph. The memory for this array is allocated by this routine.
      !> It can be freed by calling Free().
      integer(idx_t), intent(out) :: r_adjncy(*)
      !> `METIS_OK` Indicates that the function returned normally.
      !> `METIS_ERROR_INPUT` Indicates an input error.
      !> `METIS_ERROR_MEMORY` Indicates that it could not allocate the required memory.
      !> `METIS_ERROR` Indicates some other type of error.
      integer(idx_t) :: stat
   end function METIS_MeshToDual

   !> This function creates a graph corresponding to (almost) the nodal of a
   !> finite element mesh. In the nodal graph, each node is connected to the
   !> nodes corresponding to the union of nodes present in all the elements
   !> in which that node belongs.
   function METIS_MeshToNodal(ne, nn, eptr, eind, numflag, r_xadj, r_adjncy) &
                              bind(C, name="METIS_MeshToNodal") result(stat)
      import :: idx_t
      !> ne is the number of elements in the mesh.
      integer(idx_t), intent(in) :: ne
      !> nn is the number of nodes in the mesh.
      integer(idx_t), intent(in) :: nn
      !> eptr is an array of size ne+1 used to mark the start and end
      !> locations in the nind array.
      integer(idx_t), intent(in) :: eptr(ne+1)
      !> eind is an array that stores for each element the set of node IDs
      !> (indices) that it is made off. The length of this array is equal
      !> to the total number of nodes over all the mesh elements.
      integer(idx_t), intent(in) :: eind(*)
      !> numflag is either 0 or 1 indicating if the numbering of the nodes
      !> starts from 0 or 1, respectively. The same numbering is used for the
      !> returned graph as well.
      integer(idx_t), intent(in) :: numflag
      !> r_xadj indicates where the adjacency list of each vertex is stored
      !> in r_adjncy. The memory for this array is allocated by this routine.
      !> It can be freed by calling METIS_free().
      integer(idx_t), intent(out) :: r_xadj(ne+1)
      !> r_adjncy stores the adjacency list of each vertex in the generated
      !> dual graph. The memory for this array is allocated by this routine.
      !> It can be freed by calling METIS_free().
      integer(idx_t), intent(out) :: r_adjncy(*)
      !> `METIS_OK` - Indicates that the function returned normally
      !> `METIS_ERROR_INPUT` - Indicates an input error.
      !> `METIS_ERROR_MEMORY` - Indicates that it could not allocate the required memory.
      !> `METIS_ERROR` - Indicates some other type of error.
      integer(idx_t) :: stat
   end function METIS_MeshToNodal

   !> This function is used to partition a mesh into k parts based on a partitioning of the mesh’s dual graph
   function METIS_PartMeshNodal(ne, nn, eptr, eind, vwgt, vsize, nparts, tpwgts, &
                                options, objval, epart, npart) &
                                bind(C, name="METIS_PartMeshNodal") result(stat)
      import :: idx_t, real_t, METIS_NOPTIONS
      !> The number of elements in the mesh.
      integer(idx_t), intent(in) :: ne
      !> The number of nodes in the mesh.
      integer(idx_t), intent(in) :: nn
      !> An array of size ne+1 used to mark the start and end
      integer(idx_t), intent(in) :: eptr(ne+1)
      !> An array that stores for each element the set of node IDs
      integer(idx_t), intent(in) :: eind(*)
      !> An array of size nn specifying the weights of the nodes.
      !> A `NULL` value can be passed to indicate that all nodes have an equal weight.
      integer(idx_t), intent(in), optional :: vwgt(nn)
      !> An array of size nn specifying the size of the nodes that is used for computing the total
      !> communication volume as described in Section 5.7.
      !> A `NULL` value can be passed when the objective is cut or when all nodes have an equal size
      integer(idx_t), intent(in), optional :: vsize(nn)
      !> The number of parts to partition the mesh.
      integer(idx_t), intent(in) :: nparts
      !> This is an array of size nparts that specifies the desired weight for each partition.
      !> The target partition weight for the ith partition is specified at `tpwgts[i]`
      !> (the numbering for the partitions starts from `0`).
      !> The sum of the `tpwgts[]` entries must be `1.0`.
      !> A `NULL` value can be passed to indicate that the graph should be equally divided among the partitions.
      real(real_t), intent(in), optional :: tpwgts(nparts)
      !> This is the array of options as described in Section 5.4. The following options are valid:
      !> `METIS_OPTION_PTYPE`, `METIS_OPTION_OBJTYPE`, `METIS_OPTION_CTYPE`,
      !> `METIS_OPTION_IPTYPE`, `METIS_OPTION_RTYPE`, `METIS_OPTION_NCUTS`,
      !> `METIS_OPTION_NITER`, `METIS_OPTION_SEED`, `METIS_OPTION_UFACTOR`,
      !> `METIS_OPTION_NUMBERING`, `METIS_OPTION_DBGLVL`
      integer(idx_t), intent(in), optional :: options(METIS_NOPTIONS)
      !> Upon successful completion, this variable stores either the edgecut or
      !> the total communication volume of the nodal graph’s partitioning.
      integer(idx_t), intent(out) :: objval
      !> This is a vector of size ne that upon successful completion stores the partition vector for the elements
      !> of the mesh. The numbering of this vector starts from either 0 or 1, depending on the value of
      !> `options[METIS_OPTION_NUMBERING]`.
      integer(idx_t), intent(out) :: epart(ne)
      !> This is a vector of size nn that upon successful completion stores the partition vector for the nodes
      !> of the mesh. The numbering of this vector starts from either 0 or 1, depending on the value of
      !> `options[METIS_OPTION_NUMBERING]`.
      integer(idx_t), intent(out) :: npart(nn)
      !> `METIS_OK` Indicates that the function returned normally.
      !> `METIS_ERROR_INPUT` Indicates an input error.
      !> `METIS_ERROR_MEMORY` Indicates that it could not allocate the required memory.
      !> `METIS_ERROR` Indicates some other type of error.
      integer(idx_t) :: stat
   end function METIS_PartMeshNodal

   !> This function is used to partition a mesh into k parts based on a partitioning of the mesh’s dual graph.
   function METIS_PartMeshDual(ne, nn, eptr, eind, vwgt, vsize, ncommon, nparts, &
                               tpwgts, options, objval, epart, npart) &
                               bind(C, name="METIS_PartMeshDual") result(stat)
      import :: idx_t, real_t, METIS_NOPTIONS
      !> The number of elements in the mesh.
      integer(idx_t), intent(in) :: ne
      !> The number of nodes in the mesh.
      integer(idx_t), intent(in) :: nn
      !> The pair of arrays storing the mesh as described in Section 5.6.
      integer(idx_t), intent(in) :: eptr(ne+1)
      integer(idx_t), intent(in) :: eind(*)
      !> An array of size `ne` specifying the weights of the elements.
      !> A `NULL` value can be passed to indicate that all elements have an equal weight.
      integer(idx_t), intent(in), optional :: vwgt(ne)
      !> An array of size `ne` specifying the size of the elements that is used for computing the total
      !> communication volume as described in Section 5.7.
      !>  A `NULL` value can be passed when the objective is cut or when all elements have an equal size
      integer(idx_t), intent(in), optional :: vsize(ne)
      !> Specifies the number of common nodes that two elements must have in order to put an edge between
      !> them in the dual graph. Given two elements e1 and e2, containing n1 and n2 nodes, respectively,
      !> then an edge will connect the vertices in the dual graph corresponding to e1 and e2 if the number of
      !> common nodes between them is greater than or equal to `min(ncommon, n1 − 1, n2 − 1)`.
      !> The default value is 1, indicating that two elements will be connected via an edge as long as they
      !> share one node. However, this will tend to create too many edges (increasing the memory and time
      !> requirements of the partitioning). The user should select higher values that are better suited for the
      !> element types of the mesh that wants to partition. For example, for tetrahedron meshes, ncommon
      !> should be 3, which creates an edge between two tets when they share a triangular face (i.e., 3 nodes).
      integer(idx_t), intent(in) :: ncommon
      !> The number of parts to partition the mesh.
      integer(idx_t), intent(in) :: nparts
      !> This is an array of size nparts that specifies the desired weight for each partition.
      !> The target partition weight for the ith partition is specified at `tpwgts[i]`
      !> (the numbering for the partitions starts from 0). The sum of the `tpwgts[]` entries must be 1.0.
      !> A `NULL` value can be passed to indicate that the graph should be equally divided among the partitions.
      real(real_t), intent(in), optional :: tpwgts(nparts)
      !> This is the array of options as described in Section 5.4. The following options are valid:
      !> `METIS_OPTION_PTYPE`, `METIS_OPTION_OBJTYPE`, `METIS_OPTION_CTYPE`,
      !> `METIS_OPTION_IPTYPE`, `METIS_OPTION_RTYPE`, `METIS_OPTION_NCUTS`,
      !> `METIS_OPTION_NITER`, `METIS_OPTION_SEED`, `METIS_OPTION_UFACTOR`,
      !> `METIS_OPTION_NUMBERING`, `METIS_OPTION_DBGLVL`
      integer(idx_t), intent(in), optional :: options(METIS_NOPTIONS)
      !> Upon successful completion, this variable stores either the edgecut
      !> or the total communication volume of the dual graph’s partitioning.
      integer(idx_t), intent(out) :: objval
      !> This is a vector of size `ne` that upon successful completion stores the partition vector for the elements
      !> of the mesh. The numbering of this vector starts from either 0 or 1, depending on the value of
      !> `options[METIS OPTION NUMBERING]`.
      integer(idx_t), intent(out) :: epart(ne)
      !> This is a vector of size `nn` that upon successful completion stores the partition vector for the nodes
      !> of the mesh. The numbering of this vector starts from either 0 or 1, depending on the value of
      !> `options[METIS OPTION NUMBERING]`.
      integer(idx_t), intent(out) :: npart(nn)
      !> `METIS_OK` Indicates that the function returned normally.
      !> `METIS_ERROR_INPUT` Indicates an input error.
      !> `METIS_ERROR_MEMORY` Indicates that it could not allocate the required memory.
      !> `METIS ERROR` Indicates some other type of error.
      integer(idx_t) :: stat
   end function METIS_PartMeshDual

   !> This function is the entry point for the multilevel nested dissection
   !> ordering code. At each bisection, a node-separator is computed using
   !> a node-based refinement approach.
   function METIS_NodeND(nvtxs, xadj, adjncy, vwgt, options, perm, iperm) &
                         bind(C, name="METIS_NodeND") result(stat)
      import :: idx_t, METIS_NOPTIONS
      !> The number of vertices in the graph.
      integer(idx_t), intent(in) :: nvtxs
      !> Of length `nvtxs+1` marking the start of the adjancy list of each vertex in adjncy.
      integer(idx_t), intent(in) :: xadj(nvtxs+1)
      !> Stores the adjacency lists of the vertices. The adjnacy
      !> list of a vertex should not contain the vertex itself.
      integer(idx_t), intent(in) :: adjncy(*)
      !> An array of size `nvtxs` storing the weight of each
      !> vertex. If vwgt is `NULL`, then the vertices are considered
      !> to have unit weight.
      integer(idx_t), intent(in) :: vwgt(nvtxs)
      !> An array of size `nvtxs` specifying the weights of the vertices.
      integer(idx_t), intent(in) :: options(METIS_NOPTIONS)
      !> An array of size nvtxs such that if `A` and `A'` are the original and
      !> permuted matrices, then `A'[i] = A[perm[i]]`.
      integer(idx_t), intent(out) :: perm(nvtxs)
      !> An array of size nvtxs such that if `A` and `A'` are the original and
      !> permuted matrices, then `A[i] = A'[iperm[i]]`.
      integer(idx_t), intent(out) :: iperm(nvtxs)
      !> `METIS_OK` - Indicates that the function returned normally.
      !> `METIS_ERROR_INPUT` - Indicates an input error.
      !> `METIS_ERROR_MEMORY` - Indicates that it could not allocate the required memory.
      !> `METIS_ERROR` - Indicates some other type of error.
      integer(idx_t) :: stat
   end function METIS_NodeND

   !> Frees the memory that was allocated by either the `METIS_MeshToDual` or the
   !> `METIS_MeshToNodal` routines for returning the dual or nodal graph of a mesh.
   function METIS_Free(ptr) bind(C, name="METIS_Free") result(stat)
      import :: c_ptr, idx_t
      !> ptr points to the memory that was previously allocated by METIS.
      type(c_ptr), value :: ptr
      !> `METIS_OK` - Indicates that the function returned normally.
      integer(idx_t) :: stat
   end function METIS_Free


   !> Initializes the options array into its default values.
   function METIS_SetDefaultOptions(options) bind(C, name="METIS_SetDefaultOptions") result(stat)
      import :: idx_t, METIS_NOPTIONS
      !> The array of options that will be initialized. It’s size should be at least `METIS_NOPTIONS`
      integer(idx_t), intent(out) :: options(METIS_NOPTIONS)
      !> `METIS_OK` Indicates that the function returned normally.
      integer(idx_t) :: stat
   end function METIS_SetDefaultOptions

   ! These functions are used by ParMETIS */
   ! function NodeNDP(nvtxs, xadj, adjncy, vwgt, npes, options, perm, iperm, sizes) &
   !                  bind(C, name="METIS_NodeNDP") result(stat)
   !    integer(idx_t), intent(in) :: nvtxs, xadj(nvtxs+1), adjncy(*), vwgt(nvtxs), npes, options(*)
   !    integer(idx_t), intent(out) :: perm(nvtxs), iperm(nvtxs), sizes(npes)
   !    integer(idx_t) :: stat
   ! end function NodeNDP

   ! function ComputeVertexSeparator(nvtxs, xadj, adjncy, vwgt, options, sepsize, part) &
   !                                 bind(C, name="METIS_ComputeVertexSeparator") result(stat)
   !    integer(idx_t), intent(in) :: nvtxs, xadj(nvtxs+1), adjncy(*), vwgt(nvtxs), options(*)
   !    integer(idx_t), intent(out) :: sepsize, part(nvtxs)
   !    integer(idx_t) :: stat
   ! end function ComputeVertexSeparator

   ! function NodeRefine(nvtxs, xadj, vwgt, adjncy, where, hmarker, ubfactor) &
   !                     bind(C, name="METIS_NodeRefine") result(stat)
   !    integer(idx_t), intent(in) :: nvtxs, xadj(nvtxs+1), vwgt(nvtxs), adjncy(*), where(nvtxs), hmarker(nvtxs)
   !    real(kind=real_t), intent(in) :: ubfactor
   !    integer(idx_t) :: stat
   ! end function NodeRefine

   ! These functions are used by DGL */
   ! function CacheFriendlyReordering(nvtxs, xadj, adjncy, part, old2new)
   !                                  bind(C, name="METIS_CacheFriendlyReordering") result(stat)
   !    integer(idx_t), intent(in) :: nvtxs, xadj(nvtxs+1), adjncy(*), part(nvtxs)
   !    integer(idx_t), intent(out) :: old2new(nvtxs)
   !    integer(idx_t) :: stat
   ! end function CacheFriendlyReordering

   end interface

end module metis_c_interface
