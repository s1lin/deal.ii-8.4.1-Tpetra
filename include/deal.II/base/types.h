// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2015 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef dealii__types_h
#define dealii__types_h


#include <deal.II/base/config.h>
#include <cstddef>
//MPI:
#  ifdef DEAL_II_WITH_MPI
#	 include <Tpetra_MpiPlatform.hpp>
# else
#    include <Tpetra_SerialPlatform.hpp>
#    include <Teuchos_DefaultComm.hpp>
#  endif

#include <Tpetra_Map_decl.hpp>
#include <Tpetra_Vector_decl.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Vector_decl.hpp>
#include <Tpetra_Version.hpp>
#include <Tpetra_CrsMatrix_decl.hpp>
#include <Tpetra_CrsGraph_decl.hpp>
#include <Tpetra_RowMatrix_decl.hpp>
#include <Tpetra_MultiVector_decl.hpp>
#include <Tpetra_Operator.hpp>
#include <Tpetra_Import_decl.hpp>
#include <Tpetra_CombineMode.hpp>

#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>//tefos
#include <Teuchos_OrdinalTraits.hpp>

#include <MueLu.hpp>
#include <MueLu_TpetraOperator.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>
#include <MueLu_Utilities.hpp>

#include <BelosIteration.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosConfigDefs.hpp>
#include <BelosSolverManager.hpp>
#include <BelosSolverFactory.hpp>

#include <Amesos.h>



DEAL_II_NAMESPACE_OPEN

/**
 * A namespace in which we declare typedefs for types used in deal.II, as well
 * as special values for these types.
 */
namespace types
{
  /**
   * The type used to denote subdomain_ids of cells.
   *
   * See the
   * @ref GlossSubdomainId "glossary"
   * for more information.
   *
   * There is a special value, numbers::invalid_subdomain_id that is used to
   * indicate an invalid value of this type.
   */
  typedef unsigned int subdomain_id;

  /**
   * The type used for global indices of vertices.
   */
  typedef unsigned long long int global_vertex_index;

  /**
   * An identifier that denotes the MPI type associated with
   * types::global_vertex_index.
   */
#  define DEAL_II_VERTEX_INDEX_MPI_TYPE MPI_UNSIGNED_LONG_LONG

#ifdef DEAL_II_WITH_64BIT_INDICES
  /**
   * The type used for global indices of degrees of freedom. While in
   * sequential computations the 4 billion indices of 32-bit unsigned integers
   * is plenty, parallel computations using the
   * parallel::distributed::Triangulation class can overflow this number and
   * we need a bigger index space.
   *
   * The data type always indicates an unsigned integer type.
   *
   * See the
   * @ref GlobalDoFIndex
   * page for guidance on when this type should or should not be used.
   */
  // TODO: we should check that unsigned long long int
  // has the same size as uint64_t
  typedef unsigned long long int global_dof_index;

  /**
   * An identifier that denotes the MPI type associated with
   * types::global_dof_index.
   */
#  define DEAL_II_DOF_INDEX_MPI_TYPE MPI_UNSIGNED_LONG_LONG
#else
  /**
   * The type used for global indices of degrees of freedom. While in
   * sequential computations the 4 billion indices of 32-bit unsigned integers
   * is plenty, parallel computations using the
   * parallel::distributed::Triangulation class can overflow this number and
   * we need a bigger index space.
   *
   * The data type always indicates an unsigned integer type.
   */
  typedef unsigned int global_dof_index;

  /**
   * An identifier that denotes the MPI type associated with
   * types::global_dof_index.
   */
#  define DEAL_II_DOF_INDEX_MPI_TYPE MPI_UNSIGNED
#endif

  /**
   * The type used to denote boundary indicators associated with every piece
   * of the boundary and, in the case of meshes that describe manifolds in
   * higher dimensions, associated with every cell.
   *
   * There is a special value, numbers::internal_face_boundary_id that is used
   * to indicate an invalid value of this type and that is used as the
   * boundary indicator for faces that are in the interior of the domain and
   * therefore not part of any addressable boundary component.
   *
   * @see
   * @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  typedef unsigned char boundary_id;

  /**
   * The type used to denote manifold indicators associated with every object
   * of the mesh.
   *
   * There is a special value, numbers::flat_manifold_id that is used to
   * indicate the standard cartesian manifold.
   *
   * @see
   * @ref GlossManifoldIndicator "Glossary entry on manifold indicators"
   */
  typedef unsigned int manifold_id;

  /**
   * The type used to denote material indicators associated with every cell.
   *
   * There is a special value, numbers::invalid_material_id that is used to
   * indicate an invalid value of this type.
   */
  typedef unsigned char material_id;
}

namespace TrilinosWrappers
{
  namespace types
  {
      typedef int int_type;
      typedef double                                      SC;
      typedef int                                         LO;
      typedef int                                         GO;
      typedef KokkosClassic::DefaultNode::DefaultNodeType NT;


      typedef Tpetra::Map<SC, LO, GO>                     map_type;
      typedef Tpetra::Vector<SC, LO, GO, NT>              vector_type;
      typedef Tpetra::RowMatrix<SC, LO, GO, NT>           row_matrix_type;
      typedef Tpetra::Export<SC, LO, GO>                  export_type;
      typedef Tpetra::Import<SC, LO, GO>                  import_type;
      typedef Tpetra::CrsMatrix<SC, LO, GO, NT>           crs_matrix_type;
      typedef Tpetra::CrsGraph<SC, LO, GO>                crs_graph_type;
      typedef Tpetra::MultiVector<SC, LO, GO, NT>         multi_vector_type;
      typedef Tpetra::Operator<SC, LO, GO, NT>            operator_type;

      typedef MueLu::MLParameterListInterpreter<SC, LO, GO, NT> ML_ParameterListInterpreter;

      typedef Teuchos::RCP<Teuchos::Comm<int>> comm_type;

      typedef Tpetra::MpiPlatform<int>    Tpetra_MpiComm;
      typedef Tpetra::SerialPlatform<int> Tpetra_SerialComm;

      typedef Xpetra::CrsMatrix<SC, LO, GO, NT> xcrs_matrix_type;
      typedef Xpetra::Matrix<SC, LO, GO, NT> xmatrix_type;
      typedef Xpetra::CrsMatrixWrap<SC, LO, GO, NT> xcrs_matrix_wrap;
  }

}


// this part of the namespace numbers got moved to the bottom types.h file,
// because otherwise we get a circular inclusion of config.h, types.h, and
// numbers.h
namespace numbers
{
  /**
   * Representation of the largest number that can be put into an unsigned
   * integer. This value is widely used throughout the library as a marker for
   * an invalid unsigned integer value, such as an invalid array index, an
   * invalid array size, and the like.
   */
  static const unsigned int
  invalid_unsigned_int = static_cast<unsigned int> (-1);

  /**
   * Representation of the largest number that can be put into a size_type.
   * This value is used throughout the library as a marker for an invalid
   * size_type value, such as an invalid array index, an invalid array size,
   * and the like. Invalid_size_type is equivalent to invalid_dof_index.
   */
  const types::global_dof_index
  invalid_size_type = static_cast<types::global_dof_index> (-1);

  /**
   * An invalid value for indices of degrees of freedom.
   */
  const types::global_dof_index invalid_dof_index = static_cast<types::global_dof_index>(-1);

  /**
   * Invalid material_id which we need in several places as a default value.
   * We assume that all material_ids lie in the range [0,
   * invalid_material_id).
   */
  const types::material_id invalid_material_id = static_cast<types::material_id>(-1);

  /**
   * Invalid boundary_id which we need in several places as a default value.
   * We assume that all valid boundary_ids lie in the range [0,
   * invalid_boundary_id).
   *
   * @see
   * @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  const types::boundary_id invalid_boundary_id = static_cast<types::boundary_id>(-1);

  /**
   * A boundary indicator number that we reserve for internal faces.  We
   * assume that all valid boundary_ids lie in the range [0,
   * internal_face_boundary_id).
   *
   * This is an indicator that is used internally (by the library) to
   * differentiate between faces that lie at the boundary of the domain and
   * faces that lie in the interior of the domain. You should never try to
   * assign this boundary indicator to anything in user code.
   *
   * @see
   * @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  const types::boundary_id internal_face_boundary_id = static_cast<types::boundary_id>(-1);

  /**
   * Invalid manifold_id which we need in several places as a default value.
   * We assume that all valid manifold_ids lie in the range [0,
   * invalid_manifold_id).
   *
   * @see
   * @ref GlossManifoldIndicator "Glossary entry on manifold indicators"
   */
  const types::manifold_id invalid_manifold_id = static_cast<types::manifold_id>(-1);

  /**
   * A manifold_id we reserve for the default flat Cartesian manifold.
   *
   * @see
   * @ref GlossManifoldIndicator "Glossary entry on manifold indicators"
   */
  const types::manifold_id flat_manifold_id = static_cast<types::manifold_id>(-1);

  /**
   * A special id for an invalid subdomain id. This value may not be used as a
   * valid id but is used, for example, for default arguments to indicate a
   * subdomain id that is not to be used.
   *
   * See the
   * @ref GlossSubdomainId "glossary"
   * for more information.
   */
  const types::subdomain_id invalid_subdomain_id = static_cast<types::subdomain_id>(-1);

  /**
   * The subdomain id assigned to a cell whose true subdomain id we don't
   * know, for example because it resides on a different processor on a mesh
   * that is kept distributed on many processors. Such cells are called
   * "artificial".
   *
   * See the glossary entries on
   * @ref GlossSubdomainId "subdomain ids"
   * and
   * @ref GlossArtificialCell "artificial cells"
   * as well as the
   * @ref distributed
   * module for more information.
   */
  const types::subdomain_id artificial_subdomain_id = static_cast<types::subdomain_id>(-2);
}


DEAL_II_NAMESPACE_CLOSE

#endif
