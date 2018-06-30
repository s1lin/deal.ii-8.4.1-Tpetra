// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2015 by the deal.II authors
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

#include <deal.II/lac/trilinos_vector.h>

#ifdef DEAL_II_WITH_TRILINOS

#  include <deal.II/lac/trilinos_sparse_matrix.h>
#  include <deal.II/lac/trilinos_block_vector.h>
#  include <Teuchos_RCPDecl.hpp>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#  include <Tpetra_Vector_decl.hpp>
#  include <Tpetra_Import_decl.hpp>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

#  include <cmath>



DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{
  namespace
  {
    // define a helper function that queries the size of an map_type object
    // by calling either the 32- or 64-bit function necessary, and returns the
    // result in the correct data type so that we can use it in calling other
    // Tpetra member functions that are overloaded by index type
    int n_global_elements (const map_type &map)
    {
      return map.getGlobalNumElements();
    }
    // define a helper function that queries the pointer to internal array
    // containing list of global IDs assigned to the calling processor
    // by calling either the 32- or 64-bit function necessary, and returns the
    // result in the correct data type so that we can use it in calling other
    // Tpetra member functions that are overloaded by index type
    int *my_global_elements(const map_type &map)
    {
      return static_cast<int *>(map.getNodeNumElements());//?
    }
    // define a helper function that queries the global vector length of an
    // vector_type object  by calling either the 32- or 64-bit
    // function necessary.
    int global_length(const vector_type &vector)
    {
      return vector.getGlobalLength();
    }
  }

  namespace MPI
  {


    Vector::Vector ()
    {
      last_action = Zero;
      RCP<const map_type> tempMap = rcp(new map_type(0,0, Utilities::Trilinos::comm_self()));
      vector.reset(new vector_type(tempMap));
    }



    Vector::Vector (const map_type &parallel_partitioning)
    {
      reinit (parallel_partitioning);
    }



    Vector::Vector (const IndexSet &parallel_partitioning,
                    const MPI_Comm &communicator)
    {
      reinit (parallel_partitioning, communicator);
    }



    Vector::Vector (const Vector &v)
      :
      VectorBase()
    {
      last_action = Zero;
      vector.reset (new vector_type(*v.vector));
      has_ghosts = v.has_ghosts;
    }



#ifdef DEAL_II_WITH_CXX11
    Vector::Vector (Vector &&v)
    {
      // initialize a minimal, valid object and swap
      last_action = Zero;
      vector.reset(new vector_type(map_type(0,0,0,Utilities::Trilinos::comm_self())));

      swap(v);
    }
#endif



    Vector::Vector (const map_type &input_map,
                    const VectorBase &v)
      :
      VectorBase()
    {
      AssertThrow (n_global_elements(input_map) == n_global_elements(v.vector->getMap().get()),
                   ExcDimensionMismatch (n_global_elements(input_map),
                                         n_global_elements(v.vector->getMap())));

      last_action = Zero;

      if (input_map.isSameAs(v.vector->getMap().get()) == true)
        vector.reset (new vector_type(*v.vector));
      else
        {
          RCP<const map_type> tempMap = rcp(new map_type(input_map));
          vector.reset (new vector_type(tempMap));
          reinit (v, false, true);
        }
    }



    Vector::Vector (const IndexSet   &parallel_partitioner,
                    const VectorBase &v,
                    const MPI_Comm   &communicator)
      :
      VectorBase()
    {
      AssertThrow (parallel_partitioner.size() ==
                   static_cast<size_type>(n_global_elements(v.vector->getMap())),
                   ExcDimensionMismatch (parallel_partitioner.size(),
                                         n_global_elements(v.vector->getMap())));

      last_action = Zero;

      vector.reset (new vector_type
                    (parallel_partitioner.make_trilinos_map(communicator,
                                                            true)));
      reinit (v, false, true);
    }

    Vector::Vector (const IndexSet &local,
                    const IndexSet &ghost,
                    const MPI_Comm &communicator)
      :
      VectorBase()
    {
      IndexSet parallel_partitioning = local;
      parallel_partitioning.add_indices(ghost);
      reinit(parallel_partitioning, communicator);
    }



    Vector::~Vector ()
    {}



    void
    Vector::reinit (const map_type &input_map,
                    const bool        omit_zeroing_entries)
    {
      nonlocal_vector.reset();

      if (vector->getMap().SameAs(input_map)==false)
        vector.reset (new vector_type(input_map));
      else if (omit_zeroing_entries == false)
        {
          vector->putScalar(0.);
        }

      has_ghosts = vector->getMap().UniqueGIDs()==false;
      last_action = Zero;
    }



    void
    Vector::reinit (const IndexSet &parallel_partitioner,
                    const MPI_Comm &communicator,
                    const bool      omit_zeroing_entries)
    {
      nonlocal_vector.reset();

      map_type map = parallel_partitioner.make_trilinos_map (communicator,
                                                               true);
      reinit (map, omit_zeroing_entries);
    }



    void
    Vector::reinit (const VectorBase &v,
                    const bool        omit_zeroing_entries,
                    const bool        allow_different_maps)
    {
      nonlocal_vector.reset();

      // In case we do not allow to have different maps, this call means that
      // we have to reset the vector. So clear the vector, initialize our map
      // with the map in v, and generate the vector.
      if (allow_different_maps == false)
        {
          if (vector->getMap().SameAs(v.vector->getMap()) == false)
            {
              vector.reset (new vector_type(v.vector->getMap()));
              has_ghosts = v.has_ghosts;
              last_action = Zero;
            }
          else if (omit_zeroing_entries == false)
            {
              // old and new vectors
              // have exactly the
              // same map, i.e. size
              // and parallel
              // distribution
              int ierr;
              ierr = vector->GlobalAssemble (last_action);
              (void)ierr;
              Assert (ierr == 0, ExcTrilinosError(ierr));

              vector->putScalar(0.0);

              last_action = Zero;
            }
        }

      // Otherwise, we have to check that the two vectors are already of the
      // same size, create an object for the data exchange and then insert all
      // the data. The first assertion is only a check whether the user knows
      // what she is doing.
      else
        {
          Assert (omit_zeroing_entries == false,
                  ExcMessage ("It is not possible to exchange data with the "
                              "option 'omit_zeroing_entries' set, which would not write "
                              "elements."));

          AssertThrow (size() == v.size(),
                       ExcDimensionMismatch (size(), v.size()));

          import_type data_exchange (vector->getMap(), v.vector->getMap());

          const int ierr = vector->doImport(*v.vector, data_exchange, Insert);
          AssertThrow (ierr == 0, ExcTrilinosError(ierr));

          last_action = Insert;
        }

    }



    void
    Vector::reinit (const BlockVector &v,
                    const bool         import_data)
    {
      nonlocal_vector.reset();

      // In case we do not allow to have different maps, this call means that
      // we have to reset the vector. So clear the vector, initialize our map
      // with the map in v, and generate the vector.
      if (v.n_blocks() == 0)
        return;

      // create a vector that holds all the elements contained in the block
      // vector. need to manually create an map_type.
      size_type n_elements = 0, added_elements = 0, block_offset = 0;
      for (size_type block=0; block<v.n_blocks(); ++block)
        n_elements += v.block(block).local_size();
      std::vector<TrilinosWrappers::types::int_type> global_ids (n_elements, -1);
      for (size_type block=0; block<v.n_blocks(); ++block)
        {
          TrilinosWrappers::types::int_type *glob_elements =
            my_global_elements(v.block(block).vector_partitioner());
          for (size_type i=0; i<v.block(block).local_size(); ++i)
            global_ids[added_elements++] = glob_elements[i] + block_offset;
          block_offset += v.block(block).size();
        }

      Assert (n_elements == added_elements, ExcInternalError());


      RCP<const map_type> new_map = rcp(new map_type(v.size(), n_elements, &global_ids[0], 0,
                                                     v.block(0).vector_partitioner().getComm()));

      std_cxx11::shared_ptr<vector_type> actual_vec;
      if ( import_data == true )
        actual_vec.reset (new vector_type (new_map));
      else
        {
          vector.reset (new vector_type (new_map));
          actual_vec = vector;
        }

      TrilinosScalar *entries = (*actual_vec)[0];
      block_offset = 0;
      for (size_type block=0; block<v.n_blocks(); ++block)
        {
          v.block(block).trilinos_vector().ExtractCopy (entries, 0);
          entries += v.block(block).local_size();
        }

      if (import_data == true)
        {
          AssertThrow (static_cast<size_type>(global_length(*actual_vec))
                       == v.size(),
                       ExcDimensionMismatch (global_length(*actual_vec),
                                             v.size()));

          import_type data_exchange (vector->getMap(), actual_vec->getMap());

          vector->doImport(*actual_vec, data_exchange, Insert);
          
          last_action = Insert;
        }

    }


    void Vector::reinit(const IndexSet &locally_owned_entries,
                        const IndexSet &ghost_entries,
                        const MPI_Comm &communicator,
                        const bool      vector_writable)
    {
      nonlocal_vector.reset();
      if (vector_writable == false)
        {
          IndexSet parallel_partitioning = locally_owned_entries;
          parallel_partitioning.add_indices(ghost_entries);
          reinit(parallel_partitioning, communicator);
        }
      else
        {
          map_type map = locally_owned_entries.make_trilinos_map (communicator,
                                                                    true);
          Assert (map.isOneToOne(),
                  ExcMessage("A writable vector must not have ghost entries in "
                             "its parallel partitioning"));
          reinit (map);

          IndexSet nonlocal_entries(ghost_entries);
          nonlocal_entries.subtract_set(locally_owned_entries);
          if (Utilities::MPI::n_mpi_processes(communicator) > 1)
            {
              map_type nonlocal_map =
                nonlocal_entries.make_trilinos_map(communicator, true);
              nonlocal_vector.reset(new Tpetra_MultiVector(nonlocal_map, 1));
            }
        }
    }


    Vector &
    Vector::operator = (const Vector &v)
    {
      // distinguish three cases. First case: both vectors have the same
      // layout (just need to copy the local data, not reset the memory and
      // the underlying map_type). The third case means that we have to
      // rebuild the calling vector.
      if (vector->getMap().SameAs(v.vector->getMap()))
        {
          *vector = *v.vector;
          if (v.nonlocal_vector.get() != 0)
            nonlocal_vector.reset(new Tpetra_MultiVector(v.nonlocal_vector->getMap(), 1));
          last_action = Zero;
        }
      // Second case: vectors have the same global
      // size, but different parallel layouts (and
      // one of them a one-to-one mapping). Then we
      // can call the import/export functionality.
      else if (size() == v.size() &&
               (v.vector->getMap().UniqueGIDs() || vector->getMap().UniqueGIDs()))
        {
          reinit (v, false, true);
        }
      // Third case: Vectors do not have the same
      // size.
      else
        {
          vector.reset (new vector_type(*v.vector));
          last_action = Zero;
          has_ghosts = v.has_ghosts;
        }

      if (v.nonlocal_vector.get() != 0)
        nonlocal_vector.reset(new Tpetra_MultiVector(v.nonlocal_vector->getMap(), 1));

      return *this;
    }



#ifdef DEAL_II_WITH_CXX11
    Vector &Vector::operator= (Vector &&v)
    {
      swap(v);
      return *this;
    }
#endif



    Vector &
    Vector::operator = (const TrilinosWrappers::Vector &v)
    {
      nonlocal_vector.reset();

      Assert (size() == v.size(), ExcDimensionMismatch(size(), v.size()));

      import_type data_exchange (vector->getMap(), v.vector->getMap());
      vector->doImport(*v.vector, data_exchange, Insert);

      last_action = Insert;

      return *this;
    }



    void
    Vector::import_nonlocal_data_for_fe (const TrilinosWrappers::SparseMatrix &m,
                                         const Vector                         &v)
    {
      Assert (m.trilinos_matrix().isFillComplete() == true,
              ExcMessage ("Matrix is not compressed. "
                          "Cannot find exchange information!"));
      Assert (v.vector->getMap().UniqueGIDs() == true,
              ExcMessage ("The input vector has overlapping data, "
                          "which is not allowed."));

      if (vector->getMap().SameAs(m.trilinos_matrix().getColMap()) == false)
        {
          vector.reset (new vector_type(
                          m.trilinos_matrix().getColMap()
                        ));
        }

      import_type data_exchange (vector->getMap(), v.vector->getMap());
      const int ierr = vector->doImport(*v.vector, data_exchange, Insert);

      AssertThrow (ierr == 0, ExcTrilinosError(ierr));

      last_action = Insert;
    }

  } /* end of namespace MPI */




  Vector::Vector ()
  {
    last_action = Zero;
    map_type map (0, 0, Utilities::Trilinos::comm_self());
    vector.reset (new vector_type(map));
  }



  Vector::Vector (const size_type n)
  {
    last_action = Zero;
    map_type map ((TrilinosWrappers::types::int_type)n, 0, Utilities::Trilinos::comm_self());
    vector.reset (new vector_type (map));
  }



  Vector::Vector (const map_type &input_map)
  {
    last_action = Zero;
    map_type map (n_global_elements(input_map),
                         input_map.IndexBase(),
                         input_map.Comm());
    vector.reset (new vector_type(map));
  }



  Vector::Vector (const IndexSet &partitioning,
                  const MPI_Comm &communicator)
  {
    last_action = Zero;
    map_type map (static_cast<TrilinosWrappers::types::int_type>(partitioning.size()),
                         0,
#ifdef DEAL_II_WITH_MPI
                         Tpetra_MpiComm(communicator));
#else
                         Tpetra_SerialComm());
    (void)communicator;
#endif
    vector.reset (new vector_type(map));
  }



  Vector::Vector (const VectorBase &v)
  {
    last_action = Zero;
    map_type map (n_global_elements(v.vector->getMap()),
                         v.vector->getMap().IndexBase(),
                         v.vector->getMap().Comm());
    vector.reset (new vector_type(map));

    if (vector->getMap().SameAs(v.vector->getMap()) == true)
      {
        const int ierr = vector->update(1.0, *v.vector, 0.0);
        AssertThrow (ierr == 0, ExcTrilinosError(ierr));
      }
    else
      reinit (v, false, true);

  }



  void
  Vector::reinit (const size_type n,
                  const bool      omit_zeroing_entries)
  {
    if (size() != n)
      {
        map_type map ((TrilinosWrappers::types::int_type)n, 0,
                             Utilities::Trilinos::comm_self());
        vector.reset (new vector_type (map));
      }
    else if (omit_zeroing_entries == false)
      {
        int ierr;
        ierr = vector->GlobalAssemble(last_action);
        (void)ierr;
        Assert (ierr == 0, ExcTrilinosError(ierr));

        ierr = vector->putScalar(0.0);
        Assert (ierr == 0, ExcTrilinosError(ierr));
      }

    last_action = Zero;
  }



  void
  Vector::reinit (const map_type &input_map,
                  const bool        omit_zeroing_entries)
  {
    if (n_global_elements(vector->getMap()) != n_global_elements(input_map))
      {
        map_type map (n_global_elements(input_map),
                             input_map.IndexBase(),
                             input_map.Comm());
        vector.reset (new vector_type (map));
      }
    else if (omit_zeroing_entries == false)
      {
        int ierr;
        ierr = vector->GlobalAssemble(last_action);
        (void)ierr;
        Assert (ierr == 0, ExcTrilinosError(ierr));

        ierr = vector->putScalar(0.0);
        Assert (ierr == 0, ExcTrilinosError(ierr));
      }

    last_action = Zero;
  }



  void
  Vector::reinit (const IndexSet &partitioning,
                  const MPI_Comm &communicator,
                  const bool      omit_zeroing_entries)
  {
    if (n_global_elements(vector->getMap()) !=
        static_cast<TrilinosWrappers::types::int_type>(partitioning.size()))
      {
        map_type map (static_cast<TrilinosWrappers::types::int_type>(partitioning.size()),
                             0,
#ifdef DEAL_II_WITH_MPI
                             Tpetra_MpiComm(communicator));
#else
                             Tpetra_SerialComm());
        (void)communicator;
#endif
        vector.reset (new vector_type(map));
      }
    else if (omit_zeroing_entries == false)
      {
        int ierr;
        ierr = vector->GlobalAssemble(last_action);
        (void)ierr;
        Assert (ierr == 0, ExcTrilinosError(ierr));

        ierr = vector->putScalar(0.0);
        Assert (ierr == 0, ExcTrilinosError(ierr));
      }

    last_action = Zero;
  }



  void
  Vector::reinit (const VectorBase &v,
                  const bool        omit_zeroing_entries,
                  const bool        allow_different_maps)
  {
    // In case we do not allow to
    // have different maps, this
    // call means that we have to
    // reset the vector. So clear
    // the vector, initialize our
    // map with the map in v, and
    // generate the vector.
    (void)omit_zeroing_entries;
    if (allow_different_maps == false)
      {
        if (local_range() != v.local_range())
          {
            map_type map (global_length(*(v.vector)),
                                 v.vector->getMap().IndexBase(),
                                 v.vector->Comm());
            vector.reset (new vector_type(map));
          }
        else
          {
            int ierr;
            Assert (vector->getMap().SameAs(v.vector->getMap()) == true,
                    ExcMessage ("The Tpetra maps in the assignment operator ="
                                " do not match, even though the local_range "
                                " seems to be the same. Check vector setup!"));

            ierr = vector->GlobalAssemble(last_action);
            (void)ierr;
            Assert (ierr == 0, ExcTrilinosError(ierr));

            ierr = vector->putScalar(0.0);
            Assert (ierr == 0, ExcTrilinosError(ierr));
          }
        last_action = Zero;
      }

    // Otherwise, we have to check
    // that the two vectors are
    // already of the same size,
    // create an object for the data
    // exchange and then insert all
    // the data.
    else
      {
        Assert (omit_zeroing_entries == false,
                ExcMessage ("It is not possible to exchange data with the "
                            "option 'omit_zeroing_entries' set, which would not write "
                            "elements."));

        AssertThrow (size() == v.size(),
                     ExcDimensionMismatch (size(), v.size()));

        import_type data_exchange (vector->getMap(), v.vector->getMap());
 
        vector->doImport(*v.vector, data_exchange, Insert);

        last_action = Insert;
      }

  }



  Vector &
  Vector::operator = (const MPI::Vector &v)
  {
    if (size() != v.size())
      {
        map_type map (n_global_elements(v.vector->getMap()),
                             v.vector->getMap().IndexBase(),
                             v.vector->Comm());

        vector.reset (new vector_type(map));
      }

    reinit (v, false, true);
    return *this;
  }



  Vector &
  Vector::operator = (const Vector &v)
  {
    if (size() != v.size())
      {
        map_type map (n_global_elements(v.vector->getMap()),
                             v.vector->getMap().IndexBase(),
                             v.vector->Comm());
        vector.reset (new vector_type(map));
      }

     vector->update(1.0, *v.vector, 0.0);

     return *this;
  }

}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_TRILINOS
