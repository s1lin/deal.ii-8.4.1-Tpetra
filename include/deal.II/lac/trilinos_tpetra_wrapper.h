// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by Shilei Lin
//
// Some comments are adopting from Tpetra Example/Turtorials.
//
// ---------------------------------------------------------------------
#ifndef dealii__trilinos_tpetra_wrapper_h
#define dealii__trilinos_tpetra_wrapper_h


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


using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::arcp;
using Teuchos::ArrayRCP;
using Teuchos::ArrayView;
using Teuchos::Time;
using Teuchos::TimeMonitor;
using Teuchos::tuple;

//^(?![ \t]*//).*Epetra
//Epetra    		=	Tpetra
//-------------------------------------------

//General
//Comm() 			= 	getComm()
//SameAs()			=	isSameAs()
//-------------------------------------------

//Map:
//MinAllGID()    	= 	getMinAllGlobalIndex ()
//MaxAllGID()    	= 	getMaxAllGlobalIndex ()
//MinMyGID() 		= 	getMinGlobalIndex ()
//MaxMyGID() 		= 	getMaxGlobalIndex ()
//MinLID() 			= 	getMinLocalIndex()
//MaxLID() 			= 	getMaxLocalIndex ()
//NumMyElements() 	= 	getNodeNumElements ()
//NumGlobalElements() = getGlobalNumElements ()
//IndexBase() 		= 	getIndexBase()
//DistributedGlobal() = isDistributed()
//LinearMap() 		 = 	isContiguous ()
//LID() 			=   getLocalElement() *****
//GID() 			=   getGlobalElement() *****
//MyGlobalElements()= 	getNodeNumElements () *****
//MyGID()			=	isNodeGlobalElement()
//MyLID()			=	isNodeLocalElement()
//IsOneToOne()      =   isOneToOne()
//-------------------------------------------

//MPI
//NumProc()			=	getSize ()
//MyPID()			=	getRank	()
//-------------------------------------------

//CrsGraph
//RangeMap()		=	getRangeMap()
//RowMap()			=   getRowMap()
//DomainMap()		=	getDomainMap()
//GlobalLength()    =   getGlobalLength()
//ColMap()			=	getColMap()
//graph()			=	getCrsGraph()
//FillComplete()	=	fillComplete()
//Filled()			=	isFillComplete()
//IndexBase()		=	getIndexBase ()
//LRID()			=
//GRID()			=
//MaxNumIndices()   =   getNodeMaxNumRowEntries()
//NumMyRows () 		=	getNodeNumRows ()
//NumMyIndices()    =   getNumEntriesInLocalRow()
//NumGlobalCols ()  =	getGlobalNumCols()
//NumGlobalIndices() = getNumEntriesInGlobalRow ()
//GlobalMaxNumIndices() = getGlobalMaxNumRowEntries()
//ReplaceGlobalValues() = replaceGlobalValues()
//-------------------------------------------

//CrsMatrix
//RangeMap()		=	getRangeMap()
//RowMap()			=   getRowMap()
//DomainMap()		=	getDomainMap()
//GlobalLength()    =   getGlobalLength()
//ColMap()			=	getColMap()
//OperatorRangeMap()  =	getRangeMap()
//OperatorDomainMap() =	getDomainMap()
//FillComplete()	=	fillComplete()
//Filled()			=	isFillComplete()
//IndexBase()		=	getIndexBase ()
//LRID()			=
//GRID()			=
//GlobalAssmeble()  =
//MaxNumEntries()	=	getNodeMaxNumRowEntries()
//NumMyRows () 		=	getNodeNumRows ()
//NumGlobalCols ()  =	getGlobalNumCols()
//ReplaceGlobalValues() = replaceGlobalValues()
//NumMyEntries()	=	getNumEntriesInLocalRow ()
//ExtractGlobalRowCopy() =
//-------------------------------------------

//Operator
//OperatorRangeMap()  =	getRangeMap()
//OperatorDomainMap() =	getDomainMap()
//Apply()			  = Apply()
//ApplyInverse()	  = Apply()
//SetUseTranspose()   =

//-------------------------------------------

//Solver:
//reset()			=	setProblem()
//TrueResidual()	=	achievedTol()
//-------------------------------------------

//Vector
//PutScalar()		=	putScalar()
//Map()				=	getMap()
//ExtractView()		=
//Norm1() = norm1()
//Norm2() = norm2()
//MinValue() = 
//MaxValue() = 
//MeanValue() = meanValue()
//Update() = update()
//PutScalar() = putScalar()
//NumGlobalRows()   =

//-------------------------------------------
// MueLu uses Teuchos::RCP which is Trilinos version of std::shared_ptr.

typedef double                                      SC;
typedef int                                         LO;
typedef int                                         GO;
typedef KokkosClassic::DefaultNode::DefaultNodeType NT;


typedef Tpetra::Map<SC, LO, GO>                     map_type;
typedef Tpetra::Vector<SC, LO,	GO, NT>             vector_type;
typedef Tpetra::RowMatrix<SC, LO, GO, NT>           row_matrix_type;
typedef Tpetra::Export<SC, LO, GO>                  export_type;
typedef Tpetra::Import<SC, LO, GO>                  import_type;
typedef Tpetra::CrsMatrix<SC, LO, GO, NT>           crs_matrix_type;
typedef Tpetra::CrsGraph<SC, LO, GO>                crs_graph_type;
typedef Tpetra::MultiVector<SC, LO, GO, NT>         multi_vector_type;
typedef Tpetra::Operator<SC, LO, GO, NT>            operator_type;

//"reference-counted pointer."
//typedef RCP<const map_type>                         rcprowMap_type;
//typedef RCP<const map_type>                         rcpcolMap_type;
//typedef RCP<const crs_matrix_type>                  rcpcrsMatrix_type;
typedef MueLu::MLParameterListInterpreter<SC, LO, GO, NT> ML_ParameterListInterpreter;

typedef Teuchos::RCP<Teuchos::Comm<int>> comm_type;

typedef Tpetra::MpiPlatform<int>    Tpetra_MpiComm;
typedef Tpetra::SerialPlatform<int> Tpetra_SerialComm;

typedef Xpetra::CrsMatrix<SC, LO, GO, NT> xcrs_matrix_type;
typedef Xpetra::Matrix<SC, LO, GO, NT> xmatrix_type;
typedef Xpetra::CrsMatrixWrap<SC, LO, GO, NT> xcrs_matrix_wrap;

#endif


