// ************************************************************************* //
//                            avtOpenFOAMFileFormat.h                           //
// ************************************************************************* //

//THIS READER IS FOR DOUBLE PRECISION ONLY
//MUST DO A REOPEN UPON CHANGING TIMESTEPS TO
//REPOPULATE THE VECTORS/SCALARS LISTS

#ifndef AVT_OpenFOAM_FILE_FORMAT_H
#define AVT_OpenFOAM_FILE_FORMAT_H

#include <avtMTMDFileFormat.h>

#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <ostream>

#include <vtkUnstructuredGrid.h>
#include <vtkIntArray.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkPoints.h>
#include <vtkCellData.h>
#include <vtkVertex.h>
#include <vtkHexahedron.h>
#include <vtkWedge.h>
#include <vtkPyramid.h>
#include <vtkTetra.h>
#include <vtkConvexPointSet.h>
#include <vtkTriangle.h>
#include <vtkQuad.h>
#include <vtkPolygon.h>
#include <vtkUnstructuredGrid.h>
#include <vtkObjectFactory.h>
#include <vtkDirectory.h>


// ****************************************************************************
//  Class: avtOpenFOAMFileFormat
//
//  Purpose:
//      Reads in OpenFOAM files as a plugin to VisIt.
//
//  Programmer: root -- generated by xml2avt
//  Creation:   Wed Jun 7 16:01:15 PST 2006
//
// ****************************************************************************
typedef struct
{
  int faceIndex;
  bool neighborFace;
}face;

class avtOpenFOAMFileFormat : public avtMTMDFileFormat
{
  public:
  avtOpenFOAMFileFormat(const char *);
  virtual ~avtOpenFOAMFileFormat() {;};
  virtual int GetNTimesteps(void);
    virtual void           GetTimes(std::vector<double> &times);
  virtual const char  *GetType(void)   { return "OpenFOAM"; };
  virtual void  FreeUpResources(void); 
  virtual vtkDataSet  *GetMesh(int, int, const char *);
  virtual vtkDataArray *GetVar(int, int, const char *);
  virtual vtkDataArray *GetVectorVar(int, int, const char *);
  virtual bool          HasInvariantMetaData(void) const { return false; }; //Number of variables dynamic
  virtual bool          HasInvariantSIL(void) const      { return false; };  //Number of Domains dynamic
  virtual void ActivateTimestep(int);


  protected:
  virtual void           PopulateDatabaseMetaData(avtDatabaseMetaData *, int);

  //Members
  bool CreateFaces;
  bool FirstVar;
  bool FirstVectorVar;
  std::string Path;
  std::string PathPrefix;
  int NumberOfTimeSteps;
  int StartFace;
  int NFaces;
  double * Steps;
  vtkFloatArray * TempData;
  std::vector< std::string > PolyMeshPointsDir;
  std::vector< std::string > PolyMeshFacesDir;
  vtkPoints * Points;
  vtkIdType NumFaces;
  vtkIdType NumPoints;
  vtkIntArray * FaceOwner;
  vtkIntArray * FaceNeighbor;
  vtkIdType NumCells;
  std::vector< std::vector<int> > FacePoints;
  std::vector< std::vector<int> > FacesOwnerCell;
  std::vector< std::vector<int> > FacesNeighborCell;
  std::vector< std::vector<face> > FacesOfCell;
  int NumBlocks;
  int NumBoundaries;
  int NumPointZones;
  int NumFaceZones;
  int NumCellZones;
  //create temporary vectors
  std::vector< std::string > BoundaryNames;
  std::vector< std::string > PointZoneNames;
  std::vector< std::string > FaceZoneNames;
  std::vector< std::string > CellZoneNames;
  std::vector< std::string > ScalarNames;
  std::vector< std::string > VectorNames;

  //Methods
  double ControlDictDataParser(std::string);  //Parser ControlDict Entries
  void ReadControlDict();  //Read the ControlDict File
  void ReadFacesFile(std::string);  //Read the faces into a vector
  void GetPoints(int);  //Read the Points File
  void ReadOwnerFile(std::string);  //read the owner faces into a vector
  void ReadNeighborFile(std::string);  //read the neighbor faces into a vector
  void CombineOwnerNeigbor();  //Create a vector of cell faces
  vtkUnstructuredGrid * MakeInternalMesh();  //calls the functions to create internal mesh
  void PopulatePolyMeshDirArrays();  //Creates a vector that tells you at what time step the points and faces file reside
  std::string GetDataType(std::string, std::string);  //Parses the files to quickly find out what type of variables you have for the meta data
  vtkFloatArray * GetInternalVariableAtTimestep( std::string, int);  //Returns the values a requested variable for the internal mesh
  vtkFloatArray * GetBoundaryVariableAtTimestep(int, std::string, int);  //Returns the values a requested variable for the boundary meshed
  std::vector< std::string > GatherBlocks(std::string, int);  //creates a vector of all the blocks in a region
  vtkUnstructuredGrid * GetBoundaryMesh(int, int);  //returns a requested boundary mesh
  vtkUnstructuredGrid * GetPointZoneMesh(int, int);  //returns a requested point zone mesh
  vtkUnstructuredGrid * GetCellZoneMesh(int, int);  //returns a requested cell zone mesh
  vtkUnstructuredGrid * GetFaceZoneMesh(int, int);  //returns a requested face zone mesh
};


#endif
