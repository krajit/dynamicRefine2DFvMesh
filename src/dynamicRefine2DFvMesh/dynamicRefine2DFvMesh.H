/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM-v1806.

Class
    Foam::dynamicRefine2DFvMesh

SourceFiles
    dynamicRefine2DFvMesh.C

Author
    Luca Cornolti - written for foam-extend

Modified for OpenFOAM-v1806 by
    Ajit Kumar, Shiv Nadar University, India
    


Description
    A fvMesh with built-in refinement for 2D configuration based on dynamicRefineFvMesh.
    Determines which cells to refine/unrefine and does all in update().


\*---------------------------------------------------------------------------*/

#ifndef dynamicRefine2DFvMesh_H
#define dynamicRefine2DFvMesh_H

#include "dynamicFvMesh.H"
#include "hexRef2D.H"
#include "PackedBoolList.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class dynamicRefine2DFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicRefine2DFvMesh
    : public dynamicFvMesh
{
  protected:
    //- Mesh cutting engine
    hexRef2D meshCutter_;

    //- Dump cellLevel for postprocessing
    Switch dumpLevel_;

    //- Fluxes to map
    List<Pair<word>> correctFluxes_;

    //- Number of refinement/unrefinement steps done so far.
    label nRefinementIterations_;

    //- Protected cells (usually since not hexes)
    PackedBoolList protectedCell_;

    // Private Member Functions

    //- Count set/unset elements in packedlist.
    static label count(const PackedBoolList &, const unsigned int);

    //- Calculate cells that cannot be refined since would trigger
    //  refinement of protectedCell_ (since 2:1 refinement cascade)
    void calculateProtectedCells(PackedBoolList &unrefineableCell) const;

    //- Read the projection parameters from dictionary
    void readDict();

    //- Refine cells. Update mesh and fields.
    autoPtr<mapPolyMesh> refine(const labelList &);

    //- Unrefine cells. Gets passed in centre points of cells to combine.
    autoPtr<mapPolyMesh> unrefine(const labelList &);

    // Selection of cells to un/refine

    //- Calculates approximate value for refinement level so
    //  we don't go above maxCell
    scalar getRefineLevel(
        const label maxCells,
        const label maxRefinement,
        const scalar refineLevel,
        const scalarField &) const;

    //- Get per cell max of connected point
    scalarField maxPointField(const scalarField &) const;

    //- Get point min of connected cell
    scalarField minCellField(const volScalarField &) const;

    scalarField cellToPoint(const scalarField &vFld) const;

    scalarField error(
        const scalarField &fld,
        const scalar minLevel,
        const scalar maxLevel) const;

    //- Select candidate cells for refinement
    virtual void selectRefineCandidates(
        const scalar lowerRefineLevel,
        const scalar upperRefineLevel,
        const scalarField &vFld,
        PackedBoolList &candidateCell) const;

    //- Subset candidate cells for refinement
    virtual labelList selectRefineCells(
        const label maxCells,
        const label maxRefinement,
        const PackedBoolList &candidateCell) const;

    //- Select edges that can be unrefined (2D case).
    virtual labelList selectUnrefineEdges(
        const scalar unrefineLevel,
        const PackedBoolList &markedCell,
        const scalarField &pFld) const;

    //- Extend markedCell with cell-face-cell.
    void extendMarkedCells(PackedBoolList &markedCell) const;

  private:
    //- Disallow default bitwise copy construct
    dynamicRefine2DFvMesh(const dynamicRefine2DFvMesh &);

    //- Disallow default bitwise assignment
    void operator=(const dynamicRefine2DFvMesh &);

  public:
    //- Runtime type information
    TypeName("dynamicRefine2DFvMesh");

    // Constructors

    //- Construct from IOobject
    explicit dynamicRefine2DFvMesh(const IOobject &io);

    // Destructor

    virtual ~dynamicRefine2DFvMesh();

    // Member Functions

    //- Direct access to the refinement engine
    const hexRef2D &meshCutter() const
    {
        return meshCutter_;
    }

    //- Cells which should not be refined/unrefined
    const PackedBoolList &protectedCell() const
    {
        return protectedCell_;
    }

    //- Cells which should not be refined/unrefined
    PackedBoolList &protectedCell()
    {
        return protectedCell_;
    }

    //- Update the mesh for both mesh motion and topology change
    virtual bool update();

    // Writing

    //- Write using given format, version and compression
    virtual bool writeObject(
        IOstreamOption streamOpt,
        const bool writeOnProc
    ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
