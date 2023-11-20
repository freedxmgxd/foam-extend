/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     5.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Description
    Selects a section of mesh based on a cellSet.

    The utility sub-sets the mesh to choose only a part of interest. Check
    the setSet/cellSet utilities to see how to select cells based on various.

    The mesh will subset all points, faces and cells needed to make a sub-mesh
    but will not preserve attached boundary types.

\*---------------------------------------------------------------------------*/

#include "fvMeshSubset.H"
#include "argList.H"
#include "cellSet.H"
#include "IOobjectList.H"
#include "volFields.H"

#include "labelIOField.H"
#include "scalarIOField.H"
#include "vectorIOField.H"
#include "sphericalTensorIOField.H"
#include "symmTensorIOField.H"
#include "tensorIOField.H"
#include "lagrangianFieldDecomposer.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


template<class Type>
void subsetVolFields
(
    const fvMeshSubset& meshSubset,
    const wordList& fieldNames,
    PtrList<GeometricField<Type, fvPatchField, volMesh> >& subFields
)
{
    const fvMesh& baseMesh = meshSubset.baseMesh();

    forAll (fieldNames, i)
    {
        const word& fieldName = fieldNames[i];

        Info<< "Subsetting volume field " << fieldName << endl;

        GeometricField<Type, fvPatchField, volMesh> fld
        (
            IOobject
            (
                fieldName,
                baseMesh.time().timeName(),
                baseMesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            baseMesh
        );

        subFields.set(i, meshSubset.interpolate(fld));
    }
}


template<class Type>
void subsetSurfaceFields
(
    const fvMeshSubset& meshSubset,
    const wordList& fieldNames,
    PtrList<GeometricField<Type, fvsPatchField, surfaceMesh> >& subFields
)
{
    const fvMesh& baseMesh = meshSubset.baseMesh();

    forAll (fieldNames, i)
    {
        const word& fieldName = fieldNames[i];

        Info<< "Subsetting surface field " << fieldName << endl;

        GeometricField<Type, fvsPatchField, surfaceMesh> fld
        (
            IOobject
            (
                fieldName,
                baseMesh.time().timeName(),
                baseMesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            baseMesh
        );

        subFields.set(i, meshSubset.interpolate(fld));
    }
}


template<class Type>
void subsetPointFields
(
    const fvMeshSubset& meshSubset,
    const pointMesh& pMesh,
    const wordList& fieldNames,
    PtrList<GeometricField<Type, pointPatchField, pointMesh> >& subFields
)
{
    const fvMesh& baseMesh = meshSubset.baseMesh();

    forAll (fieldNames, i)
    {
        const word& fieldName = fieldNames[i];

        Info<< "Subsetting point field " << fieldName << endl;

        GeometricField<Type, pointPatchField, pointMesh> fld
        (
            IOobject
            (
                fieldName,
                baseMesh.time().timeName(),
                baseMesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            pMesh
        );

        subFields.set(i, meshSubset.interpolate(fld));
    }
}


// Main program:

int main(int argc, char *argv[])
{
    argList::validArgs.append("set");
    argList::validOptions.insert("patch", "patch name");
    argList::validOptions.insert("overwrite", "");

#   include "setRootCase.H"
#   include "createTime.H"
    runTime.functionObjects().off();
#   include "createMesh.H"
    const word oldInstance = mesh.pointsInstance();

    word setName(args.additionalArgs()[0]);
    bool overwrite = args.optionFound("overwrite");


    Info<< "Reading cell set from " << setName << endl << endl;

    // Create mesh subsetting engine
    fvMeshSubset meshSubset
    (
        IOobject
        (
            "set",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );


    label patchI = -1;

    if (args.optionFound("patch"))
    {
        word patchName(args.option("patch"));

        patchI = mesh.boundaryMesh().findPatchID(patchName);

        if (patchI == -1)
        {
            FatalErrorIn(args.executable()) << "Illegal patch " << patchName
                << nl << "Valid patches are " << mesh.boundaryMesh().names()
                << exit(FatalError);
        }

        Info<< "Adding exposed internal faces to patch " << patchName << endl
            << endl;
    }
    else
    {
        Info<< "Adding exposed internal faces to a patch called"
            << " \"oldInternalFaces\" (created if necessary)" << endl
            << endl;
    }


    cellSet currentSet(mesh, setName);

    meshSubset.setLargeCellSubset(currentSet, patchI, true);

    IOobjectList objects(mesh, runTime.timeName());


    // Read vol fields and subset
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~
    Info<< "Volume fields" << nl << endl;

    wordList scalarNames(objects.names(volScalarField::typeName));
    PtrList<volScalarField> scalarFlds(scalarNames.size());
    subsetVolFields(meshSubset, scalarNames, scalarFlds);

    wordList vectorNames(objects.names(volVectorField::typeName));
    PtrList<volVectorField> vectorFlds(vectorNames.size());
    subsetVolFields(meshSubset, vectorNames, vectorFlds);

    wordList sphericalTensorNames
    (
        objects.names(volSphericalTensorField::typeName)
    );
    PtrList<volSphericalTensorField> sphericalTensorFlds
    (
        sphericalTensorNames.size()
    );
    subsetVolFields(meshSubset, sphericalTensorNames, sphericalTensorFlds);

    wordList symmTensorNames(objects.names(volSymmTensorField::typeName));
    PtrList<volSymmTensorField> symmTensorFlds(symmTensorNames.size());
    subsetVolFields(meshSubset, symmTensorNames, symmTensorFlds);

    wordList tensorNames(objects.names(volTensorField::typeName));
    PtrList<volTensorField> tensorFlds(tensorNames.size());
    subsetVolFields(meshSubset, tensorNames, tensorFlds);


    // Read surface fields and subset
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Info<< "Surface fields" << nl << endl;

    wordList surfScalarNames(objects.names(surfaceScalarField::typeName));
    PtrList<surfaceScalarField> surfScalarFlds(surfScalarNames.size());
    subsetSurfaceFields(meshSubset, surfScalarNames, surfScalarFlds);

    wordList surfVectorNames(objects.names(surfaceVectorField::typeName));
    PtrList<surfaceVectorField> surfVectorFlds(surfVectorNames.size());
    subsetSurfaceFields(meshSubset, surfVectorNames, surfVectorFlds);

    wordList surfSphericalTensorNames
    (
        objects.names(surfaceSphericalTensorField::typeName)
    );
    PtrList<surfaceSphericalTensorField> surfSphericalTensorFlds
    (
        surfSphericalTensorNames.size()
    );
    subsetSurfaceFields
    (
        meshSubset,
        surfSphericalTensorNames,
        surfSphericalTensorFlds
    );

    wordList surfSymmTensorNames
    (
        objects.names(surfaceSymmTensorField::typeName)
    );
    PtrList<surfaceSymmTensorField> surfSymmTensorFlds
    (
        surfSymmTensorNames.size()
    );
    subsetSurfaceFields(meshSubset, surfSymmTensorNames, surfSymmTensorFlds);

    wordList surfTensorNames(objects.names(surfaceTensorField::typeName));
    PtrList<surfaceTensorField> surfTensorFlds(surfTensorNames.size());
    subsetSurfaceFields(meshSubset, surfTensorNames, surfTensorFlds);


    // Read point fields and subset
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Info<< "Point fields" << nl << endl;

    const pointMesh& pMesh = pointMesh::New(mesh);

    wordList pointScalarNames(objects.names(pointScalarField::typeName));
    PtrList<pointScalarField> pointScalarFlds(pointScalarNames.size());
    subsetPointFields(meshSubset, pMesh, pointScalarNames, pointScalarFlds);

    wordList pointVectorNames(objects.names(pointVectorField::typeName));
    PtrList<pointVectorField> pointVectorFlds(pointVectorNames.size());
    subsetPointFields(meshSubset, pMesh, pointVectorNames, pointVectorFlds);

    wordList pointSphericalTensorNames
    (
        objects.names(pointSphericalTensorField::typeName)
    );
    PtrList<pointSphericalTensorField> pointSphericalTensorFlds
    (
        pointSphericalTensorNames.size()
    );
    subsetPointFields
    (
        meshSubset,
        pMesh,
        pointSphericalTensorNames,
        pointSphericalTensorFlds
    );

    wordList pointSymmTensorNames
    (
        objects.names(pointSymmTensorField::typeName)
    );
    PtrList<pointSymmTensorField> pointSymmTensorFlds
    (
        pointSymmTensorNames.size()
    );
    subsetPointFields
    (
        meshSubset,
        pMesh,
        pointSymmTensorNames,
        pointSymmTensorFlds
    );

    wordList pointTensorNames(objects.names(pointTensorField::typeName));
    PtrList<pointTensorField> pointTensorFlds(pointTensorNames.size());
    subsetPointFields(meshSubset, pMesh, pointTensorNames, pointTensorFlds);

    // Lagrangian fields
    fileNameList cloudDirs
    (
        readDir(runTime.timePath()/cloud::prefix, fileName::DIRECTORY)
    );
    Info<< "READ " << cloudDirs.size() << " clouds" << endl;
    // Particles
    PtrList<Cloud<indexedParticle> > lagrangianPositions(cloudDirs.size());
    // Particles per cell
    PtrList< List<SLList<indexedParticle*>*> > cellParticles(cloudDirs.size());

    PtrList<PtrList<labelIOField> > lagrangianLabelFields(cloudDirs.size());
    PtrList<PtrList<scalarIOField> > lagrangianScalarFields(cloudDirs.size());
    PtrList<PtrList<vectorIOField> > lagrangianVectorFields(cloudDirs.size());
    PtrList<PtrList<sphericalTensorIOField> > lagrangianSphericalTensorFields
    (
        cloudDirs.size()
    );
    PtrList<PtrList<symmTensorIOField> > lagrangianSymmTensorFields
    (
        cloudDirs.size()
    );
    PtrList<PtrList<tensorIOField> > lagrangianTensorFields
    (
        cloudDirs.size()
    );

    label cloudI = 0;

    forAll (cloudDirs, i)
    {
        IOobjectList sprayObjs
        (
            mesh,
            runTime.timeName(),
            cloud::prefix/cloudDirs[i]
        );

        IOobject* positionsPtr = sprayObjs.lookup("positions");

        if (positionsPtr)
        {
            // Read lagrangian particles
            // ~~~~~~~~~~~~~~~~~~~~~~~~~

            Info<< "Identified lagrangian data set: " << cloudDirs[i] << endl;

            lagrangianPositions.set
            (
                cloudI,
                new Cloud<indexedParticle>
                (
                    mesh,
                    cloudDirs[i],
                    false
                )
            );


            // Sort particles per cell
            // ~~~~~~~~~~~~~~~~~~~~~~~

            cellParticles.set
            (
                cloudI,
                new List<SLList<indexedParticle*>*>
                (
                    mesh.nCells(),
                    static_cast<SLList<indexedParticle*>*>(nullptr)
                )
            );

            label i = 0;

            forAllIter
            (
                Cloud<indexedParticle>,
                lagrangianPositions[cloudI],
                iter
            )
            {
                iter().index() = i++;

                label celli = iter().cell();

                // Check
                if (celli < 0 || celli >= mesh.nCells())
                {
                    FatalErrorIn(args.executable())
                        << "Illegal cell number " << celli
                        << " for particle with index " << iter().index()
                        << " at position " << iter().position() << nl
                        << "Cell number should be between 0 and "
                        << mesh.nCells()-1 << nl
                        << "On this mesh the particle should be in cell "
                        << mesh.findCell(iter().position())
                        << exit(FatalError);
                }

                if (!cellParticles[cloudI][celli])
                {
                    cellParticles[cloudI][celli] =
                        new SLList<indexedParticle*>();
                }

                cellParticles[cloudI][celli]->append(&iter());
            }

            // Read fields
            // ~~~~~~~~~~~

            IOobjectList lagrangianObjects
            (
                mesh,
                runTime.timeName(),
                cloud::prefix/cloudDirs[cloudI]
            );

            lagrangianFieldDecomposer::readFields
            (
                cloudI,
                lagrangianObjects,
                lagrangianLabelFields
            );

            lagrangianFieldDecomposer::readFields
            (
                cloudI,
                lagrangianObjects,
                lagrangianScalarFields
            );

            lagrangianFieldDecomposer::readFields
            (
                cloudI,
                lagrangianObjects,
                lagrangianVectorFields
            );

            lagrangianFieldDecomposer::readFields
            (
                cloudI,
                lagrangianObjects,
                lagrangianSphericalTensorFields
            );

            lagrangianFieldDecomposer::readFields
            (
                cloudI,
                lagrangianObjects,
                lagrangianSymmTensorFields
            );

            lagrangianFieldDecomposer::readFields
            (
                cloudI,
                lagrangianObjects,
                lagrangianTensorFields
            );

            cloudI++;
        }
    }

    lagrangianPositions.setSize(cloudI);
    cellParticles.setSize(cloudI);
    lagrangianLabelFields.setSize(cloudI);
    lagrangianScalarFields.setSize(cloudI);
    lagrangianVectorFields.setSize(cloudI);
    lagrangianSphericalTensorFields.setSize(cloudI);
    lagrangianSymmTensorFields.setSize(cloudI);
    lagrangianTensorFields.setSize(cloudI);
    

    // Write mesh and fields to new time
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (!overwrite)
    {
        runTime++;
    }
    else
    {
        meshSubset.subMesh().setInstance(oldInstance);
    }

    Info<< "Writing subsetted mesh and fields to time " << runTime.timeName()
        << endl;
    meshSubset.subMesh().write();


    // Write out mesh mapping information
    labelIOList pointMap
    (
        IOobject
        (
            "pointMap",
            meshSubset.subMesh().dbDir(),
            polyMesh::meshSubDir,
            meshSubset.subMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        meshSubset.pointMap()
    );
    pointMap.write();

    labelIOList faceMap
    (
        IOobject
        (
            "faceMap",
            meshSubset.subMesh().dbDir(),
            polyMesh::meshSubDir,
            meshSubset.subMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        meshSubset.faceMap()
    );
    faceMap.write();

    labelIOList cellMap
    (
        IOobject
        (
            "cellMap",
            meshSubset.subMesh().dbDir(),
            polyMesh::meshSubDir,
            meshSubset.subMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        meshSubset.cellMap()
    );
    cellMap.write();

    labelIOList patchMap
    (
        IOobject
        (
            "patchMap",
            meshSubset.subMesh().dbDir(),
            polyMesh::meshSubDir,
            meshSubset.subMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        meshSubset.patchMap()
    );
    patchMap.write();


    // Subsetting adds 'subset' prefix. Rename field to be like original.
    forAll (scalarFlds, i)
    {
        scalarFlds[i].rename(scalarNames[i]);

        scalarFlds[i].write();
    }
    forAll (vectorFlds, i)
    {
        vectorFlds[i].rename(vectorNames[i]);

        vectorFlds[i].write();
    }
    forAll (sphericalTensorFlds, i)
    {
        sphericalTensorFlds[i].rename(sphericalTensorNames[i]);

        sphericalTensorFlds[i].write();
    }
    forAll (symmTensorFlds, i)
    {
        symmTensorFlds[i].rename(symmTensorNames[i]);

        symmTensorFlds[i].write();
    }
    forAll (tensorFlds, i)
    {
        tensorFlds[i].rename(tensorNames[i]);

        tensorFlds[i].write();
    }

    // Surface ones.
    forAll (surfScalarFlds, i)
    {
        surfScalarFlds[i].rename(surfScalarNames[i]);

        surfScalarFlds[i].write();
    }
    forAll (surfVectorFlds, i)
    {
        surfVectorFlds[i].rename(surfVectorNames[i]);

        surfVectorFlds[i].write();
    }
    forAll (surfSphericalTensorFlds, i)
    {
        surfSphericalTensorFlds[i].rename(surfSphericalTensorNames[i]);

        surfSphericalTensorFlds[i].write();
    }
    forAll (surfSymmTensorFlds, i)
    {
        surfSymmTensorFlds[i].rename(surfSymmTensorNames[i]);

        surfSymmTensorFlds[i].write();
    }
    forAll (surfTensorNames, i)
    {
        surfTensorFlds[i].rename(surfTensorNames[i]);

        surfTensorFlds[i].write();
    }

    // Point ones
    forAll (pointScalarFlds, i)
    {
        pointScalarFlds[i].rename(pointScalarNames[i]);

        pointScalarFlds[i].write();
    }
    forAll (pointVectorFlds, i)
    {
        pointVectorFlds[i].rename(pointVectorNames[i]);

        pointVectorFlds[i].write();
    }
    forAll (pointSphericalTensorFlds, i)
    {
        pointSphericalTensorFlds[i].rename(pointSphericalTensorNames[i]);

        pointSphericalTensorFlds[i].write();
    }
    forAll (pointSymmTensorFlds, i)
    {
        pointSymmTensorFlds[i].rename(pointSymmTensorNames[i]);

        pointSymmTensorFlds[i].write();
    }
    forAll (pointTensorNames, i)
    {
        pointTensorFlds[i].rename(pointTensorNames[i]);

        pointTensorFlds[i].write();
    }


    // If there is lagrangian data write it out
    Info<< "Lagrangian data" << nl << endl;
    
    forAll (lagrangianPositions, cloudI)
    {
        if (lagrangianPositions[cloudI].size())
        {
            lagrangianFieldDecomposer fieldDecomposer
            (
                mesh,
                meshSubset.subMesh(),
                cellMap,
                cloudDirs[cloudI],
                lagrangianPositions[cloudI],
                cellParticles[cloudI]
            );

            // Lagrangian fields
            if
            (
                lagrangianLabelFields[cloudI].size()
             || lagrangianScalarFields[cloudI].size()
             || lagrangianVectorFields[cloudI].size()
             || lagrangianSphericalTensorFields[cloudI].size()
             || lagrangianSymmTensorFields[cloudI].size()
             || lagrangianTensorFields[cloudI].size()
            )
            {
                fieldDecomposer.decomposeFields
                (
                    cloudDirs[cloudI],
                    lagrangianLabelFields[cloudI]
                );
                fieldDecomposer.decomposeFields
                (
                    cloudDirs[cloudI],
                    lagrangianScalarFields[cloudI]
                );
                fieldDecomposer.decomposeFields
                (
                    cloudDirs[cloudI],
                    lagrangianVectorFields[cloudI]
                );
                fieldDecomposer.decomposeFields
                (
                    cloudDirs[cloudI],
                    lagrangianSphericalTensorFields[cloudI]
                );
                fieldDecomposer.decomposeFields
                (
                    cloudDirs[cloudI],
                    lagrangianSymmTensorFields[cloudI]
                );
                fieldDecomposer.decomposeFields
                (
                    cloudDirs[cloudI],
                    lagrangianTensorFields[cloudI]
                );
            }
        }
    }

    Info << nl << "End" << endl;

    return 0;
}


// ************************************************************************* //
