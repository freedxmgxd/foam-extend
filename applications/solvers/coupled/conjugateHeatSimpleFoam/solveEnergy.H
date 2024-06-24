{
    while (simpleSolid.correctNonOrthogonal())
    {
        coupledFvScalarMatrix TEqns(2);

        // Fluid side
        rhoCpPhi = rho*Cp*phi;

        fvScalarMatrix* TFluidEqn = new fvScalarMatrix
        (
            fvm::div(rhoCpPhi, T)
          - fvm::laplacian(kappaEff, T)
         ==
            radiation->Ru()
          - fvm::Sp(4.0*radiation->Rp()*pow3(T), T)
          + 3.0*radiation->Rp()*pow4(T)
        );

        if (thermalModelFluid.active())
        {
            fvScalarMatrix& TFEqn = *TFluidEqn;
            TFEqn -= thermalModelFluid.S();
        }

        TFluidEqn->relax();

        // Solid side
        fvScalarMatrix* TSolidEqn = new fvScalarMatrix
        (
          - fvm::laplacian(kSolidf, Tsolid, "laplacian(k,T)")
        );

        if (solidThermalModel.active())
        {
            fvScalarMatrix& TSEqn = *TSolidEqn;
            TSEqn -= solidThermalModel.source().S();
        }

        TSolidEqn->relax();

        // Add fluid equation
        TEqns.set(0, TFluidEqn);

        // Add solid equation
        TEqns.set(1, TSolidEqn);

        TEqns.solve();
    }
}
