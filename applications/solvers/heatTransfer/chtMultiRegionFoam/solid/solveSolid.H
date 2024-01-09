{
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        tmp<fvScalarMatrix> TEqn
        (
            fvm::ddt(rho*cp, T)
          - fvm::laplacian(Kappa, T)
        );
        TEqn.ref().relax();
        solve(TEqn);
    }

    Info<< "Min/max T:" << min(T) << ' ' << max(T) << endl;
}
