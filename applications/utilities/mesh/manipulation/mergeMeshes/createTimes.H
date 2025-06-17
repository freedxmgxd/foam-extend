    Info<< nl << "Create Times" << endl;

    Time runTimeMaster
    (
        Time::controlDictName,
        rootDirMaster,
        caseDirMaster,
        "system",
        "constant",
        false  // Disable function objects
    );

    Time runTimeToAdd
    (
        Time::controlDictName,
        rootDirToAdd,
        caseDirToAdd,
        "system",
        "constant",
        false  // Disable function objects
    );
