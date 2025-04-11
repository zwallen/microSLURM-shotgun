This directory contains a `.yml` file that can be used for building an environment with `conda` with all the necessary packages for running the `microSLURM_shotgun`. To create the environment (titled `microSLURM_shotgun`) use the following command:
```
$ conda env create -f microSLURM_shotgun.yml
```
and once the environment is built, activate it with
```
$ conda activate microSLURM_shotgun
```
and make sure you move the `resources` directory found in this directory to the `microSLURM_shotgun/bin` directory as these reference files are needed when running BBDuk, and that the MetaPhlAn databases are downloaded in the `microSLURM_shotgun/lib/python3.7/site-packages/metaphlan/metaphlan_databases/` directory. If MetaPhlAn databases are not downloaded, then run
```
$ metaphlan --install
```

If this environment is built, then `conda activate microSLURM_shotgun` should suffice for the `-p` parameter in the `microSLURM_shotgun` scripts.
