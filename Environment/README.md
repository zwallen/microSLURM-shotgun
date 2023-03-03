This directory contains a `.yml` file that can be used for building an environment with `conda` with all the necessary packages for running the `microSLURM_metagenome`. To create the environment (titled `microSLURM_metagenome`) use the following command:
```
$ conda env create -f microSLURM_metagenome.yml
```
and once the environment is built, activate it with
```
$ conda activate microSLURM_metagenome
```
and make sure you move the `resources` directory found in this directory to the `microSLURM_metagenome/bin` directory as these reference files are needed when running BBDuk, and that the MetaPhlAn databases are downloaded in the `microSLURM_metagenome/lib/python3.7/site-packages/metaphlan/metaphlan_databases/` directory. If MetaPhlAn databases are not downloaded, then run
```
$ metaphlan --install
```

If this environment is built, then `conda activate microSLURM_metagenome` should suffice for the `-p` parameter in the `microSLURM_metagenome` scripts.
