# Instructions

## Required

```bash
conda install conda-build anaconda-client
```

## Building and pushing to https://anaconda.org/williamgilpin

```bash
conda build . --no-anaconda-upload
PACKAGE_OUTPUT=`conda build . --output`
anaconda login
anaconda upload --user williamgilpin $PACKAGE_OUTPUT
conda build purge
anaconda logout
```

## Install

```
conda install -c williamgilpin pypdb
```

## Additional Info
https://docs.anaconda.com/anaconda-cloud/user-guide/tasks/work-with-packages
