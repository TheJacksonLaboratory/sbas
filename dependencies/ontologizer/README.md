# Running Ontologizer in a Jupyter Notebook session (no sudo)

You can use `Ontologizer` in an environment with no sudo privileges by installing the programme by source. Follow the steps below to install `Ontologizer` in a JupyterLab Notebook Session where you have no sudo:

## 1. Activate the base conda environment

```bash
conda init
```

and close the terminal for your changes to take effect.

```bash
conda activate base
```

## 2. Install `Java` and the `bash kernel` for Jupyter via conda install

```bash
conda install java-jdk bash_kernel -y
```

## 3. Follow the Ontologizer instructions to get the .jar file

Read instructions from Ontologizer docs [here](http://ontologizer.de).
Download the ontologizer.jar file from [here](http://ontologizer.de/cmdline/Ontologizer.jar) by copying link.

```bash
wget <ontologizer-link>
```

This will now be accessible only from the currend directory.

## 4. Test that you have Ontologizer installed by typing the --help

```bash
# or the full path if you are in a different directory in JupyterLab /mnt/shared/gcp-user/session_data/Ontologizer.jar
java -jar Ontologizer.jar --help
```
