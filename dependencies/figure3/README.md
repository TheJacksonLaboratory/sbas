# Starting in a new Jupyter Lab session in CloudOS
For this step of the analysis, start a fairly large instance if running from scratch, say about an r2.8xlarge (8 cpus with 61 GB RAM).
After your Notebook Session has initilased in CloudOS, open a terminal type the following command:


Postpone opening an `.ipynb` file until you have completed all teh following steps, to make sure your Notebook environment is the updated one, after the installations have been successfully completed.


```bash
git clone https://github.com/TheJacksonLaboratory/sbas.git
```

Follow the instructions to configure the git config, your GitHub associated email and name.

# `cd` into the `sbas` directory

```bash

cd sbas

```

# Initialise your environment

Open a terminal and type the following command:

```bash

conda init

```

This will prepare your environment to be able to use `conda`. After the command has been executed, close this terminal and open a fresh one.
You might be prompted to do so by `conda` as well.


# Update the `base` conda environment

```bash

cd sbas
conda env update --name base --file dependencies/figure3/requirements.txt
conda activate base

```

# Ready!

Your environment is ready now, you can open your  `.ipynb` file and start working on your Jupyter Notebook.
