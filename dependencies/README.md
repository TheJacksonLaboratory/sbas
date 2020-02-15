# Starting in a new Jupyter Lab session in CloudOS
For this paper, start a fairly large instance if running from scratch, say about an r2.8xlarge (8 cpus with 61 GB RAM).
After your Notebook Session has initilased in CloudOS, open a terminal type the following command:

```bash
git clone https://github.com/TheJacksonLaboratory/lifebitCloudOSDRE.git

```

Follow the instructions to configure the git config, your GitHub associated email and name.

# `cd` into the `lifebitCloudOSDRE` directory

```bash

cd lifebitCloudOSDRE
git checkout cristina-jupyter

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

cd lifebitCloudOSDRE
conda env update --name base --file dependencies/requirements.txt 
conda activate base
```

#  Install `TheJacksonLaboratory/yarn`

Along with some other dependencies:


```bash

Rscript dependencies/install.R 

```


# Ready!

Your environment is ready now, you can start working on the Jupyter Notebook.
