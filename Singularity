BootStrap: debootstrap
DistType "debian"
MirrorURL: http://us.archive.ubuntu.com/ubuntu/
OSVersion: xenial

%labels

    AUTHOR Thomas Cokelaer
    # For version 1.0.0, the container needs 2.5Go

%post

    apt-get install -y wget
    apt-get install -y bzip2
    apt-get install -y vim

    # install anaconda
    if [ ! -d /usr/local/anaconda ]; then
        #wget https://repo.continuum.io/miniconda/Miniconda3-4.3.14-Linux-x86_64.sh\
        # for now, we use 4.2.12 to have python3.5 by default so no need to
        # create a new env saving space in the process. The reason for using 3.5
        # is inherent to the packages used at the moment.
        wget https://repo.continuum.io/miniconda/Miniconda3-4.2.12-Linux-x86_64.sh\
           -O ~/anaconda.sh && \
        bash ~/anaconda.sh -b -p /usr/local/anaconda && \
        rm ~/anaconda.sh
    fi

    # set anaconda path
    export PATH=$PATH:/usr/local/anaconda/bin
    conda update conda

    conda config --add channels r
    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels bioconda

    # The main packages for sequana:
    conda install --file https://raw.githubusercontent.com/CancerRxGene/gdsctools/master/requirements.txt

    # Let us save some space
    conda clean --packages -y

    # Sequana source code
    pip install gdsctools==0.20.1

    conda clean --all -y # next requires lots of space
    rm -rf /usr/local/anaconda/pkgs

    if [ ! -d /data ]; then mkdir /data; fi
    if [ ! -d /scripts ]; then mkdir /scripts; fi
    if [ ! -d /scratch ]; then mkdir /scratch; fi
    if [ ! -d /mounting ]; then mkdir /mounting; fi
    if [ ! -d /pasteur ]; then mkdir /pasteur; fi

%environment
    export PATH=$PATH:/usr/local/anaconda/bin
    echo "backend:agg" > matplotlibrc

