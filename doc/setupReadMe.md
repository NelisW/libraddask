# Setup Procedure for libRadtranDask 

## Create the server account and software

1. Create a user account on the server from which libRadtran is served. This is a standard Linux admin task.  For the present case the account name is dgriffith.

1. Download libRadtran from `http://www.libradtran.org/doku.php?id=download` and follow the installation instructions. The instructions are repeated below.

    In this case version 2.0.3 is downloaded and installed.

1. Unzip the compressed tar file with 

    tar -xvf  libradtran-2.0.3.tar.gz

1. Ensure that you have Python 2.7.XXX

1. Compile the distribution

        cd libRadtran-2.0.3
        ./configure
        make

1. Test the program
    make check (make sure to use GNU make)




## Creating the environment

Create environment with 

    conda create --name mort python=3.7 anaconda

## Additional packages

Search for packages here: https://anaconda.org/conda-forge/

conda install  dill easygui ipyparallel easygui pyephem xarray

`paramiko` and `sphinx` are required in morticia, but installed with the Anaconda distribution.  `xray` is now renamed to `xarray` in Python. All the code has been updated to this effect; the package `xray` should no longer appear in morticia.

The packages required by morticia, what must be user installed are:

- dill Serialize all of python (almost)
- easygui EasyGUI is a module for very simple, very easy GUI programming in Python.
- ipyparallel  Interactive Parallel Computing with IPython
- pyephem Compute positions of the planets and stars
- xarray N-D labeled arrays and datasets in Python.

The  following packages may not be found on the standard anaconda source channels.
In this case use the install format (replace the channel with any other appropriate channel if this does not work):  

    conda install -c conda-forge easygui

- easygui
- pyephem


## Make Morticia known

The morticia library is not installed in the Python site-packages tree, so Python must know where to find it. To tell Python where to find this library create a file with the filename `morticia.pth` in the `site-packages` folder of the `mort` environment. On my PC the environment is here: 

    C:\ProgramData\Anaconda3\envs\mort\lib\site-packages\morticia.pth

On Linux it is in a location similar to this

    ~/anaconda2/envs/mordevpy37/lib/python3.7/site-packages/morticia.pth

The file must have only one line, and this line must be the morticia library location.    In my case the library was cloned from the repository into the following folder, and this must be the contents of the file:

    W:\Morticia\OSS-gitlab\morticia

or on the nimbus server (146.64.246.94) the contents must be

    /home/dgriffith/GitHub/morticia



## Prepare for Jupyter

Make the environment visible in a Jupyter notebook"

    conda install notebook ipykernel
    ipython kernel install

See [here])https://github.com/NelisW/ComputationalRadiometry/blob/master/00-Installing-Python-pyradi.ipynb) for more detail.


## Code Changes

1. Change all instances of `xray` to `xarray`

1. Change all print statement lines to print function lines.

1. Changed all np.linspace() to use integer numbers in third parameter

1. All local imports must be `from . import xxx`. This rippled through to hundreds of changes.

1. Updated Python 2 `file` type tests to Python 3 `io.IOBase` tests.

1. Numerous cases where integers were expected and floats were provided (Python 3 is stricter and does not do automatic conversions)





## Test Status


The following files are all the *.py and *.ipynb files.
Files status is as follows:
1. 'v'  no error, no warning, no further action required.
1. '+'  no error, with NEP18 warning, fix warning later.
1. '!'  with errors and other non-NEP18 warnings, more work required.
1. '.'  not yet fully tested, may be edited and could be working.

### *.ipynb files

    ! sensor/01a-EyeCTFExamples.ipynb [issue with xarray m attribute]
    + sensor/02a-OpticsMTFTutorials.ipynb
    + sensor/03a-LensTutorial.ipynb
    + sensor/04a-FocalPlaneArraySQE.ipynb
    + sensor/05a-Imager-Lens-and-Camera.ipynb


    + rad/01a-SpectralChannelsAndFilters.ipynb
    ! rad/02a-libRadtranTutorial.ipynb [major redesign of the code required]
    v rad/03a-libRadtranPyExplore.ipynb
    ! rad/04a-Radiant-Env-Maps-RadEnv.ipynb [libradtran cluster not available to test]
    ! rad/05a-libRadtran-ipyparallel.ipynb [libradtran cluster not available to test]
    ! rad/06a-Mitsuba-Introduction.ipynb [mitsuba not available to test]

    + tools/01a-xrayExercises.ipynb


# Work on nimbus (146.64.246.94) for user `dgriffith`

Home directory `~` is `/home/dgriffith`

## Bring nimbus to Python 3 status

Installed a python 3 environment for morticia: mordevpy37, with the morticia dependencies listed above.

Changed the folder name  `~/GitHub/morticia` to `~/GitHub/morticia2` to keep the old copy available.

Cloned the python 3 version of morticia into folder `~/GitHub/morticia`, from `https://github.com/NelisW/morticia`

To see where the `site-packages` folder on your anaconda installation are use this code:

    from distutils.sysconfig import get_python_lib
    print(get_python_lib())

Added `morticia.pth` to the prepared environment's `site-packages` folder, e.g.,  `~/anaconda2/envs/mordevpy37/lib/python3.7/site-packages/morticia.pth`. The contents of this file is the single line. 

    /home/dgriffith/GitHub/morticia

## libRadTran server

1. On the PC activate the Python 3 morticia environment 

1. Start up the VPN (if used) and find the IP number allocated to the VPN ethernet0. On Windows uUse `ipconfig` and look for something like the following (the Ethernet adapter connection number may be different, look for the DNS suffix that says `csir.co.za`). This should be an IP address on the CSIR network `146.64.xxx.xxx`:

        Ethernet adapter Local Area Connection 2:

        Connection-specific DNS Suffix  . : csir.co.za
        Link-local IPv6 Address . . . . . : fe80::8992:152a:c1f0:1260%27
        IPv4 Address. . . . . . . . . . . : 146.64.202.118
        Subnet Mask . . . . . . . . . . . : 255.255.0.0
        Default Gateway . . . . . . . . . :

    To find the IP address on Linux use `ip addr show`.


1. On the PC start the Dask scheduler, using the VPN IP address

      dask-scheduler --host 146.64.202.118

1. On nimbus activate the Python 3 morticia environment 

1. On nimbus cd to the libRadTran bin folder  `/home/dgriffith/libRadtran/libRadtran-2.0.2/bin/`

1. On nimbus start the Dask scheworker process, using the VPN IP address and the Dask port number

    dask-worker tcp://146.64.202.118:8786

1. Open the notebook and set the scheduler using the VPN host address


https://docs.dask.org/en/latest/  
https://docs.dask.org/en/latest/futures.html  
https://github.com/dask/dask-tutorial/blob/master/05_distributed.ipynb

if Remote:
    serverusername = 'tcp://146.64.202.118:8786'
    client = Client(serverusername)


#     cluster = SSHCluster(
#         ["146.64.246.94", "146.64.246.94", "146.64.246.94", "146.64.246.94"],
#         connect_options={"known_hosts": None},
#         worker_options={"nthreads": 2},
#         scheduler_options={"port": 0, "dashboard_address": ":8797"}
#         )

else:
    # Setup a local cluster.
    # By default this sets up 1 worker per core    
    client = Client() # use local host

#client.get_versions(check=True)

client



dask-scheduler --host 146.64.246.94 --port 8786 --bokeh-port <open-port>

dask-worker --host 172.16.4.30 146.64.246.94:8786 --worker-port 8786
    
 https://stackoverflow.com/questions/51354166/using-dask-in-ec2-instances-throws-couldnt-gather-1-keys/51355260   