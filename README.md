<p align="center">
  <img src="header_short_course.png" />
</p>

# A warm welcome to the `fdaPDE` world!

`fdaPDE` is a library for *physics-informed statistical learning of spatial and functional data*, at the intersection of statistics and numerical analysis, designed for data located over *complex multidimensional domains*,
 ranging from irregular planar regions and curved surfaces to linear networks and volumes. The use of Partial Differential Equations (PDEs) allows incorporating information derived from the *physics* of 
 the problem under study into the *statistical modeling*, making `fdaPDE` an extremely flexible tool for the analysis of complex data.

`fdaPDE` offers a wide range of modeling capabilities -- including regression, nonparametric density estimation, functional data analysis, and more -- for data located over a spatial domain, possibly evolving over time.

# Getting Started

`fdaPDE` is a C++ library with an interface to **R**, one of the most widely used languages for data analysis. 
It runs on all major platforms (**Linux**, **macOS**, and **Windows**) and has been available on CRAN for nearly 10 years.
We are currently developing a new version of the library, available [here](https://github.com/fdaPDE/fdaPDE-R), 
which introduces major improvements but also comes with stricter system requirements. 
Specifically, it must be built using **GCC version 15**, which is not readily available in some popular Linux distributions (e.g., **Ubuntu**) and is difficult to install on **macOS**.
We are actively working to relax these requirements to make the build process easier. Meanwhile, we provide a **Docker** image 
with a fully configured environment that includes the development version of the package.

---

## Recommended Setup

* **macOS** and **Linux** users: we strongly recommend using the **Docker** image, which already contains the development version of the `fdaPDE` package along with other useful packages.

* **Windows** users: 

  * If your machine supports **R version $\geq$ 4.5.0** (for reference, R version 4.5.0 was released on April 11, 2025), install it along with `RTools45` (download available [here](https://cran.r-project.org/bin/windows/Rtools/rtools45/rtools.html)). 
    This is the preferred approach, if compatible. Once R version $\geq$ 4.5.0 and `RTools45` are installed, install the library from GitHub. In the R console, type:

    ```
    install.packages("devtools")
    devtools::install_github("https://github.com/fdaPDE/fdaPDE-R")
    ``` 
    
    Additionally, we recommend installing the following R packages:

    ```
    install.packages(c("fields", "ggmap", "ggplot2", "latex2exp", "leafsync", 
                        "mapview", "patchwork", "R6", "Rcpp", "RcppEigen",
                        "RTriangle", "raster", "sf", "viridis"))

    ```

  * Alternatively, use the **Docker** image (this requires Windows Subsystem for Linux -- see instructions below). 
    


---

## Data 

The data for this short course can be downloaded from [this link](https://polimi365-my.sharepoint.com/:f:/g/personal/10539238_polimi_it/Ev1ynJ1q9DdNvmMZ0HUhafABeAjj1_Rezl7MqU3Jbdbieg?e=DyVDTc). 
Please store the downloaded content inside the folder associated with this repository. You can either clone this GitHub repository or download it as a `.zip` file and extract it.
Your working directory should look as follows:

```
📁 IWSM25_short_course
    📁 data
    📁 scripts
    📁 utils
    📁 vignettes
```

The `vignettes/` folder contains the vignettes that will be presented during the course, while the `scripts/` folder includes the same examples in plain R script format.

## Docker Image

If **Docker** is not yet installed on your machine, please refer to the **Installation** section below. To pull the **Docker** image, run the following command in a terminal:

```
docker pull aldoclemente/fdapde-docker:rstudio
```

To run a container, execute (replacing `/path/to/IWSM25_short_course` with the **absolute path** to your local copy of the course material):

```bash
docker run --rm -d -p 8787:8787 -v /path/to/IWSM25_short_course:/home/user/IWSM25_short_course --name rstudio -e PASSWORD=password aldoclemente/fdapde-docker:rstudio
```

This will launch an `RStudio Server` instance inside **Docker**. You can then access `RStudio` in your browser at [http://localhost:8787](http://localhost:8787).
We suggest using **Google Chrome** as browser (we encountered some small rendering issues using **Firefox**).   

Log in with:

* **Username**: `user`
* **Password**: `password`


Inside `RStudio`, set your working directory by running:

```
setwd("IWSM25_short_course/")
```

You have read/write permissions in this directory, so you can modify scripts, save results, etc. All changes are reflected on your local machine. You can safely close and reconnect to [http://localhost:8787](http://localhost:8787) at any time.

To stop the container when you're done:

```
docker stop rstudio
```

---

## Docker Installation

### Windows Users

#### Step 1: Enable Virtualization

1. Restart your computer and enter the **BIOS setup** (usually via F2, F10, Del, or Esc).
2. Enable **virtualization** (e.g., Intel VT-x or AMD-V).
3. Save and exit the BIOS.

#### Step 2: Enable WSL and Virtual Machine Platform

1. Open the **"Turn Windows features on or off"** menu (search via Start or Cortana).
2. Enable:

   * **Windows Subsystem for Linux**
   * **Virtual Machine Platform**

#### Step 3: Install Ubuntu via WSL

1. Open **Command Prompt** or **PowerShell** as Administrator.

2. Run:
    ```
    wsl --install
    ```

    By default, the **Ubuntu** distribution gets installed in Windows Subsystem for Linux (WSL2)

3. Follow the instructions to set up your user account.


#### Step 4: Install Docker Desktop

1. Download and install **Docker Desktop** from [this link](https://docs.docker.com/desktop/setup/install/windows-install/).
   
    ❗ The latest version of Docker Desktop may not be compatible with your machine. 
        If you encounter any issues, please refer to [this page](https://docs.docker.com/desktop/release-notes/#400) 
        and install a version of **Docker Desktop** that is compatible with your system.

2. Log in or create a **Docker Hub** account;
3. Go to **Settings > Resources > WSL Integration** and enable integration for Ubuntu;
4. Consider increasing the resources allocated to **Docker Desktop**. Specifically, you can increase the 
   maximum number of CPUs and the amount of RAM in the **Advanced** settings. 
5. Restart your machine.

#### Step 5: Test Docker

1. Open your WSL terminal by running from **PowerShell**:

    ```
    wsl
    ```

2. Run:

    ```
    docker --version
    ```

3. You can now pull the **Docker** image as described above, and run the **Docker** container inside WSL.
---

### macOS Users

1. Download and install **Docker Desktop** from [this link](https://docs.docker.com/desktop/setup/install/mac-install/).

    ❗ The latest version of Docker Desktop may not be compatible with your machine. 
        If you encounter any issues, please refer to [this page](https://docs.docker.com/desktop/release-notes/#400) 
        and install a version of **Docker Desktop** that is compatible with your system.


2. Open **Docker Desktop** and grant the necessary permissions;
3. Log in or create a **Docker Hub** account;
4. Consider increasing the resources allocated to **Docker Desktop**. Specifically, you can increase the 
   maximum number of CPUs and the amount of RAM in the **Advanced** settings. 
5. Verify the **Docker** installation:

    ```
    docker --version
    ```
6. You can now pull the **Docker** image and run the container as described above.

---

### Ubuntu Users

1. Install **Docker**:

    ```
    sudo apt-get update
    sudo apt-get install -y docker.io
    ```

2. Check that **Docker** is running:

    ```
    systemctl status docker
    ```

    If **Docker** is inactive, start the service:

    ```
    sudo systemctl start docker
    ```

3. To avoid typing `sudo` whenever you run the `docker` command, add your username to the **Docker** group:

    ```
    sudo usermod -aG docker $USER
    ```

4. Restart your machine.

5. Log in to **Docker Hub**:

    ```
    docker login
    ```

6. You can now pull the **Docker** image and run the container as described above.

---

## Final Notes

* Ensure **Docker** is running before pulling or running any image.
* If you encounter issues, consult the [official Docker documentation](https://docs.docker.com/) for your platform.


**<h2 style="text-align: center;">You are ready to use the brand-new version of *fdaPDE* ! 😄 </h2>**
